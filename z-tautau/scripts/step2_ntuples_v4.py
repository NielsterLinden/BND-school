#!/usr/bin/env python
"""Step 2 (v4) -- v4 skims -> flat ntuples of the lepton channels (docs/10-v4-plan.md sections 2-3, 10).

    python scripts/step2_ntuples_v4.py                          # every sample with v4 skims, all channels
    python scripts/step2_ntuples_v4.py --samples DY_NLO --channels mutau --max-files 2

Channels: mutau (SingleMuon, IsoMu24), etau (SingleElectron, Ele27), emu (MuonEG, cross triggers) and the
trigger-efficiency side samples emu_mu (SingleMuon-triggered e mu events) and emu_el (SingleElectron-
triggered). Per event with a selected pair: MET filters, the trigger of the channel, the pair choice
(leptons.select_ltau / select_emu). The nominal tau energy scale is applied to genuine tau_h in the simulation
(propagated to the MET) before the pair choice. Everything else (isolation, charge, vetoes, m_T, D_zeta, the
tau_h working point) is stored and cut in analysis_v4.regions, so the sidebands and the systematic variations
that move events across the cuts are all available downstream.

Output: $BND_TAUTAU_CACHE/ntuples_v4/<sample>_<channel>.root (tree `ntuple`) + .meta.json (cutflow, generator sums).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import multiprocessing as mp
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import config, corrections, io, leptons, mass, samples, skim  # noqa: E402

CHANNELS = ["mutau", "etau", "emu", "emu_mu", "emu_el"]
STREAM_OF = {"mutau": ("SingleMuon", skim.CAT_MT), "etau": ("SingleElectron", skim.CAT_ET), "emu": ("MuonEG", skim.CAT_EM),
             "emu_mu": ("SingleMuon", skim.CAT_EM_MU), "emu_el": ("SingleElectron", skim.CAT_EM_EL)}
TRIGGERS_OF = {"mutau": config.MU_TRIGGERS, "etau": config.EL_TRIGGERS, "emu": config.EMU_TRIGGERS,
               "emu_mu": config.MU_TRIGGERS, "emu_el": config.EL_TRIGGERS}
READ_COMMON = (["run", "luminosityBlock", "event", "PV_npvsGood", "MET_pt", "MET_phi", "MET_covXX", "MET_covXY", "MET_covYY",
                "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY", "PuppiMET_pt", "skim_cat"]
               + config.MET_FILTERS + config.MU_TRIGGERS + config.EL_TRIGGERS + config.EMU_TRIGGERS
               + ["nTau"] + [f"Tau_{f}" for f in leptons.TAU_FIELDS]
               + ["nMuon"] + [f"Muon_{f}" for f in leptons.MU_FIELDS]
               + ["nElectron"] + [f"Electron_{f}" for f in leptons.EL_FIELDS]
               + ["nJet", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_btagDeepFlavB",
                  "nTrigObj", "TrigObj_pt", "TrigObj_eta", "TrigObj_phi", "TrigObj_id", "TrigObj_filterBits"])
READ_DATA = ["era"]
READ_MC = ["Tau_genPartFlav", "Muon_genPartFlav", "Electron_genPartFlav", "Jet_hadronFlavour", "genWeight", "Pileup_nTrueInt",
           "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Dn", "GenMET_pt", "gen_lhe_flavour",
           "gen_mll_lhe", "gen_fid", "gen_n_vistau", "gen_vis1_pt", "gen_vis2_pt", "LHE_Njets", "LHE_NpNLO", "LHE_Vpt",
           "gen_top_pt", "gen_antitop_pt"]
VECTORS = {"LHEScaleWeight": 9, "PSWeight": 4, "LHEPdfWeight": 103}


def _np(x, fill=0.0, dtype=np.float32):
    return ak.to_numpy(ak.fill_none(x, fill)).astype(dtype)


def _lepton_fields(out, prefix, lep, is_el: bool):
    out[f"{prefix}_pt"] = _np(lep.pt)
    out[f"{prefix}_eta"] = _np(lep.eta)
    out[f"{prefix}_phi"] = _np(lep.phi)
    out[f"{prefix}_mass"] = _np(lep.mass)
    out[f"{prefix}_charge"] = _np(lep.charge, 0, np.int8)
    out[f"{prefix}_iso"] = _np(lep.pfRelIso03_all if is_el else lep.pfRelIso04_all, 99.0)
    out[f"{prefix}_genflav"] = _np(lep.genflav, -1, np.int8)
    if is_el:
        out[f"{prefix}_sceta"] = _np(lep.sceta)
        e = out[f"{prefix}_pt"] * np.cosh(out[f"{prefix}_eta"])
        for v in ("dEscaleUp", "dEscaleDown", "dEsigmaUp", "dEsigmaDown"):
            out[f"{prefix}_{v}"] = (_np(lep[v]) / np.maximum(e, 1e-6)).astype(np.float32)   # relative shift
        out[f"{prefix}_wp80"] = _np(lep.mvaFall17V2noIso_WP80, False, bool)
        out[f"{prefix}_match_ele27"] = _np(lep.match_ele27, False, bool)
        out[f"{prefix}_match_emu12"] = _np(lep.match_emu12, False, bool)
        out[f"{prefix}_match_emu23"] = _np(lep.match_emu23, False, bool)
    else:
        out[f"{prefix}_match_iso24"] = _np(lep.match_iso24, False, bool)
        out[f"{prefix}_match_trkiso8"] = _np(lep.match_trkiso8, False, bool)
        out[f"{prefix}_match_trkiso23"] = _np(lep.match_trkiso23, False, bool)


def _tau_fields(out, t):
    out["t_pt"] = _np(t.pt)
    out["t_pt_raw"] = _np(t.pt_raw)
    out["t_eta"] = _np(t.eta)
    out["t_phi"] = _np(t.phi)
    out["t_mass"] = _np(t.mass)
    out["t_charge"] = _np(t.charge, 0, np.int8)
    out["t_dm"] = _np(t.decayMode, -1, np.int8)
    out["t_vsjet"] = _np(t.idDeepTau2017v2p1VSjet, 0, np.uint8)
    out["t_vse"] = _np(t.idDeepTau2017v2p1VSe, 0, np.uint8)
    out["t_vsmu"] = _np(t.idDeepTau2017v2p1VSmu, 0, np.uint8)
    out["t_rawvsjet"] = _np(t.rawDeepTau2017v2p1VSjet)
    out["t_genflav"] = _np(t.genflav, -1, np.int8)


def _event_fields(out, ev, is_mc):
    out.update(run=_np(ev.run, 0, np.uint32), lumi=_np(ev.luminosityBlock, 0, np.uint32), event=_np(ev.event, 0, np.uint64),
               npv=_np(ev.PV_npvsGood, 0, np.int16), skim_cat=_np(ev.skim_cat, 0, np.uint8),
               met_covxx=_np(ev.MET_covXX), met_covxy=_np(ev.MET_covXY), met_covyy=_np(ev.MET_covYY),
               met_uncl_dx=_np(ev.MET_MetUnclustEnUpDeltaX), met_uncl_dy=_np(ev.MET_MetUnclustEnUpDeltaY),
               puppimet_pt=_np(ev.PuppiMET_pt))
    for name in config.MU_TRIGGERS + config.EL_TRIGGERS + config.EMU_TRIGGERS:
        out[name] = _np(ev[name], False, bool) if name in ak.fields(ev) else np.zeros(len(ev), bool)
    if is_mc:
        for b in ("genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Dn",
                  "GenMET_pt", "gen_mll_lhe", "gen_vis1_pt", "gen_vis2_pt"):
            out[b] = _np(ev[b])
        out["gen_lhe_flavour"] = _np(ev.gen_lhe_flavour, 0, np.uint8)
        out["gen_fid"] = _np(ev.gen_fid, False, bool)
        out["gen_n_vistau"] = _np(ev.gen_n_vistau, 0, np.int8)
        for b, name in (("LHE_Njets", "lhe_njets"), ("LHE_NpNLO", "lhe_npnlo")):
            out[name] = _np(ev[b], 0, np.uint8) if b in ak.fields(ev) else np.zeros(len(ev), np.uint8)
        out["lhe_vpt"] = _np(ev.LHE_Vpt) if "LHE_Vpt" in ak.fields(ev) else np.zeros(len(ev), np.float32)
        for b in ("gen_top_pt", "gen_antitop_pt"):
            out[b] = _np(ev[b]) if b in ak.fields(ev) else np.zeros(len(ev), np.float32)
        for b, size in VECTORS.items():
            if b in ak.fields(ev):
                out[b] = ak.to_numpy(ak.fill_none(ak.pad_none(ev[b], size, axis=1, clip=True), 1.0)).astype(np.float32)
    else:
        out["era"] = _np(ev.era, 0, np.uint8)


def _met_after_tes(ev, objs):
    """MET with the nominal tau energy-scale shift of the given tau records propagated."""
    metx = _np(ev.MET_pt).astype(np.float64) * np.cos(_np(ev.MET_phi))
    mety = _np(ev.MET_pt).astype(np.float64) * np.sin(_np(ev.MET_phi))
    for t in objs:
        d = _np(t.pt) - _np(t.pt_raw)
        metx -= d * np.cos(_np(t.phi))
        mety -= d * np.sin(_np(t.phi))
    return metx, mety


def _kinematics(out, p1, p2, leptonic, metx, mety):
    ms = mass.all_masses(p1, p2, metx, mety, out["met_covxx"].astype(np.float64), out["met_covxy"].astype(np.float64),
                         out["met_covyy"].astype(np.float64), leptonic=leptonic)
    for k, v in ms.items():
        out[k] = v.astype(np.float32)
    out["met_x"], out["met_y"] = metx.astype(np.float32), mety.astype(np.float32)
    out["dr"] = leptons.delta_r(p1[1], p1[2], p2[1], p2[2]).astype(np.float32)
    out["dzeta"] = mass.dzeta(p1[0], p1[2], p2[0], p2[2], metx, mety).astype(np.float32)


def process_chunk(ev, is_mc, tes, channel, cut, stream_bit):
    n = len(ev)
    cut["skim"] += n
    ev = ak.with_field(ev, ak.zip({"pt": ev.TrigObj_pt, "eta": ev.TrigObj_eta, "phi": ev.TrigObj_phi,
                                   "id": ev.TrigObj_id, "filterBits": ev.TrigObj_filterBits}), "TrigObj")
    sel = io.met_filters(ev) & ((ak.to_numpy(ev.skim_cat) & stream_bit) > 0)
    sel &= io.trigger_or(ev, TRIGGERS_OF[channel])
    cut["trigger"] += int(sel.sum())
    mu = leptons.muons(ev, is_mc)
    el = leptons.electrons(ev, is_mc)
    if channel in ("mutau", "etau"):
        is_mu = channel == "mutau"
        taus = leptons.ltau_taus(ev, is_mc, tes, channel)
        lep, tau, has, n_mu, n_el = leptons.select_ltau(ev, mu if is_mu else el, taus, is_mu)
        sel &= has
        cut["pair"] += int(sel.sum())
        if not sel.any():
            return None
        ev, lep, tau = ev[sel], lep[sel], tau[sel]
        out = {}
        _event_fields(out, ev, is_mc)
        _lepton_fields(out, "l", lep, is_el=not is_mu)
        _tau_fields(out, tau)
        other_mu, other_el = leptons.count_veto(mu[sel], el[sel])
        out["n_veto_mu"] = (n_mu[sel] if is_mu else other_mu).astype(np.int8)
        out["n_veto_el"] = (other_el if is_mu else n_el[sel]).astype(np.int8)
        out["os"] = out["l_charge"] * out["t_charge"] < 0
        metx, mety = _met_after_tes(ev, [tau])
        p1 = tuple(out[f"l_{k}"].astype(np.float64) for k in ("pt", "eta", "phi", "mass"))
        p2 = tuple(out[f"t_{k}"].astype(np.float64) for k in ("pt", "eta", "phi", "mass"))
        _kinematics(out, p1, p2, (True, False), metx, mety)
        out.update(leptons.jets_v4(ev, [lep.eta, tau.eta], [lep.phi, tau.phi], 2))
        return out
    # e mu (signal channel or trigger side samples)
    m, e, has, n_mu, n_el = leptons.select_emu(ev, mu, el, channel)
    sel &= has
    cut["pair"] += int(sel.sum())
    if not sel.any():
        return None
    ev, m, e = ev[sel], m[sel], e[sel]
    out = {}
    _event_fields(out, ev, is_mc)
    _lepton_fields(out, "mu", m, is_el=False)
    _lepton_fields(out, "el", e, is_el=True)
    out["n_veto_mu"] = n_mu[sel].astype(np.int8)
    out["n_veto_el"] = n_el[sel].astype(np.int8)
    out["os"] = out["mu_charge"] * out["el_charge"] < 0
    # a tau_h candidate in the event (ttbar fake-factor determination region of the paper): the best one
    taus = leptons.ltau_taus(ev, is_mc, tes, "mutau")
    far = (leptons.delta_r(taus.eta, taus.phi, ak.fill_none(m.eta, 99.0), ak.fill_none(m.phi, 0.0)) > 0.5) & \
          (leptons.delta_r(taus.eta, taus.phi, ak.fill_none(e.eta, 99.0), ak.fill_none(e.phi, 0.0)) > 0.5)
    tc = taus[far]
    best = ak.firsts(tc[ak.argmax(tc.rawDeepTau2017v2p1VSjet, axis=1, keepdims=True)])
    out["has_tau"] = _np(~ak.is_none(best), False, bool)
    _tau_fields(out, best)
    metx, mety = _met_after_tes(ev, [])
    p1 = tuple(out[f"mu_{k}"].astype(np.float64) for k in ("pt", "eta", "phi", "mass"))
    p2 = tuple(out[f"el_{k}"].astype(np.float64) for k in ("pt", "eta", "phi", "mass"))
    _kinematics(out, p1, p2, (True, True), metx, mety)
    out.update(leptons.jets_v4(ev, [m.eta, e.eta], [m.phi, e.phi], 2))
    return out


def process_file(args):
    path, key, channels, max_events = args
    s = samples.SAMPLES[key]
    is_mc = s["is_mc"]
    tes = corrections.tes_nominal() if is_mc else None
    cuts = {ch: {"skim": 0, "trigger": 0, "pair": 0} for ch in channels}
    parts = {ch: [] for ch in channels}
    fh = uproot.open(path)
    if "Events" not in [k.split(";")[0] for k in fh.keys(recursive=False)]:
        return parts, cuts
    tree = fh["Events"]
    keys = set(tree.keys())
    branches = READ_COMMON + (READ_MC if is_mc else READ_DATA)
    if is_mc:
        branches += [b for b in VECTORS if b in keys]
    branches = [b for b in branches if b in keys]
    for ev in tree.iterate(branches, step_size="120 MB", entry_stop=max_events):
        for f in config.MET_FILTERS:
            if f not in ak.fields(ev):
                ev = ak.with_field(ev, np.ones(len(ev), dtype=bool), f)
        for ch in channels:
            res = process_chunk(ev, is_mc, tes, ch, cuts[ch], STREAM_OF[ch][1])
            if res is not None:
                parts[ch].append(res)
    return parts, cuts


def concat(parts):
    if not parts:
        return {}
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def channels_of(key: str, requested):
    s = samples.SAMPLES[key]
    chans = [ch for ch in CHANNELS if s["is_mc"] or STREAM_OF[ch][0] == s["stream"]]
    return [ch for ch in chans if ch in requested]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--channels", nargs="*", default=CHANNELS)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--max-events", type=int, default=None, help="per file (testing)")
    ap.add_argument("--workers", type=int, default=config.N_WORKERS)
    ap.add_argument("--skim-dir", type=Path, default=config.SKIM_DIR_V4)
    ap.add_argument("--out-dir", type=Path, default=config.NTUPLE_DIR_V4)
    args = ap.parse_args()
    keys = args.samples or [k for k in samples.V4_SKIM_KEYS if skim.skim_files(k, args.skim_dir)]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for key in keys:
        files = skim.skim_files(key, args.skim_dir)[: args.max_files]
        chans = channels_of(key, args.channels)
        if not files or not chans:
            print(f"[ntuple] {key}: no skims or no channel, skipped")
            continue
        t0 = time.time()
        with mp.get_context("spawn").Pool(args.workers) as pool:
            results = pool.map(process_file, [(str(p), key, chans, args.max_events) for p in files], chunksize=1)
        gens = skim.load_gensums(key, args.skim_dir) if samples.SAMPLES[key]["is_mc"] else None
        for ch in chans:
            parts = [p for r, _ in results for p in r[ch]]
            cut = {}
            for _, c in results:
                for k, v in c[ch].items():
                    cut[k] = cut.get(k, 0) + v
            data = concat(parts)
            prov = {"sample": key, "channel": ch, "n_skim_files": len(files), "cutflow": cut,
                    "n_events": int(len(data.get("run", [])))}
            if gens is not None:
                prov["gensums"] = {k: (v.tolist() if np.ndim(v) else float(v)) for k, v in gens.items()}
            out = args.out_dir / f"{key}_{ch}.root"
            tmp = out.with_name(out.name + ".tmp")
            with uproot.recreate(tmp, compression=uproot.ZSTD(5)) as f:
                if data:
                    f.mktree("ntuple", data)
            tmp.rename(out)
            (args.out_dir / f"{key}_{ch}.meta.json").write_text(json.dumps(prov, indent=1))
            print(f"[ntuple] {key} {ch}: {prov['n_events']:,} events, cutflow {cut} -> {out}", flush=True)
        print(f"[ntuple] {key}: {len(files)} files in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
