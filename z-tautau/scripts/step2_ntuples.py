#!/usr/bin/env python
"""Step 2 -- skims -> flat analysis ntuples, one per sample (docs/02-selection.md, docs/04-ditau-mass.md).

    python scripts/step2_ntuples.py                      # all samples that have skims
    python scripts/step2_ntuples.py --samples DY_NLO --max-files 2

Per event with a selected tau_h tau_h pair (any charge, any isolation above VVVLoose):
MET filters, pair choice (objects.best_pair), both legs trigger-matched, extra-lepton veto. The nominal
tau energy scale is applied to genuine tau_h in the simulation (and propagated to the MET) *before* the
pair choice. The ntuple stores the two taus, MET and its covariance, jets, the three di-tau masses
(m_vis, m_col, m_tt), and in the simulation the generator weights, pileup, prefiring and generator
truth. The final tau pT > 40 GeV cut, the isolation/charge regions and all weights are applied in
step 3/4, so systematic variations can move events across them.

Output: $BND_TAUTAU_CACHE/ntuples_v1/<sample>.root (tree `ntuple`) + <sample>.meta.json (cutflow, generator sums).
These files are the laptop bundle: everything downstream of this step reads only them.
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

from ztautau import config, corrections, io, mass, objects, samples, skim  # noqa: E402

READ = (["run", "luminosityBlock", "event", "PV_npvsGood", "MET_pt", "MET_phi", "MET_covXX", "MET_covXY",
         "MET_covYY", "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY", "PuppiMET_pt", "PuppiMET_phi"]
        + config.MET_FILTERS + ["nTau"] + [f"Tau_{f}" for f in objects.TAU_FIELDS]
        + ["nMuon", "Muon_pt", "Muon_eta", "Muon_dxy", "Muon_dz", "Muon_mediumId", "Muon_pfRelIso04_all",
           "nElectron", "Electron_pt", "Electron_eta", "Electron_dxy", "Electron_dz", "Electron_mvaFall17V2noIso_WP90",
           "Electron_pfRelIso03_all", "Electron_convVeto", "Electron_lostHits",
           "nJet", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_btagDeepFlavB",
           "nTrigObj", "TrigObj_pt", "TrigObj_eta", "TrigObj_phi", "TrigObj_id", "TrigObj_filterBits"])
READ_DATA = ["era"]
READ_MC = ["Tau_genPartFlav", "genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up",
           "L1PreFiringWeight_Dn", "GenMET_pt", "gen_lhe_flavour", "gen_mll_lhe", "gen_fid", "gen_n_vistau",
           "gen_vis1_pt", "gen_vis2_pt", "LHE_Njets", "LHE_NpNLO", "LHE_Vpt"]
VECTORS = {"LHEScaleWeight": 9, "PSWeight": 4, "LHEPdfWeight": 103}


def _np(x, fill=0.0, dtype=np.float32):
    return ak.to_numpy(ak.fill_none(x, fill)).astype(dtype)


def process_chunk(ev, is_mc, tes, cut):
    n = len(ev)
    cut["skim"] += n
    ev = ak.with_field(ev, ak.zip({"pt": ev.TrigObj_pt, "eta": ev.TrigObj_eta, "phi": ev.TrigObj_phi,
                                   "id": ev.TrigObj_id, "filterBits": ev.TrigObj_filterBits}), "TrigObj")
    taus = objects.tau_collection(ev, is_mc, tes)
    cands = objects.tau_candidates(taus)
    t1, t2, has = objects.best_pair(cands)
    sel = io.met_filters(ev)
    cut["met_filters"] += int(sel.sum())
    sel &= has
    cut["pair"] += int(sel.sum())
    sel &= objects.trigger_matched(ev, t1) & objects.trigger_matched(ev, t2)
    cut["trigger_match"] += int(sel.sum())
    sel &= objects.extra_lepton_veto(ev)
    cut["lepton_veto"] += int(sel.sum())
    if not sel.any():
        return None
    ev, t1, t2 = ev[sel], t1[sel], t2[sel]
    out = {"run": _np(ev.run, 0, np.uint32), "lumi": _np(ev.luminosityBlock, 0, np.uint32),
           "event": _np(ev.event, 0, np.uint64), "npv": _np(ev.PV_npvsGood, 0, np.int16)}
    for i, t in ((1, t1), (2, t2)):
        out[f"t{i}_pt"] = _np(t.pt)
        out[f"t{i}_pt_raw"] = _np(t.pt_raw)
        out[f"t{i}_eta"] = _np(t.eta)
        out[f"t{i}_phi"] = _np(t.phi)
        out[f"t{i}_mass"] = _np(t.mass)
        out[f"t{i}_dm"] = _np(t.decayMode, -1, np.int8)
        out[f"t{i}_charge"] = _np(t.charge, 0, np.int8)
        out[f"t{i}_vsjet"] = _np(t.idDeepTau2017v2p1VSjet, 0, np.uint8)
        out[f"t{i}_vse"] = _np(t.idDeepTau2017v2p1VSe, 0, np.uint8)
        out[f"t{i}_vsmu"] = _np(t.idDeepTau2017v2p1VSmu, 0, np.uint8)
        out[f"t{i}_rawvsjet"] = _np(t.rawDeepTau2017v2p1VSjet)
        out[f"t{i}_genflav"] = _np(t.genflav, -1, np.int8)
    out["os"] = (out["t1_charge"] * out["t2_charge"] < 0)
    # MET with the nominal tau energy scale propagated (only differs from NanoAOD in the simulation)
    metx = _np(ev.MET_pt).astype(np.float64) * np.cos(_np(ev.MET_phi))
    mety = _np(ev.MET_pt).astype(np.float64) * np.sin(_np(ev.MET_phi))
    for i in (1, 2):
        d = out[f"t{i}_pt"] - out[f"t{i}_pt_raw"]
        metx -= d * np.cos(out[f"t{i}_phi"])
        mety -= d * np.sin(out[f"t{i}_phi"])
    out.update(met_x=metx.astype(np.float32), met_y=mety.astype(np.float32),
               met_covxx=_np(ev.MET_covXX), met_covxy=_np(ev.MET_covXY), met_covyy=_np(ev.MET_covYY),
               met_uncl_dx=_np(ev.MET_MetUnclustEnUpDeltaX), met_uncl_dy=_np(ev.MET_MetUnclustEnUpDeltaY),
               puppimet_pt=_np(ev.PuppiMET_pt))
    nj, nb, j1 = objects.jets(ev, t1, t2)
    out.update(njets=nj.astype(np.int8), nbjets=nb.astype(np.int8), jet1_pt=j1.astype(np.float32))
    tt1 = tuple(out[f"t1_{k}"].astype(np.float64) for k in ("pt", "eta", "phi", "mass"))
    tt2 = tuple(out[f"t2_{k}"].astype(np.float64) for k in ("pt", "eta", "phi", "mass"))
    ms = mass.all_masses(tt1, tt2, metx, mety, out["met_covxx"].astype(np.float64),
                         out["met_covxy"].astype(np.float64), out["met_covyy"].astype(np.float64))
    for k, v in ms.items():
        out[k] = v.astype(np.float32)
    out["dr_tt"] = objects.delta_r(out["t1_eta"], out["t1_phi"], out["t2_eta"], out["t2_phi"]).astype(np.float32)
    if is_mc:
        for b in ("genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up",
                  "L1PreFiringWeight_Dn", "GenMET_pt", "gen_mll_lhe", "gen_vis1_pt", "gen_vis2_pt"):
            out[b] = _np(ev[b])
        out["gen_lhe_flavour"] = _np(ev.gen_lhe_flavour, 0, np.uint8)
        out["gen_fid"] = _np(ev.gen_fid, False, bool)
        out["gen_n_vistau"] = _np(ev.gen_n_vistau, 0, np.int8)
        # jet-binned stitching (analysis.dy_norm) needs the NLO parton multiplicity; old skims lack it
        out["lhe_njets"] = _np(ev.LHE_Njets, 0, np.uint8) if "LHE_Njets" in ak.fields(ev) else np.zeros(len(ev), np.uint8)
        out["lhe_npnlo"] = _np(ev.LHE_NpNLO, 0, np.uint8) if "LHE_NpNLO" in ak.fields(ev) else np.zeros(len(ev), np.uint8)
        out["lhe_vpt"] = _np(ev.LHE_Vpt, 0.0) if "LHE_Vpt" in ak.fields(ev) else np.zeros(len(ev), np.float32)
        for b, size in VECTORS.items():
            if b in ak.fields(ev):
                arr = ak.fill_none(ak.pad_none(ev[b], size, axis=1, clip=True), 1.0)
                out[b] = ak.to_numpy(arr).astype(np.float32)
    else:
        out["era"] = _np(ev.era, 0, np.uint8)
    return out


def process_file(args):
    path, key, max_events = args
    s = samples.SAMPLES[key]
    is_mc = s["is_mc"]
    tes = corrections.tes_nominal() if is_mc else None
    cut = {"skim": 0, "met_filters": 0, "pair": 0, "trigger_match": 0, "lepton_veto": 0}
    parts = []
    fh = uproot.open(path)
    if "Events" not in [k.split(";")[0] for k in fh.keys(recursive=False)]:
        return parts, cut                  # no event of this file passed the skim (only GenSums)
    tree = fh["Events"]
    keys = set(tree.keys())
    branches = READ + (READ_MC if is_mc else READ_DATA)
    if is_mc:
        branches += [b for b in VECTORS if b in keys]
    branches = [b for b in branches if b in keys]
    for ev in tree.iterate(branches, step_size="150 MB", entry_stop=max_events):
        for f in config.MET_FILTERS:
            if f not in ak.fields(ev):
                ev = ak.with_field(ev, np.ones(len(ev), dtype=bool), f)
        res = process_chunk(ev, is_mc, tes, cut)
        if res is not None:
            parts.append(res)
    return parts, cut


def concat(parts):
    if not parts:
        return {}
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--max-events", type=int, default=None, help="per file (testing)")
    ap.add_argument("--workers", type=int, default=config.N_WORKERS)
    ap.add_argument("--out-dir", type=Path, default=config.NTUPLE_DIR)
    args = ap.parse_args()
    keys = args.samples or [k for k in samples.SAMPLES if skim.skim_files(k)]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for key in keys:
        files = skim.skim_files(key)[: args.max_files]
        if not files:
            print(f"[ntuple] {key}: no skims, skipped")
            continue
        t0 = time.time()
        # spawn, not fork: after the first sample the parent holds uproot/threading locks, and forked
        # workers deadlock on them
        with mp.get_context("spawn").Pool(args.workers) as pool:
            results = pool.map(process_file, [(str(p), key, args.max_events) for p in files], chunksize=1)
        parts = [p for r, _ in results for p in r]
        cut = {}
        for _, c in results:
            for k, v in c.items():
                cut[k] = cut.get(k, 0) + v
        data = concat(parts)
        prov = {"sample": key, "n_skim_files": len(files), "n_parent_files": len(samples.sources(key)),
                "cutflow": cut, "n_events": int(len(data.get("run", []))), "config": {
                    "tau_pt_ntuple": config.TAU_PT_NTUPLE, "tau_eta": config.TAU_ETA_MAX,
                    "vse_bit": config.TAU_VSE_BIT, "vsmu_bit": config.TAU_VSMU_BIT}}
        if samples.SAMPLES[key]["is_mc"]:
            g = skim.load_gensums(key)
            prov["gensums"] = {k: (v.tolist() if np.ndim(v) else float(v)) for k, v in g.items()}
        out = args.out_dir / f"{key}.root"
        tmp = out.with_name(out.name + ".tmp")
        with uproot.recreate(tmp, compression=uproot.ZSTD(5)) as f:
            if data:
                f.mktree("ntuple", data)
        tmp.rename(out)
        (args.out_dir / f"{key}.meta.json").write_text(json.dumps(prov, indent=1))
        print(f"[ntuple] {key}: {len(files)} files, {prov['n_events']:,} events, cutflow {cut} "
              f"({time.time() - t0:.0f}s) -> {out}", flush=True)


if __name__ == "__main__":
    main()
