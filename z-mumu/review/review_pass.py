#!/usr/bin/env python
"""Review pass over the v2 skims: e-mu OS/SS with MC flavour composition, and PU/MET studies.

Outputs (scratchpad/review_pass.pkl): dict of numpy histograms keyed
  emu|<sample>|<os/ss>|<pp/np/all>|<var>          var in pt_el, met, njet, npv, mass, njet_clean
  emu_flav|<sample>|<os/ss>                        (flav_mu, flav_el) weighted counts, 16x16
  sr|<sample>|ntrue_npv                            2D (100, 60) genWeight*norm*prefiring, no PU weight (MC)
  sr|<sample>|ntrue_jet_met                        3D (100, 3, 30): njet_clean 0/1/>=2 x MET
  sr|<sample>|ntrue_jet_puppimet                   3D (100, 3, 30)
  sr|Data|npv, sr|Data|npv_jet_met (60,3,30), sr|Data|npv_jet_puppimet, sr|Data|era_npv (2,60)
"""
from __future__ import annotations
import json, os, pickle, sys, glob
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import awkward as ak
import uproot

ROOT = Path("/project/atlas/Users/nterlind/BND-school/z-mumu")
sys.path.insert(0, str(ROOT))
from zmumu import config, objects, regions, samples, pileup, tnp, histograms as H

OUT = ROOT / "output" / "v2"
SKIM = Path(os.environ["BND_SKIM_DIR"])
BR = ["run", "event", "PV_npvsGood", "MET_pt", "PuppiMET_pt", "HLT_IsoMu24", "HLT_IsoTkMu24", "Muon_*", "TrigObj_*",
      "Flag_*", "FsrPhoton_*", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Electron_*", "genWeight", "Pileup_nTrueInt",
      "L1PreFiringWeight_Nom", "skim_cat", "skim_prescale", "gen_lhe_flavour"]
MET_E = np.linspace(0, 150, 31)
NPV_E = np.arange(0, 61, 1.0)
NTRUE_E = np.arange(0, 101, 1.0)
VARS = {"pt_el": np.linspace(20, 200, 37), "met": MET_E, "njet": np.arange(0, 9, 1.0), "njet_clean": np.arange(0, 9, 1.0),
        "npv": np.linspace(0, 50, 51), "mass": np.linspace(60, 120, 13), "pt1": np.linspace(20, 200, 37)}
FLAVOUR_SAMPLE = {13: "DYmumu", 11: "DYee", 15: "DYtautau"}


def emu_both(ev, masks, matched, clean):
    tight, anti = masks["tight"], masks["anti"]
    el = regions.electron_mask(ev)
    n_t = ak.to_numpy(ak.sum(tight, axis=1)); n_a = ak.to_numpy(ak.sum(anti, axis=1)); n_e = ak.to_numpy(ak.sum(el, axis=1))
    idx = np.nonzero((n_t == 1) & (n_a == 0) & (n_e == 1) & clean)[0]
    if idx.size == 0:
        return None
    sub = ev[idx]; t = tight[idx]; e = el[idx]
    f = lambda a: ak.to_numpy(ak.flatten(a))
    mu_pt, mu_eta, mu_phi, mu_m, mu_q = f(sub.Muon_pt[t]), f(sub.Muon_eta[t]), f(sub.Muon_phi[t]), f(sub.Muon_mass[t]), f(sub.Muon_charge[t])
    mu_match = f(matched[idx][t])
    el_pt, el_eta, el_phi, el_m, el_q = f(sub.Electron_pt[e]), f(sub.Electron_eta[e]), f(sub.Electron_phi[e]), f(sub.Electron_mass[e]), f(sub.Electron_charge[e])
    px1, py1, pz1, e1 = objects.p4(mu_pt, mu_eta, mu_phi, mu_m)
    px2, py2, pz2, e2 = objects.p4(el_pt, el_eta, el_phi, el_m)
    mass = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
    base = (mu_pt > config.MU_PT_LEAD) & mu_match & (mass > config.MASS_LO) & (mass < config.MASS_HI)
    out = {}
    d = {"idx": idx, "pt1": mu_pt, "eta1": mu_eta, "pt_el": el_pt, "eta_el": el_eta, "mass": mass, "mu_eta": mu_eta, "mu_phi": mu_phi,
         "el_eta": el_eta, "el_phi": el_phi}
    if "Muon_genPartFlav" in ak.fields(ev):
        d["flav1"] = f(sub.Muon_genPartFlav[t]); d["flav_el"] = f(sub.Electron_genPartFlav[e])
    for name, sel in (("os", base & (mu_q * el_q < 0)), ("ss", base & (mu_q * el_q > 0))):
        out[name] = {k: v[sel] for k, v in d.items()}
    return out


def add(out, key, h):
    out[key] = out.get(key, 0) + h


def process(args):
    key, path = args
    s = samples.SAMPLES[key]; is_mc = s["is_mc"]
    out = {}
    pu = pileup.PileupWeights.from_json(OUT / "tnp" / "pileup_weights.json") if is_mc else None
    sf = tnp.ScaleFactors(OUT / "tnp" / "tnp_result.json") if is_mc else None
    norm = 1.0
    if is_mc:
        g = json.load(open(OUT / "gensums.json"))[key]
        norm = s["xsec_pb"] * config.LUMI_PB_NORMTAG / float(g["sumw"]) if s["xsec_pb"] else 1.0
    try:
        f = uproot.open(path)
        if "Events" not in f:
            return out
        for ev in f["Events"].iterate(filter_name=BR, step_size="150 MB", library="ak"):
            keep = regions.fired(ev) & regions.event_clean(ev)
            ev = ev[keep]
            if len(ev) == 0:
                continue
            n = len(ev)
            npv = np.asarray(ev.PV_npvsGood, dtype=float); met = np.asarray(ev.MET_pt, dtype=float)
            pmet = np.asarray(ev.PuppiMET_pt, dtype=float)
            masks = regions.muon_masks(ev); matched = regions.trigger_match(ev)
            jet = (ev.Jet_pt > 30) & (abs(ev.Jet_eta) < 2.4) & (ev.Jet_jetId >= 2)
            njet = np.asarray(ak.sum(jet, axis=1), dtype=float)
            # lepton-cleaned jets: dR > 0.4 to any tight/anti muon or selected electron
            lep_eta = ak.concatenate([ev.Muon_eta[masks["tight"] | masks["anti"]], ev.Electron_eta[regions.electron_mask(ev)]], axis=1)
            lep_phi = ak.concatenate([ev.Muon_phi[masks["tight"] | masks["anti"]], ev.Electron_phi[regions.electron_mask(ev)]], axis=1)
            j, l = ak.unzip(ak.cartesian([ak.zip({"eta": ev.Jet_eta, "phi": ev.Jet_phi}), ak.zip({"eta": lep_eta, "phi": lep_phi})], axis=1, nested=True))
            near = ak.any(objects.delta_r(j.eta, j.phi, l.eta, l.phi) < 0.4, axis=2)
            njet_clean = np.asarray(ak.sum(jet & ~near, axis=1), dtype=float)
            if is_mc:
                gen = np.asarray(ev.genWeight, dtype=float) * norm
                pref = np.asarray(ev.L1PreFiringWeight_Nom, dtype=float)
                ntrue = np.asarray(ev.Pileup_nTrueInt, dtype=float)
                w_nopu = gen * pref
                w_full = w_nopu * pu(ntrue)
                flav = np.asarray(ev.gen_lhe_flavour) if s.get("split_lhe") else None
            else:
                w_full = np.asarray(ev.skim_prescale, dtype=float); w_nopu = w_full; ntrue = None; flav = None
            def sample_of(idx):
                if flav is None:
                    return np.full(len(idx), s["fit_sample"] or key, dtype=object)
                return np.array([FLAVOUR_SAMPLE.get(int(x), "DYother") for x in flav[idx]], dtype=object)
            ones = np.ones(n, dtype=bool)
            # ---------------- SR: pileup / MET study
            sr = regions.dimuon_regions(ev, masks, matched, ones)["SR"]
            if sr is not None and len(sr["idx"]):
                idx = sr["idx"]
                jc = np.clip(njet_clean[idx], 0, 2)
                if is_mc:
                    pp = regions.is_prompt(sr["flav1"]) & regions.is_prompt(sr["flav2"])
                    w = np.where(pp, w_nopu[idx], 0.0)
                    smp = sample_of(idx)
                    for sname in np.unique(smp):
                        m = smp == sname
                        add(out, f"sr|{sname}|ntrue_npv", np.histogram2d(ntrue[idx][m], npv[idx][m], bins=[NTRUE_E, NPV_E], weights=w[m])[0])
                        add(out, f"sr|{sname}|ntrue_jet_met", np.histogramdd(np.stack([ntrue[idx][m], jc[m], met[idx][m]], 1), bins=[NTRUE_E, [0, 1, 2, 3], MET_E], weights=w[m])[0])
                        add(out, f"sr|{sname}|ntrue_jet_puppimet", np.histogramdd(np.stack([ntrue[idx][m], jc[m], pmet[idx][m]], 1), bins=[NTRUE_E, [0, 1, 2, 3], MET_E], weights=w[m])[0])
                        add(out, f"sr|{sname}|ntrue_npv_zpt", np.histogramdd(np.stack([ntrue[idx][m], np.clip(sr["zpt"][m], 0, 199.9)], 1), bins=[NTRUE_E, np.linspace(0, 200, 41)], weights=w[m])[0])
                else:
                    w = w_full[idx]
                    add(out, "sr|Data|npv", np.histogram(npv[idx], bins=NPV_E, weights=w)[0])
                    add(out, "sr|Data|npv_jet_met", np.histogramdd(np.stack([npv[idx], jc, met[idx]], 1), bins=[NPV_E, [0, 1, 2, 3], MET_E], weights=w)[0])
                    add(out, "sr|Data|npv_jet_puppimet", np.histogramdd(np.stack([npv[idx], jc, pmet[idx]], 1), bins=[NPV_E, [0, 1, 2, 3], MET_E], weights=w)[0])
                    add(out, "sr|Data|npv_zpt", np.histogram2d(npv[idx], np.clip(sr["zpt"], 0, 199.9), bins=[NPV_E, np.linspace(0, 200, 41)], weights=w)[0])
                    era = (np.asarray(ev.run)[idx] >= 281613).astype(int)
                    add(out, "sr|Data|era_npv", np.histogram2d(era, npv[idx], bins=[[-0.5, 0.5, 1.5], NPV_E], weights=w)[0])
            # ---------------- e-mu region OS/SS
            emu = emu_both(ev, masks, matched, ones)
            if emu is None:
                continue
            for charge in ("os", "ss"):
                d = emu[charge]; idx = d["idx"]
                if len(idx) == 0:
                    continue
                vals = {"pt_el": d["pt_el"], "met": met[idx], "njet": njet[idx], "njet_clean": njet_clean[idx], "npv": npv[idx],
                        "mass": d["mass"], "pt1": d["pt1"]}
                if is_mc:
                    w = w_full[idx] * sf.single_muon_weights(d["pt1"], d["eta1"])["nominal"]
                    pp = regions.is_prompt(d["flav1"]) & regions.is_prompt(d["flav_el"])
                    smp = sample_of(idx)
                    for sname in np.unique(smp):
                        m = smp == sname
                        for cat, sel in (("pp", pp & m), ("np", ~pp & m)):
                            for var, e in VARS.items():
                                add(out, f"emu|{sname}|{charge}|{cat}|{var}", np.histogram(vals[var][sel], bins=e, weights=w[sel])[0])
                                add(out, f"emu|{sname}|{charge}|{cat}|{var}|w2", np.histogram(vals[var][sel], bins=e, weights=w[sel] ** 2)[0])
                        fe = np.arange(0, 24, 1.0)
                        add(out, f"emu_flav|{sname}|{charge}", np.histogram2d(np.clip(d["flav1"][m], 0, 22), np.clip(d["flav_el"][m], 0, 22), bins=[fe, fe], weights=w[m])[0])
                        # non-prompt electron: pt_el x met for the shape
                        add(out, f"emu|{sname}|{charge}|np|ptel_met", np.histogram2d(d["pt_el"][m & ~pp], met[idx][m & ~pp], bins=[VARS["pt_el"], MET_E], weights=w[m & ~pp])[0])
                else:
                    w = w_full[idx]
                    for var, e in VARS.items():
                        add(out, f"emu|Data|{charge}|all|{var}", np.histogram(vals[var], bins=e, weights=w)[0])
    except Exception as exc:  # keep going, report
        out["error"] = f"{key} {path}: {exc!r}"
    return out


def merge(tot, part):
    for k, v in part.items():
        if k == "error":
            tot.setdefault("errors", []).append(v); continue
        tot[k] = tot.get(k, 0) + v


if __name__ == "__main__":
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    keys = sys.argv[2:] or list(samples.SAMPLES)
    tasks = []
    for key in keys:
        files = sorted(glob.glob(str(SKIM / key / "*.root")))
        tasks += [(key, f) for f in files]
    print(f"{len(tasks)} files, {workers} workers", flush=True)
    tot = {}
    with Pool(workers) as pool:
        for i, part in enumerate(pool.imap_unordered(process, tasks)):
            merge(tot, part)
            if i % 25 == 0:
                print(f"  {i}/{len(tasks)}", flush=True)
    here = Path(__file__).parent
    with open(here / "review_pass.pkl", "wb") as fh:
        pickle.dump(tot, fh)
    print("errors:", tot.get("errors", []))
    print("done")
