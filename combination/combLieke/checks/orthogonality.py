#!/usr/bin/env python
"""Event-level orthogonality of the three fitted selections, measured on 2016 G+H data.

    source ../../../fitting/setup.sh
    python checks/orthogonality.py [--workers 8] [--max-files N]

An event can only be in two channels if it passes both offline selections, so each pair is
tested on the data sample of one channel, with the other channel's offline selection applied on
top (its trigger is dropped, which can only make the count larger: every number is an upper bound).

* mumu & ee      SingleMuon v2 skim ($BND_SKIM_DIR): the exact mumu SR of z-mumu/zmumu/regions.py,
                 then the z-ee selection (exactly two electrons pT > 20, |eta| < 2.5, cutBased >= 3,
                 opposite sign, 60 < m_ee < 120).
* tautau & mumu  Tau skim of z-tautau (all events before the tau pair selection): the z-tautau extra
* tautau & ee    lepton veto (z-tautau/ztautau/objects.extra_lepton_veto), then >= 2 muons passing the
                 mumu signal-muon definition (no charge or mass requirement), or the z-ee selection.
                 The tautau fitted events are a subset of the veto-passing events (the veto is applied
                 when the ntuples are built, scripts/step2_ntuples.py). Events that survive are looked
                 up by (run, lumi, event) in the z-tautau data ntuples and counted if they are in the
                 tautau signal region (z-tautau/ztautau/analysis.regions, "SR": OS, both legs Tight).

Writes checks/orthogonality.json.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import awkward as ak
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SKIM_MUMU = Path(os.environ.get("BND_SKIM_DIR", "/data/atlas/users/nterlind/BND-school-cache/skims_v2"))
SKIM_TAUTAU = Path(os.environ.get("BND_TAUTAU_CACHE", "/data/atlas/users/sjankovy/BND-school-cache/ztautau")) / "skims_v1"
ERAS_MUMU = ["data_2016G", "data_2016H"]
ERAS_TAUTAU = ["data_2016G", "data_2016H"]

#: z-ee/z-ee.ipynb, process_sample
EE_PT, EE_ETA, EE_ID = 20.0, 2.5, 3


def _p4_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2):
    px = pt1 * np.cos(phi1) + pt2 * np.cos(phi2)
    py = pt1 * np.sin(phi1) + pt2 * np.sin(phi2)
    pz = pt1 * np.sinh(eta1) + pt2 * np.sinh(eta2)
    e = np.sqrt((pt1 * np.cosh(eta1)) ** 2 + m1 ** 2) + np.sqrt((pt2 * np.cosh(eta2)) ** 2 + m2 ** 2)
    return np.sqrt(np.maximum(e ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))


def ee_good(ev):
    return (ev.Electron_pt > EE_PT) & (abs(ev.Electron_eta) < EE_ETA) & (ev.Electron_cutBased >= EE_ID)


def ee_selected(ev):
    """z-ee offline selection without HLT_Ele27_WPTight_Gsf: exactly two good electrons, OS, 60-120."""
    good = ee_good(ev)
    two = ak.to_numpy(ak.sum(good, axis=1) == 2)
    out = np.zeros(len(ev), dtype=bool)
    if not two.any():
        return out
    sub, g = ev[two], good[two]
    pt, eta, phi = sub.Electron_pt[g], sub.Electron_eta[g], sub.Electron_phi[g]
    mass, q = sub.Electron_mass[g], sub.Electron_charge[g]
    col = lambda a, i: ak.to_numpy(a[:, i]).astype(np.float64)
    m = _p4_mass(col(pt, 0), col(eta, 0), col(phi, 0), col(mass, 0), col(pt, 1), col(eta, 1), col(phi, 1), col(mass, 1))
    ok = (col(q, 0) * col(q, 1) < 0) & (m > 60.0) & (m < 120.0)
    out[np.nonzero(two)[0]] = ok
    return out


def mumu_file(path):
    import uproot
    sys.path.insert(0, str(REPO / "z-mumu"))
    from zmumu import regions
    res = dict(files=1, events=0, mumu_sr=0, mumu_sr_and_ee_selected=0, mumu_sr_and_2ee_good=0, overlap_ids=[])
    with uproot.open(path) as f:
        if "Events" not in f:
            return res
        branches = ["run", "luminosityBlock", "event", "PV_npvsGood", "HLT_IsoMu24", "HLT_IsoTkMu24", "Muon_*",
                    "TrigObj_*", "Flag_*", "FsrPhoton_*", "Electron_*"]
        for ev in f["Events"].iterate(filter_name=branches, step_size="200 MB", library="ak"):
            res["events"] += len(ev)
            masks = regions.muon_masks(ev)
            sr = regions.dimuon_regions(ev, masks, regions.trigger_match(ev), regions.event_clean(ev))["SR"]
            if sr is None or len(sr["idx"]) == 0:
                continue
            idx = sr["idx"]
            res["mumu_sr"] += len(idx)
            sub = ev[idx]
            res["mumu_sr_and_2ee_good"] += int(np.sum(ak.to_numpy(ak.sum(ee_good(sub), axis=1) >= 2)))
            both = ee_selected(sub)
            res["mumu_sr_and_ee_selected"] += int(both.sum())
            for r, l, e in zip(ak.to_numpy(sub.run[both]), ak.to_numpy(sub.luminosityBlock[both]), ak.to_numpy(sub.event[both])):
                res["overlap_ids"].append([int(r), int(l), int(e)])
    return res


def tautau_file(path):
    import uproot
    sys.path.insert(0, str(REPO / "z-tautau"))
    from ztautau import objects
    sys.path.insert(0, str(REPO / "z-mumu"))
    from zmumu import config as mcfg
    res = dict(files=1, events=0, veto_pass=0, veto_pass_2mumu_tight=0, veto_pass_2ee_good=0, veto_pass_ee_selected=0,
               ids_mumu=[], ids_ee=[])
    with uproot.open(path) as f:
        if "Events" not in f:
            return res
        branches = ["run", "luminosityBlock", "event", "Muon_*", "Electron_*"]
        for ev in f["Events"].iterate(filter_name=branches, step_size="200 MB", library="ak"):
            res["events"] += len(ev)
            veto = objects.extra_lepton_veto(ev)
            res["veto_pass"] += int(veto.sum())
            sub = ev[veto]
            # the mumu signal-muon definition (z-mumu/zmumu/regions.muon_masks, "tight")
            tight = ((sub.Muon_pt > 20.0) & (abs(sub.Muon_eta) < 2.4) & (abs(sub.Muon_dxy) < mcfg.MU_DXY_MAX)
                     & (abs(sub.Muon_dz) < mcfg.MU_DZ_MAX) & sub.Muon_tightId & (sub.Muon_pfRelIso04_all < 0.15))
            two_mu = ak.to_numpy(ak.sum(tight, axis=1) >= 2)
            two_e = ak.to_numpy(ak.sum(ee_good(sub), axis=1) >= 2)
            res["veto_pass_2mumu_tight"] += int(two_mu.sum())
            res["veto_pass_2ee_good"] += int(two_e.sum())
            ee_sel = ee_selected(sub)
            res["veto_pass_ee_selected"] += int(ee_sel.sum())
            for key, m in (("ids_mumu", two_mu), ("ids_ee", ee_sel)):
                for r, l, e in zip(ak.to_numpy(sub.run[m]), ak.to_numpy(sub.luminosityBlock[m]), ak.to_numpy(sub.event[m])):
                    res[key].append([int(r), int(l), int(e)])
    return res


def tautau_signal_region_ids():
    """(run, lumi, event) of every data event in the tautau signal region, from the z-tautau ntuples."""
    import uproot
    sys.path.insert(0, str(REPO / "z-tautau"))
    from ztautau import analysis
    ids, n = set(), 0
    for era in ERAS_TAUTAU:
        path = SKIM_TAUTAU.parent / "ntuples_v1" / f"{era}.root"
        d = uproot.open(path)["ntuple"].arrays(["run", "lumi", "event", "t1_pt", "t2_pt", "t1_vsjet", "t2_vsjet", "os"],
                                               library="np")
        sr = analysis.regions(d)["SR"]
        n += int(sr.sum())
        ids.update(zip(d["run"][sr].tolist(), d["lumi"][sr].tolist(), d["event"][sr].tolist()))
    return ids, n


def run(tasks, fn, workers):
    total = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fn, str(t)) for t in tasks]
        for n, fut in enumerate(as_completed(futures), 1):
            part = fut.result()
            for k, v in part.items():
                total[k] = total.get(k, [] if isinstance(v, list) else 0) + v
            if n % 20 == 0 or n == len(futures):
                print(f"  {fn.__name__}: {n}/{len(futures)} files", flush=True)
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=None)
    args = ap.parse_args()

    mumu_files = [p for era in ERAS_MUMU for p in sorted((SKIM_MUMU / era).glob("*.root"))][: args.max_files]
    tau_files = [p for era in ERAS_TAUTAU for p in sorted((SKIM_TAUTAU / era).glob("*.root"))][: args.max_files]
    print(f"mumu skim: {len(mumu_files)} files; tautau skim: {len(tau_files)} files", flush=True)
    tau = run(tau_files, tautau_file, args.workers)
    mumu = run(mumu_files, mumu_file, args.workers)
    sr_ids, n_tau_sr = tautau_signal_region_ids()
    in_sr = lambda rows: [r for r in rows if tuple(r) in sr_ids]

    n_mumu_fit = 10378567           # mumu_SR__Data integral, z-mumu/fit/fitinputs/zmumu.root
    n_ee_fit = 6320097              # Data integral, z-ee/Zee_fit.tar.gz
    out = {
        "date": "2026-09-16",
        "complete": args.max_files is None,
        "method": __doc__.split("Writes")[0].strip(),
        "mumu_and_ee": {
            "sample": "SingleMuon v2 skim, Run2016G+H, golden JSON",
            "files": mumu["files"], "events_read": mumu["events"],
            "mumu_sr_events": mumu["mumu_sr"], "mumu_sr_events_in_fit": n_mumu_fit,
            "mumu_sr_with_two_or_more_ee_electrons": mumu["mumu_sr_and_2ee_good"],
            "overlap_upper_bound": mumu["mumu_sr_and_ee_selected"],
            "fraction_of_mumu": mumu["mumu_sr_and_ee_selected"] / n_mumu_fit,
            "fraction_of_ee": mumu["mumu_sr_and_ee_selected"] / n_ee_fit,
            "overlap_event_ids": mumu["overlap_ids"],
        },
        "tautau_and_mumu_ee": {
            "sample": "z-tautau Tau skim (skims_v1), Run2016G+H",
            "files": tau["files"], "events_read": tau["events"], "events_passing_tautau_lepton_veto": tau["veto_pass"],
            "veto_pass_with_two_or_more_mumu_signal_muons": tau["veto_pass_2mumu_tight"],
            "veto_pass_with_two_or_more_ee_signal_electrons": tau["veto_pass_2ee_good"],
            "veto_pass_and_ee_selected": tau["veto_pass_ee_selected"],
            "tautau_signal_region_events": n_tau_sr,
            "tautau_sr_and_two_mumu_signal_muons": len(in_sr(tau["ids_mumu"])),
            "tautau_sr_and_ee_selected": len(in_sr(tau["ids_ee"])),
            "fraction_of_tautau_sr": len(in_sr(tau["ids_ee"])) / max(n_tau_sr, 1),
            "fraction_of_ee": len(in_sr(tau["ids_ee"])) / 6320097,
            "overlap_event_ids_ee": in_sr(tau["ids_ee"]),
            "overlap_event_ids_mumu": in_sr(tau["ids_mumu"]),
        },
    }
    if out["complete"] and mumu["mumu_sr"] != n_mumu_fit:
        out["mumu_and_ee"]["warning"] = f"SR count {mumu['mumu_sr']} differs from the fit input {n_mumu_fit}"
    (HERE / "orthogonality.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if not kk.endswith("ids")} if isinstance(v, dict) else v
                      for k, v in out.items() if k != "method"}, indent=1))


if __name__ == "__main__":
    main()
