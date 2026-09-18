#!/usr/bin/env python
"""Freeze the real-data pictures of the Z -> ee selection clips (section 3): one selected event, the electron p_T
spectrum, the m_ee spectra for the five electron-ID levels, a cutflow and two electron candidates for the ID clip.

    source fitting/setup.sh && python presentation/data/extract_zee_selection.py [--json PATH] [--check-only] [--files N]

Reads (read-only) the z-ee team's NanoAOD cache of the SingleElectron primary dataset
(/data/atlas/users/lvdurenw/nano_cache/Run2016{G,H}, records 30529 / 30562, the files z-ee.ipynb processed), the first
N files of each period (default 5 + 5), certified lumi sections only (datasets/GRL/GRL.txt, runs 278820-284044).
The selection is z-ee's (z-ee/handoff.md "Selection and method", z-ee.ipynb): HLT_Ele27_WPTight_Gsf; electrons with
p_T > 20 GeV, |eta| < 2.5, cutBased >= 3 (Medium); exactly two, opposite charge; 60 < m_ee < 120 GeV.
The spectra are a sample of the data (shapes and fractions); the totals on screen are the full-dataset numbers of the
fit (zee_fit.json: 6,320,097 in 60-120 GeV) and of the notebook (6,560,718 before the mass window).
Writes presentation/data/zee_selection.json.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import REPO, Checker, finalize, load_json, provenance, standard_args  # noqa: E402

CACHE = Path("/data/atlas/users/lvdurenw/nano_cache")
GRL = REPO / "datasets" / "GRL" / "GRL.txt"
TRIGGER = "HLT_Ele27_WPTight_Gsf"
BRANCHES = ["run", "luminosityBlock", "event", TRIGGER, "Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass",
            "Electron_charge", "Electron_cutBased", "Electron_pfRelIso03_all", "Electron_sieie", "Electron_hoe",
            "Electron_deltaEtaSC"]
WP_NAMES = ["none", "veto", "loose", "medium", "tight"]      # cutBased >= 0 .. 4
M_EDGES = np.arange(60.0, 120.01, 2.0)
PT_EDGES = np.arange(0.0, 100.01, 5.0)
N_DEFAULT = 5


def lumi_mask(runs, lumis, grl):
    mask = np.zeros(len(runs), dtype=bool)
    for r in np.unique(runs):
        rng = grl.get(int(r))
        if not rng:
            continue
        sel = runs == r
        ls = lumis[sel]
        ok = np.zeros(len(ls), dtype=bool)
        for lo, hi in rng:
            ok |= (ls >= lo) & (ls <= hi)
        mask[sel] = ok
    return mask


def build(n_files=N_DEFAULT):
    import awkward as ak
    import uproot
    import vector
    vector.register_awkward()
    grl = {int(k): v for k, v in json.load(open(GRL)).items() if 278820 <= int(k) <= 284044}
    files = []
    for period in ("Run2016G", "Run2016H"):
        files += sorted(glob.glob(str(CACHE / period / "*.root")))[:n_files]

    cut = {"events_certified": 0, "trigger": 0, "two_medium": 0, "opposite_charge": 0, "window": 0}
    os_h = {w: np.zeros(len(M_EDGES) - 1) for w in WP_NAMES}
    ss_h = {w: np.zeros(len(M_EDGES) - 1) for w in WP_NAMES}
    pt_in, pt_out = np.zeros(len(PT_EDGES) - 1), np.zeros(len(PT_EDGES) - 1)
    lead_in, sub_in = np.zeros(len(PT_EDGES) - 1), np.zeros(len(PT_EDGES) - 1)
    event, fake = None, None

    for fpath in files:
        a = uproot.open(fpath)["Events"].arrays(BRANCHES)
        a = a[lumi_mask(ak.to_numpy(a.run), ak.to_numpy(a.luminosityBlock), grl)]
        cut["events_certified"] += len(a)
        a = a[a[TRIGGER]]
        cut["trigger"] += len(a)
        ele = ak.zip({"pt": a.Electron_pt, "eta": a.Electron_eta, "phi": a.Electron_phi, "mass": a.Electron_mass},
                     with_name="Momentum4D")

        def pairs(wp, ptcut):
            sel = (a.Electron_cutBased >= wp) & (a.Electron_pt > ptcut) & (abs(a.Electron_eta) < 2.5)
            two = ak.sum(sel, axis=1) == 2
            e = ele[sel][two]
            q = a.Electron_charge[sel][two]
            os_ = ak.to_numpy(q[:, 0] * q[:, 1] < 0)
            m = ak.to_numpy((e[:, 0] + e[:, 1]).mass)
            return two, sel, e, os_, m

        for i, w in enumerate(WP_NAMES):
            _, _, _, os_, m = pairs(i, 20.0)
            os_h[w] += np.histogram(m[os_], M_EDGES)[0]
            ss_h[w] += np.histogram(m[~os_], M_EDGES)[0]

        two, sel, e, os_, m = pairs(3, 20.0)
        cut["two_medium"] += int(ak.sum(two))
        cut["opposite_charge"] += int(os_.sum())
        win = os_ & (m > 60) & (m < 120)
        cut["window"] += int(win.sum())

        # p_T: Medium, |eta| < 2.5, exactly two, opposite charge, no p_T cut (NanoAOD stores p_T > 5 GeV)
        two0, sel0, e0, os0, m0 = pairs(3, 0.0)
        pts = ak.to_numpy(e0.pt)[os0]
        w0 = (m0[os0] > 60) & (m0[os0] < 120)
        pt_in += np.histogram(pts[w0].ravel(), PT_EDGES)[0]
        pt_out += np.histogram(pts[~w0].ravel(), PT_EDGES)[0]
        lead_in += np.histogram(pts[w0].max(axis=1), PT_EDGES)[0]
        sub_in += np.histogram(pts[w0].min(axis=1), PT_EDGES)[0]

        if event is None:
            idx = np.where(two)[0]
            pt2, eta2, phi2 = ak.to_numpy(e.pt), ak.to_numpy(e.eta), ak.to_numpy(e.phi)
            dphi = np.abs((phi2[:, 0] - phi2[:, 1] + np.pi) % (2 * np.pi) - np.pi)
            good = (win & (np.abs(m - 91.19) < 1.5) & np.all(np.abs(eta2) < 1.2, axis=1) & np.all((pt2 > 38) & (pt2 < 50), axis=1)
                    & (dphi > 2.4) & (dphi < 2.9) & (np.abs(eta2[:, 0] - eta2[:, 1]) > 0.6))
            if good.any():
                k = int(np.where(good)[0][0])
                j = int(idx[k])
                s = ak.to_numpy(sel[j])
                recs = []
                for n in np.where(s)[0]:
                    recs.append({"pt": float(a.Electron_pt[j][n]), "eta": float(a.Electron_eta[j][n]),
                                 "phi": float(a.Electron_phi[j][n]), "charge": int(a.Electron_charge[j][n]),
                                 "cutBased": int(a.Electron_cutBased[j][n]), "relIso": float(a.Electron_pfRelIso03_all[j][n]),
                                 "sieie": float(a.Electron_sieie[j][n]), "hoe": float(a.Electron_hoe[j][n])})
                recs.sort(key=lambda r: -r["pt"])
                event = {"run": int(a.run[j]), "lumi": int(a.luminosityBlock[j]), "event": int(a.event[j]),
                         "m_ee": float(m[k]), "electrons": recs, "file": Path(fpath).name}
        if fake is None:
            # a real electron candidate that fails every ID level, looks like a jet: wide shower, hadronic energy, not isolated
            cb, ptv, etav = a.Electron_cutBased, a.Electron_pt, a.Electron_eta
            c = (cb == 0) & (ptv > 25) & (abs(etav) < 1.2) & (a.Electron_pfRelIso03_all > 0.5) & (a.Electron_hoe > 0.15) \
                & (a.Electron_sieie > 0.012)
            ev_i = np.where(ak.to_numpy(ak.any(c, axis=1)))[0]
            if len(ev_i):
                j = int(ev_i[0])
                n = int(np.where(ak.to_numpy(c[j]))[0][0])
                fake = {"run": int(a.run[j]), "event": int(a.event[j]), "pt": float(ptv[j][n]), "eta": float(etav[j][n]),
                        "phi": float(a.Electron_phi[j][n]), "charge": int(a.Electron_charge[j][n]), "cutBased": 0,
                        "relIso": float(a.Electron_pfRelIso03_all[j][n]), "sieie": float(a.Electron_sieie[j][n]),
                        "hoe": float(a.Electron_hoe[j][n])}
        print(f"  {Path(fpath).name}: cumulative window {cut['window']}", flush=True)

    peak = {w: float(os_h[w].sum()) for w in WP_NAMES}
    ss = {w: float(ss_h[w].sum()) for w in WP_NAMES}
    wp = {w: {"os": os_h[w].tolist(), "ss": ss_h[w].tolist(), "n_os": peak[w], "n_ss": ss[w],
              "eff_rel_none": peak[w] / peak["none"], "ss_rel_none": ss[w] / ss["none"], "ss_over_os": ss[w] / peak[w]}
          for w in WP_NAMES}
    fit = load_json(Path(__file__).resolve().parent / "zee_fit.json")
    return {
        "provenance": provenance("extract_zee_selection.py", [CACHE / "Run2016G", CACHE / "Run2016H", GRL],
                                 dataset="CMS 2016 Open Data, SingleElectron Run2016G+H, NanoAODv9 (records 30529, 30562)",
                                 files=[Path(f).name for f in files], n_files=len(files),
                                 selection="z-ee/handoff.md: HLT_Ele27_WPTight_Gsf; p_T > 20 GeV, |eta| < 2.5, cutBased >= 3; "
                                           "exactly two, opposite charge; 60 < m_ee < 120 GeV",
                                 note="a sample of the data for shapes and fractions; the full-dataset totals are "
                                      "totals.window (the fit's data) and totals.selected (z-ee.ipynb cell 13)"),
        "trigger": TRIGGER,
        "cuts": {"pt_min": 20.0, "abs_eta_max": 2.5, "id": "cutBased >= 3 (Medium)", "window": [60.0, 120.0],
                 "trigger_pt": 27.0, "tracker_abs_eta": 2.5, "ecal_barrel_abs_eta": 1.479, "ecal_endcap_abs_eta": 3.0},
        "totals": {"window": int(sum(fit["data"])), "selected": 3108273 + 3452445,
                   "selected_note": "Run2016G 3,108,273 + Run2016H 3,452,445 (z-ee.ipynb cell 13 output, z-ee/handoff.md)"},
        "cutflow_sample": cut,
        "event": event,
        "fake_candidate": fake,
        "m_edges": M_EDGES.tolist(),
        "wp_names": WP_NAMES,
        "wp": wp,
        "pt_edges": PT_EDGES.tolist(),
        "pt": {"in_window": pt_in.tolist(), "outside": pt_out.tolist(), "lead_in_window": lead_in.tolist(),
               "sub_in_window": sub_in.tolist(),
               "note": "electrons of opposite-charge Medium pairs (|eta| < 2.5, exactly two, no p_T cut) after the trigger; "
                       "in_window = pairs with 60 < m_ee < 120 GeV"},
    }


def verify(d, ck: Checker):
    ck.check_true("an event was chosen", d["event"] is not None)
    ev = d["event"]
    ck.check_true("event passes the selection", all(e["cutBased"] >= 3 and e["pt"] > 20 and abs(e["eta"]) < 2.5 for e in ev["electrons"])
                  and len(ev["electrons"]) == 2 and ev["electrons"][0]["charge"] * ev["electrons"][1]["charge"] < 0
                  and 60 < ev["m_ee"] < 120, source="z-ee/handoff.md")
    ck.check_true("a failing candidate was chosen", d["fake_candidate"] is not None)
    ck.check("totals: 6,320,097 in the window (fit data)", d["totals"]["window"], 6320097, source="zee_fit.json")
    ck.check("totals: 6,560,718 selected (notebook)", d["totals"]["selected"], 6560718, source="z-ee/handoff.md")
    w = d["wp"]
    ck.check_true("tighter ID keeps fewer pairs", all(w[a]["n_os"] > w[b]["n_os"] for a, b in zip(d["wp_names"], d["wp_names"][1:])))
    ck.check_true("tighter ID removes same-sign pairs faster than the peak",
                  w["tight"]["ss_rel_none"] < w["tight"]["eff_rel_none"] and w["medium"]["ss_over_os"] < w["none"]["ss_over_os"])
    c = d["cutflow_sample"]
    ck.check("sample: window / opposite-charge pairs ~ full data (6,320,097 / 6,560,718 = 0.963)",
             c["window"] / c["opposite_charge"], 0.9633, tol=0.01, source="z-ee.ipynb, zee_fit.json")
    ck.check("medium OS window == wp['medium'] OS", w["medium"]["n_os"], c["window"], tol=0.5)
    pt_in = np.array(d["pt"]["in_window"])
    ck.check_true("Z electrons peak at 40-45 GeV (Jacobian peak near m_Z/2)", int(np.argmax(pt_in)) == 8)
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__), "zee_selection.json")
    ap.add_argument("--files", type=int, default=N_DEFAULT, help="files per period (default 5)")
    args = ap.parse_args()
    finalize(args, lambda: build(args.files), verify, "zee_selection")


if __name__ == "__main__":
    main()
