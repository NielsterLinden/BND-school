#!/usr/bin/env python
"""Freeze the Z -> tautau signal-region m_tt stack: data, the 11 simulated samples and the fake template in the 14 fit bins,
inclusive and per BDT category.

    source setup.sh && python presentation/data/extract_ztautau_sr_stack.py [--json PATH] [--check-only]

Reads (read-only, uproot): z-tautau/fit/fitinputs/ztautau.root (TH1D tautau_SR<k>__<sample>, the TRExFitter inputs) and
z-tautau/output/data/yields.json (the unclamped per-region yields of RESULTS.md). Writes presentation/data/ztautau_sr_stack.json.
Anchors: output/RESULTS.md prefit table, docs/02 cutflow (21 160), handoff.md.

The fit inputs have negative bins clamped to 0 (commit "Fix negative bins"), so their sums exceed the RESULTS.md yields by
the clamped content (WJets in SR2: -78.8 -> 0; DYlowmass in SR1: -4.6 -> 0; a few single bins elsewhere): the histograms
are frozen verbatim, the `totals` block carries the RESULTS.md numbers, and `totals_fitinputs` the histogram sums.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATASET_TAUTAU, VERSION_TAUTAU, ZTAUTAU, Checker, finalize, load_json, provenance,  # noqa: E402
                             standard_args)

ROOT_IN = ZTAUTAU / "fit" / "fitinputs" / "ztautau.root"
YIELDS = ZTAUTAU / "output" / "data" / "yields.json"
EDGES = [0, 40, 60, 70, 80, 90, 100, 110, 120, 130, 150, 175, 200, 250, 350]
REGIONS = ["tautau_SR0", "tautau_SR1", "tautau_SR2"]
SAMPLES = ["DYtautau", "DYtautau_nonfid", "DYee", "DYmumu", "DYlowmass", "WJets", "TTbar", "SingleTop", "WW", "WZ", "ZZ", "Fakes"]
GROUPS = {"rest": ["DYee", "DYmumu", "DYlowmass", "WJets", "TTbar", "SingleTop", "WW", "WZ", "ZZ"],
          "DYtautau_nonfid": ["DYtautau_nonfid"], "DYtautau": ["DYtautau"]}
STACK_ORDER = ["Fakes", "rest", "DYtautau_nonfid", "DYtautau"]


def build():
    import uproot
    y = load_json(YIELDS)
    hist, var = {}, {}
    with uproot.open(ROOT_IN) as f:
        for r in REGIONS:
            for s in SAMPLES + ["Data"]:
                h = f[f"{r}__{s}"]
                assert np.allclose(h.axis().edges(), EDGES), (r, s, h.axis().edges())
                hist[(r, s)] = h.values().astype(float)
                var[(r, s)] = h.variances().astype(float)
    data_reg = {r: [int(v) for v in hist[(r, "Data")]] for r in REGIONS}
    data_inc = np.sum([hist[(r, "Data")] for r in REGIONS], axis=0)
    assert np.all(data_inc == np.round(data_inc))
    mc = {}
    for s in SAMPLES:
        if s == "Fakes":
            continue
        mc[s] = {"inclusive": np.sum([hist[(r, s)] for r in REGIONS], axis=0).tolist(),
                 "per_region": {r: hist[(r, s)].tolist() for r in REGIONS},
                 "var_per_region": {r: var[(r, s)].tolist() for r in REGIONS},
                 "total_fitinputs": float(sum(hist[(r, s)].sum() for r in REGIONS)),
                 "total": y["yields"][s]["value"], "total_stat": y["yields"][s].get("stat"),
                 "per_region_total": {r: y["regions"][r][s]["value"] for r in REGIONS}}
    fakes = {"inclusive": np.sum([hist[(r, "Fakes")] for r in REGIONS], axis=0).tolist(),
             "per_region": {r: hist[(r, "Fakes")].tolist() for r in REGIONS},
             "var_per_region": {r: var[(r, "Fakes")].tolist() for r in REGIONS},
             "total": y["yields"]["Fakes"]["value"], "total_stat": y["yields"]["Fakes"]["stat"],
             "total_fitinputs": float(sum(hist[(r, "Fakes")].sum() for r in REGIONS)),
             "per_region_total": {r: y["regions"][r]["Fakes"]["value"] for r in REGIONS}}
    groups_counts = {}
    for g, members in GROUPS.items():
        per_region = {r: np.sum([hist[(r, s)] for s in members], axis=0) for r in REGIONS}
        inc = np.sum([per_region[r] for r in REGIONS], axis=0)
        neg = [[r, int(i)] for r in REGIONS for i in np.where(per_region[r] < 0)[0]]
        groups_counts[g] = {"inclusive": inc.tolist(), "per_region": {r: per_region[r].tolist() for r in REGIONS},
                            "inclusive_clamped": np.clip(inc, 0, None).tolist(), "negative_bins": neg,
                            "total": float(sum(y["yields"][s]["value"] for s in members)),
                            "per_region_total": {r: float(sum(y["regions"][r][s]["value"] for s in members)) for r in REGIONS}}
    mc_no_fakes = float(sum(y["yields"][s]["value"] for s in SAMPLES if s != "Fakes"))
    n_fakes = y["yields"]["Fakes"]["value"]
    pred = y["yields"]["Total"]["value"]
    pred_fi = float(sum(hist[(r, s)].sum() for r in REGIONS for s in SAMPLES))
    mc_fi = float(sum(hist[(r, s)].sum() for r in REGIONS for s in SAMPLES if s != "Fakes"))
    out = {
        "provenance": provenance("extract_ztautau_sr_stack.py", [ROOT_IN, YIELDS], dataset=DATASET_TAUTAU, version=VERSION_TAUTAU,
                                 histograms="TH1D <region>__<sample> of fit/fitinputs/ztautau.root (uproot), 14 variable-width m_tt bins, overflow in the last bin, "
                                            "negative bins clamped to 0 at fit-input writing; no rebinning",
                                 selection="OS, both tau_h DeepTau VSjet Tight, pT > 40 GeV, |eta| < 2.1, trigger-matched, extra-lepton veto; "
                                           "regions = BDT categories score < 0.55 / 0.55-0.90 / > 0.90",
                                 anchors=["z-tautau/output/RESULTS.md prefit table", "z-tautau/docs/02-selection.md cutflow", "z-tautau/handoff.md"],
                                 clamping_note="totals = output/data/yields.json (= RESULTS.md, unclamped, WJets SR2 = -78.8); "
                                               "totals_fitinputs = sums of the clamped histograms (what the stack draws)"),
        "edges": EDGES, "regions": REGIONS, "samples": SAMPLES, "groups": GROUPS, "stack_order": STACK_ORDER,
        "region_labels": {"tautau_SR0": "BDT < 0.55 (fake dominated)", "tautau_SR1": "0.55 < BDT < 0.90", "tautau_SR2": "BDT > 0.90 (signal dominated)"},
        "data": {"total": int(data_inc.sum()), "inclusive": [int(v) for v in data_inc], "per_region": data_reg,
                 "per_region_total": {r: int(sum(v)) for r, v in data_reg.items()}},
        "mc": mc,
        "fakes": fakes,
        "groups_counts": groups_counts,
        "totals": {"data": int(data_inc.sum()), "mc_no_fakes": mc_no_fakes, "fakes": n_fakes, "pred": pred,
                   "data_over_mc": float(data_inc.sum() / mc_no_fakes), "data_over_pred": float(data_inc.sum() / pred),
                   "fakes_fraction": float(n_fakes / data_inc.sum()),
                   "per_region": {r: {"data": int(sum(data_reg[r])), "pred": float(sum(y["regions"][r][s]["value"] for s in SAMPLES)),
                                      "fakes": y["regions"][r]["Fakes"]["value"]} for r in REGIONS},
                   "source": "output/data/yields.json = RESULTS.md prefit table (unclamped)"},
        "totals_fitinputs": {"mc_no_fakes": mc_fi, "fakes": fakes["total_fitinputs"], "pred": pred_fi,
                             "data_over_mc": float(data_inc.sum() / mc_fi), "data_over_pred": float(data_inc.sum() / pred_fi),
                             "per_region_pred": {r: float(sum(hist[(r, s)].sum() for s in SAMPLES)) for r in REGIONS},
                             "source": "sums of the clamped fit-input histograms"},
    }
    return out


def verify(d, ck: Checker):
    R = "z-tautau/output/RESULTS.md"
    t = d["totals"]
    ck.check("data total 21160", d["data"]["total"], 21160, 0, R + " / docs/02 cutflow")
    ck.check("data inclusive bins sum", sum(d["data"]["inclusive"]), 21160, 0, "consistency")
    ck.check("per-region data 15742 / 2704 / 2714", [d["data"]["per_region_total"][r] for r in d["regions"]], [15742, 2704, 2714], 0, R)
    ck.check("fakes total 13592 +- 1", t["fakes"], 13592, 1, R)
    ck.check("DYtautau total 3759 +- 1", d["mc"]["DYtautau"]["total"], 3759, 1, R)
    ck.check("DYtautau_nonfid total 2337 +- 1", d["mc"]["DYtautau_nonfid"]["total"], 2337, 1, R)
    ck.check("pred total 20718 +- 1", t["pred"], 20718, 1, R)
    ck.check("pred == mc_no_fakes + fakes", t["pred"], t["mc_no_fakes"] + t["fakes"], 1e-6, "consistency")
    ck.check("data / simulation (no fakes) 2.97", t["data_over_mc"], 2.97, 0.005, "briefs/ztautau.md section 3: 21160 / (20718 - 13592)")
    ck.check("data / pred 1.021", t["data_over_pred"], 1.021, 0.0005, R + " / brief")
    ck.check("fakes fraction 0.642", t["fakes_fraction"], 0.642, 0.0005, "docs/02 (64 %)")
    ck.check("per-region pred 15716 / 2556 / 2446", [t["per_region"][r]["pred"] for r in d["regions"]], [15716, 2556, 2446], 1, R)
    ck.check("per-region fakes 12351 / 1044 / 197", [t["per_region"][r]["fakes"] for r in d["regions"]], [12351, 1044, 197], 0.5, R)
    ck.check("per-region DYtautau 486 / 1131 / 2141", [d["mc"]["DYtautau"]["per_region_total"][r] for r in d["regions"]], [486, 1131, 2141], 0.5, R)
    for s in ("TTbar", "WJets", "DYee", "SingleTop", "WW", "WZ", "ZZ", "DYlowmass", "DYmumu"):
        want = {"TTbar": 410, "WJets": 358, "DYee": 110, "SingleTop": 47, "WW": 43, "WZ": 30, "ZZ": 20, "DYlowmass": 11, "DYmumu": 0}[s]
        ck.check(f"{s} total {want}", d["mc"][s]["total"], want, 0.5, R)
    # fit-input histograms vs the unclamped yields: clamping can only add, and only by the known negative content
    fi = d["totals_fitinputs"]
    ck.check_true("fit-input sums >= yields (clamping only adds)", all(d["mc"][s]["total_fitinputs"] >= d["mc"][s]["total"] - 1e-6 for s in d["mc"])
                  and d["fakes"]["total_fitinputs"] >= d["fakes"]["total"] - 1e-6, detail="")
    ck.check("fit-input pred sum 20805 (= 20718 + 86.7 clamped)", fi["pred"], 20805.07, 0.1, "fit/fitinputs/ztautau.root")
    ck.check("fit-input Fakes sum == yields.json Fakes", fi["fakes"], t["fakes"], 0.1, "no negative fake bins")
    ck.check("fit-input DYtautau sum 3760.3 (yields 3759.1: one clamped bin)", d["mc"]["DYtautau"]["total_fitinputs"], 3759.14, 1.5, "clamped bins < 1.5 events")
    ck.check("fit-input pred per region 15717.8 / 2562.1 / 2525.2", [fi["per_region_pred"][r] for r in d["regions"]], [15717.77, 2562.11, 2525.18], 0.05, "fit inputs")
    ck.check("fit-input pred vs RESULTS.md 15716 / 2556 / 2446", [fi["per_region_pred"][r] for r in d["regions"]], [15716, 2556, 2446], 1,
             R + " -- expected to differ by the clamped WJets (-78.8 in SR2) and DYlowmass (-4.6 in SR1)", soft=True)
    ck.check_true("no negative bins left in the fit-input groups", all(not g["negative_bins"] for g in d["groups_counts"].values()),
                  detail=str({k: g["negative_bins"] for k, g in d["groups_counts"].items()}))
    ck.check_true("groups cover all non-fake samples", sorted(sum(d["groups"].values(), [])) == sorted(s for s in d["samples"] if s != "Fakes"), "consistency")
    ck.check("stack sum == fit-input pred", np.sum([d["groups_counts"][g]["inclusive"] for g in ("rest", "DYtautau_nonfid", "DYtautau")], axis=0).sum() + sum(d["fakes"]["inclusive"]),
             fi["pred"], 1e-6, "consistency", rel=True)
    ck.check("edges", d["edges"], EDGES, 0, "z-tautau/ztautau/config.py FIT_BINS")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_sr_stack.json")
    args = ap.parse_args()
    finalize(args, build, verify, "ztautau_sr_stack")


if __name__ == "__main__":
    main()
