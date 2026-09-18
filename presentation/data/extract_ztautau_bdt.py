#!/usr/bin/env python
"""Freeze the Z -> tautau k-fold BDT: the signal-region score distribution (data, fakes, simulation), the same-sign closure in
the score, the three categories with their prefit yields, the held-out AUC and the feature importance.

    source fitting/setup.sh && python presentation/data/extract_ztautau_bdt.py [--json PATH] [--check-only]

Reads (read-only): z-tautau/output/data/bdt.json (scripts/step3b_bdt.py: sr_score, ss_closure_score, training), z-tautau/fit/bdt_info.json,
z-tautau/output/data/yields.json (the per-category yields of RESULTS.md, with the per-category C_OS/SS).
Writes presentation/data/ztautau_bdt.json. Anchors: output/RESULTS.md (categories, AUC), docs/09-bdt.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATASET_TAUTAU, VERSION_TAUTAU, ZTAUTAU, Checker, finalize, load_json, provenance,  # noqa: E402
                             standard_args)

BDT_JSON = ZTAUTAU / "output" / "data" / "bdt.json"
BDT_INFO = ZTAUTAU / "fit" / "bdt_info.json"
YIELDS = ZTAUTAU / "output" / "data" / "yields.json"
REGIONS = ["tautau_SR0", "tautau_SR1", "tautau_SR2"]
REST = ["DYee", "DYmumu", "DYlowmass", "WJets", "TTbar", "SingleTop", "WW", "WZ", "ZZ"]
FEATURE_LABELS = {"t1_pt": "pT(tau1)", "t2_pt": "pT(tau2)", "pt_ratio": "pT(tau2)/pT(tau1)", "t1_abseta": "|eta(tau1)|", "t2_abseta": "|eta(tau2)|",
                  "dr_tt": "dR(tau1, tau2)", "dphi_tt": "dphi(tau1, tau2)", "met": "pT_miss", "met_sig": "MET significance", "pt_vis": "pT(tautau)",
                  "dphi_met_tt": "dphi(pT_miss, tautau)", "pt_tt": "pT(tautau + pT_miss)", "njets": "N_jets", "jet1_pt": "pT(jet1)", "t1_dm": "DM(tau1)", "t2_dm": "DM(tau2)"}


def build():
    b = load_json(BDT_JSON)
    info = load_json(BDT_INFO)
    y = load_json(YIELDS)
    tr = b["training"]
    sr = b["sr_score"]
    mc = {s: [float(v) for v in sr["mc"][s]] for s in sr["mc"]}
    groups = {"DYtautau": mc["DYtautau"], "DYtautau_nonfid": mc["DYtautau_nonfid"], "rest": np.sum([mc[s] for s in REST], axis=0).tolist()}
    fakes = [float(v) for v in sr["fakes"]]
    cats = {}
    for k, r in enumerate(REGIONS):
        yr = y["regions"][r]
        rest = float(sum(yr[s]["value"] for s in REST))
        cats[r] = {"DYtautau": yr["DYtautau"]["value"], "fakes": yr["Fakes"]["value"], "DYtautau_nonfid": yr["DYtautau_nonfid"]["value"], "rest": rest,
                   "data": int(yr["Data"]["value"]), "pred": yr["DYtautau"]["value"] + yr["Fakes"]["value"] + yr["DYtautau_nonfid"]["value"] + rest,
                   "score_lo": b["category_edges"][k], "score_hi": b["category_edges"][k + 1],
                   "fakes_inclusive_C": b["category_yields_prefit"][str(k)]["fakes"]}
    ss = b["ss_closure_score"]
    ss_ratio = (np.asarray(ss["obs"]) - np.asarray(ss["obs_mc"])) / np.asarray(ss["pred"])
    out = {
        "provenance": provenance("extract_ztautau_bdt.py", [BDT_JSON, BDT_INFO, YIELDS], dataset=DATASET_TAUTAU, version=VERSION_TAUTAU,
                                 method=("XGBoost, k = 5 folds by event number, signal = fiducial Z -> tautau simulation in the SR, background = AR data x FF; "
                                         "every event scored by the model that never saw it; 16 mass-agnostic inputs (config.BDT_FEATURES)"),
                                 anchors=["z-tautau/output/RESULTS.md (categories, AUC)", "z-tautau/docs/09-bdt.md"],
                                 categories_note=("`categories` = output/data/yields.json = RESULTS.md (fake yields with C_OS/SS per category); "
                                                  "`fakes_inclusive_C` = bdt.json category_yields_prefit (step 3b, inclusive C: 12408 / 975 / 186)")),
        "sr_score": {"edges": sr["edges"], "data": [int(v) for v in sr["data"]], "fakes": fakes, "mc": mc, "groups": groups,
                     "pred": (np.asarray(fakes) + np.sum([mc[s] for s in mc], axis=0)).tolist(),
                     "note": "prefit signal-region score, 20 bins of 0.05; fakes = AR data x FF x C (inclusive C at step 3b), MC genuine-tau1 subtracted"},
        "ss_closure_score": {**ss, "ratio_mcsub": ss_ratio.tolist(),
                             "note": "same-sign closure in the score: SS_T observed (obs, minus obs_mc) vs FF x SS_L predicted (pred, MC subtracted)"},
        "category_edges": b["category_edges"],
        "categories": cats,
        "auc_test_mean": float(np.mean(tr["auc_test"])),
        "auc_test": tr["auc_test"], "auc_train": tr["auc_train"],
        "importance_sorted": [[k, v] for k, v in sorted(tr["importance"].items(), key=lambda kv: -kv[1])],
        "features": tr["features"], "feature_labels": FEATURE_LABELS,
        "params": tr["params"], "k": tr["k"], "n_sig": tr["n_sig"], "n_bkg": tr["n_bkg"],
        "stat_only_sensitivity": b["stat_only_sensitivity"],
        "s_over_b_sr2": cats["tautau_SR2"]["DYtautau"] / (cats["tautau_SR2"]["pred"] - cats["tautau_SR2"]["DYtautau"]),
    }
    assert info["training"]["auc_test"] == tr["auc_test"], "fit/bdt_info.json and output/data/bdt.json disagree"
    return out


def verify(d, ck: Checker):
    R = "z-tautau/output/RESULTS.md"
    c = d["categories"]
    ck.check("category data 15742 / 2704 / 2714", [c[r]["data"] for r in REGIONS], [15742, 2704, 2714], 0, R)
    ck.check("category DYtautau 486 / 1131 / 2141 +- 1", [c[r]["DYtautau"] for r in REGIONS], [486, 1131, 2141], 1, R)
    ck.check("category fakes 12351 / 1044 / 197 +- 1", [c[r]["fakes"] for r in REGIONS], [12351, 1044, 197], 1, R)
    ck.check("category DYtautau_nonfid 1972 / 240 / 124", [c[r]["DYtautau_nonfid"] for r in REGIONS], [1972, 240, 124], 1, R)
    ck.check("category pred 15716 / 2556 / 2446", [c[r]["pred"] for r in REGIONS], [15716, 2556, 2446], 1, R)
    ck.check("AUC mean 0.966 +- 0.001", d["auc_test_mean"], 0.966, 0.001, R)
    ck.check("AUC folds 0.967, 0.966, 0.968, 0.964, 0.965", d["auc_test"], [0.967, 0.966, 0.968, 0.964, 0.965], 5e-4, R)
    ck.check("category edges 0 / 0.55 / 0.9 / 1", d["category_edges"], [0, 0.55, 0.9, 1], 0, "config.BDT_CATEGORY_EDGES")
    ck.check("score data sum == 21160", sum(d["sr_score"]["data"]), 21160, 0, R)
    ck.check("score data per category (from the 20 bins) == 15742 / 2704 / 2714", [sum(d["sr_score"]["data"][:11]), sum(d["sr_score"]["data"][11:18]), sum(d["sr_score"]["data"][18:])],
             [15742, 2704, 2714], 0, "bins 0-10 (< 0.55), 11-17 (0.55-0.90), 18-19 (>= 0.90)")
    ck.check("score DYtautau sum == 3759", sum(d["sr_score"]["groups"]["DYtautau"]), 3759, 1, R)
    ck.check("score fakes sum == 13569 (inclusive C at step 3b)", sum(d["sr_score"]["fakes"]), 13592, 30, R + " (13592 with the per-category C)", soft=True)
    ck.check("score fakes per category == bdt.json fakes_inclusive_C", [sum(d["sr_score"]["fakes"][:11]), sum(d["sr_score"]["fakes"][11:18]), sum(d["sr_score"]["fakes"][18:])],
             [c[r]["fakes_inclusive_C"] for r in REGIONS], 1e-6, "consistency", rel=True)
    ck.check_true("importance: dR(tau1,tau2) first", d["importance_sorted"][0][0] == "dr_tt", "docs/09", d["importance_sorted"][0][0])
    ck.check("importance sums to 1", sum(v for _, v in d["importance_sorted"]), 1.0, 1e-6, "xgboost")
    ck.check_true("16 features", len(d["features"]) == 16, "config.BDT_FEATURES")
    ck.check("S/B in SR2 ~ 5.3 (handoff) / 10.9 (docs/09)?", d["s_over_b_sr2"], 5.3, 0.3, "handoff.md quotes 5.3 (with WJets = -79 the ratio is 7.0)", soft=True)
    ck.check("stat-only sensitivity 2.01 % -> 1.93 %", [100 * d["stat_only_sensitivity"]["inclusive"], 100 * d["stat_only_sensitivity"]["categories"]], [2.01, 1.93], 0.005, R)
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_bdt.json")
    args = ap.parse_args()
    finalize(args, build, verify, "ztautau_bdt")


if __name__ == "__main__":
    main()
