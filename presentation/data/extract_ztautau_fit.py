#!/usr/bin/env python
"""Freeze the Z -> tautau fit result (mu_Z, sigma_fid, sigma(60-120), grouped impacts, NPs, ranking, the 3 x 14-bin
category plots pre/post-fit).

    source fitting/setup.sh && python presentation/data/extract_ztautau_fit.py [--json PATH] [--check-only]

Reads (read-only): z-tautau/fit/results/ztautau_fit_result.json, z-tautau/output/results.json (fit.mcsub),
z-tautau/fit/results/ztautau/Plots/tautau_SR{0,1,2}_{prefit,postfit}.yaml (TRExFitter plot dumps, git-ignored, reached
through the worktree link), z-tautau/fit/ztautau.config (NP -> group). Writes presentation/data/ztautau_fit.json.
Anchors: z-tautau/handoff.md "Result (v3, nominal)" + the impact table, output/RESULTS.md.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATASET_TAUTAU, LUMI_PB, VERSION_TAUTAU, ZTAUTAU, Checker, finalize, load_json,  # noqa: E402
                             provenance, standard_args)

FIT_JSON = ZTAUTAU / "fit" / "results" / "ztautau_fit_result.json"
RES_JSON = ZTAUTAU / "output" / "results.json"
PLOTS = ZTAUTAU / "fit" / "results" / "ztautau" / "Plots"
CONFIG = ZTAUTAU / "fit" / "ztautau.config"
REGIONS = ["tautau_SR0", "tautau_SR1", "tautau_SR2"]
EDGES = [0, 40, 60, 70, 80, 90, 100, 110, 120, 130, 150, 175, 200, 250, 350]
SAMPLES = ["DYtautau", "DYtautau_nonfid", "DYee", "DYmumu", "DYlowmass", "WJets", "TTbar", "SingleTop", "WW", "WZ", "ZZ", "Fakes"]
# TRExFitter sample titles (z-tautau/scripts/step5_fit.py) -> docs/CONVENTIONS.md names
TITLE_TO_SAMPLE = {"Z#rightarrow#tau#tau (fiducial)": "DYtautau", "Z/#gamma*#rightarrow#tau#tau (non-fid.)": "DYtautau_nonfid",
                   "Z#rightarrowee": "DYee", "Z#rightarrow#mu#mu": "DYmumu", "Z/#gamma*#rightarrowll (m<50)": "DYlowmass",
                   "W+jets": "WJets", "t#bar{t}": "TTbar", "single t": "SingleTop", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ",
                   "jet#rightarrow#tau_{h} (FF)": "Fakes"}
NPS_SHOWN = ["TauID_DM0", "TauID_DM1", "TauID_DM10", "TauID_DM11", "FakeOSSS_tautau_c0", "FakeClosure_tautau_c2_hi", "TauTrigger_DM1", "Lumi"]


def read_plot_yaml(path):
    import yaml
    with open(path) as fh:
        y = yaml.safe_load(fh)
    samples = {TITLE_TO_SAMPLE[s["Name"]]: [float(v) for v in s["Yield"]] for s in y["Samples"]}
    tot = y["Total"][0]
    return {"samples": samples, "total": [float(v) for v in tot["Yield"]],
            "unc_up": [float(v) for v in tot["UncertaintyUp"]], "unc_down": [float(v) for v in tot["UncertaintyDown"]],
            "data": [int(v) for v in y["Data"][0]["Yield"]], "edges": [float(v) for v in y["Figure"][0]["BinEdges"]]}


def np_groups():
    """NP name -> group as the fit spells it (Category, or SubCategory for the 'Tau' category), from fit/ztautau.config."""
    groups, cur, cat, sub = {}, None, None, None
    for line in open(CONFIG):
        m = re.match(r'\s*Systematic:\s*"([^"]+)"', line)
        if m:
            if cur:
                groups[cur] = sub or cat
            cur, cat, sub = m.group(1), None, None
        m = re.match(r'\s*Category:\s*"([^"]+)"', line)
        if m:
            cat = m.group(1)
        m = re.match(r'\s*SubCategory:\s*"([^"]+)"', line)
        if m:
            sub = m.group(1)
    if cur:
        groups[cur] = sub or cat
    return groups


def build():
    fit = load_json(FIT_JSON)
    res = load_json(RES_JSON)
    mc = res["fit"]["mcsub"]
    groups = np_groups()
    plots = {"prefit": {}, "postfit": {}}
    for r in REGIONS:
        for tag in ("prefit", "postfit"):
            p = read_plot_yaml(PLOTS / f"{r}_{tag}.yaml")
            assert p["edges"] == EDGES, (r, tag, p["edges"])
            plots[tag][r] = {"samples": {s: p["samples"][s] for s in SAMPLES}, "total": p["total"], "unc_up": p["unc_up"],
                             "unc_down": p["unc_down"], "data": p["data"]}
    nps = []
    for name, (pull, up, down) in fit["nps"].items():
        if name.startswith("gamma_"):
            continue
        nps.append({"name": name, "pull": pull, "constraint": 0.5 * (up + down), "constr_up": up, "constr_down": down,
                    "group": groups.get(name, "Background normalisation" if name.startswith("XS_") else "other")})
    order = list(groups)
    nps.sort(key=lambda r: order.index(r["name"]) if r["name"] in order else 999)
    gi = dict(mc["grouped_impact"])
    impact_order = sorted((g for g in gi if g != "FullSyst"), key=lambda g: -gi[g])
    sf, s6 = mc["sigma_fid_pb"], mc["sigma_60_120_pb"]
    out = {
        "provenance": provenance("extract_ztautau_fit.py", [FIT_JSON, RES_JSON, CONFIG] + [PLOTS / f"{r}_{t}.yaml" for r in REGIONS for t in ("prefit", "postfit")],
                                 dataset=DATASET_TAUTAU, version=VERSION_TAUTAU,
                                 fit_tool="TRExFitter v1.8.0, profile-likelihood fit of m_tt (MET-likelihood di-tau mass) in 3 BDT categories x 14 bins, "
                                          "tautau_SR0 fitted above 110 GeV only (DropBins 1-7), 39 NPs + 42 gammas",
                                 anchors=["z-tautau/handoff.md (Result v3 nominal; impact table)", "z-tautau/output/RESULTS.md",
                                          "z-tautau/docs/08-fit-and-results.md"],
                                 units="cross sections in pb; grouped impacts and uncertainties on mu_Z are relative fractions (not %)",
                                 title_to_sample=TITLE_TO_SAMPLE),
        "poi": {"name": fit["poi"], "value": fit["poi_value"], "err_up": fit["poi_err_up"], "err_down": fit["poi_err_down"],
                "err_sym": fit["poi_err"], "stat": mc["mu_stat"], "syst": mc["mu_syst"],
                "stat_only_fit": dict(fit["poi_stat_only"]),
                "reference": "aMC@NLO DYJetsToLL_M-50 fiducial template, sigma(Z/gamma* -> ll, m > 50) = 6077.22/3 pb per flavour"},
        "asimov": dict(fit["poi_asimov"]),
        "gof_p": fit["gof"]["gof_probability"],
        "grouped_impact": gi,
        "impact_order": impact_order,
        "stat_impact": mc["mu_stat"],
        "sigma_fid": {"value": sf["value"], "stat": sf["stat"], "syst": sf["syst"], "err_up": sf["err_up"], "err_down": sf["err_down"],
                      "pred": sf["prediction"], "unit": "pb",
                      "definition": "pp -> Z/gamma* -> tau tau, both tau hadronic, visible tau pT > 40 GeV, |eta| < 2.1, 60 < m_LHE < 120 GeV"},
        "sigma_60_120": {"value": s6["value"], "stat": s6["stat"], "syst": s6["syst"], "err_up": s6["err_up"], "err_down": s6["err_down"],
                         "acc": s6["acceptance"], "pred": s6["prediction"], "A": res["prediction"]["A"],
                         "A_rel_unc": fit["sigma_60_120_pb"]["acceptance_rel_unc"], "unit": "pb",
                         "definition": "sigma(Z/gamma* -> tau tau, 60 < m_LHE < 120 GeV) = mu_Z x sigma^pred(60-120); err_up/down = syst (+) stat (MINOS); "
                                       "acc = the acceptance uncertainty (3.7 %), quoted separately, not printed"},
        "prediction": dict(res["prediction"]),
        "lumi_pb": fit["lumi_pb"],
        "nps": nps,
        "nps_shown": NPS_SHOWN,
        "ranking": mc["ranking"][:12],
        "ranking_note": "output/results.json fit.mcsub.ranking, top 12 by post-fit impact; impact_up/down are the shifts of mu_Z",
        "edges": EDGES,
        "regions": plots,
        "region_labels": {"tautau_SR0": "BDT < 0.55 (fake dominated)", "tautau_SR1": "0.55 < BDT < 0.90", "tautau_SR2": "BDT > 0.90 (signal dominated)"},
        "dropped_bins": {"tautau_SR0": {"m_tt_max": 110.0, "bins": [0, 1, 2, 3, 4, 5, 6]}},
        "plots_note": "prefit = mu_Z = 1 templates (Plots/*_prefit.yaml, negative bins already clamped to 0 in the fit inputs); "
                      "postfit = best-fit NPs and mu_Z; total uncertainty band from the same yaml",
    }
    return out


def verify(d, ck: Checker):
    H = "z-tautau/handoff.md (Result v3)"
    R = "z-tautau/output/RESULTS.md"
    p, sf, s6 = d["poi"], d["sigma_fid"], d["sigma_60_120"]
    ck.check("mu_Z 1.07074", p["value"], 1.07074, 5e-6, R)
    ck.check("mu_Z err_up 0.114103", p["err_up"], 0.114103, 5e-7, R)
    ck.check("mu_Z err_down 0.0999551", p["err_down"], 0.0999551, 5e-8, R)
    ck.check("mu_Z rounded 1.071 +0.114 -0.100", [p["value"], p["err_up"], p["err_down"]], [1.071, 0.114, 0.100], 5e-4, H)
    ck.check("mu_Z stat 0.0210107", p["stat"], 0.0210107, 5e-8, R)
    ck.check("mu_Z syst 0.10495", p["syst"], 0.10495, 5e-6, R + " (stat 0.021, syst 0.105)")
    ck.check("Asimov +0.106 -0.093", [d["asimov"]["err_up"], d["asimov"]["err_down"]], [0.106, 0.093], 5e-4, H)
    ck.check("GoF p 0.25", d["gof_p"], 0.25, 5e-3, R + " (handoff table quotes 0.22 from an earlier pass)")
    ck.check("sigma_fid 4.8189", sf["value"], 4.8189, 5e-5, "results.json fit.mcsub")
    ck.check("sigma_fid 4.82 +- 0.09 +- 0.47", [sf["value"], sf["stat"], sf["syst"]], [4.82, 0.09, 0.47], 5e-3, H)
    ck.check("sigma_fid stat 0.0946 syst 0.4723", [sf["stat"], sf["syst"]], [0.0946, 0.4723], 5e-5, "results.json")
    ck.check("sigma_fid pred 4.5005", sf["pred"], 4.5005, 5e-5, R)
    ck.check("sigma(60-120) 2082.46", s6["value"], 2082.46, 5e-3, "results.json")
    ck.check("sigma(60-120) 2082 +222 -194", [s6["value"], s6["err_up"], s6["err_down"]], [2082, 222, 194], 0.5, H)
    ck.check("sigma(60-120) stat 40.86 syst 204.11 acc 76.30", [s6["stat"], s6["syst"], s6["acc"]], [40.86, 204.11, 76.30], 5e-3, "results.json")
    ck.check("sigma(60-120) exact err_up 221.92 err_down 194.40", [s6["err_up"], s6["err_down"]], [221.92, 194.40], 5e-3, "results.json")
    ck.check("prediction 1944.88 pb", s6["pred"], 1944.88, 5e-3, R + " (handoff: 1944.9)")
    ck.check("prediction 1944.9 pb rounded", s6["pred"], 1944.9, 0.05, H)
    ck.check("A 0.002314", s6["A"], 0.002314, 5e-7, H)
    ck.check("lumi_pb", d["lumi_pb"], LUMI_PB, 1e-6, "docs/CONVENTIONS.md")
    gi = d["grouped_impact"]
    for name, want in (("FullSyst", 0.1050), ("Tau ID", 0.0844), ("Fakes", 0.0650), ("Gammas", 0.0389), ("Tau trigger", 0.0302),
                       ("Background normalisation", 0.0257), ("Tau energy scale", 0.0147), ("Signal modelling", 0.0125), ("MET", 0.0081),
                       ("Luminosity", 0.0070), ("Pileup", 0.0055), ("L1 prefiring", 0.0018)):
        ck.check(f"grouped impact {name} {want}", gi[name], want, 5e-5 if name not in ("Tau ID", "Fakes") else 1e-3, R)
    ck.check("grouped impact Tau ID 8.4 % / Fakes 6.5 %", [100 * gi["Tau ID"], 100 * gi["Fakes"]], [8.4, 6.5], 0.1, H)
    ck.check("stat impact 0.021", d["stat_impact"], 0.021, 5e-4, H + " (data statistics 2.1 %)")
    ck.check_true("impact_order starts Tau ID, Fakes, Gammas, Tau trigger", d["impact_order"][:4] == ["Tau ID", "Fakes", "Gammas", "Tau trigger"],
                  detail=str(d["impact_order"][:4]))
    nps = {r["name"]: r for r in d["nps"]}
    # handoff.md says "39 NPs + 42 gammas"; the fit result lists 40 non-gamma names (every FakeClosure/FakeOSSS/Tau*/XS_*
    # parameter of docs/07) and 35 gammas: a documentation count, so soft
    ck.check("non-gamma NPs (handoff says 39; the workspace has 40)", len(d["nps"]), 39, tol=1,
             source="z-tautau/handoff.md (39 NPs + 42 gammas), docs/07 table", soft=True)
    ck.check_true("every nps_shown is a fitted NP", all(n in nps for n in d["nps_shown"]), detail=str([n for n in d["nps_shown"] if n not in nps]))
    ck.check_true("every NP has a group", all(r["group"] != "other" for r in d["nps"]), detail=str([r["name"] for r in d["nps"] if r["group"] == "other"]))
    for name, pull, constr in (("FakeOSSS_tautau_c0", 0.316, 0.681), ("TauID_DM0", 0.137, 0.843), ("TauID_DM10", -0.312, 0.904),
                               ("FakeClosure_tautau_c2_hi", 1.544, 0.828), ("Lumi", 0.008, 0.991)):
        ck.check(f"pull {name} {pull:+.3f} ({constr:.3f})", [nps[name]["pull"], nps[name]["constraint"]], [pull, constr], 5e-4, "fit/results/ztautau_fit_result.json")
    ck.check("FakeOSSS_tautau_c0 pulled +0.63 sigma?", nps["FakeOSSS_tautau_c0"]["pull"], 0.63, 0.05, "docs/08 quotes +0.63 (earlier pass); stored 0.316", soft=True)
    ck.check_true("ranking top entry FakeOSSS_tautau_c0, then TauID_DM0, TauID_DM1", [r["name"] for r in d["ranking"][:3]] == ["FakeOSSS_tautau_c0", "TauID_DM0", "TauID_DM1"],
                  detail=str([r["name"] for r in d["ranking"][:3]]))
    # region plots vs RESULTS.md prefit table (the yaml / fit inputs have negative bins clamped to 0: WJets in SR2 is -78.8 in the
    # unclamped yields.json / RESULTS.md, 0 in the fit inputs; DYlowmass in SR1 -4.6 -> 0; a few single bins elsewhere)
    pre, post = d["regions"]["prefit"], d["regions"]["postfit"]
    ck.check("data per region 15742 / 2704 / 2714", [sum(pre[r]["data"]) for r in REGIONS], [15742, 2704, 2714], 0, R)
    ck.check("data total 21160", sum(sum(pre[r]["data"]) for r in REGIONS), 21160, 0, R)
    ck.check("prefit DYtautau per region 486 / 1131 / 2141", [sum(pre[r]["samples"]["DYtautau"]) for r in REGIONS], [486, 1131, 2141], 1.6, R + " (SR0 487.5 in the yaml: clamped bins)")
    ck.check("prefit Fakes per region 12351 / 1044 / 197", [sum(pre[r]["samples"]["Fakes"]) for r in REGIONS], [12351, 1044, 197], 0.5, R)
    ck.check("prefit DYtautau_nonfid per region 1972 / 240 / 124", [sum(pre[r]["samples"]["DYtautau_nonfid"]) for r in REGIONS], [1972, 240, 124], 1.0, R)
    ck.check("prefit total per region 15716 / 2556 / 2446 (RESULTS.md, unclamped)", [sum(pre[r]["total"]) for r in REGIONS], [15716, 2556, 2446], 1.0,
             R + " -- the yaml totals are 15718 / 2562 / 2525: negative bins clamped in the fit inputs (WJets SR2 -78.8 -> 0)", soft=True)
    ck.check("prefit total per region == yaml sample sum", [sum(pre[r]["total"]) for r in REGIONS], [sum(sum(v) for v in pre[r]["samples"].values()) for r in REGIONS], 1e-6, "consistency", rel=True)
    ck.check("prefit total per region 15718 / 2562 / 2525 (clamped fit inputs)", [sum(pre[r]["total"]) for r in REGIONS], [15717.77, 2562.11, 2525.18], 0.05, "fit/fitinputs/ztautau.root sums")
    ck.check("postfit DYtautau / prefit DYtautau ~ mu_Z", sum(sum(post[r]["samples"]["DYtautau"]) for r in REGIONS) / sum(sum(pre[r]["samples"]["DYtautau"]) for r in REGIONS),
             d["poi"]["value"], 0.02, "consistency (NP shifts)")
    ck.check_true("postfit total within 1.5 % of data per region", all(abs(sum(post[r]["total"]) / sum(post[r]["data"]) - 1) < 0.015 for r in REGIONS),
                  detail=str([round(sum(post[r]["total"]) / sum(post[r]["data"]), 4) for r in REGIONS]))
    ck.check("edges", d["edges"], EDGES, 0, "z-tautau/ztautau/config.py FIT_BINS")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_fit.json")
    args = ap.parse_args()
    finalize(args, build, verify, "ztautau_fit")


if __name__ == "__main__":
    main()
