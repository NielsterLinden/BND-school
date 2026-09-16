#!/usr/bin/env python
"""Freeze the Z -> mumu fit result (mu_Z, sigma_fid, sigma(60-120), NPs, ranking, 12-bin SR plots).

    source setup.sh && python presentation/data/extract_zmumu_fit.py [--json PATH] [--check-only]

Reads (read-only): z-mumu/fit/results/zmumu_fit_result.json, z-mumu/output/v2/results_v2.json,
z-mumu/fit/results/zmumu/Plots/mumu_SR_{prefit,postfit}.yaml (TRExFitter plot dumps),
z-mumu/output/v2/histograms.pkl and fakes.json (for the 5 GeV rebin cross-checks).
Writes presentation/data/zmumu_fit.json. Anchors: z-mumu/handoff.md:11-13, output/v2/RESULTS_v2.md.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (FIT_RESULTS, LUMI_PB, V2, ZMUMU, Checker, add_zmumu_path, finalize, load_json,  # noqa: E402
                             provenance, rebin, standard_args)

add_zmumu_path()

FIT_JSON = FIT_RESULTS / "zmumu_fit_result.json"
RES_JSON = V2 / "results_v2.json"
PLOTS = FIT_RESULTS / "zmumu" / "Plots"
HIST_PKL = V2 / "histograms.pkl"
FAKES_JSON = V2 / "fakes.json"

# TRExFitter sample titles (z-mumu/scripts/v2_5_fit.py:163-173) -> fitting/CONVENTIONS.md names
TITLE_TO_SAMPLE = {"Z/#gamma* #rightarrow #mu#mu": "DYmumu", "Z/#gamma* #rightarrow #tau#tau": "DYtautau",
                   "Z/#gamma* #rightarrow ee": "DYee", "t#bar{t}": "TTbar", "tW": "SingleTop", "WW": "WW", "WZ": "WZ",
                   "ZZ": "ZZ", "W+jets (jet #rightarrow e)": "WJets", "Non-prompt (fake factor)": "Fakes"}
# NP -> systematic group (z-mumu/scripts/v2_5_fit.py:43-44, 178-200)
GROUP = {"Lumi": "Luminosity", "MuonReco": "Muon efficiency", "MuonID": "Muon efficiency", "MuonIso": "Muon efficiency",
         "MuonTrigger": "Muon efficiency", "Pileup": "Pileup", "L1Prefiring": "L1 prefiring", "MuonScale": "Muon momentum",
         "MuonRes": "Muon momentum", "PDF": "Signal modelling", "AlphaS": "Signal modelling", "QCDScale": "Signal modelling",
         "PS_ISR": "Signal modelling", "PS_FSR": "Signal modelling", "SigModel": "Signal modelling",
         "FakeStat_mumu": "Fakes", "FakeMethod_mumu": "Fakes", "ElectronEff_mumu": "Electron efficiency"}
NP_ORDER = ["Lumi", "MuonReco", "MuonID", "MuonIso", "MuonTrigger", "Pileup", "L1Prefiring", "MuonScale", "MuonRes",
            "PDF", "AlphaS", "QCDScale", "PS_ISR", "PS_FSR", "SigModel", "XS_TTbar", "XS_SingleTop", "XS_WW", "XS_WZ",
            "XS_ZZ", "XS_DYtautau", "FakeStat_mumu", "FakeMethod_mumu", "ElectronEff_mumu"]


def read_plot_yaml(path):
    import yaml
    with open(path) as fh:
        y = yaml.safe_load(fh)
    samples = {TITLE_TO_SAMPLE[s["Name"]]: [float(v) for v in s["Yield"]] for s in y["Samples"]}
    tot = y["Total"][0]
    return {"samples": samples, "total": [float(v) for v in tot["Yield"]],
            "unc_up": [float(v) for v in tot["UncertaintyUp"]], "unc_down": [float(v) for v in tot["UncertaintyDown"]],
            "data": [int(v) for v in y["Data"][0]["Yield"]], "edges": [float(v) for v in y["Figure"][0]["BinEdges"]]}


def build():
    fit = load_json(FIT_JSON)
    res = load_json(RES_JSON)
    pre, post = read_plot_yaml(PLOTS / "mumu_SR_prefit.yaml"), read_plot_yaml(PLOTS / "mumu_SR_postfit.yaml")
    assert pre["edges"] == post["edges"] and pre["data"] == post["data"]
    data = np.array(pre["data"], dtype=float)
    ratio_pre = (data / np.array(pre["total"])).tolist()
    ratio_post = (data / np.array(post["total"])).tolist()
    ratio_err = (np.sqrt(data) / np.array(pre["total"])).tolist()

    nps = []
    gammas = []
    for name, (pull, up, down) in fit["nps"].items():
        if name.startswith("gamma_"):
            gammas.append({"name": name, "bin": int(name.rsplit("_", 1)[1]), "value": pull, "err_up": up, "err_down": down})
        else:
            grp = GROUP.get(name, "Background normalisation" if name.startswith("XS_") else "other")
            nps.append({"name": name, "pull": pull, "constr_up": up, "constr_down": down, "constraint": 0.5 * (up + down), "group": grp})
    order = {n: i for i, n in enumerate(NP_ORDER)}
    nps.sort(key=lambda r: order.get(r["name"], 999))
    gammas.sort(key=lambda r: r["bin"])

    meta = fit["meta"]
    pred_60_120 = fit["sigma_fid_pred_pb"] / fit["A_60_120"]
    out = {
        "provenance": provenance("extract_zmumu_fit.py", [FIT_JSON, RES_JSON, PLOTS / "mumu_SR_prefit.yaml", PLOTS / "mumu_SR_postfit.yaml", HIST_PKL, FAKES_JSON],
                                 fit_tool="TRExFitter v1.8.0, profile-likelihood fit of m(mumu) in 12 x 5 GeV bins, MINOS on every parameter",
                                 anchors=["z-mumu/handoff.md:11-13", "z-mumu/output/v2/RESULTS_v2.md:7-37, 55-70", "fitting/CONVENTIONS.md section 6"],
                                 units="cross sections in pb; grouped impacts and uncertainties on mu_Z are relative fractions (not %)",
                                 title_to_sample=TITLE_TO_SAMPLE),
        "poi": {"name": fit["poi"], "value": fit["mu"], "err_up": fit["mu_err_up"], "err_down": fit["mu_err_down"], "err_sym": fit["mu_err"],
                "stat_only": fit["mu_stat_only_fit"], "stat_from_cov": fit["mu_stat_from_cov"], "syst_total_incl_lumi": fit["mu_syst"],
                "syst_excl_lumi": fit["mu_syst_excl_lumi"], "lumi_profiled": fit["mu_lumi_profiled"], "lumi_external": fit["mu_lumi"],
                "reference": "aMC@NLO DYJetsToLL_M-50, sigma_fid^pred below"},
        "sigma_fid": {"value": fit["sigma_fid_pb"], "stat": fit["sigma_fid_stat_pb"], "syst": fit["sigma_fid_syst_pb"], "lumi": fit["sigma_fid_lumi_pb"],
                      "total": fit["sigma_fid_tot_pb"], "pred": fit["sigma_fid_pred_pb"], "unit": "pb",
                      "definition": "pp -> Z/gamma* -> mu mu; dressed muons (dR < 0.1), pT > 26/20 GeV, |eta| < 2.4, 60 < m < 120 GeV"},
        "sigma_60_120": {"value": fit["sigma_60_120_pb"], "stat": fit["sigma_60_120_stat_pb"], "syst": fit["sigma_60_120_syst_pb"],
                         "lumi": fit["sigma_60_120_lumi_pb"], "acc": fit["sigma_60_120_acc_pb"], "total": fit["sigma_60_120_tot_pb"],
                         "A": fit["A_60_120"], "A_rel_unc": fit["acceptance"]["A_rel_unc"], "pred": pred_60_120, "unit": "pb",
                         "definition": "sigma(Z/gamma* -> mu mu, 60 < m_LHE < 120 GeV) = sigma_fid / A; pred = sigma_fid^pred / A (fitting/CONVENTIONS.md section 6)"},
        "sigma_m50": {"value": fit["sigma_m50_pb"], "total": fit["sigma_m50_tot_pb"], "A": fit["A_m50"], "pred": fit["meta"]["dy_xsec_pb"] / 3.0, "unit": "pb"},
        "C": meta["C_factor"],
        "counting": dict(meta["counting"]),
        "gof_p": fit["gof"]["gof_probability"],
        "n_fid_expected": meta["n_fid_expected"],
        "lumi_pb": meta["lumi_pb"],
        "grouped_impacts": dict(fit["grouped_impacts_mu"]),
        "nps": nps,
        "gammas": gammas,
        "ranking": fit["ranking"],
        "sr_12bin": {"edges": pre["edges"], "bin_width_gev": meta["sr_bin_width_gev"], "samples": list(pre["samples"]),
                     "prefit": pre["samples"], "postfit": post["samples"],
                     "total_prefit": {"yield": pre["total"], "unc_up": pre["unc_up"], "unc_down": pre["unc_down"]},
                     "total_postfit": {"yield": post["total"], "unc_up": post["unc_up"], "unc_down": post["unc_down"]},
                     "data": pre["data"], "ratio_prefit": ratio_pre, "ratio_postfit": ratio_post, "ratio_err": ratio_err,
                     "note": "prefit = mu_Z = 1 templates (TRExFitter Plots/mumu_SR_prefit.yaml); postfit = best-fit NPs and mu_Z"},
        "stability": [{k: v for k, v in s.items()} for s in res["stability"]],
        "np_count_note": ("23 nuisance parameters are constrained in the signal-region fit (+ 12 gammas); fit/zmumu.config declares 25 "
                          "Systematic blocks, of which XS_WJets and ElectronEff_mumu act on the validation region mumu_CRemu only; "
                          "handoff.md:18 quotes 24."),
    }
    return out


def verify(d, ck: Checker):
    fit = load_json(FIT_JSON)
    res = load_json(RES_JSON)
    H = "z-mumu/handoff.md:11-13"
    R = "z-mumu/output/v2/RESULTS_v2.md"
    p, sf, s6 = d["poi"], d["sigma_fid"], d["sigma_60_120"]
    ck.check("mu_Z", p["value"], 0.988134, 5e-7, R + ":9 / fit result")
    ck.check("mu_Z err_up", p["err_up"], 0.0158981, 5e-8, R + ":9")
    ck.check("mu_Z err_down", p["err_down"], 0.0155715, 5e-8, R + ":9")
    ck.check("mu_Z rounded 0.988 +- 0.016", [p["value"], p["err_sym"]], [0.988, 0.016], 5e-4, H)
    ck.check("mu_Z stat-only 0.031%", 100 * p["stat_only"], 0.031, 5e-4, R + ":37")
    ck.check("sigma_fid 790.1 pb", sf["value"], 790.1, 0.05, H)
    ck.check("sigma_fid stat 0.2", sf["stat"], 0.2, 0.05, H)
    ck.check("sigma_fid syst 8.5", sf["syst"], 8.5, 0.05, H)
    ck.check("sigma_fid lumi 9.6", sf["lumi"], 9.6, 0.05, H)
    ck.check("sigma_fid exact 790.078 +- 0.248 +- 8.473 +- 9.595", [sf["value"], sf["stat"], sf["syst"], sf["lumi"]],
             [790.078, 0.248, 8.473, 9.595], 5e-4, "fit/results/zmumu_fit_result.json")
    ck.check("sigma_fid^pred 799.566 pb", sf["pred"], 799.566, 5e-4, "z-mumu/handoff.md:77")
    ck.check("sigma(60-120) 1931 +- 33 pb", [s6["value"], s6["total"]], [1931.0, 33.0], 0.5, H)
    ck.check("sigma(60-120) exact 1930.746 +- 32.921", [s6["value"], s6["total"]], [1930.746, 32.921], 5e-4, "fit result")
    ck.check("A_60_120 0.409209", s6["A"], 0.409209, 5e-7, "z-mumu/handoff.md:78")
    ck.check("sigma^pred(60-120) 1953.9 pb", s6["pred"], 1953.9, 0.05, "fitting/CONVENTIONS.md section 6 / handoff.md:79")
    ck.check("C factor 0.7914", d["C"], 0.7914, 5e-5, R + ":17")
    ck.check("C factor exact 0.791415", d["C"], 0.791415, 5e-7, "fit result meta")
    ck.check("counting cross-check 794.7 pb", d["counting"]["sigma_fid_pb"], 794.7, 0.05, R + ":19")
    ck.check("counting n_obs", d["counting"]["n_obs"], 10378567, 0, "z-mumu/handoff.md:35")
    ck.check("GoF p 0.792273", d["gof_p"], 0.792273, 5e-7, R + ":9")
    ck.check("lumi_pb", d["lumi_pb"], LUMI_PB, 1e-6, "fitting/CONVENTIONS.md")
    gi = d["grouped_impacts"]
    for name, want in (("FullSyst", 1.573), ("Luminosity", 1.163), ("Muon efficiency", 0.869), ("L1 prefiring", 0.506),
                       ("Signal modelling", 0.491), ("Gammas", 0.382), ("Muon momentum", 0.367), ("Pileup", 0.134),
                       ("Background normalisation", 0.081), ("Fakes", 0.044), ("Electron efficiency", 0.0)):
        ck.check(f"grouped impact {name} {want}%", 100 * gi[name], want, 5e-4, R + ":26-37")
    nps = {r["name"]: r for r in d["nps"]}
    for name, pull, constr in (("SigModel", 0.66, 0.14), ("MuonScale", -0.48, 0.14), ("MuonRes", 0.40, 0.12), ("PDF", 0.56, 0.93),
                               ("QCDScale", 0.76, 0.79), ("Lumi", 0.08, 1.00), ("Pileup", 0.06, 0.82), ("L1Prefiring", -0.06, 1.00)):
        ck.check(f"pull {name} {pull:+.2f} ({constr:.2f})", [nps[name]["pull"], nps[name]["constraint"]], [pull, constr], 5e-3, R + ":55-70")
    ck.check_true("every NP has a group", all(r["group"] != "other" for r in d["nps"]), detail=str([r["name"] for r in d["nps"] if r["group"] == "other"]))
    # the config declares 25 Systematic blocks; XS_WJets and ElectronEff_mumu act on the VALIDATION region mumu_CRemu only,
    # so the SR fit constrains 23 (handoff.md:18 quotes "24": a documentation count, reported, not adjusted here)
    ck.check_true("23 constrained NPs in the fit result + 12 gammas", len(d["nps"]) == 23 and len(d["gammas"]) == 12,
                  "fit/results/zmumu_fit_result.json nps", f"{len(d['nps'])} NPs, {len(d['gammas'])} gammas")
    n_cfg = sum(1 for line in open(ZMUMU / "fit" / "zmumu.config") if line.startswith("Systematic:"))
    ck.check("25 Systematic blocks in fit/zmumu.config (23 in the SR + XS_WJets, ElectronEff_mumu on the validation region)", n_cfg, 25, 0, "z-mumu/fit/zmumu.config")
    # 12-bin plots vs the fit tables, the report ratio and the 1 GeV histograms
    sr = d["sr_12bin"]
    tab = {TITLE_TO_SAMPLE.get(s.get("Sample"), s.get("Sample")): s["Yield"] for s in fit["table_prefit"][0]["Samples"] if "Sample" in s}
    for smp in sr["samples"]:
        ck.check(f"YAML prefit sum == table_prefit {smp}", sum(sr["prefit"][smp]), tab[smp], 1e-9, "results_v2.json /fit/table_prefit", rel=True)
    ck.check("YAML prefit total == table_prefit Total", sum(sr["total_prefit"]["yield"]), tab["Total"], 1e-9, "results_v2.json", rel=True)
    ck.check("data 12 bins sum", sum(sr["data"]), 10378567, 0, "z-mumu/handoff.md:35")
    ck.check("ratio_prefit == lineshape/data_over_pred_prefit_5gev", sr["ratio_prefit"], res["lineshape"]["data_over_pred_prefit_5gev"], 1e-9,
             "results_v2.json /lineshape", rel=True)
    ck.check("ratio_prefit range 0.9644 .. 1.0133 (plan table quotes 0.965 .. 1.013)", [min(sr["ratio_prefit"]), max(sr["ratio_prefit"])], [0.9644, 1.0133], 5e-5,
             "results_v2.json /lineshape/data_over_pred_prefit_5gev")
    with open(HIST_PKL, "rb") as fh:
        hall = pickle.load(fh)
    for smp in ("DYmumu", "DYtautau", "TTbar", "SingleTop", "WW", "WZ", "ZZ", "Data"):
        k = f"{smp}|SR|mass_fit|nominal"
        ck.check(f"YAML prefit {smp} == 5 GeV rebin of histograms.pkl", sr["prefit"][smp] if smp != "Data" else sr["data"], rebin(hall[k], 5), 1e-9,
                 "output/v2/histograms.pkl", rel=True)
    f = np.array(load_json(FAKES_JSON)["templates"]["nominal"])
    ck.check("YAML prefit Fakes == 5 GeV rebin of fakes.json template", sr["prefit"]["Fakes"], np.clip(rebin(f, 10), 0, None), 1e-9, "output/v2/fakes.json", rel=True)
    ck.check("postfit DYmumu total / prefit == mu_Z (within NP shifts)", sum(sr["postfit"]["DYmumu"]) / sum(sr["prefit"]["DYmumu"]), d["poi"]["value"], 0.01, "consistency")
    ck.check_true("ratio_postfit within 0.3% of 1", max(abs(r - 1) for r in sr["ratio_postfit"]) < 0.003, detail=f"max |r-1| = {max(abs(r - 1) for r in sr['ratio_postfit']):.4f}")
    ck.check("stability nominal mu", d["stability"][0]["mu"], 0.988134, 5e-7, "fit/results/stability.json")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "zmumu_fit.json")
    args = ap.parse_args()
    finalize(args, build, verify, "zmumu_fit")


if __name__ == "__main__":
    main()
