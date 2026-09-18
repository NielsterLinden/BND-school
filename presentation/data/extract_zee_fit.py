#!/usr/bin/env python
"""Freeze the Z -> ee template fit (z-ee delivery of 16 Sep 14:44, unchanged by the freeze) for section 3.

    source fitting/setup.sh && python presentation/data/extract_zee_fit.py [--json PATH] [--check-only]

Reads (read-only) z-ee/Zee_fit.tar.gz, unpacked into presentation/work/extract/zee_fit/:
  Plots/ee_SR_prefit.yaml, Plots/ee_SR_postfit.yaml    TRExFitter per-bin yields (what its plots show)
  Histograms/Zee_fit_histos.root                       nominal + the +-1 sigma templates the fit used
  Fits/Zee_fit.txt                                     post-fit nuisance parameters, gammas, mu_Z
  Fits/Zee_fit_errDecomp_mu_Z.txt, Fits/Zee_fit_group_errDecomp_mu_Z.txt
  reference_cross_section.txt                          sigma(Z/gamma* -> ee, 60 < m_LHE < 120) of the aMC@NLO sample
and, for the result frame, combination/combLieke/config/references.json (published CMS / ATLAS values) and
presentation/data/ztautau_reference.json (the prediction band the tau tau chapter already shows, anchor rule 06 B5).
Writes presentation/data/zee_fit.json. Anchors: docs/FREEZE.md "Z -> ee" (mu = 0.942 +- 0.015, 1840.8 +- 29.9 pb),
z-ee/fit.config (samples, systematics).
"""

from __future__ import annotations

import argparse
import sys
import tarfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATA_DIR, LUMI_PB, REPO, WORK, Checker, finalize, load_json, provenance,  # noqa: E402
                             standard_args)

TARBALL = REPO / "z-ee" / "Zee_fit.tar.gz"
UNPACK = WORK / "zee_fit"
FIT = UNPACK / "Zee_fit"
REFS = REPO / "combination" / "combLieke" / "config" / "references.json"
TT_REF = DATA_DIR / "ztautau_reference.json"
DATASET_EE = "CMS 2016 Open Data, SingleElectron Run2016G+H, NanoAODv9 (records 30529, 30562)"

# TRExFitter titles (z-ee/fit.config) -> sample keys; fit.config names in the histogram file
TITLE_TO_SAMPLE = {"DY #rightarrow ee": "DYee", "DY #rightarrow #tau#tau": "DYtautau", "Diboson": "Diboson",
                   "t#bar{t}": "TTbar", "W+jets": "WJets"}
FILE_NAME = {"DYee": "DYee", "DYtautau": "DY_tautau", "Diboson": "Diboson", "TTbar": "ttbar", "WJets": "Wjets"}
STACK_ORDER = ["WJets", "TTbar", "Diboson", "DYtautau", "DYee"]            # bottom-up
NPS = ["Lumi", "Pileup", "L1Prefiring", "ElectronRECO", "ElectronID", "PDF", "QCDScale", "PS_FSR", "ElectronES",
       "ElectronER", "XS_DYtautau", "XS_TTbar", "XS_VV", "XS_WJets"]
# names on screen (06 B5: descriptive, never the fit-internal code)
NP_LABEL = {"Lumi": "luminosity", "Pileup": "pile-up", "L1Prefiring": "L1 prefiring", "ElectronRECO": "electron reco",
            "ElectronID": "electron ID", "PDF": "PDF", "QCDScale": "QCD scale", "PS_FSR": "FSR (shower)",
            "ElectronES": "electron energy scale", "ElectronER": "electron energy res.",
            "XS_DYtautau": "σ(Z→ττ)", "XS_TTbar": "σ(tt̄)", "XS_VV": "σ(VV)", "XS_WJets": "σ(W+jets)"}


def unpack():
    if not (FIT / "Fits" / "Zee_fit.txt").exists():
        UNPACK.mkdir(parents=True, exist_ok=True)
        with tarfile.open(TARBALL) as tf:
            tf.extractall(UNPACK)


def read_plot_yaml(path):
    import yaml
    with open(path) as fh:
        y = yaml.safe_load(fh)
    samples = {TITLE_TO_SAMPLE[s["Name"]]: [float(v) for v in s["Yield"]] for s in y["Samples"]}
    tot = y["Total"][0]
    return {"samples": samples, "total": [float(v) for v in tot["Yield"]],
            "unc_up": [float(v) for v in tot["UncertaintyUp"]], "unc_down": [float(v) for v in tot["UncertaintyDown"]],
            "data": [int(v) for v in y["Data"][0]["Yield"]], "edges": [float(v) for v in y["Figure"][0]["BinEdges"]]}


def read_fit_txt(path):
    nps, in_np = {}, False
    for line in open(path):
        s = line.split()
        if line.startswith("NUISANCE_PARAMETERS"):
            in_np = True
            continue
        if line.startswith("CORRELATION_MATRIX"):
            break
        if in_np and len(s) == 4:
            nps[s[0]] = (float(s[1]), float(s[2]), float(s[3]))
    return nps


def read_decomp(path):
    out = {}
    for line in open(path):
        if not line.strip():
            continue
        name, rest = (line.split("(")[0].rsplit(None, 1)) if "(" in line else (line.split()[0], line.split()[1])
        out[name.strip()] = float(rest)
    return out


def read_err(path):
    out = {}
    for line in open(path):
        s = line.split()
        if len(s) >= 2:
            out[s[0]] = float(s[1])
    return out


def build():
    import uproot
    unpack()
    pre = read_plot_yaml(FIT / "Plots" / "ee_SR_prefit.yaml")
    post = read_plot_yaml(FIT / "Plots" / "ee_SR_postfit.yaml")
    f = uproot.open(FIT / "Histograms" / "Zee_fit_histos.root")
    edges = f["ee_SR/Data_2016/nominal/ee_SR_Data_2016"].axis().edges().tolist()
    data_root = f["ee_SR/Data_2016/nominal/ee_SR_Data_2016"].values().tolist()
    nominal, templates = {}, {}
    for s in STACK_ORDER:
        fn = FILE_NAME[s]
        nominal[s] = f[f"ee_SR/{fn}/nominal/ee_SR_{fn}"].values().tolist()
        templates[s] = {}
        for np_ in NPS:
            key_up = f"ee_SR/{fn}/{np_}/ee_SR_{fn}_{np_}_Up"
            if key_up in f:
                templates[s][np_] = {"up": f[key_up].values().tolist(),
                                     "down": f[f"ee_SR/{fn}/{np_}/ee_SR_{fn}_{np_}_Down"].values().tolist()}
    fit = read_fit_txt(FIT / "Fits" / "Zee_fit.txt")
    mu, mu_up, mu_down = fit["mu_Z"]
    ref = float(open(FIT / "reference_cross_section.txt").read().split()[0])
    err = read_err(FIT / "Fits" / "Zee_fit_errDecomp_mu_Z.txt")
    groups = read_decomp(FIT / "Fits" / "Zee_fit_group_errDecomp_mu_Z.txt")

    nps = [{"name": n, "label": NP_LABEL[n], "pull": fit[n][0], "constraint": 0.5 * (fit[n][1] - fit[n][2])} for n in NPS]
    gammas = sorted(((int(k.rsplit("_", 1)[1]), v[0], v[1]) for k, v in fit.items() if k.startswith("gamma_stat")))

    # total prediction with one parameter at +-1 (others nominal): what the one-at-a-time clip shows
    tot_nom = np.sum([nominal[s] for s in STACK_ORDER], axis=0)
    one_at_a_time = {}
    for np_ in NPS:
        up = np.sum([templates[s][np_]["up"] if np_ in templates[s] else nominal[s] for s in STACK_ORDER], axis=0)
        dn = np.sum([templates[s][np_]["down"] if np_ in templates[s] else nominal[s] for s in STACK_ORDER], axis=0)
        one_at_a_time[np_] = {"total_rel_up": float(up.sum() / tot_nom.sum() - 1.0),
                              "total_rel_down": float(dn.sum() / tot_nom.sum() - 1.0),
                              "bin_rel_up_min": float((up / tot_nom - 1).min()), "bin_rel_up_max": float((up / tot_nom - 1).max())}

    # the model at the post-fit point (piecewise-linear templates x gammas) vs TRExFitter's own post-fit yields
    def model(mu_, theta, gam):
        out = np.zeros(len(tot_nom))
        for s in STACK_ORDER:
            v = np.array(nominal[s], dtype=float)
            for np_, t in theta.items():
                if np_ not in templates[s]:
                    continue
                d = np.array(templates[s][np_]["up"] if t >= 0 else templates[s][np_]["down"]) - np.array(nominal[s])
                v = v + abs(t) * d
            out += v * (mu_ if s == "DYee" else 1.0)
        return out * gam
    theta_hat = {n: fit[n][0] for n in NPS}
    gam = np.array([g[1] for g in gammas])
    m_post = model(mu, theta_hat, gam)
    post_tot = np.array(post["total"])

    n_data = int(sum(pre["data"]))
    b_pre = float(sum(np.sum(pre["samples"][s]) for s in STACK_ORDER if s != "DYee"))
    s_pre = float(np.sum(pre["samples"]["DYee"]))
    eff_acc = s_pre / (ref * LUMI_PB)
    sigma_count = (n_data - b_pre) / (eff_acc * LUMI_PB)

    refs = load_json(REFS)["combined"]
    tt = load_json(TT_REF)
    pub = {}
    for r in refs:
        e = r["errors"]
        pub[r["label"].lower()] = {"label": r["label"], "value": r["value"], "stat": e["stat"], "syst": e["syst"],
                                   "lumi": e["lumi"], "total": float(np.sqrt(e["stat"] ** 2 + e["syst"] ** 2 + e["lumi"] ** 2)),
                                   "window": r["window"], "detail": r["detail"], "citation": r["citation"]}
    th = tt["theory"]
    theory = {"value": ref, "unc_up": ref * th["unc_rel_up"], "unc_down": ref * th["unc_rel_down"],
              "unc_rel_up": th["unc_rel_up"], "unc_rel_down": th["unc_rel_down"], "unc_source": th["unc_source"],
              "label": "aMC@NLO, NNLO norm.",
              "note": "sigma(Z/gamma* -> ee, 60 < m_LHE < 120) of the aMC@NLO DYJetsToLL_M-50 sample normalised to 6077.22 pb "
                      "(z-ee/Zee_fit reference_cross_section.txt, z-ee.ipynb cell 11); band = the relative uncertainty the "
                      "tau tau chapter uses (ztautau_reference.json theory.unc_rel_up/down)"}

    return {
        "provenance": provenance("extract_zee_fit.py", [TARBALL, REFS, TT_REF], dataset=DATASET_EE,
                                 version="z-ee delivery 16 Sep 2026 14:44 (Zee_fit.tar.gz), unchanged by the freeze (docs/FREEZE.md)",
                                 note="the channel's own TRExFitter fit; the combination's ee line (2094 +122 -114 pb) uses a "
                                      "shape/normalisation split and is not shown in this chapter (deck owner, 17 Sep 2026)"),
        "edges": edges,
        "data": pre["data"],
        "stack_order": STACK_ORDER,
        "prefit": {"samples": pre["samples"], "total": pre["total"], "unc_up": pre["unc_up"], "unc_down": pre["unc_down"],
                   "yields": {s: float(np.sum(pre["samples"][s])) for s in STACK_ORDER}},
        "postfit": {"samples": post["samples"], "total": post["total"], "unc_up": post["unc_up"],
                    "yields": {s: float(np.sum(post["samples"][s])) for s in STACK_ORDER}},
        "nominal_root": nominal,
        "data_root": data_root,
        "templates": templates,
        "template_note": "+-1 sigma templates exactly as the fit used them (Histograms/Zee_fit_histos.root, after TRExFitter's "
                         "symmetrisation and smoothing); a scene may interpolate linearly between nominal and a template "
                         "to draw an intermediate parameter value (display only)",
        "one_at_a_time": one_at_a_time,
        "nps": nps,
        "gammas": [{"bin": b, "value": v, "err": e} for b, v, e in gammas],
        "mu": {"value": mu, "err_up": mu_up, "err_down": -mu_down, "stat": err["STAT_ERROR"], "syst": err["SYST_ERROR"],
               "mcstat": err["MCSTAT_ERROR"], "totsyst": err["TOTSYST_ERROR"]},
        "sigma": {"value": mu * ref, "err": 0.5 * (mu_up - mu_down) * ref, "stat": err["STAT_ERROR"] * ref,
                  "syst": err["TOTSYST_ERROR"] * ref, "ref": ref, "window": "60-120"},
        "groups_rel": groups,
        "counting": {"n_data": n_data, "n_bkg": b_pre, "n_sig_pred": s_pre, "eff_acc": eff_acc, "lumi_pb": LUMI_PB,
                     "sigma": sigma_count,
                     "note": "prefit: sigma = (N - B) / (eps A L), eps A = N_DYee^pred / (sigma_ref L); no nuisance parameters"},
        "model_check": {"max_rel_dev_model_vs_postfit": float(np.max(np.abs(m_post / post_tot - 1.0))),
                        "note": "piecewise-linear templates x gammas at the post-fit point vs Plots/ee_SR_postfit.yaml"},
        "theory": theory,
        "published": pub,
    }


def verify(d, ck: Checker):
    pre, post = d["prefit"], d["postfit"]
    ck.check("30 bins, 60-120 GeV, 2 GeV", d["edges"], np.arange(60, 121, 2), source="z-ee/z-ee.ipynb cell 12")
    ck.check("N data = 6,320,097 (TRExFitter yields table)", sum(d["data"]), 6320097, source="Tables/Yields.txt")
    ck.check("data yaml == data histogram", d["data"], d["data_root"], source="Histograms/Zee_fit_histos.root")
    for s, want in (("DYee", 6468985.9), ("DYtautau", 6681.1), ("Diboson", 13232.4), ("TTbar", 21616.6), ("WJets", 2826.6)):
        ck.check(f"prefit {s}", pre["yields"][s], want, tol=0.1, source="Plots/ee_SR.png legend, Tables/Yields.txt")
        ck.check(f"prefit {s} yaml == nominal template", pre["samples"][s], d["nominal_root"][s], tol=1e-6, rel=True,
                 source="Histograms/Zee_fit_histos.root")
    for s, want in (("DYee", 6275350.4), ("DYtautau", 6221.7), ("Diboson", 13616.4), ("TTbar", 22678.7), ("WJets", 2264.9)):
        ck.check(f"postfit {s}", post["yields"][s], want, tol=0.1, source="Plots/ee_SR_postFit.png legend")
    ck.check("prefit total 6,513,342.6", sum(pre["total"]), 6513342.6, tol=0.1, source="Plots/ee_SR.png")
    ck.check("postfit total 6,320,132.1", sum(post["total"]), 6320132.1, tol=0.1, source="Plots/ee_SR_postFit.png")
    ck.check("mu_Z = 0.942 +- 0.015", [round(d["mu"]["value"], 3), round(d["mu"]["err_up"], 3), round(d["mu"]["err_down"], 3)],
             [0.942, 0.015, 0.015], tol=1e-9, source="docs/FREEZE.md Z -> ee, Fits/Zee_fit.txt")
    ck.check("reference 1954.10 pb", d["sigma"]["ref"], 1954.1032472811223, tol=1e-6, source="reference_cross_section.txt")
    ck.check("sigma = 1840.8 +- 29.9 pb", [d["sigma"]["value"], d["sigma"]["err"]], [1840.8, 29.9], tol=0.05,
             source="docs/FREEZE.md Z -> ee")
    ck.check("stat (+) syst (+) MC stat == total", np.sqrt(d["mu"]["stat"] ** 2 + d["mu"]["syst"] ** 2 + d["mu"]["mcstat"] ** 2),
             d["mu"]["err_up"], tol=2e-5,
             source="Fits/Zee_fit_errDecomp_mu_Z.txt")
    ck.check("luminosity group 1.15 %", d["groups_rel"]["Luminosity"], 0.0114574, tol=1e-7, source="group_errDecomp")
    names = [n["name"] for n in d["nps"]]
    ck.check_true("14 nuisance parameters as in z-ee/fit.config", len(names) == 14, source="z-ee/fit.config")
    pulls = {n["name"]: n for n in d["nps"]}
    ck.check("pulls ElectronID / L1Prefiring / QCDScale / Pileup", [pulls[k]["pull"] for k in ("ElectronID", "L1Prefiring", "QCDScale", "Pileup")],
             [0.594403, 2.02122, -1.97745, -1.57884], tol=1e-6, source="Fits/Zee_fit.txt, Pulls/All/NuisPar.png")
    ck.check("ElectronID constraint 0.080", pulls["ElectronID"]["constraint"], 0.079595, tol=1e-6, source="Fits/Zee_fit.txt")
    ck.check("30 gammas", len(d["gammas"]), 30, source="Fits/Zee_fit.txt")
    oat = d["one_at_a_time"]
    ck.check("luminosity template +-1.2 %", [oat["Lumi"]["total_rel_up"], oat["Lumi"]["total_rel_down"]], [0.012, -0.012],
             tol=1e-6, source="z-ee/fit.config Lumi OverallUp 0.012")
    ck.check("electron ID template +-5.8 % (the ECAL-gap artefact, combLieke/README.md)", oat["ElectronID"]["total_rel_up"],
             0.0583, tol=0.0005, source="combination/combLieke/README.md 'The ee channel' (5.9 %)")
    ck.check_true("energy scale / resolution templates have no effect (not animated)",
                  abs(oat["ElectronES"]["total_rel_up"]) < 1e-4 and abs(oat["ElectronER"]["total_rel_up"]) < 1e-4,
                  source="Histograms/Zee_fit_histos.root")
    ck.check("model(mu, theta, gamma) reproduces the post-fit yields", d["model_check"]["max_rel_dev_model_vs_postfit"], 0.0,
             tol=0.01, source="Plots/ee_SR_postfit.yaml", soft=True)
    c = d["counting"]
    ck.check("counting: B = 44,356.7 (prefit backgrounds)", c["n_bkg"], 44356.7, tol=0.1, source="Tables/Yields.txt")
    ck.check("counting: eps A = 0.2019", round(c["eff_acc"], 4), 0.2019, tol=1e-9, source="computed")
    ck.check("counting: sigma = 1895.7 pb", round(c["sigma"], 1), 1895.7, tol=1e-9, source="computed")
    p = d["published"]
    ck.check("CMS 1952 +- 4 +- 18 +- 45", [p["cms"]["value"], p["cms"]["stat"], p["cms"]["syst"], p["cms"]["lumi"]],
             [1952, 4, 18, 45], source="combLieke/config/references.json (arXiv:2408.03744 Table 13)")
    ck.check("ATLAS 1981 +- 7 +- 38 +- 42", [p["atlas"]["value"], p["atlas"]["stat"], p["atlas"]["syst"], p["atlas"]["lumi"]],
             [1981, 7, 38, 42], source="combLieke/config/references.json (arXiv:1603.09222 Table 3)")
    tt = load_json(TT_REF)
    ck.check("published totals == the tau tau chapter's", [p["cms"]["total"], p["atlas"]["total"]],
             [tt["cms"]["total"], tt["atlas"]["total"]], tol=1e-9, source="ztautau_reference.json")
    ck.check("prediction band +15 -21 pb", [round(d["theory"]["unc_up"]), round(d["theory"]["unc_down"])], [15, 21],
             source="arXiv:2408.03744 Table 5 (relative), ztautau_reference.json")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__), "zee_fit.json")
    finalize(ap.parse_args(), build, verify, "zee_fit")


if __name__ == "__main__":
    main()
