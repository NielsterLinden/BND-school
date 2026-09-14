#!/usr/bin/env python
"""v2 step 5 -- TRExFitter inputs, config, fit, and the cross section.

    python scripts/v2_5_fit.py [--no-fit] [--emu-control]

Reads output/v2/histograms.pkl, output/v2/fakes.json, output/v2/gensums.json and the
acceptance numbers; writes fit/fitinputs/zmumu.root, fit/zmumu.config, runs trex-fitter
(h, w, f, dp, r, i and a stat-only fit) into fit/results/zmumu/, and writes
fit/results/zmumu_fit_result.json with mu, its uncertainty breakdown and the cross sections.
See docs/14-fit-and-systematics.md and fitting/CONVENTIONS.md.
"""

from __future__ import annotations

import argparse
import json
import pickle
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ROOT))

import hist
import numpy as np

from fitting import run_trex, trexconfig as tc, trexhist
from zmumu import config, fakes, histograms as H, samples, weights

OUT = config.OUTPUT_DIR / "v2"
FIT = config.REPO_DIR / "fit"
JOB = "zmumu"
MC_SAMPLES = ["DYmumu", "DYtautau", "DYee", "TTbar", "SingleTop", "WW", "WZ", "ZZ"]
ALL_MC_NP = ["Pileup", "L1Prefiring", "MuonID", "MuonIso", "MuonTrigger", "MuonScale", "MuonRes"]
CATEGORY = {"Pileup": "Pileup", "L1Prefiring": "L1 prefiring", "MuonID": "Muon efficiency", "MuonIso": "Muon efficiency",
            "MuonTrigger": "Muon efficiency", "MuonScale": "Muon momentum", "MuonRes": "Muon momentum"}


def to_hist(values, variances, edges):
    h = hist.Hist(hist.axis.Regular(len(edges) - 1, edges[0], edges[-1], name="x"), storage=hist.storage.Weight())
    v = h.view()
    v.value[...] = values
    v.variance[...] = variances
    return h


def rebin(values, factor):
    return values.reshape(-1, factor).sum(axis=1)


def build_inputs(hists_all, fakes_res, gens, regions_used):
    """{histogram name: hist.Hist} following fitting/CONVENTIONS.md, plus meta."""
    out, meta = {}, {}
    for region in regions_used:
        e = H.edges(region, "mass_fit")
        key = lambda s, v="nominal": f"{s}|{region}|mass_fit|{v}"
        # data
        out[trexhist.hname(region, "Data")] = to_hist(hists_all[key("Data")], hists_all[key("Data") + "|w2"], e)
        for smp in MC_SAMPLES:
            if key(smp) not in hists_all or np.sum(hists_all[key(smp)]) <= 0:
                continue
            out[trexhist.hname(region, smp)] = to_hist(hists_all[key(smp)], hists_all[key(smp) + "|w2"], e)
            for np_ in ALL_MC_NP:
                for d in ("Up", "Down"):
                    k = key(smp, np_ + d)
                    if k in hists_all:
                        out[trexhist.hname(region, smp, np_, d)] = to_hist(hists_all[k], hists_all[k + "|w2"], e)
        # signal theory variations (DYmumu only): envelopes from the weight members
        nom = hists_all.get(key("DYmumu"))
        if nom is not None:
            var_nom = hists_all[key("DYmumu") + "|w2"]
            th = {k: hists_all.get(key("DYmumu", f"theory_{k}")) for k in ("pdf", "scale", "ps")}
            if th["pdf"] is not None:
                up, down = weights.pdf_envelope(th["pdf"], nom)
                out[trexhist.hname(region, "DYmumu", "PDF", "Up")] = to_hist(up, var_nom, e)
                out[trexhist.hname(region, "DYmumu", "PDF", "Down")] = to_hist(down, var_nom, e)
                au, ad = weights.alphas_variation(th["pdf"], nom)
                out[trexhist.hname(region, "DYmumu", "AlphaS", "Up")] = to_hist(au, var_nom, e)
                out[trexhist.hname(region, "DYmumu", "AlphaS", "Down")] = to_hist(ad, var_nom, e)
            if th["scale"] is not None:
                up, down = weights.scale_envelope(th["scale"], nom)
                out[trexhist.hname(region, "DYmumu", "QCDScale", "Up")] = to_hist(up, var_nom, e)
                out[trexhist.hname(region, "DYmumu", "QCDScale", "Down")] = to_hist(down, var_nom, e)
            if th["ps"] is not None:
                for name, (u, d) in weights.ps_variations(th["ps"]).items():
                    out[trexhist.hname(region, "DYmumu", name, "Up")] = to_hist(u, var_nom, e)
                    out[trexhist.hname(region, "DYmumu", name, "Down")] = to_hist(d, var_nom, e)
            # SigModel: powheg template normalised to the same fiducial yield as the NLO nominal
            pw = hists_all.get(f"DYmumu_powheg|{region}|mass_fit|nominal")
            if pw is not None and "DY_powheg" in gens and "DY_NLO" in gens:
                g_nlo, g_pw = gens["DY_NLO"], gens["DY_powheg"]
                n_fid_expected = weights.LUMI_PB * config.DY_XSEC_PB * float(g_nlo["sumw_fid_dressed"]) / float(g_nlo["sumw"])
                scale = n_fid_expected / float(g_pw["sumw_fid_dressed"])
                out[trexhist.hname(region, "DYmumu", "SigModel", "Up")] = to_hist(
                    pw * scale, hists_all[f"DYmumu_powheg|{region}|mass_fit|nominal|w2"] * scale ** 2, e)
                if region == "SR":
                    meta["sigmodel_powheg_over_nlo"] = float(np.sum(pw * scale) / np.sum(nom))
        # fakes (SR only; the CRemu prompt subtraction is not part of the fake estimate)
        if region == "SR" and fakes_res is not None:
            f = np.array(fakes_res["templates"]["nominal"]); fv = np.array(fakes_res["templates"]["nominal_var"])
            factor = len(f) // (len(e) - 1)
            nom_f = np.clip(rebin(f, factor), 0, None); var_f = rebin(fv, factor)
            out[trexhist.hname(region, "Fakes")] = to_hist(nom_f, var_f, e)
            up = np.clip(rebin(np.array(fakes_res["templates"]["ff_stat_up"]), factor), 0, None)
            dn = np.clip(rebin(np.array(fakes_res["templates"]["ff_stat_down"]), factor), 0, None)
            out[trexhist.hname(region, "Fakes", "FakeStat_mumu", "Up")] = to_hist(up, var_f, e)
            out[trexhist.hname(region, "Fakes", "FakeStat_mumu", "Down")] = to_hist(dn, var_f, e)
            rel = min(max(fakes_res["yields"]["method_rel_unc"], 0.3), 1.0)
            out[trexhist.hname(region, "Fakes", "FakeMethod_mumu", "Up")] = to_hist(nom_f * (1 + rel), var_f, e)
            out[trexhist.hname(region, "Fakes", "FakeMethod_mumu", "Down")] = to_hist(nom_f * (1 - rel), var_f, e)
            meta["fakes_method_rel_unc"] = rel
    return out, meta


def make_config(regions_used, present, emu_control):
    mc = [s for s in MC_SAMPLES if trexhist.hname("SR", s) in present]
    blocks = [
        tc.job(JOB, ExperimentLabel="CMS Open Data", Label="Z #rightarrow #mu#mu", CmeLabel="13 TeV",
               LumiLabel="16.4 fb^{-1}", POI="mu_Z",
               ReadFrom="HIST", HistoPath="fitinputs", HistoFile=JOB, OutputDir="results", MCstatThreshold=0,
               UseGammaPulls=True, SystControlPlots=True, SystCategoryTables=True, SystLarge=0.5,
               SuppressNegativeBinWarnings=True, GetChi2=True, DoSummaryPlot=True, DoTables=True, RankingMaxNP=25,
               RankingPlot="ALL", ImageFormat="png", PlotOptions="YIELDS", RatioYmin=0.95, RatioYmax=1.05,
               RatioYminPostFit=0.98, RatioYmaxPostFit=1.02, DebugLevel=1),
        tc.fit("fit", FitType="SPLUSB", FitRegion="CRSR", UseMinos="mu_Z", NumCPU=4),
        tc.region("mumu_SR", Type="SIGNAL", HistoName="mumu_SR", VariableTitle="m_{#mu#mu} [GeV]",
                  Label="#mu#mu, 60 < m < 120 GeV", ShortLabel="SR", BinWidth=2, LogScale=True),
    ]
    if "CRemu" in regions_used:
        blocks.append(tc.region("mumu_CRemu", Type="CONTROL" if emu_control else "VALIDATION", HistoName="mumu_CRemu",
                                VariableTitle="m_{e#mu} [GeV]", Label="e#mu (OS)", ShortLabel="CR e#mu", BinWidth=5))
    blocks.append(tc.sample("Data", Type="DATA", Title="Data", HistoNameSuff="__Data"))
    colours = {"DYmumu": 800, "DYtautau": 616, "DYee": 400, "TTbar": 632, "SingleTop": 634, "WW": 432, "WZ": 433, "ZZ": 434}
    titles = {"DYmumu": "Z/#gamma* #rightarrow #mu#mu", "DYtautau": "Z/#gamma* #rightarrow #tau#tau",
              "DYee": "Z/#gamma* #rightarrow ee", "TTbar": "t#bar{t}", "SingleTop": "tW", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ"}
    for s in mc:
        regs = [r for r in ("mumu_SR", "mumu_CRemu") if trexhist.hname(r.split("_", 1)[1] if r != "mumu_SR" else "SR", s) in present]
        blocks.append(tc.sample(s, Type="SIGNAL" if s == "DYmumu" else "BACKGROUND", Title=titles[s], FillColor=colours[s],
                                LineColor=1, HistoNameSuff=f"__{s}", Regions=",".join(regs)))
    if trexhist.hname("SR", "Fakes") in present:
        blocks.append(tc.sample("Fakes", Type="BACKGROUND", Title="Non-prompt (fake factor)", FillColor=920, LineColor=1,
                                HistoNameSuff="__Fakes", Regions="mumu_SR"))
    blocks.append(tc.normfactor("mu_Z", Title="#mu_{Z}", Nominal=1, Min=0.5, Max=1.5, Samples="DYmumu"))
    if emu_control:
        blocks.append(tc.normfactor("mu_top", Title="#mu_{top}", Nominal=1, Min=0.0, Max=3.0, Samples="TTbar,SingleTop"))
    mc_list = ",".join(mc)
    blocks.append(tc.overall_syst("Lumi", 0.012, -0.012, mc_list, "Luminosity", title="Luminosity"))
    blocks.append(tc.overall_syst("MuonReco", 0.004, -0.004, mc_list, "Muon efficiency", title="Muon reconstruction"))
    for s, unc in samples.XSEC_UNC.items():
        if s in mc and not (emu_control and s in ("TTbar", "SingleTop")):
            blocks.append(tc.overall_syst(f"XS_{s}", unc, -unc, s, "Background normalisation", title=f"{titles[s]} cross section"))
    for np_ in ALL_MC_NP:
        smoothing = 40 if np_ in ("MuonScale", "MuonRes") else None
        blocks.append(tc.histo_syst(np_, mc_list, CATEGORY[np_], title=np_, smoothing=smoothing))
    for np_, title in (("PDF", "PDF (NNPDF3.1)"), ("AlphaS", "#alpha_{s}"), ("QCDScale", "QCD scales"),
                       ("PS_ISR", "Parton shower ISR"), ("PS_FSR", "Parton shower FSR")):
        if trexhist.hname("SR", "DYmumu", np_, "Up") in present:
            blocks.append(tc.histo_syst(np_, "DYmumu", "Signal modelling", title=title))
    if trexhist.hname("SR", "DYmumu", "SigModel", "Up") in present:
        blocks.append(tc.histo_syst("SigModel", "DYmumu", "Signal modelling", title="Generator (powheg vs aMC@NLO)",
                                    one_sided=True, smoothing=40))
    if trexhist.hname("SR", "Fakes") in present:
        blocks.append(tc.histo_syst("FakeStat_mumu", "Fakes", "Fakes", title="Fake factor stat.", regions="mumu_SR"))
        blocks.append(tc.histo_syst("FakeMethod_mumu", "Fakes", "Fakes", title="Fake factor method", regions="mumu_SR"))
    if "CRemu" in regions_used:
        blocks.append(tc.overall_syst("ElectronEff_mumu", 0.03, -0.03, mc_list, "Electron efficiency",
                                      title="Electron ID/iso (e#mu region)", regions="mumu_CRemu"))
    header = ("Z -> mumu cross section, CMS Open Data 2016 G+H, TRExFitter v1.8.0\n"
              "generated by scripts/v2_5_fit.py -- run from z-mumu/fit/:  trex-fitter <actions> zmumu.config")
    return tc.render(blocks, header=header)


def cross_sections(res, gens, acc):
    """sigma_fid, sigma(60-120), sigma(m>50) from mu and the acceptance numbers."""
    g = gens["DY_NLO"]
    frac_fid = float(g["sumw_fid_dressed"]) / float(g["sumw"])
    sigma_fid_pred = config.DY_XSEC_PB * frac_fid
    a_m50 = 3.0 * frac_fid
    a_60120 = float(g["sumw_fid_dressed"]) / float(g["sumw_lhe_mumu_60_120"])
    mu, up, down = res["poi_value"], res["poi_err_up"], res["poi_err_down"]
    tot = 0.5 * (up + down)
    grouped = res.get("grouped_impact", {})
    full_syst = grouped.get("FullSyst", 0.0)
    stat_cov = float(np.sqrt(max(tot ** 2 - full_syst ** 2, 0.0)))
    stat_fit = res.get("poi_stat_only", {}).get("err", stat_cov)
    lumi = grouped.get("Luminosity", 0.0)
    out = {"mu": mu, "mu_err_up": up, "mu_err_down": down, "mu_err": tot, "mu_stat_only_fit": stat_fit,
           "mu_stat_from_cov": stat_cov, "mu_syst": full_syst, "mu_lumi": lumi,
           "mu_syst_excl_lumi": float(np.sqrt(max(full_syst ** 2 - lumi ** 2, 0.0))),
           "grouped_impacts_mu": grouped, "sigma_fid_pred_pb": sigma_fid_pred, "A_m50": a_m50, "A_60_120": a_60120,
           "acceptance": acc}
    out["sigma_fid_pb"] = mu * sigma_fid_pred
    out["sigma_fid_stat_pb"] = stat_fit * sigma_fid_pred
    out["sigma_fid_syst_pb"] = out["mu_syst_excl_lumi"] * sigma_fid_pred
    out["sigma_fid_lumi_pb"] = lumi * sigma_fid_pred
    out["sigma_fid_tot_pb"] = tot * sigma_fid_pred
    a_rel = acc.get("A_rel_unc", 0.0)
    for name, a in (("60_120", a_60120), ("m50", a_m50)):
        s = out["sigma_fid_pb"] / a
        out[f"sigma_{name}_pb"] = s
        out[f"sigma_{name}_stat_pb"] = out["sigma_fid_stat_pb"] / a
        out[f"sigma_{name}_syst_pb"] = out["sigma_fid_syst_pb"] / a
        out[f"sigma_{name}_lumi_pb"] = out["sigma_fid_lumi_pb"] / a
        out[f"sigma_{name}_acc_pb"] = s * a_rel
        out[f"sigma_{name}_tot_pb"] = float(np.sqrt((out["sigma_fid_tot_pb"] / a) ** 2 + (s * a_rel) ** 2))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-fit", action="store_true", help="write inputs and config only")
    ap.add_argument("--emu-control", action="store_true", help="use the e-mu region as a control region with mu_top")
    ap.add_argument("--no-emu", action="store_true")
    args = ap.parse_args()

    hists_all = pickle.load(open(OUT / "histograms.pkl", "rb"))
    gens = json.load(open(OUT / "gensums.json"))
    fakes_res = json.load(open(OUT / "fakes.json")) if (OUT / "fakes.json").exists() else None
    regions_used = ["SR"] + ([] if args.no_emu or "Data|CRemu|mass_fit|nominal" not in hists_all else ["CRemu"])
    inputs, meta = build_inputs(hists_all, fakes_res, gens, regions_used)
    # rename regions to the channel-prefixed convention
    inputs = {n.replace("SR__", "mumu_SR__", 1).replace("CRemu__", "mumu_CRemu__", 1): h for n, h in inputs.items()}
    present = {n.replace("mumu_", "", 1) if n.startswith("mumu_") else n for n in inputs}
    acc = {}
    acc_pkl = config.DATA_DIR / "mc_acceptance_nlo_result.pkl"
    if acc_pkl.exists():
        r = pickle.load(open(acc_pkl, "rb"))
        acc = {k: float(r[k]) for k in ("A", "A_stat", "A_pdf_rel", "A_alphas_rel", "A_scale_rel", "A_lhe_mumu_60_120") if k in r}
        acc["A_rel_unc"] = float(np.sqrt((acc["A_stat"] / acc["A"]) ** 2 + acc["A_pdf_rel"] ** 2 + acc["A_alphas_rel"] ** 2 + acc["A_scale_rel"] ** 2))
    g = gens["DY_NLO"]
    meta.update(lumi_pb=weights.LUMI_PB, dy_xsec_pb=config.DY_XSEC_PB,
                sigma_fid_pred_pb=config.DY_XSEC_PB * float(g["sumw_fid_dressed"]) / float(g["sumw"]),
                A_m50=3.0 * float(g["sumw_fid_dressed"]) / float(g["sumw"]),
                A_60_120=float(g["sumw_fid_dressed"]) / float(g["sumw_lhe_mumu_60_120"]),
                n_fid_expected=weights.LUMI_PB * config.DY_XSEC_PB * float(g["sumw_fid_dressed"]) / float(g["sumw"]),
                acceptance=acc, regions=regions_used, samples=sorted({n.split("__")[1] for n in inputs}))
    sr_sig = inputs["mumu_SR__DYmumu"].values().sum()
    meta["C_factor"] = float(sr_sig / meta["n_fid_expected"])
    FIT.mkdir(parents=True, exist_ok=True)
    rep = trexhist.write_fitinputs(FIT / "fitinputs" / f"{JOB}.root", inputs, meta=meta)
    trexhist.check_fitinputs(FIT / "fitinputs" / f"{JOB}.root", list(inputs))
    print(f"[fit] {rep['n_hists']} histograms written; clipped negative bins in: {list(rep['clipped'])}")
    print(f"[fit] C = {meta['C_factor']:.4f}, sigma_fid^pred = {meta['sigma_fid_pred_pb']:.1f} pb, A(60-120) = {meta['A_60_120']:.4f}, A(m>50) = {meta['A_m50']:.4f}")
    (FIT / f"{JOB}.config").write_text(make_config(regions_used, present, args.emu_control))
    # counting cross-check
    data = inputs["mumu_SR__Data"].values().sum()
    bkg = sum(inputs[n].values().sum() for n in inputs if n.startswith("mumu_SR__") and n.count("__") == 1 and n not in ("mumu_SR__Data", "mumu_SR__DYmumu"))
    meta["counting"] = {"n_obs": float(data), "n_bkg": float(bkg), "sigma_fid_pb": float((data - bkg) / (weights.LUMI_PB * meta["C_factor"]))}
    print(f"[fit] counting: N_obs {data:,.0f}, N_bkg {bkg:,.0f} -> sigma_fid = {meta['counting']['sigma_fid_pb']:.1f} pb")
    if args.no_fit:
        return
    logs = FIT / "results" / JOB / "logs"
    if (FIT / "results" / JOB).exists():
        shutil.rmtree(FIT / "results" / JOB)
    for actions in ("h", "w", "f", "dp", "r", "i"):
        print(f"[fit] trex-fitter {actions}", flush=True)
        run_trex.run(f"{JOB}.config", actions, cwd=FIT, log=f"results/{JOB}/logs/{actions}.log")
    run_trex.run(f"{JOB}.config", "wf", options="StatOnly=TRUE:Suffix=_statOnly", cwd=FIT, log=f"results/{JOB}/logs/statonly.log")
    res = run_trex.summarise(FIT / "results" / JOB, JOB, poi="mu_Z")
    res["gof"] = run_trex.parse_gof(logs / "f.log")
    xs = cross_sections(res, gens, acc)
    result = {"channel": "zmumu", "poi": "mu_Z", **xs, "meta": meta, "nps": res["nps"],
              "ranking": res.get("ranking", [])[:25], "err_decomp": res.get("err_decomp"), "gof": res["gof"],
              "table_postfit": res.get("table_postfit"), "table_prefit": res.get("table_prefit")}
    with open(FIT / "results" / f"{JOB}_fit_result.json", "w") as fh:
        json.dump(result, fh, indent=1, default=lambda o: float(o))
    print(f"[fit] mu_Z = {xs['mu']:.4f} +{xs['mu_err_up']:.4f} -{xs['mu_err_down']:.4f} "
          f"(stat {xs['mu_stat_only_fit']:.4f}, syst {xs['mu_syst']:.4f} incl. lumi {xs['mu_lumi']:.4f}); GoF p = {res['gof']['gof_probability']}")
    print(f"[fit] sigma_fid = {xs['sigma_fid_pb']:.1f} +- {xs['sigma_fid_stat_pb']:.1f} (stat) +- {xs['sigma_fid_syst_pb']:.1f} (syst) +- {xs['sigma_fid_lumi_pb']:.1f} (lumi) pb")
    print(f"[fit] sigma(60<m<120) = {xs['sigma_60_120_pb']:.0f} +- {xs['sigma_60_120_tot_pb']:.0f} pb;  sigma(m>50) = {xs['sigma_m50_pb']:.0f} +- {xs['sigma_m50_tot_pb']:.0f} pb")


if __name__ == "__main__":
    main()
