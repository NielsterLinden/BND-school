#!/usr/bin/env python
"""Step 5 -- binned profile-likelihood fit with TRExFitter v1.8.0 (docs/08-fit-and-results.md).

    python scripts/step5_fit.py [--ff-variant mcsub|nosub] [--skip-ranking] [--config-only]

One region per BDT category (config.REGIONS), m_tt in every region, POI mu_Z on the fiducial signal.
Writes fit/ztautau[_nosub].config (committed) from the histograms of step 4 using the shared
fitting/trexconfig.py, then runs in fit/:
    trex-fitter h / w / f / d / p / i / r        (read, workspace, fit, pre/post-fit plots, grouped
                                                   impacts, ranking)
    trex-fitter w f  StatOnly=TRUE:Suffix=_statOnly
    trex-fitter f    FitBlind=TRUE:Suffix=_asimov   (expected result on the Asimov data set)
and parses everything with fitting/run_trex.py into fit/results/ztautau[_nosub]_fit_result.json.

The binary is the build of the TRExFitter-v1.8.0 submodule; a worktree without its own build uses
the main checkout's (config.TREX_FALLBACK_HOME) -- the same binary the z-mumu channel uses.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fitting import run_trex, trexconfig as tc  # noqa: E402
from ztautau import analysis, config, samples  # noqa: E402

TITLES = {"DYtautau": "Z#rightarrow#tau#tau (fiducial)", "DYtautau_nonfid": "Z/#gamma*#rightarrow#tau#tau (non-fid.)",
          "DYee": "Z#rightarrowee", "DYmumu": "Z#rightarrow#mu#mu", "DYlowmass": "Z/#gamma*#rightarrowll (m<50)",
          "WJets": "W+jets", "TTbar": "t#bar{t}", "SingleTop": "single t", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ",
          "Fakes": "jet#rightarrow#tau_{h} (FF)"}
COLORS = {"DYtautau": 796, "DYtautau_nonfid": 800, "DYee": 857, "DYmumu": 862, "DYlowmass": 795, "WJets": 632, "TTbar": 616, "SingleTop": 617,
          "WW": 418, "WZ": 419, "ZZ": 420, "Fakes": 606}


def trex_environment():
    if shutil.which("trex-fitter"):
        return
    home = Path(os.environ.get("TREXFITTER_HOME", ""))
    if not (home / "build/bin/trex-fitter").exists():
        home = config.TREX_FALLBACK_HOME
    if not (home / "build/bin/trex-fitter").exists():
        raise RuntimeError("no TRExFitter build: run fitting/build_trexfitter.sh (after source setup.sh)")
    os.environ["TREXFITTER_HOME"] = str(home)
    os.environ["PATH"] = f"{home / 'build/bin'}:{os.environ['PATH']}"
    os.environ["LD_LIBRARY_PATH"] = f"{home / 'build/lib'}:{os.environ.get('LD_LIBRARY_PATH', '')}"


def build_config(job: str, histo_file: str, meta: dict) -> str:
    mc = [samples.SIGNAL] + samples.FIT_BACKGROUNDS
    mc_str = ",".join(mc)
    dy_str = f"{samples.SIGNAL},{samples.SIGNAL_NONFID}"
    blocks = [
        tc.job(job, Label="Z #rightarrow #tau_{h}#tau_{h}", CmeLabel="13 TeV", LumiLabel="16.4 fb^{-1}",
               ExperimentLabel="CMS", PlotLabel="Open Data", POI="mu_Z", ReadFrom="HIST", HistoPath="fitinputs",
               HistoFile=histo_file, OutputDir="results", MCstatThreshold=0.0, UseGammaPulls=True, DebugLevel=1,
               ImageFormat="png", SystControlPlots=True, DoSummaryPlot=True, DoTables=True,
               DoSignalRegionsPlot=False, DoPieChartPlot=True, RankingMaxNP=20, RankingPlot="SYSTS",
               HistoChecks="NOCRASH", SystPruningShape=0.001, SystPruningNorm=0.001, GetChi2="TRUE",
               SystCategoryTables=True, RatioYmax=1.5, RatioYmin=0.5, POIPrecision=3, SummaryPlotYmin=1,
               LegendNColumns=2),
        tc.fit("fit", FitType="SPLUSB", FitRegion="CRSR", UseMinos="mu_Z", NumCPU=4),
    ]
    for region, label in zip(meta["regions"], meta.get("region_labels", meta["regions"])):
        blocks.append(tc.region(region, Type="SIGNAL", HistoName=region, VariableTitle="m_{#tau#tau} [GeV]",
                                Label=label, ShortLabel=region.replace("tautau_", "")))
    blocks.append(tc.sample("Data", Type="DATA", Title="Data", HistoNameSuff="__Data"))
    smoothed = meta.get("smoothed", {})
    for s in mc + ["Fakes"]:
        blocks.append(tc.sample(s, Type="SIGNAL" if s == samples.SIGNAL else "BACKGROUND", Title=TITLES[s],
                                HistoNameSuff=f"__{s}", FillColor=COLORS[s], LineColor=1,
                                UseMCstat=s not in smoothed))
    blocks.append(tc.normfactor("mu_Z", Title="#mu_{Z}", Nominal=1, Min=0, Max=3, Samples=samples.SIGNAL))
    # normalisation uncertainties
    blocks.append(tc.overall_syst("Lumi", config.LUMI_REL_UNC, -config.LUMI_REL_UNC, mc_str, "Luminosity",
                                  title="Luminosity"))
    blocks.append(tc.overall_syst("XS_DYll", 0.05, -0.05, "DYee,DYmumu", "Background normalisation",
                                  title="#sigma(Z#rightarrowee/#mu#mu)"))
    for s in (samples.SIGNAL_NONFID, "DYlowmass", "WJets", "TTbar", "SingleTop", "WW", "WZ", "ZZ"):
        u = samples.XSEC_UNC[s]
        blocks.append(tc.overall_syst(f"XS_{s}", u, -u, s, "Background normalisation", title=f"#sigma({TITLES[s]})"))
    for s, info in smoothed.items():
        u = info["sr_stat_rel"]
        blocks.append(tc.overall_syst(f"MCStatNorm_{s}_tautau", u, -u, s, "Background normalisation",
                                      title=f"{TITLES[s]} MC stat. (norm.)"))
    # experimental shape/normalisation variations of the simulation
    blocks.append(tc.histo_syst("Pileup", mc_str, "Pileup", title="Pileup"))
    blocks.append(tc.histo_syst("L1Prefiring", mc_str, "L1 prefiring", title="L1 prefiring"))
    for dm in config.TAU_DMS:
        blocks.append(tc.histo_syst(f"TauID_DM{dm}", mc_str, "Tau", title=f"#tau_{{h}} ID DM{dm}", subcategory="Tau ID"))
        blocks.append(tc.histo_syst(f"TauTrigger_DM{dm}", mc_str, "Tau", title=f"#tau_{{h}} trigger DM{dm}",
                                    subcategory="Tau trigger"))
        blocks.append(tc.histo_syst(f"TauES_DM{dm}", mc_str, "Tau", title=f"#tau_{{h}} energy scale DM{dm}",
                                    smoothing=40, subcategory="Tau energy scale"))
    blocks.append(tc.histo_syst("TauFakeEle", mc_str, "Tau", title="e#rightarrow#tau_{h} ID", subcategory="Tau ID"))
    blocks.append(tc.histo_syst("TauFakeMu", mc_str, "Tau", title="#mu#rightarrow#tau_{h} ID", subcategory="Tau ID"))
    blocks.append(tc.histo_syst("MET_Unclustered", mc_str, "MET", title="MET unclustered energy", smoothing=40))
    # signal modelling: fiducial-normalised on the signal (only C varies), the same multipliers on the
    # non-fiducial part (its ratio to the fiducial yield varies)
    for s, t in (("QCDScale", "QCD scales"), ("PDF", "PDF"), ("PS_ISR", "PS ISR"), ("PS_FSR", "PS FSR")):
        blocks.append(tc.histo_syst(s, dy_str, "Signal modelling", title=t))
    if meta.get("has_sigmodel"):
        # fiducial C factor with LO madgraph instead of aMC@NLO: normalisation only, one-sided, this channel
        # only (the z-mumu `SigModel` is a powheg-vs-aMC@NLO shape and must not be correlated with it)
        blocks.append(tc.block("Systematic", "SigModel_tautau", Title="generator (LO vs NLO), fiducial C", Type="HISTO",
                               HistoNameSufUp="__SigModel_tautauUp", Samples=samples.SIGNAL, Symmetrisation="ONESIDED",
                               DropShapeIn="all", Category="Signal modelling"))
    # fakes
    for name, info in meta.get("osss_nps", {"FakeOSSS_tautau": {"region": None}}).items():
        blocks.append(tc.histo_syst(name, "Fakes", "Fakes", title="FF OS/SS extrapolation " + (info["region"] or "").replace("tautau_", ""),
                                    regions=info["region"]))
    for name, info in meta.get("closure_nps", {}).items():
        cat, tag = name.split("_")[-2:]
        blocks.append(tc.histo_syst(name, "Fakes", "Fakes", title=f"FF non-closure {info['region'].replace('tautau_', '')} "
                                    f"{'m<' if tag == 'lo' else 'm>'}{config.FF_CLOSURE_MASS_SPLIT:.0f}", regions=info["region"]))
    header = (f"Z -> tau_h tau_h, generated by z-tautau/scripts/step5_fit.py (do not edit by hand); tau_h working point {config.TAU_WP}\n"
              f"fake-factor variant: {meta['variant']}; fit variable {meta['fit_variable']} in {len(meta['regions'])} BDT "
              f"categories; TRExFitter v1.8.0")
    return tc.render(blocks, header=header)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    nominal_variant = "mcsub" if config.FF_SUBTRACT_MC else "nosub"
    ap.add_argument("--ff-variant", choices=["mcsub", "nosub"], default=nominal_variant)
    ap.add_argument("--skip-ranking", action="store_true")
    ap.add_argument("--config-only", action="store_true")
    args = ap.parse_args()
    suffix = "" if args.ff_variant == nominal_variant else f"_{args.ff_variant}"
    job = f"{config.JOB}{suffix}"
    fitdir = config.FIT_DIR
    meta = json.loads(Path(f"{fitdir}/fitinputs/{config.JOB}{suffix}.root.meta.json").read_text())
    cfg = fitdir / f"{job}.config"
    cfg.write_text(build_config(job, f"{config.JOB}{suffix}", meta))
    print(f"config -> {cfg}")
    if args.config_only:
        return
    trex_environment()
    logs = f"results/{job}/logs"
    for actions in ("h", "w", "f", "d", "p", "i") + (() if args.skip_ranking else ("r",)):
        run_trex.run(cfg.name, actions, cwd=fitdir, log=f"{logs}/{actions}.log")
    run_trex.run(cfg.name, "wf", options="StatOnly=TRUE:Suffix=_statOnly", cwd=fitdir, log=f"{logs}/statonly.log")
    run_trex.run(cfg.name, "wf", options="FitBlind=TRUE:Suffix=_asimov", cwd=fitdir, log=f"{logs}/asimov.log")

    res = run_trex.summarise(fitdir / "results" / job, job, poi="mu_Z")
    res["gof"] = run_trex.parse_gof(fitdir / logs / "f.log")
    asimov = fitdir / "results" / job / "Fits" / f"{job}_asimov.txt"
    if asimov.exists():
        a = run_trex.parse_fit_txt(asimov)["nps"].get("mu_Z")
        res["poi_asimov"] = {"value": a[0], "err_up": a[1], "err_down": a[2]} if a else None
    sig = analysis.signal_prediction(samples.DY_INCLUSIVE)
    mu, up, dn = res["poi_value"], res["poi_err_up"], res["poi_err_down"]
    stat = res.get("poi_stat_only", {}).get("err")
    res["prediction"] = sig
    res["sigma_fid_pb"] = {"value": mu * sig["sigma_fid_pb"], "up": up * sig["sigma_fid_pb"], "down": dn * sig["sigma_fid_pb"],
                           "stat": stat * sig["sigma_fid_pb"] if stat else None}
    s60 = sig["sigma_tautau_60_120_pb"]
    a_unc = float(sum(v ** 2 for v in sig["A_unc"].values()) ** 0.5)
    res["sigma_60_120_pb"] = {"value": mu * s60, "up": up * s60, "down": dn * s60, "stat": stat * s60 if stat else None,
                              "acceptance_rel_unc": a_unc,
                              "note": "sigma(pp -> Z/gamma* -> tautau, 60 < m_LHE < 120 GeV) = mu x prediction; the "
                                      "theory variations in the fit move only C, the acceptance uncertainty is separate"}
    res["lumi_pb"] = config.LUMI_PB
    res["ff_variant"] = args.ff_variant
    out = fitdir / "results" / f"{job}_fit_result.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    print(f"mu_Z = {mu:.4f} +{up:.4f} -{dn:.4f}  (stat {stat})")
    print(f"sigma_fid = {res['sigma_fid_pb']['value']:.3f} pb,  sigma(60-120) = {res['sigma_60_120_pb']['value']:.1f} pb")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
