#!/usr/bin/env python
"""Step 5 (v4) -- the four-channel profile-likelihood fit with TRExFitter v1.8.0 (docs/10-v4-plan.md section 8).

    python scripts/step5_fit_v4.py [--skip-ranking] [--config-only] [--channels tautau mutau etau emu]

Regions tautau_SR0/1/2, mutau_SR, etau_SR, emu_SR, emu_CRtt; POI mu_Z on the Z/gamma* -> tautau (60-120 GeV)
templates of every channel; free NormFactors TauIDSF_DM0/1/10/11 (nominal = TauPOG value) on the templates with
one genuine tau_h of that decay mode and their products (Expression) on the tau_h tau_h templates with two;
mu_ttbar on ttbar. The systematic registry of the fit inputs (step 4) becomes the Systematic blocks.
Writes fit_v4/ztautau_v4.config and runs h/w/f/d/p/i(/r) + stat-only + Asimov into fit_v4/results/.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fitting import run_trex, trexconfig as tc  # noqa: E402
from ztautau import analysis_v4 as an, config  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step5_fit import trex_environment  # noqa: E402

COLORS = {"DYtautau": 796, "DYtautau_out": 800, "DYee": 857, "DYmumu": 862, "DYlowmass": 795, "WJets": 632, "TTbar": 616, "SingleTop": 617,
          "WW": 418, "WZ": 419, "ZZ": 420, "Fakes": 606}
REGION_LABELS = {"tautau_SR0": "#tau_{h}#tau_{h}, BDT < 0.55", "tautau_SR1": "#tau_{h}#tau_{h}, 0.55 < BDT < 0.90", "tautau_SR2": "#tau_{h}#tau_{h}, BDT > 0.90",
                 "mutau_SR": "#mu#tau_{h}", "etau_SR": "e#tau_{h}", "emu_SR": "e#mu", "emu_CRtt": "e#mu t#bar{t} CR"}
for _ch, _lab in (("mutau", "#mu#tau_{h}"), ("etau", "e#tau_{h}")):
    for _dm in config.TAU_DMS:
        REGION_LABELS[f"{_ch}_SR_dm{_dm}"] = f"{_lab}, DM {_dm}"


def build_config(job: str, meta: dict, channels, fix_tauid: bool = False) -> str:
    regions = [r for r in meta["regions"] if r.split("_")[0] in channels]
    samples_ = {t: info for t, info in meta["samples"].items() if any(r in regions for r in info["regions"])}
    titles = meta["titles"]
    blocks = [
        tc.job(job, Label="Z #rightarrow #tau#tau (#tau_{h}#tau_{h}, #mu#tau_{h}, e#tau_{h}, e#mu)", CmeLabel="13 TeV", LumiLabel="16.4 fb^{-1}",
               ExperimentLabel="CMS", PlotLabel="Open Data", POI="mu_Z", ReadFrom="HIST", HistoPath="fitinputs",
               HistoFile=job, OutputDir="results", MCstatThreshold=0.0, UseGammaPulls=True, DebugLevel=1,
               ImageFormat="png", SystControlPlots=False, DoSummaryPlot=True, DoTables=True, DoSignalRegionsPlot=False,
               DoPieChartPlot=True, RankingMaxNP=25, RankingPlot="SYSTS", HistoChecks="NOCRASH", SystPruningShape=0.001,
               SystPruningNorm=0.001, GetChi2="TRUE", SystCategoryTables=True, RatioYmax=1.5, RatioYmin=0.5, POIPrecision=3,
               SummaryPlotYmin=1, LegendNColumns=2, PlotOptions="NOSIG,NOXERR", SummaryPlotRegions=",".join(regions)),
        tc.fit("fit", FitType="SPLUSB", FitRegion="CRSR", UseMinos="mu_Z", NumCPU=6),
    ]
    edges = list(meta["bins"]["tautau_SR0"]) if "tautau_SR0" in regions else []
    drop = [i + 1 for i in range(len(edges) - 1) if edges[i + 1] <= meta["sideband_mtt_min"]]
    for r in regions:
        blocks.append(tc.region(r, Type="CONTROL" if r.endswith("CRtt") else "SIGNAL", HistoName=r, VariableTitle="m_{#tau#tau} [GeV]",
                                Label=REGION_LABELS.get(r, r), ShortLabel=r, DropBins=",".join(map(str, drop)) if r == "tautau_SR0" and drop else None))
    blocks.append(tc.sample("Data", Type="DATA", Title="Data", HistoNameSuff="__Data"))
    smoothed = set(meta["yields_extra"].get("wjets_smoothed", {}).get("samples", []))
    for t, info in sorted(samples_.items(), key=lambda kv: (not kv[1]["is_signal"], kv[0])):
        base, key = info["base"], info["key"]
        blocks.append(tc.sample(t, Type="SIGNAL" if info["is_signal"] else "BACKGROUND", Title=f"{titles[base]} [{key}]", Group=titles[base],
                                HistoNameSuff=f"__{t}", FillColor=COLORS[base], LineColor=1, Regions=",".join(r for r in info["regions"] if r in regions),
                                UseMCstat=t not in smoothed))
    blocks.append(tc.sample("Fakes", Type="BACKGROUND", Title=titles["Fakes"], Group=titles["Fakes"], HistoNameSuff="__Fakes", FillColor=COLORS["Fakes"], LineColor=1))
    # normalisation factors
    sig = [t for t, i in samples_.items() if i["is_signal"]]
    blocks.append(tc.normfactor("mu_Z", Title="#mu_{Z}", Nominal=1, Min=0, Max=3, Samples=",".join(sig)))
    ttbar = [t for t in samples_ if t.startswith("TTbar_tDM")]
    if ttbar:
        blocks.append(tc.normfactor("mu_ttbar", Title="#mu_{t#bar{t}}", Nominal=1, Min=0, Max=3, Samples=",".join(ttbar)))
    pog = meta["tau_id_sf_pog"]
    single = {str(dm): [t for t, i in samples_.items() if i["key"] == str(dm)] for dm in config.TAU_DMS}
    pairs = sorted({i["key"] for i in samples_.values() if "_" in i["key"]})
    needed = {dm for pk in pairs for dm in pk.split("_")}
    for dm, ts in single.items():
        if not ts:
            if dm in needed:
                raise RuntimeError(f"TauIDSF_DM{dm} is needed by a product but has no single-tau template")
            continue                      # a NormFactor without samples would apply to *every* sample (e mu alone)
        nom = pog[dm][0]
        blocks.append(tc.normfactor(f"TauIDSF_DM{dm}", Title=f"#tau_{{h}} ID SF DM{dm}", Nominal=round(nom, 4), Min=0.5, Max=1.5,
                                    Samples=",".join(ts), Category="Tau ID (fitted)", Constant=True if fix_tauid else None))
    for pk in pairs:
        a, b = pk.split("_")
        ts = [t for t, i in samples_.items() if i["key"] == pk]
        blocks.append(tc.normfactor(f"TauIDSF_{a}x{b}", Title=f"#tau_{{h}} ID SF DM{a} #times DM{b}", Samples=",".join(ts),
                                    Expression=f"TauIDSF_DM{a}*TauIDSF_DM{b}:TauIDSF_DM{a},TauIDSF_DM{b}" if a != b else f"TauIDSF_DM{a}*TauIDSF_DM{a}:TauIDSF_DM{a}"))
    # systematics from the registry
    for name, s in meta["systs"].items():
        regs = [r for r in s["regions"] if r in regions]
        smp = [x for x in s["samples"] if x in samples_ or x == "Fakes"]
        if not regs or not smp:
            continue
        if s["type"] == "OVERALL":
            blocks.append(tc.overall_syst(name, s["up"], s["down"], ",".join(smp), s["category"], title=s["title"], regions=",".join(regs), subcategory=s["subcategory"]))
        else:
            blocks.append(tc.histo_syst(name, ",".join(smp), s["category"], title=s["title"], regions=",".join(regs), smoothing=s["smoothing"], subcategory=s["subcategory"]))
    header = (f"Z -> tautau v4 (four channels), generated by z-tautau/scripts/step5_fit_v4.py (do not edit by hand); tau_h working point {config.TAU_WP}\n"
              f"tau_h ID scale factors free per decay mode (nominal = TauPOG), tau energy scale prior {100 * meta['tes_prior']:.0f}%; TRExFitter v1.8.0")
    return tc.render(blocks, header=header)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-ranking", action="store_true")
    ap.add_argument("--config-only", action="store_true")
    ap.add_argument("--channels", nargs="*", default=config.CHANNELS)
    ap.add_argument("--job", default=config.JOB_V4)
    ap.add_argument("--fix-tauid", action="store_true", help="tau_h ID scale factors fixed to the TauPOG values (single-channel cross-checks)")
    args = ap.parse_args()
    job = args.job
    fitdir = config.FIT_DIR_V4
    # a single channel reads its own exported file (step 4b) when it exists, so its config, fit inputs and
    # workspace are self-contained for the MultiFit combination; otherwise the combined file
    histo_file = config.JOB_V4
    if len(args.channels) == 1 and (fitdir / "fitinputs" / f"{config.JOB_V4}_{args.channels[0]}.root").exists():
        histo_file = f"{config.JOB_V4}_{args.channels[0]}"
    meta = json.loads(Path(f"{fitdir}/fitinputs/{histo_file}.root.meta.json").read_text())
    cfg = fitdir / f"{job}.config"
    txt = build_config(job, meta, args.channels, args.fix_tauid).replace(f'HistoFile: "{job}"', f'HistoFile: "{histo_file}"')
    txt = re.sub(r'Expression: "([^"]+)"', r"Expression: \1", txt)       # TRExFitter splits the value on ':' before unquoting
    cfg.write_text(txt)
    print(f"config -> {cfg}  ({len(meta['samples'])} templates, {len(meta['systs'])} systematics)")
    if args.config_only:
        return
    trex_environment()
    logs = f"results/{job}/logs"
    for actions in ("h", "w", "f", "d", "p", "i") + (() if args.skip_ranking else ("r",)):
        run_trex.run(cfg.name, actions, cwd=fitdir, log=f"{logs}/{actions}.log")
    # keep the workspace under the conventional name as well (TRExFitter adds "allBinsFitRegions" with DropBins)
    ws = fitdir / "results" / job / "RooStats"
    for p in ws.glob(f"{job}_allBinsFitRegions_combined_{job}_model.root"):
        shutil.copy(p, ws / f"{job}_combined_{job}_model.root")
    run_trex.run(cfg.name, "wf", options="StatOnly=TRUE:Suffix=_statOnly", cwd=fitdir, log=f"{logs}/statonly.log")
    run_trex.run(cfg.name, "wf", options="FitBlind=TRUE:Suffix=_asimov", cwd=fitdir, log=f"{logs}/asimov.log")
    res = run_trex.summarise(fitdir / "results" / job, job, poi="mu_Z")
    res["gof"] = run_trex.parse_gof(fitdir / logs / "f.log")
    asimov = fitdir / "results" / job / "Fits" / f"{job}_asimov.txt"
    if asimov.exists():
        a = run_trex.parse_fit_txt(asimov)["nps"].get("mu_Z")
        res["poi_asimov"] = {"value": a[0], "err_up": a[1], "err_down": a[2]} if a else None
    fit_txt = fitdir / "results" / job / "Fits" / f"{job}.txt"
    nps = run_trex.parse_fit_txt(fit_txt)["nps"] if fit_txt.exists() else {}
    res["tau_id_sf"] = {f"DM{dm}": {"value": nps[f"TauIDSF_DM{dm}"][0], "err_up": nps[f"TauIDSF_DM{dm}"][1], "err_down": nps[f"TauIDSF_DM{dm}"][2],
                                    "pog": meta["tau_id_sf_pog"][str(dm)]} for dm in config.TAU_DMS if f"TauIDSF_DM{dm}" in nps}
    res["tau_es"] = {f"DM{dm}": {"pull": nps[f"TauES_DM{dm}"][0], "constraint": nps[f"TauES_DM{dm}"][1], "prior_pct": 100 * meta["tes_prior"]}
                     for dm in config.TAU_DMS if f"TauES_DM{dm}" in nps}
    res["mu_ttbar"] = nps.get("mu_ttbar")
    sig = meta["signal_prediction"]
    mu, up, dn = res["poi_value"], res["poi_err_up"], res["poi_err_down"]
    stat = res.get("poi_stat_only", {}).get("err")
    s60 = sig["sigma_tautau_60_120_pb"]
    res["prediction"] = sig
    res["sigma_60_120_pb"] = {"value": mu * s60, "up": up * s60, "down": dn * s60, "stat": stat * s60 if stat else None,
                              "note": "sigma(pp -> Z/gamma* -> tautau, 60 < m_LHE < 120 GeV) = mu_Z x aMC@NLO prediction normalised to 6077.22 pb (m > 50)"}
    res["lumi_pb"] = config.LUMI_PB
    res["channels"] = args.channels
    res["tau_id_fixed"] = args.fix_tauid
    out = fitdir / "results" / f"{job}_fit_result.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    print(f"mu_Z = {mu:.4f} +{up:.4f} -{dn:.4f}  (stat {stat});  sigma(60-120) = {res['sigma_60_120_pb']['value']:.1f} pb")
    for dm, v in res["tau_id_sf"].items():
        print(f"  TauIDSF {dm}: {v['value']:.3f} +{v['err_up']:.3f} -{v['err_down']:.3f}  (POG {v['pog'][0]:.3f} +- {v['pog'][1]:.3f})")
    for dm, v in res["tau_es"].items():
        print(f"  TauES {dm}: pull {v['pull']:+.2f}, constraint {v['constraint']:.2f} x 3% = {100 * meta['tes_prior'] * v['constraint']:.2f}%")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
