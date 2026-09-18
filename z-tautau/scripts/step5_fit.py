#!/usr/bin/env python
"""Step 5 (v4) -- the four-channel profile-likelihood fit with TRExFitter v1.8.0 (docs/10-v4-plan.md section 8).

    python scripts/step5_fit.py [--skip-ranking] [--config-only] [--channels tautau mutau etau emu]

Regions tautau_SR0/1/2, mutau_SR, etau_SR, emu_SR, emu_CRtt; POI mu_Z on the Z/gamma* -> tautau (60-120 GeV)
templates of every channel; free NormFactors TauIDSF_DM0/1/10/11 (nominal = TauPOG value) on the templates with
one genuine tau_h of that decay mode and their products (Expression) on the tau_h tau_h templates with two;
mu_ttbar on ttbar. The systematic registry of the fit inputs (step 4) becomes the Systematic blocks. Every fit
of the measurement's model (four channels, one scale factor per decay mode) also carries `TauIDpT_tautau`, the
spread between that model and the pT-split one (`pt_model`): run the `ztautau_flatsf` and `ztautau_ptsplit`
jobs first (run_all.py and condor/orchestrator.sh do).
Writes fit_v4/ztautau_v4.config and runs h/w/f/d/p/i(/r) + stat-only + Asimov into fit_v4/results/.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import uproot  # noqa: E402

from fitting import run_trex, trexconfig as tc  # noqa: E402
from ztautau import analysis_v4 as an, config  # noqa: E402



def trex_environment():
    """Put the TRExFitter v1.8.0 build on PATH ($TREXFITTER_HOME, else the shared build next to the main checkout)."""
    if shutil.which("trex-fitter"):
        return
    home = Path(os.environ.get("TREXFITTER_HOME", ""))
    if not (home / "build/bin/trex-fitter").exists():
        home = config.TREX_FALLBACK_HOME
    if not (home / "build/bin/trex-fitter").exists():
        raise RuntimeError("no TRExFitter build: run fitting/build_trexfitter.sh (after source fitting/setup.sh)")
    os.environ["TREXFITTER_HOME"] = str(home)
    os.environ["PATH"] = f"{home / 'build/bin'}:{os.environ['PATH']}"
    os.environ["LD_LIBRARY_PATH"] = f"{home / 'build/lib'}:{os.environ.get('LD_LIBRARY_PATH', '')}"

COLORS = {"DYtautau": 796, "DYtautau_out": 800, "DYee": 857, "DYmumu": 862, "DYlowmass": 795, "WJets": 632, "TTbar": 616, "SingleTop": 617,
          "WW": 418, "WZ": 419, "ZZ": 420, "Fakes": 606}
REGION_LABELS = {"tautau_SR0": "#tau_{h}#tau_{h}, BDT < 0.55", "tautau_SR1": "#tau_{h}#tau_{h}, 0.55 < BDT < 0.90", "tautau_SR2": "#tau_{h}#tau_{h}, BDT > 0.90",
                 "mutau_SR": "#mu#tau_{h}", "etau_SR": "e#tau_{h}", "emu_SR": "e#mu", "emu_CRtt": "e#mu t#bar{t} CR"}
for _ch, _lab in (("mutau", "#mu#tau_{h}"), ("etau", "e#tau_{h}")):
    for _dm in config.TAU_DMS:
        REGION_LABELS[f"{_ch}_SR_dm{_dm}"] = f"{_lab}, DM {_dm}"


def empty_bins(histo_path: Path, regions) -> dict:
    """{region: [1-based bin indices]} of the bins with no data and no prediction at all.

    Such a bin carries no information, but with `MCstatThreshold: 0` its MC-statistics gamma is
    unconstrained and sits at zero, which makes the per-bin offset log(0): HESSE is then forced
    positive-definite and MINOS wanders into a region where every parameter is NaN. That is what made the
    Asimov fit of the first v4 round fail (review/REVIEW_v4.md finding 1) and what the combination had to work
    around in its own copy of these inputs. They are dropped from the fit instead.
    """
    out = {}
    with uproot.open(histo_path) as f:
        keys = [k.split(";")[0] for k in f.keys()]
        for r in regions:
            nominal = [k for k in keys if k.startswith(f"{r}__") and "__" not in k[len(r) + 2:]]
            if not nominal:
                continue
            total = None
            for k in nominal:
                v = np.asarray(f[k].values(), dtype=float)
                total = v.copy() if total is None else total + v
            # step 4 writes 1e-6 placeholders for templates a region does not really have
            out[r] = [i + 1 for i, v in enumerate(total) if abs(v) < 1e-3]
    return out


def select_regions(meta, channels, region_set="nominal"):
    """Regions of a fit: the channels asked for, in the region set asked for. 'nominal' is the measurement
    (one l tau_h region per decay mode); 'ptsplit' replaces the l tau_h regions by their pT(tau_h) split
    (review/REVIEW_v4.md finding 4) and keeps the other channels as they are."""
    sets = meta.get("region_sets", {})
    out = []
    for r in meta["regions"]:
        ch = r.split("_")[0]
        if ch not in channels:
            continue
        rset = sets.get(r, "nominal")
        if rset == region_set or (region_set != "nominal" and rset == "nominal" and ch not in ("mutau", "etau")):
            out.append(r)
    return out


TAUID_PT_NP = "TauIDpT_tautau"


def pt_model(fitdir: Path) -> dict:
    """Size of `TauIDpT_tautau`, the modelling uncertainty of the one-scale-factor-per-decay-mode assumption.

    The tau_h tau_h / (l tau_h)^2 lever assumes the tau_h ID scale factor is flat in pT from 30 GeV upwards.
    `ztautau_ptsplit` gives the l tau_h regions below 40 GeV their own scale factors and moves mu_Z by 6 %,
    1.8 times the uncertainty of the fit without this parameter (review/REVIEW_v4_RESPONSE.md section 6). The split is
    not the nominal model, so the relative difference between the two fits (`ztautau_flatsf`: this model
    without the parameter; `ztautau_ptsplit`) is carried as one OVERALL parameter on the signal. It is
    degenerate with mu_Z by construction: the central value stays, the uncertainty grows.
    """
    mu = {}
    for j in ("flatsf", "ptsplit"):
        path = fitdir / "results" / f"{config.JOB_V4}_{j}_fit_result.json"
        if not path.exists():
            raise SystemExit(f"{path} is missing: run the {config.JOB_V4}_flatsf (--no-tauid-pt) and "
                             f"{config.JOB_V4}_ptsplit jobs before a fit that carries {TAUID_PT_NP}")
        mu[j] = json.loads(path.read_text())["poi_value"]
    return {"name": TAUID_PT_NP, "rel": abs(mu["ptsplit"] - mu["flatsf"]) / mu["flatsf"],
            "mu_flatsf": mu["flatsf"], "mu_ptsplit": mu["ptsplit"]}


def build_config(job: str, meta: dict, channels, fix_tauid: bool = False, region_set: str = "nominal",
                 scale_systs: dict | None = None, fit_strategy: int | None = None, num_cpu: int = 6,
                 empty: dict | None = None, pt: dict | None = None) -> str:
    scale_systs = scale_systs or {}
    empty = empty or {}
    regions = select_regions(meta, channels, region_set)
    lowpt = [r for r in regions if "_SRlo_" in r]
    samples_ = {t: info for t, info in meta["samples"].items() if any(r in regions for r in info["regions"])}
    titles = meta["titles"]
    blocks = [
        tc.job(job, Label="Z #rightarrow #tau#tau (#tau_{h}#tau_{h}, #mu#tau_{h}, e#tau_{h}, e#mu)", CmeLabel="13 TeV", LumiLabel="16.4 fb^{-1}",
               ExperimentLabel="CMS", PlotLabel="Open Data", POI="mu_Z", ReadFrom="HIST", HistoPath="fitinputs",
               HistoFile=job, OutputDir="results", MCstatThreshold=0.0, UseGammaPulls=True, DebugLevel=1,
               ImageFormat="png", SystControlPlots=False, DoSummaryPlot=True, DoTables=True, DoSignalRegionsPlot=False,
               # PlotOptions has no NOSIG: with it, TRExFitter leaves the signal out of the plotted total,
               # of `h_tot_postFit` and of the per-bin Plots/<region>_postfit.yaml `Total` -- so every
               # data/prediction ratio built from those files excludes Z -> tautau, which is most of the
               # prediction in the signal regions.
               DoPieChartPlot=True, RankingMaxNP=25, RankingPlot="SYSTS", HistoChecks="NOCRASH", SystPruningShape=0.001,
               SystPruningNorm=0.001, GetChi2="TRUE", SystCategoryTables=True, RatioYmax=1.5, RatioYmin=0.5, POIPrecision=3,
               SummaryPlotYmin=1, LegendNColumns=2, PlotOptions="NOXERR", SummaryPlotRegions=",".join(regions)),
        tc.fit("fit", FitType="SPLUSB", FitRegion="CRSR", UseMinos="mu_Z", NumCPU=num_cpu, FitStrategy=fit_strategy,
               MaximumNumberFCNcalls=2000000 if fit_strategy == 3 else None),      # HESSE at strategy 3 needs more than the default
    ]
    edges = list(meta["bins"]["tautau_SR0"]) if "tautau_SR0" in regions else []
    sideband = [i + 1 for i in range(len(edges) - 1) if edges[i + 1] <= meta["sideband_mtt_min"]]
    for r in regions:
        drop = sorted(set(empty.get(r, [])) | (set(sideband) if r == "tautau_SR0" else set()))
        blocks.append(tc.region(r, Type="CONTROL" if r.endswith("CRtt") else "SIGNAL", HistoName=r, VariableTitle="m_{#tau#tau} [GeV]",
                                Label=REGION_LABELS.get(r, r), ShortLabel=r, DropBins=",".join(map(str, drop)) if drop else None))
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
    hipt = [r for r in regions if r not in lowpt]
    for dm, ts in single.items():
        if not ts:
            if dm in needed:
                raise RuntimeError(f"TauIDSF_DM{dm} is needed by a product but has no single-tau template")
            continue                      # a NormFactor without samples would apply to *every* sample (e mu alone)
        nom = pog[dm][0]
        blocks.append(tc.normfactor(f"TauIDSF_DM{dm}", Title=f"#tau_{{h}} ID SF DM{dm}", Nominal=round(nom, 4), Min=0.5, Max=1.5,
                                    Samples=",".join(ts), Regions=",".join(hipt) if lowpt else None,
                                    Category="Tau ID (fitted)", Constant=True if fix_tauid else None))
        if lowpt:
            blocks.append(tc.normfactor(f"TauIDSF_DM{dm}_lowpt", Title=f"#tau_{{h}} ID SF DM{dm}, p_{{T}} < {config.LTAU_PT_SPLIT:.0f} GeV",
                                        Nominal=round(nom, 4), Min=0.5, Max=1.5, Samples=",".join(ts), Regions=",".join(lowpt),
                                        Category="Tau ID (fitted)", Constant=True if fix_tauid else None))
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
            blk = tc.overall_syst(name, s["up"], s["down"], ",".join(smp), s["category"], title=s["title"], regions=",".join(regs), subcategory=s["subcategory"])
        else:
            blk = tc.histo_syst(name, ",".join(smp), s["category"], title=s["title"], regions=",".join(regs), smoothing=s["smoothing"], subcategory=s["subcategory"])
        if name in scale_systs:
            blk[2].extend([("ScaleUp", float(scale_systs[name])), ("ScaleDown", float(scale_systs[name]))])
        blocks.append(blk)
    if pt:
        blocks.append(tc.overall_syst(pt["name"], pt["rel"], -pt["rel"], ",".join(sig), "Tau ID pT dependence",
                                      title="#tau_{h} ID SF p_{T} dependence"))
    header = (f"Z -> tautau v4 (four channels), generated by z-tautau/scripts/step5_fit.py (do not edit by hand); tau_h working point {config.TAU_WP}\n"
              f"tau_h ID scale factors free per decay mode (nominal = TauPOG), tau energy scale prior {100 * meta['tes_prior']:.0f}%; TRExFitter v1.8.0")
    return tc.render(blocks, header=header)


MINOS_RX = re.compile(r"minos status\s+(-?\d+)")


def minos_status(log_path):
    """Last 'minos status' printed by TRExFitter's FittingTool (0 = both errors found, -1 = MINOS failed)."""
    path = Path(log_path)
    if not path.exists():
        return None
    hits = MINOS_RX.findall(path.read_text())
    return int(hits[-1]) if hits else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-ranking", action="store_true")
    ap.add_argument("--skip-asimov", action="store_true", help="no expected (Asimov) fit: for the cross-check jobs")
    ap.add_argument("--config-only", action="store_true")
    ap.add_argument("--channels", nargs="*", default=config.CHANNELS)
    ap.add_argument("--job", default=config.JOB_V4)
    ap.add_argument("--fix-tauid", action="store_true", help="tau_h ID scale factors fixed to the TauPOG values (single-channel cross-checks)")
    ap.add_argument("--region-set", default="nominal", choices=("nominal", "ptsplit"),
                    help="'ptsplit': the l tau_h regions split at pT(tau_h) = 40 GeV with their own tau_h ID scale factors below it")
    ap.add_argument("--no-tauid-pt", action="store_true",
                    help=f"leave out {TAUID_PT_NP} (the job ztautau_flatsf, which with ztautau_ptsplit sets its size). It is "
                         "carried by every fit of the measurement's model: four channels, one scale factor per decay mode")
    ap.add_argument("--scale-syst", action="append", default=[], metavar="NAME=FACTOR",
                    help="scale one systematic variation (cross-check, e.g. EmuTrigger=2.0)")
    ap.add_argument("--fit-strategy", type=int, default=None, help="Minuit2 strategy of every fit (default: TRExFitter's)")
    ap.add_argument("--num-cpu", type=int, default=6, help="NumCPU of the Fit block (RooFit parallelisation)")
    ap.add_argument("--summarise-only", action="store_true",
                    help="do not fit: only re-read the existing TRExFitter outputs into the result json "
                         "(used after the ranking has been produced in parallel, see condor/)")
    args = ap.parse_args()
    job = args.job
    scale_systs = dict(kv.split("=") for kv in args.scale_syst)
    fitdir = config.FIT_DIR_V4
    # a single channel reads its own exported file (step 4b) when it exists, so its config, fit inputs and
    # workspace are self-contained for the MultiFit combination; otherwise the combined file
    histo_file = config.JOB_V4
    if len(args.channels) == 1 and (fitdir / "fitinputs" / f"{config.JOB_V4}_{args.channels[0]}.root").exists():
        histo_file = f"{config.JOB_V4}_{args.channels[0]}"
    meta = json.loads(Path(f"{fitdir}/fitinputs/{histo_file}.root.meta.json").read_text())
    cfg = fitdir / f"{job}.config"
    carries_pt = (not args.no_tauid_pt and not args.fix_tauid and args.region_set == "nominal"
                  and set(args.channels) == set(config.CHANNELS))
    args.pt = pt_model(fitdir) if carries_pt else None
    if args.pt:
        print(f"{TAUID_PT_NP}: +-{100 * args.pt['rel']:.2f} % on the signal "
              f"(mu_Z {args.pt['mu_flatsf']:.4f} flat, {args.pt['mu_ptsplit']:.4f} pT-split)")

    def write(path, strategy=None):
        txt = build_config(job, meta, args.channels, args.fix_tauid, args.region_set, scale_systs, strategy,
                           args.num_cpu, empty, args.pt)
        txt = txt.replace(f'HistoFile: "{job}"', f'HistoFile: "{histo_file}"')
        txt = re.sub(r'Expression: "([^"]+)"', r"Expression: \1", txt)   # TRExFitter splits the value on ':' before unquoting
        Path(path).write_text(txt)

    empty = empty_bins(fitdir / "fitinputs" / f"{histo_file}.root",
                       select_regions(meta, args.channels, args.region_set))
    n_empty = sum(len(v) for v in empty.values())
    if n_empty:
        print(f"dropping {n_empty} empty bin(s): " + ", ".join(f"{r}{v}" for r, v in empty.items() if v))

    if args.summarise_only:          # the fits already ran (condor/): only re-read their outputs
        summarise(args, job, fitdir, meta, scale_systs)
        return
    write(cfg, args.fit_strategy)
    # The expected (Asimov) fit is run from its own copy of the config with Minuit2 strategy 3: with the
    # default strategy MINOS failed on mu_Z on the Asimov data set ("Invalid lower error", Hessian forced
    # pos-def) and the HESSE error of that non-positive-definite matrix was quoted as the expected
    # uncertainty (review/REVIEW_v4.md finding 1). Strategy 2 cured that until TauIDpT_tautau, which is degenerate
    # with mu_Z: with it strategy 2 ends on a Hessian forced pos-def (MINOS status 1), strategy 3 is clean
    # (the combination found the same for this likelihood). If MINOS still fails, nothing is quoted.
    cfg_asimov = fitdir / f"{job}_asimov.config"
    write(cfg_asimov, 3)
    print(f"config -> {cfg}  ({len(meta['samples'])} templates, {len(meta['systs'])} systematics, "
          f"{len(select_regions(meta, args.channels, args.region_set))} regions, region set {args.region_set})")
    if args.config_only:
        return
    trex_environment()
    logs = f"results/{job}/logs"
    # Order matters: the stat-only and Asimov fits overwrite the unsuffixed correlation matrix and the
    # NP plots, so the observed fit and everything drawn from it (d, p, i, r) run *last* (finding 1).
    run_trex.run(cfg.name, "h", cwd=fitdir, log=f"{logs}/h.log")
    run_trex.run(cfg.name, "w", cwd=fitdir, log=f"{logs}/w.log")
    run_trex.run(cfg.name, "wf", options="StatOnly=TRUE:Suffix=_statOnly", cwd=fitdir, log=f"{logs}/statonly.log")
    if not args.skip_asimov:
        run_trex.run(cfg_asimov.name, "wf", options="FitBlind=TRUE:Suffix=_asimov", cwd=fitdir, log=f"{logs}/asimov.log")
    for actions in ("f", "d", "p", "i") + (() if args.skip_ranking else ("r",)):
        run_trex.run(cfg.name, actions, cwd=fitdir, log=f"{logs}/{actions}.log")
    # With DropBins TRExFitter writes two workspaces: <job>_combined_<job>_model.root is the fitted model (bins
    # dropped) and <job>_allBinsFitRegions_... the one it draws from. Never copy the second over the first (done
    # here until 17 Sep 2026): the ranking, the MultiFit and any combination of workspaces then fit the fake
    # sideband of tautau_SR0 and the empty bin as well (mu_Z 1.0172 instead of 1.0187).
    summarise(args, job, fitdir, meta, scale_systs)


def summarise(args, job, fitdir, meta, scale_systs):
    """Read everything TRExFitter wrote for `job` into fit/results/<job>_fit_result.json."""
    logs = f"results/{job}/logs"
    res = run_trex.summarise(fitdir / "results" / job, job, poi="mu_Z")
    res["gof"] = run_trex.parse_gof(fitdir / logs / "f.log")
    res["minos_status"] = minos_status(fitdir / logs / "f.log")
    asimov = fitdir / "results" / job / "Fits" / f"{job}_asimov.txt"
    st = minos_status(fitdir / logs / "asimov.log")
    res["poi_asimov"] = None
    res["asimov_minos_status"] = st
    if asimov.exists() and not args.skip_asimov:
        a = run_trex.parse_fit_txt(asimov)["nps"].get("mu_Z")
        if a and st == 0:
            res["poi_asimov"] = {"value": a[0], "err_up": a[1], "err_down": a[2], "minos_status": st}
        elif a:
            # MINOS failed: the numbers in the file are the HESSE errors of a Hessian that was forced
            # positive definite. They are not an expected uncertainty and must not be quoted.
            res["asimov_hesse_not_quoted"] = {"value": a[0], "err_up": a[1], "err_down": a[2], "minos_status": st}
            print(f"  WARNING: Asimov MINOS failed (status {st}); no expected uncertainty is quoted")
    fit_txt = fitdir / "results" / job / "Fits" / f"{job}.txt"
    nps = run_trex.parse_fit_txt(fit_txt)["nps"] if fit_txt.exists() else {}
    sf_names = [f"TauIDSF_DM{dm}" for dm in config.TAU_DMS] + [f"TauIDSF_DM{dm}_lowpt" for dm in config.TAU_DMS]
    res["tau_id_sf"] = {n.replace("TauIDSF_", ""): {"value": nps[n][0], "err_up": nps[n][1], "err_down": nps[n][2],
                                                    "pog": meta["tau_id_sf_pog"][n.replace("TauIDSF_DM", "").replace("_lowpt", "")]}
                        for n in sf_names if n in nps}
    res["tau_es"] = {f"DM{dm}": {"pull": nps[f"TauES_DM{dm}"][0], "constraint": nps[f"TauES_DM{dm}"][1], "prior_pct": 100 * meta["tes_prior"]}
                     for dm in config.TAU_DMS if f"TauES_DM{dm}" in nps}
    res["mu_ttbar"] = ({"value": nps["mu_ttbar"][0], "err_up": nps["mu_ttbar"][1], "err_down": nps["mu_ttbar"][2]}
                       if "mu_ttbar" in nps else None)
    sig = meta["signal_prediction"]
    mu, up, dn = res["poi_value"], res["poi_err_up"], res["poi_err_down"]
    stat = res.get("poi_stat_only", {}).get("err")
    s60 = sig["sigma_tautau_60_120_pb"]
    res["prediction"] = sig
    res["sigma_60_120_pb"] = {"value": mu * s60, "up": up * s60, "down": dn * s60, "stat": stat * s60 if stat else None,
                              "note": "sigma(pp -> Z/gamma* -> tautau, 60 < m_LHE < 120 GeV) = mu_Z x aMC@NLO prediction normalised to 6077.22 pb (m > 50)"}
    # The grouped impacts share nuisance parameters between categories (mu_ttbar <-> EmuTrigger <-> mu_Z),
    # so their quadrature sum over-shoots the MINOS total. Record the over-shoot instead of hiding it
    # (review/REVIEW_v4.md finding 5): the total to quote is always the MINOS one.
    groups = {k: v for k, v in res.get("grouped_impact", {}).items() if k not in ("FullSyst", "Total")}
    quad = float(sum(v ** 2 for v in groups.values()) + (stat or 0.0) ** 2) ** 0.5
    res["grouped_impact_quadrature_sum"] = quad
    res["grouped_impact_scale"] = (0.5 * (up + dn)) / quad if quad > 0 else None
    res["lumi_pb"] = config.LUMI_PB
    res["channels"] = args.channels
    res["tau_id_fixed"] = args.fix_tauid
    res["region_set"] = args.region_set
    res["scaled_systematics"] = scale_systs
    res["tauid_pt_model"] = args.pt
    res["channels"] = args.channels
    out = fitdir / "results" / f"{job}_fit_result.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    print(f"mu_Z = {mu:.4f} +{up:.4f} -{dn:.4f}  (stat {stat});  sigma(60-120) = {res['sigma_60_120_pb']['value']:.1f} pb")
    if res["poi_asimov"]:
        print(f"  expected (Asimov, MINOS): +{res['poi_asimov']['err_up']:.4f} -{res['poi_asimov']['err_down']:.4f}")
    print(f"  grouped impacts in quadrature {quad:.4f} vs MINOS total {0.5 * (up + dn):.4f} (x{res['grouped_impact_scale'] or 0:.2f})")
    for dm, v in res["tau_id_sf"].items():
        print(f"  TauIDSF {dm}: {v['value']:.3f} +{v['err_up']:.3f} -{v['err_down']:.3f}  (POG {v['pog'][0]:.3f} +- {v['pog'][1]:.3f})")
    for dm, v in res["tau_es"].items():
        print(f"  TauES {dm}: pull {v['pull']:+.2f}, constraint {v['constraint']:.2f} x 3% = {100 * meta['tes_prior'] * v['constraint']:.2f}%")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
