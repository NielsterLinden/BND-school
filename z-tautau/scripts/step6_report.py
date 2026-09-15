#!/usr/bin/env python
"""Step 6 -- collect everything into output/results.json, output/RESULTS.md and summary plots.

    python scripts/step6_report.py

Reads the fit results of both fake-factor variants (fit/results/ztautau[_nosub]_fit_result.json), the
yields, fake-factor and BDT summaries of steps 3-4 and the ntuple cutflows. Copies the key TRExFitter plots to
output/plots/fit_*.png so they are committed next to the analysis plots.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import analysis, config, samples  # noqa: E402

NOMINAL = "mcsub" if config.FF_SUBTRACT_MC else "nosub"
OTHER = "nosub" if NOMINAL == "mcsub" else "mcsub"
LABEL = {"nosub": "FF without MC subtraction", "mcsub": "FF with genuine-tau subtraction"}


def job(variant):
    return config.JOB if variant == NOMINAL else f"{config.JOB}_{variant}"


def load_json(path):
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else None


def summarise_fit(res):
    if res is None:
        return None
    mu, up, dn = res["poi_value"], res["poi_err_up"], res["poi_err_down"]
    stat = res["poi_stat_only"]["err"]
    tot = 0.5 * (up + dn)
    syst = float(np.sqrt(max(tot ** 2 - stat ** 2, 0.0)))
    sig = res["prediction"]
    out = {"mu": mu, "mu_err_up": up, "mu_err_down": dn, "mu_stat": stat, "mu_syst": syst,
           "mu_expected_asimov": res.get("poi_asimov"), "gof_probability": res.get("gof", {}).get("gof_probability"),
           "grouped_impact": res.get("grouped_impact", {})}
    for name, pred in (("sigma_fid_pb", sig["sigma_fid_pb"]), ("sigma_60_120_pb", sig["sigma_tautau_60_120_pb"])):
        out[name] = {"value": mu * pred, "stat": stat * pred, "syst": syst * pred, "err_up": up * pred,
                     "err_down": dn * pred, "prediction": pred}
    out["sigma_60_120_pb"]["acceptance"] = mu * sig["sigma_tautau_60_120_pb"] * float(
        np.sqrt(sum(v ** 2 for v in sig["A_unc"].values()) + sig["A_mc_stat"] ** 2))
    ranking = sorted(res.get("ranking", []), key=lambda r: -max(abs(r["dpoi_up_post"]), abs(r["dpoi_down_post"])))
    out["ranking"] = [{"name": r["name"], "pull": r["pull"], "constraint": 0.5 * (abs(r["err_up"]) + abs(r["err_down"])),
                       "impact_up": r["dpoi_up_post"], "impact_down": r["dpoi_down_post"]} for r in ranking[:15]]
    out["pulls"] = {k: v for k, v in res["nps"].items() if not k.startswith("gamma")}
    return out


def copy_trex_plots(variant):
    base = config.FIT_DIR / "results" / job(variant)
    tag = "" if variant == NOMINAL else f"_{variant}"
    pairs = [("NuisPar.png", f"fit_pulls{tag}.png"), ("CorrMatrix.png", f"fit_corrmatrix{tag}.png"),
             ("Rankings/RankingSysts_mu_Z_systs.png", f"fit_ranking{tag}.png"),
             ("Rankings/RankingSysts_mu_Z_Breakdown_syst.png", f"fit_ranking_breakdown{tag}.png"),
             ("Gammas.png", f"fit_gammas{tag}.png"), ("NormFactors.png", f"fit_normfactors{tag}.png"),
             ("Plots/Summary_postFit.png", f"fit_summary_postfit{tag}.png"), ("Plots/Summary.png", f"fit_summary_prefit{tag}.png")]
    for region in config.REGIONS:
        short = region.replace("tautau_", "")
        pairs += [(f"Plots/{region}_postFit.png", f"fit_{short}_postfit{tag}.png"), (f"Plots/{region}.png", f"fit_{short}_prefit{tag}.png")]
    for src, dst in pairs:
        p = base / src
        if not p.exists():
            cand = sorted(base.glob(f"**/{Path(src).name}"))
            p = cand[0] if cand else None
        if p and p.exists():
            shutil.copy(p, config.PLOT_DIR / dst)


def main():
    res = {v: load_json(config.FIT_DIR / "results" / f"{job(v)}_fit_result.json") for v in (NOMINAL, OTHER)}
    fits = {v: summarise_fit(r) for v, r in res.items()}
    yields = {v: load_json(config.DATA_DIR / f"yields{'' if v == NOMINAL else '_' + v}.json") for v in (NOMINAL, OTHER)}
    ff = load_json(config.DATA_DIR / "fakefactors.json")
    bdt = load_json(config.DATA_DIR / "bdt.json")
    sig = analysis.signal_prediction(samples.DY_INCLUSIVE)
    cutflows = {k: analysis.load(k)[1]["cutflow"] for k in samples.DATA_KEYS + [samples.DY_INCLUSIVE]}

    y = yields[NOMINAL]["yields"]
    n_obs = y["Data"]["value"]
    n_sig = y[samples.SIGNAL]["value"]
    n_bkg = sum(v["value"] for k, v in y.items() if k not in ("Data", "Total", samples.SIGNAL))
    acc_eff = n_sig / (sig["sigma_tautau_60_120_pb"] * config.LUMI_PB)
    counting = (n_obs - n_bkg) / (acc_eff * config.LUMI_PB)
    results = {
        "channel": "Z -> tau_h tau_h", "dataset": "CMS Open Data 2016 Tau Run2016G+H (UL NanoAODv9)",
        "lumi_pb": config.LUMI_PB, "lumi_rel_unc": config.LUMI_REL_UNC, "fit_variable": config.FIT_VARIABLE,
        "regions": config.REGIONS, "category_edges": config.BDT_CATEGORY_EDGES,
        "nominal_variant": NOMINAL, "variants": {v: LABEL[v] for v in (NOMINAL, OTHER)},
        "prediction": sig, "fit": fits,
        "for_combination": {"n_obs": n_obs, "n_bkg_prefit": n_bkg, "n_sig_prefit": n_sig, "acc_eff": acc_eff,
                            "acc_eff_definition": "N_SR(Z->tautau fiducial, prefit, all corrections, all categories) / (sigma_pred(Z/gamma*->tautau, 60<m_LHE<120) x L)",
                            "sigma_counting_pb": counting, "A_fid": sig["A"], "C": n_sig / (sig["sigma_fid_pb"] * config.LUMI_PB)},
        "yields_prefit": {v: yields[v]["yields"] if yields[v] else None for v in (NOMINAL, OTHER)},
        "yields_prefit_per_region": yields[NOMINAL]["regions"],
        "fake_factors": {v: {"C_OS_SS_per_njet": ff[v]["osss"]["C"], "C_OS_SS_stat": ff[v]["osss"]["stat"],
                             "C_OS_SS_inclusive": ff[v]["osss"]["inclusive"]["C"],
                             "closure_corrections": ff[v]["ff"].get("closure")} for v in ("mcsub", "nosub")},
        "closure_nps": yields[NOMINAL]["meta"].get("closure_nps"),
        "C_osss_per_category": yields[NOMINAL]["meta"].get("C_osss"),
        "sigmodel": yields[NOMINAL]["meta"].get("sigmodel"), "dy_samples": yields[NOMINAL]["meta"].get("dy_samples"),
        "bdt": {k: bdt[k] for k in ("training", "category_edges", "category_yields_prefit", "stat_only_sensitivity")} if bdt else None,
        "genuine_tau_contamination": ff["contamination"], "cutflow": cutflows,
        "stat_only_sensitivity_fixed_bkg": yields[NOMINAL].get("stat_only_sensitivity"),
    }
    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (config.OUTPUT_DIR / "results.json").write_text(json.dumps(results, indent=1, default=float))
    for v in (NOMINAL, OTHER):
        if res[v]:
            copy_trex_plots(v)

    # ------------------------------------------------------------------ RESULTS.md
    f = fits[NOMINAL]
    lines = ["# Z -> tau_h tau_h: results", "",
             f"CMS Open Data 2016 (Run2016G+H, Tau dataset), L = {config.LUMI_PB:.1f} pb^-1 (+-1.2%). "
             f"Binned profile-likelihood fit of m_tt (MET-likelihood di-tau mass) in {len(config.REGIONS)} BDT categories of the "
             f"signal region with TRExFitter v1.8.0. Generated by scripts/step6_report.py.", ""]
    for v in (NOMINAL, OTHER):
        fv = fits[v]
        if not fv:
            continue
        s60, sf = fv["sigma_60_120_pb"], fv["sigma_fid_pb"]
        lines += [f"## {LABEL[v]}{' (nominal)' if v == NOMINAL else ' (cross-check)'}", "",
                  "| quantity | value |", "|---|---|",
                  f"| mu_Z | {fv['mu']:.3f} +{fv['mu_err_up']:.3f} -{fv['mu_err_down']:.3f} (stat {fv['mu_stat']:.3f}, syst {fv['mu_syst']:.3f}) |",
                  f"| expected (Asimov) | +{fv['mu_expected_asimov']['err_up']:.3f} -{fv['mu_expected_asimov']['err_down']:.3f} |" if fv["mu_expected_asimov"] else "| expected | n/a |",
                  f"| sigma_fid (tau_h tau_h, vis pT > 40, abs(eta) < 2.1, 60 < m < 120) | {sf['value']:.3f} +- {sf['stat']:.3f} (stat) +- {sf['syst']:.3f} (syst) pb  (pred. {sf['prediction']:.3f}) |",
                  f"| sigma(pp -> Z/gamma* -> tautau, 60 < m < 120) | {s60['value']:.0f} +- {s60['stat']:.0f} (stat) +- {s60['syst']:.0f} (syst) +- {s60['acceptance']:.0f} (acc) pb  (pred. {s60['prediction']:.0f}) |",
                  f"| goodness of fit (saturated) | p = {fv['gof_probability']} |", "",
                  "Grouped impacts on mu_Z:", "", "| group | impact |", "|---|---|"]
        for k, val in sorted(fv["grouped_impact"].items(), key=lambda kv: -kv[1]):
            lines.append(f"| {k} | {val:.4f} |")
        lines.append("")
    lines += ["## Prefit yields per BDT category (nominal)", "",
              "| sample | " + " | ".join(f"{r} ({lab})" for r, lab in zip(config.REGIONS, config.REGION_LABELS)) + " | all |",
              "|---|" + "---:|" * (len(config.REGIONS) + 1)]
    for k, v in y.items():
        if k == "Total":
            continue
        per = [results["yields_prefit_per_region"][r][k]["value"] for r in config.REGIONS]
        lines.append(f"| {k} | " + " | ".join(f"{p:.0f}" for p in per) + f" | {v['value']:.0f}{' +- %.0f' % v['stat'] if 'stat' in v else ''} |")
    lines.append(f"| Total prediction | " + " | ".join(f"{sum(results['yields_prefit_per_region'][r][s]['value'] for s in y if s not in ('Data', 'Total')):.0f}" for r in config.REGIONS) + f" | {y['Total']['value']:.0f} |")
    if bdt:
        t = bdt["training"]
        lines += ["", f"BDT: held-out AUC {np.mean(t['auc_test']):.3f} (folds {', '.join(f'{a:.3f}' for a in t['auc_test'])}); "
                  f"stat-only sensitivity (fakes fixed) inclusive {100 * bdt['stat_only_sensitivity']['inclusive']:.2f}% -> "
                  f"categories {100 * bdt['stat_only_sensitivity']['categories']:.2f}%"]
    lines += ["", "## Inputs for the combination", "", "| variable | value |", "|---|---|"]
    for k, v in results["for_combination"].items():
        lines.append(f"| {k} | {v:.6g} |" if isinstance(v, float) else f"| {k} | {v} |")
    lines += ["", f"Signal prediction: sigma(Z/gamma*->tautau, m>50) = {sig['sigma_tautau_pb']:.1f} pb, 60-120: "
              f"{sig['sigma_tautau_60_120_pb']:.1f} pb, fiducial {sig['sigma_fid_pb']:.3f} pb, A = {sig['A']:.5f} "
              f"(scale {100 * sig['A_unc']['A_scale']:.1f}%, PDF {100 * sig['A_unc']['A_pdf']:.1f}%, "
              f"alpha_s {100 * sig['A_unc']['A_alphas']:.1f}%, ISR {100 * sig['A_unc']['A_isr']:.1f}%, "
              f"FSR {100 * sig['A_unc']['A_fsr']:.1f}%, MC stat {100 * sig['A_mc_stat']:.1f}%)", ""]
    (config.OUTPUT_DIR / "RESULTS.md").write_text("\n".join(lines) + "\n")

    # ------------------------------------------------------------------ plots
    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    rows = []
    for v in (NOMINAL, OTHER):
        if fits[v]:
            s = fits[v]["sigma_60_120_pb"]
            rows.append((LABEL[v], s["value"], s["stat"], np.hypot(0.5 * (s["err_up"] + s["err_down"]), s["acceptance"]),
                         s["err_down"], s["err_up"]))
    pred = sig["sigma_tautau_60_120_pb"]
    ax.axvspan(pred * (1 - 0.04), pred * (1 + 0.04), color="#ffcc66", alpha=0.5)
    ax.axvline(pred, color="#c08a00")
    for i, (lab, val, stat, tot, dn, up) in enumerate(rows):
        yv = len(rows) - i
        ax.errorbar([val], [yv], xerr=[[dn], [up]], fmt="o", color="black", capsize=4, lw=1.5)
        ax.errorbar([val], [yv], xerr=[[stat], [stat]], fmt="none", color="tab:red", capsize=0, lw=4)
        ax.text(1230, yv + 0.33, f"{lab}", fontsize=10, ha="left")
        ax.text(3270, yv + 0.33, f"{val:.0f} +{up:.0f} -{dn:.0f} pb (stat. {stat:.0f})", fontsize=10, ha="right")
    ax.text(pred + 15, 0.45, f"prediction {pred:.0f} pb", fontsize=9, color="#8a6200")
    ax.set_ylim(0.3, len(rows) + 0.9)
    ax.set_xlim(1200, 3300)
    ax.set_yticks([])
    ax.set_xlabel(r"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau)$, $60<m<120$ GeV [pb]")
    ax.set_title(r"Z$\rightarrow\tau_h\tau_h$, CMS Open Data 2016 G+H, 16.4 fb$^{-1}$ (red: stat.)", fontsize=11)
    fig.savefig(config.PLOT_DIR / "step6_summary.png", bbox_inches="tight", dpi=130)
    fig.savefig(config.PLOT_DIR / "step6_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    gi = {k: v for k, v in f["grouped_impact"].items() if k not in ("FullSyst",)}
    gi["Data statistics"] = f["mu_stat"]
    items = sorted(gi.items(), key=lambda kv: kv[1])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh([k for k, _ in items], [100 * v for _, v in items], color="#4a90d9")
    for i, (_, val) in enumerate(items):
        ax.text(100 * val + 0.2, i, f"{100 * val:.1f}%", va="center", fontsize=9)
    ax.set_xlabel(r"impact on $\mu_Z$ [%]")
    ax.set_title(f"Grouped uncertainties ({LABEL[NOMINAL]})", fontsize=11)
    fig.savefig(config.PLOT_DIR / "step6_impacts.png", bbox_inches="tight", dpi=130)
    fig.savefig(config.PLOT_DIR / "step6_impacts.pdf", bbox_inches="tight")
    plt.close(fig)
    print((config.OUTPUT_DIR / "RESULTS.md").read_text())


if __name__ == "__main__":
    main()
