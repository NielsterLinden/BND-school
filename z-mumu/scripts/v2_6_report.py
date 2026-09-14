#!/usr/bin/env python
"""v2 step 6 -- data/MC plots, results_v2.json, RESULTS_v2.md and the summary plot.

    python scripts/v2_6_report.py
"""

from __future__ import annotations

import json
import pickle
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import numpy as np

from zmumu import config, hists, plotting

OUT = config.OUTPUT_DIR / "v2"
FIT = config.REPO_DIR / "fit" / "results" / "zmumu_fit_result.json"
V1 = {"sigma_fid": 773.2, "sigma_fid_err": 11.9, "revised": 776.9, "revised_err": 14.8, "reviewer": 797.2,
      "reviewer_err": 797.2 * np.hypot(0.013, 0.012), "pred_nlo": 799.6, "pred_lo": 825.4}


def main():
    hall = pickle.load(open(OUT / "histograms.pkl", "rb"))
    fakes_res = json.load(open(OUT / "fakes.json")) if (OUT / "fakes.json").exists() else None
    made = plotting.all_plots(hall, fakes_res)
    print(f"[report] {len(made)} data/MC plots")
    fit = json.load(open(FIT)) if FIT.exists() else None
    tnp = json.load(open(OUT / "tnp" / "tnp_result.json")) if (OUT / "tnp" / "tnp_result.json").exists() else None
    mom = json.load(open(OUT / "momentum.json")) if (OUT / "momentum.json").exists() else None
    yields = {}
    for k, v in hall.items():
        smp, region, var, variation = k.split("|")[:4]
        if var == "mass_fit" and variation == "nominal" and not k.endswith("|w2"):
            yields.setdefault(region, {})[smp] = float(np.sum(v))
    if fakes_res:
        yields["SR"]["Fakes"] = fakes_res["yields"]["sr_fakes"]
    result = {"generated": str(date.today()), "lumi_pb": config.LUMI_PB_NORMTAG, "yields": yields, "fit": fit,
              "fakes": fakes_res["yields"] if fakes_res else None, "momentum": mom,
              "tnp_meta": tnp["meta"] if tnp else None, "v1_reference": V1}
    with open(OUT / "results_v2.json", "w") as fh:
        json.dump(result, fh, indent=1, default=float)
    lines = [f"# Z -> mu+ mu- cross section, v2 (MC-based, TRExFitter fit)", "",
             f"Generated {date.today()} from CMS 2016 Open Data (SingleMuon, Run2016G+H, {config.LUMI_PB_NORMTAG/1000:.3f} fb^-1, normtag).", ""]
    if fit:
        lines += ["## Result", "",
                  f"**sigma_fid(pp -> Z/gamma* -> mu mu; dressed, pT > 26/20 GeV, |eta| < 2.4, 60 < m < 120 GeV) = "
                  f"{fit['sigma_fid_pb']:.1f} +- {fit['sigma_fid_stat_pb']:.1f} (stat) +- {fit['sigma_fid_syst_pb']:.1f} (syst) +- {fit['sigma_fid_lumi_pb']:.1f} (lumi) pb**", "",
                  f"mu_Z = {fit['mu']:.4f} +{fit['mu_err_up']:.4f} -{fit['mu_err_down']:.4f}; goodness of fit p = {fit['gof']['gof_probability']}", "",
                  f"| quantity | value |", "|---|---:|",
                  f"| sigma(Z/gamma* -> mu mu, 60 < m < 120 GeV), A = {fit['A_60_120']:.4f} | {fit['sigma_60_120_pb']:.0f} +- {fit['sigma_60_120_tot_pb']:.0f} pb |",
                  f"| sigma(Z/gamma* -> mu mu, m > 50 GeV), A = {fit['A_m50']:.4f} | {fit['sigma_m50_pb']:.0f} +- {fit['sigma_m50_tot_pb']:.0f} pb |",
                  f"| C factor (reco/fiducial, all corrections) | {fit['meta']['C_factor']:.4f} |",
                  f"| NLO prediction sigma_fid (6077.22 pb x fiducial fraction) | {fit['sigma_fid_pred_pb']:.1f} pb |",
                  f"| counting cross-check (N_obs - N_bkg)/(C L) | {fit['meta']['counting']['sigma_fid_pb']:.1f} pb |",
                  f"| v1 (data-only counting) / v1 revised / reviewer | {V1['sigma_fid']} / {V1['revised']} / {V1['reviewer']} pb |", "",
                  "## Uncertainty breakdown (impact on mu_Z, from the grouped-impact fit)", "", "| group | relative |", "|---|---:|"]
        for k, v in sorted(fit["grouped_impacts_mu"].items(), key=lambda kv: -kv[1]):
            lines.append(f"| {k} | {100*v:.3f}% |")
        lines += [f"| statistical (stat-only fit) | {100*fit['mu_stat_only_fit']:.3f}% |", ""]
        lines += ["## Leading nuisance parameters (ranking)", "", "| NP | pull | constraint | +impact | -impact |", "|---|---:|---:|---:|---:|"]
        for r in fit["ranking"][:12]:
            lines.append(f"| {r['name']} | {r['pull']:+.2f} | {0.5*(r['err_up']+r['err_down']):.2f} | {100*r['dpoi_up_post']:+.3f}% | {100*r['dpoi_down_post']:+.3f}% |")
        lines.append("")
    lines += ["## Yields (mass_fit templates, prompt-prompt MC + data-driven fakes)", "", "| region | sample | events |", "|---|---|---:|"]
    for region, d in yields.items():
        for smp, n in sorted(d.items(), key=lambda kv: -kv[1]):
            lines.append(f"| {region} | {smp} | {n:,.1f} |")
    if fakes_res:
        y = fakes_res["yields"]; c = fakes_res["closure"]
        lines += ["", "## Fake factor", "", f"SR non-prompt estimate {y['sr_fakes']:.0f} +- {y['sr_fakes_stat']:.0f} (stat); method uncertainty {100*y['method_rel_unc']:.0f}% (largest: {y['method_worst']}).",
                  f"Same-sign closure: predicted {c['predicted_ss']:.0f} +- {c['predicted_ss_err']:.0f}, observed {c['observed_ss']:.0f} +- {c['observed_ss_err']:.0f}."]
    if tnp:
        m = tnp["meta"]
        lines += ["", "## Tag-and-probe", "", f"Pairs: {m['n_pairs_data']:,} (data). Trigger plateau efficiency per muon: data {m['trig_plateau_data']:.4f}, MC {m['trig_plateau_mc']:.4f}.",
                  f"Isolation SF vs pileup half-spread: {m['iso_sf_npv_spread']:.4f}.", "",
                  "| efficiency | <SF> | <stat> | <syst> |", "|---|---:|---:|---:|"]
        for eff in ("id", "iso", "antiiso"):
            lines.append(f"| {eff} | {np.mean(tnp['sf'][eff]):.4f} | {np.mean(tnp['sf_stat'][eff]):.4f} | {np.mean(tnp['sf_syst'][eff]):.4f} |")
    if mom:
        lines += ["", "## Momentum calibration (MC pT scale kappa and extra smearing per |eta| bin)", ""]
        lines.append("| " + " | ".join(f"{a:.1f}-{b:.1f}" for a, b in zip(mom["eta_edges"][:-1], mom["eta_edges"][1:])) + " |")
        lines.append("|" + "---:|" * (len(mom["kappa"])))
        lines.append("| " + " | ".join(f"{k:.5f} +- {e:.5f}, {100*s:.2f}%" for k, e, s in zip(mom["kappa"], mom["kappa_err"], mom["smear"])) + " |")
    (OUT / "RESULTS_v2.md").write_text("\n".join(lines) + "\n")
    print(f"[report] wrote {OUT / 'RESULTS_v2.md'}")
    if fit:
        import matplotlib.pyplot as plt
        rows = [("DY aMC@NLO prediction", V1["pred_nlo"], V1["pred_nlo"] * fit["acceptance"].get("A_rel_unc", 0.006), "#1f77b4"),
                ("DY madgraph LO prediction", V1["pred_lo"], V1["pred_lo"] * 0.01, "#aec7e8"),
                ("v1: data-only counting (773.2)", V1["sigma_fid"], V1["sigma_fid_err"], "grey"),
                ("v1 revised (review corrections)", V1["revised"], V1["revised_err"], "grey"),
                ("reviewer's analysis (797.2)", V1["reviewer"], V1["reviewer_err"], "dimgrey"),
                ("v2: MC-based, TRExFitter fit", fit["sigma_fid_pb"], fit["sigma_fid_tot_pb"], "black")]
        fig, ax = plt.subplots(figsize=(10, 6))
        for y, (label, value, err, colour) in enumerate(rows[::-1]):
            ax.errorbar(value, y, xerr=err, fmt="o", color=colour, capsize=4, markersize=7)
            ax.text(value, y + 0.22, f"{value:.1f} $\\pm$ {err:.1f} pb", ha="center", fontsize=10)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=11)
        ax.set_ylim(-0.6, len(rows) - 0.3); ax.set_xlabel(r"$\sigma_{fid}(pp \to Z/\gamma^* \to \mu\mu)$ [pb]")
        hists._decorate(ax, lumi_fb=config.LUMI_PB_NORMTAG / 1000)
        hists.save_fig(fig, "summary_v2.png")


if __name__ == "__main__":
    main()
