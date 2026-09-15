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
config.PLOT_DIR = OUT / "plots"          # v2 plots live next to the v2 results
FIT = config.REPO_DIR / "fit" / "results" / "zmumu_fit_result.json"
STAB = config.REPO_DIR / "fit" / "results" / "stability.json"
NOT_IN_FIT = {"DYmumu_powheg": "raw generator weights; SigModel template only", "WJets": "prompt-prompt part is a negative-weight fluctuation; e-mu regions only",
              "DYother": "LHE flavour not e/mu/tau"}
V1 = {"sigma_fid": 773.2, "sigma_fid_err": 11.9, "revised": 776.9, "revised_err": 14.8, "reviewer": 797.2,
      "reviewer_err": 797.2 * np.hypot(0.013, 0.012), "pred_nlo": 799.6, "pred_lo": 825.4}


def lineshape_plot(hall, gens, fakes_res):
    """powheg vs aMC@NLO lineshape: generator level (LHE mass, all generated mu mu events) and
    reconstructed level inside the powheg generator window, next to the pre-fit data/MC ratio.
    Documents where the SigModel template comes from (REVIEW.md F3/F4)."""
    import matplotlib.pyplot as plt
    from zmumu import histograms as H
    hn = np.array(gens["DY_NLO"]["h_lhe_mll"]); hp = np.array(gens["DY_powheg"]["h_lhe_mll"])   # 1 GeV bins, 0-200
    win = slice(60, 120)
    rb = lambda a, k: a.reshape(-1, k).sum(axis=1)
    n_lhe, p_lhe = rb(hn[win], 5) / hn[win].sum(), rb(hp[win], 5) / hp[win].sum()
    e = H.edges("SR", "mass_fit")[::5]; c = 0.5 * (e[1:] + e[:-1])
    dy = rb(hall["DYmumu|SR|mass_fit|lhe50120"], 5); pw = rb(hall["DYmumu_powheg|SR|mass_fit|nominal"], 5)
    data = rb(hall["Data|SR|mass_fit|nominal"], 5)
    mc = sum(rb(hall[f"{s}|SR|mass_fit|nominal"], 5) for s in ("DYmumu", "DYtautau", "DYee", "TTbar", "SingleTop", "WW", "WZ", "ZZ")
             if f"{s}|SR|mass_fit|nominal" in hall)
    if fakes_res:
        f = np.array(fakes_res["templates"]["nominal"]); mc = mc + np.clip(rb(f, len(f) // len(c)), 0, None)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.step(e[:-1], p_lhe / n_lhe, where="post", color="tab:red", label="powheg / aMC@NLO, generator level (m_LHE, all mu mu events)")
    ax.step(e[:-1], (pw / pw.sum()) / (dy / dy.sum()), where="post", color="tab:orange", linestyle="--",
            label="powheg / aMC@NLO, reconstructed, both 50 < m_LHE < 120 GeV")
    ax.errorbar(c, data / mc, yerr=np.sqrt(data) / mc, fmt="o", color="black", markersize=4, label=r"data / prediction, pre-fit ($\mu_Z$ = 1)")
    ax.axhline(1, color="grey", linestyle=":")
    ax.set_xlim(60, 120); ax.set_ylim(0.9, 1.06); ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("ratio (shapes normalised in 60-120 GeV)")
    ax.legend(fontsize=10, loc="lower right")
    hists._decorate(ax, lumi_fb=config.LUMI_PB_NORMTAG / 1000); hists._title(ax, "Signal lineshape: generator comparison vs data")
    hists.save_fig(fig, "sigmodel_lineshape.png")
    return {"lhe_ratio_5gev": (p_lhe / n_lhe).tolist(), "reco_ratio_5gev": ((pw / pw.sum()) / (dy / dy.sum())).tolist(),
            "data_over_pred_prefit_5gev": (data / mc).tolist()}


def main():
    hall = pickle.load(open(OUT / "histograms.pkl", "rb"))
    fakes_res = json.load(open(OUT / "fakes.json")) if (OUT / "fakes.json").exists() else None
    made = plotting.all_plots(hall, fakes_res)
    print(f"[report] {len(made)} data/MC plots")
    gens = json.load(open(OUT / "gensums.json"))
    lineshape = lineshape_plot(hall, gens, fakes_res) if "DYmumu|SR|mass_fit|lhe50120" in hall else None
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
    stab = json.load(open(STAB)) if STAB.exists() else None
    result = {"generated": str(date.today()), "lumi_pb": config.LUMI_PB_NORMTAG, "yields": yields, "fit": fit, "stability": stab,
              "lineshape": lineshape,
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
                  f"mu_Z = {fit['mu']:.4f} +{fit['mu_err_up']:.4f} -{fit['mu_err_down']:.4f}; goodness of fit p = {fit['gof']['gof_probability']}; "
                  f"signal region in {fit['meta'].get('sr_bin_width_gev', 2)} GeV bins", "",
                  "The statistical uncertainty is the data one; the effective statistical limit of the fit is the MC "
                  "statistics (gammas, see the breakdown), not the data.", "",
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
        lines += [f"| statistical (stat-only fit) | {100*fit['mu_stat_only_fit']:.3f}% |", "",
                  f"Luminosity quoted as the external {100*fit['mu_lumi']:.1f}% (profiled impact {100*fit['mu_lumi_profiled']:.3f}%).", ""]
        if stab:
            lines += ["## Stability of mu_Z against the fit configuration (scripts/v2_5_fit_variants.py)", "",
                      "| configuration | mu_Z | GoF p | notable pulls (constraint) |", "|---|---:|---:|---|"]
            for x in stab:
                lines.append(f"| {x['label']} | {x['mu']:.4f} +{x['err_up']:.4f} -{x['err_down']:.4f} | {x['gof_p']} | {'; '.join(x['notable'][:5])} |")
            mus = [x["mu"] for x in stab]
            lines += ["", f"Spread of mu_Z over the variants: {min(mus):.4f} - {max(mus):.4f} (half-spread {50*(max(mus)-min(mus)):.2f}%).", ""]
        lines += ["## Leading nuisance parameters (ranking)", "", "| NP | pull | constraint | +impact | -impact |", "|---|---:|---:|---:|---:|"]
        for r in fit["ranking"][:12]:
            lines.append(f"| {r['name']} | {r['pull']:+.2f} | {0.5*(r['err_up']+r['err_down']):.2f} | {100*r['dpoi_up_post']:+.3f}% | {100*r['dpoi_down_post']:+.3f}% |")
        lines.append("")
    lines += ["## Yields (mass_fit templates; SR/SS: prompt-prompt MC + data-driven fakes; e-mu regions: all MC)", "",
              "| region | sample | events |", "|---|---|---:|"]
    for region, d in yields.items():
        for smp, n in sorted(d.items(), key=lambda kv: -kv[1]):
            if smp in NOT_IN_FIT and not (smp == "WJets" and region in ("CRemu", "SSemu")):
                continue
            lines.append(f"| {region} | {smp} | {n:,.1f} |")
    lines += ["", "Not listed: " + "; ".join(f"`{k}` ({v})" for k, v in NOT_IN_FIT.items()) + "."]
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
