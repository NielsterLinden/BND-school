#!/usr/bin/env python
"""Step 6 (v4) -- collect the four-channel result into output_v4/results.json, output_v4/RESULTS.md and summary plots.

    python scripts/step6_report_v4.py

Reads fit_v4/results/ztautau_v4_fit_result.json (+ the per-channel fits if they were run with
step5_fit_v4.py --channels <ch> --job ztautau_v4_<ch>), output_v4/data/yields_v4.json, fakes_*.json and the
in-situ trigger efficiencies; copies the key TRExFitter plots to output_v4/plots/fit_*.png.
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

from ztautau import analysis_v4 as an, config, plotting  # noqa: E402

V3 = config.OUTPUT_DIR / "results.json"


def load_json(path):
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else None


def summarise(res):
    if res is None:
        return None
    mu, up, dn = res["poi_value"], res["poi_err_up"], res["poi_err_down"]
    stat = (res.get("poi_stat_only") or {}).get("err")
    tot = 0.5 * (up + dn)
    syst = float(np.sqrt(max(tot ** 2 - (stat or 0) ** 2, 0.0)))
    s60 = res["prediction"]["sigma_tautau_60_120_pb"]
    out = {"mu": mu, "mu_err_up": up, "mu_err_down": dn, "mu_stat": stat, "mu_syst": syst, "mu_expected_asimov": res.get("poi_asimov"),
           "gof_probability": (res.get("gof") or {}).get("gof_probability"), "grouped_impact": res.get("grouped_impact", {}),
           "sigma_60_120_pb": {"value": mu * s60, "stat": (stat or 0) * s60, "syst": syst * s60, "err_up": up * s60, "err_down": dn * s60, "prediction": s60},
           "tau_id_sf": res.get("tau_id_sf"), "tau_es": res.get("tau_es"), "mu_ttbar": res.get("mu_ttbar"), "channels": res.get("channels"),
           "tau_id_fixed": res.get("tau_id_fixed", False)}
    ranking = sorted(res.get("ranking", []), key=lambda r: -max(abs(r["dpoi_up_post"]), abs(r["dpoi_down_post"])))
    out["ranking"] = [{"name": r["name"], "pull": r["pull"], "constraint": 0.5 * (abs(r["err_up"]) + abs(r["err_down"])),
                       "impact_up": r["dpoi_up_post"], "impact_down": r["dpoi_down_post"]} for r in ranking[:20]]
    out["pulls"] = {k: v for k, v in res["nps"].items() if not k.startswith("gamma")}
    return out


def copy_trex_plots(job, tag):
    base = config.FIT_DIR_V4 / "results" / job
    pairs = [("NuisPar.png", f"fit_pulls{tag}.png"), ("CorrMatrix.png", f"fit_corrmatrix{tag}.png"),
             ("Rankings/RankingSysts_mu_Z_systs.png", f"fit_ranking{tag}.png"), ("NormFactors.png", f"fit_normfactors{tag}.png"),
             ("Plots/Summary_postFit.png", f"fit_summary_postfit{tag}.png"), ("Plots/Summary.png", f"fit_summary_prefit{tag}.png")]
    for region in config.REGIONS_V4:
        pairs += [(f"Plots/{region}_postFit.png", f"fit_{region}_postfit{tag}.png"), (f"Plots/{region}.png", f"fit_{region}_prefit{tag}.png")]
    for src, dst in pairs:
        p = base / src
        if not p.exists():
            cand = sorted(base.glob(f"**/{Path(src).name}"))
            p = cand[0] if cand else None
        if p and p.exists():
            shutil.copy(p, config.PLOT_DIR_V4 / dst)


def main():
    config.PLOT_DIR_V4.mkdir(parents=True, exist_ok=True)
    # per channel: the cross-check with the POG tau ID scale factors fixed (<job>_<ch>_fixedid) is what the
    # per-channel table shows; <job>_<ch> (tau ID free) are the workspaces of the MultiFit combination
    jobs = {"combined": config.JOB_V4, **{ch: f"{config.JOB_V4}_{ch}_fixedid" for ch in config.CHANNELS}}
    free = {ch: load_json(config.FIT_DIR_V4 / "results" / f"{config.JOB_V4}_{ch}_fit_result.json") for ch in config.CHANNELS}
    comb_mf = load_json(config.FIT_DIR_V4 / "results" / "comb_v4_fit_result.json")
    res = {k: load_json(config.FIT_DIR_V4 / "results" / f"{j}_fit_result.json") for k, j in jobs.items()}
    fits = {k: summarise(r) for k, r in res.items() if r}
    yields = load_json(config.DATA_DIR_V4 / "yields_v4.json")
    fakes = {ch: load_json(config.DATA_DIR_V4 / f"fakes_{ch}.json") for ch in ("mutau", "etau", "emu")}
    trig = load_json(an.TRIG_INSITU)
    v3 = load_json(V3)
    meta = load_json(config.FIT_DIR_V4 / "fitinputs" / f"{config.JOB_V4}.root.meta.json")
    results = {"version": config.VERSION, "channels": config.CHANNELS, "dataset": "CMS Open Data 2016 Run2016G+H (UL NanoAODv9): Tau, SingleMuon, SingleElectron, MuonEG",
               "lumi_pb": config.LUMI_PB, "signal_definition": "Z/gamma* -> tautau, 60 < m_LHE < 120 GeV, all decays",
               "fit": fits, "yields_prefit": yields, "fakes": {ch: ({k: v for k, v in f.items() if k in ("osss", "w_mt", "closure_ss", "sr_fakes", "sr_subtracted_mc", "n_data", "SB1", "SB2", "sr_multijet", "ss_region")}) for ch, f in fakes.items() if f},
               "trigger_insitu": {k: {"sf_plateau": np.asarray(v["sf"])[-3:, :].tolist(), "err_plateau": np.asarray(v["err"])[-3:, :].tolist()} for k, v in (trig or {}).items()},
               "v3_reference": (v3 or {}).get("fit", {}).get("mcsub", {}).get("sigma_60_120_pb"), "n_templates": meta and len(meta["samples"]), "n_systs": meta and len(meta["systs"]),
               "per_channel_workspaces": {ch: {"fitinputs": f"fit_v4/fitinputs/{config.JOB_V4}_{ch}.root", "config": f"fit_v4/{config.JOB_V4}_{ch}.config",
                                               "workspace": f"fit_v4/results/{config.JOB_V4}_{ch}/RooStats/{config.JOB_V4}_{ch}_combined_{config.JOB_V4}_{ch}_model.root",
                                               "mu_Z_alone_tauid_free": (r["poi_value"], r["poi_err_up"], r["poi_err_down"]) if r else None}
                                          for ch, r in free.items()},
               "multifit_check": comb_mf,
               "tau_id_sf_pog": meta and meta["tau_id_sf_pog"]}
    config.OUTPUT_DIR_V4.mkdir(parents=True, exist_ok=True)
    (config.OUTPUT_DIR_V4 / "results.json").write_text(json.dumps(results, indent=1, default=float))
    for k, j in jobs.items():
        if res.get(k):
            copy_trex_plots(j, "" if k == "combined" else f"_{k}")
    f = fits.get("combined")
    lines = [f"# Z -> tautau ({config.VERSION}): four-channel result", "",
             f"CMS Open Data 2016 (Run2016G+H), L = {config.LUMI_PB:.1f} pb^-1. Channels: tau_h tau_h (3 BDT categories), mu tau_h, e tau_h, e mu "
             f"(+ e mu ttbar control region). Simultaneous binned profile-likelihood fit of m_tt with TRExFitter v1.8.0; the tau_h ID scale factors "
             f"(per decay mode) are free parameters, the tau_h energy scale has a {100 * config.TES_PRIOR_V4:.0f}% prior. Generated by scripts/step6_report_v4.py.", ""]
    if f:
        s = f["sigma_60_120_pb"]
        lines += ["## Combined fit", "", "| quantity | value |", "|---|---|",
                  f"| mu_Z | {f['mu']:.3f} +{f['mu_err_up']:.3f} -{f['mu_err_down']:.3f} (stat {f['mu_stat']}, syst {f['mu_syst']:.3f}) |",
                  f"| sigma(pp -> Z/gamma* -> tautau, 60 < m < 120) | {s['value']:.0f} +{s['err_up']:.0f} -{s['err_down']:.0f} pb (stat {s['stat']:.0f}) (pred. {s['prediction']:.0f} pb) |",
                  f"| expected (Asimov) | +{f['mu_expected_asimov']['err_up']:.3f} -{f['mu_expected_asimov']['err_down']:.3f} |" if f["mu_expected_asimov"] else "| expected | n/a |",
                  f"| goodness of fit | p = {f['gof_probability']} |",
                  f"| mu_ttbar | {f['mu_ttbar']} |", ""]
        if f.get("tau_id_sf"):
            lines += ["tau_h ID scale factors (fitted, TauPOG prior-free):", "", "| DM | fitted | TauPOG |", "|---|---|---|"]
            for dm, v in f["tau_id_sf"].items():
                lines.append(f"| {dm} | {v['value']:.3f} +{v['err_up']:.3f} -{v['err_down']:.3f} | {v['pog'][0]:.3f} +- {v['pog'][1]:.3f} |")
            lines.append("")
        if f.get("tau_es"):
            lines += ["tau_h energy scale (pull and post-fit constraint in units of the 3% prior):", "", "| DM | pull | constraint | -> % |", "|---|---|---|---|"]
            for dm, v in f["tau_es"].items():
                lines.append(f"| {dm} | {v['pull']:+.2f} | {v['constraint']:.2f} | {v['constraint'] * v['prior_pct']:.2f}% |")
            lines.append("")
        lines += ["Grouped impacts on mu_Z:", "", "| group | impact |", "|---|---|"]
        for k, val in sorted(f["grouped_impact"].items(), key=lambda kv: -kv[1]):
            lines.append(f"| {k} | {val:.4f} |")
        lines.append("")
    per = [(ch, fits[ch]) for ch in config.CHANNELS if ch in fits]
    if per:
        lines += ["## Per-channel fits (same model, one channel each; tau_h ID scale factors fixed to the TauPOG values where marked)", "",
                  "| channel | tau_h ID SF | mu_Z | sigma(60-120) [pb] |", "|---|---|---|---|"]
        for ch, fc in per:
            s = fc["sigma_60_120_pb"]
            lines.append(f"| {ch} | {'fixed (POG)' if fc.get('tau_id_fixed') else 'free'} | {fc['mu']:.3f} +{fc['mu_err_up']:.3f} -{fc['mu_err_down']:.3f} | {s['value']:.0f} +{s['err_up']:.0f} -{s['err_down']:.0f} |")
        lines.append("")
    lines += ["## Per-channel exports for the combination", "",
              "| channel | fit inputs | config | workspace | mu_Z alone (tau ID free) |", "|---|---|---|---|---|"]
    for ch, r in free.items():
        w = f"fit_v4/results/{config.JOB_V4}_{ch}/RooStats/{config.JOB_V4}_{ch}_combined_{config.JOB_V4}_{ch}_model.root"
        ok = (config.FIT_DIR_V4.parent / w).exists()
        mu = f"{r['poi_value']:.3f} +{r['poi_err_up']:.3f} -{r['poi_err_down']:.3f}" if r else "n/a"
        lines.append(f"| {ch} | `fit_v4/fitinputs/{config.JOB_V4}_{ch}.root` | `fit_v4/{config.JOB_V4}_{ch}.config` | {'present' if ok else 'missing'} | {mu} |")
    if comb_mf:
        lines += ["", f"MultiFit of the four workspaces (`fit_v4/comb_v4.config`): mu_Z = {comb_mf['poi_value']:.3f} +{comb_mf['poi_err_up']:.3f} -{comb_mf['poi_err_down']:.3f}"
                  f" (single-file fit: {f['mu']:.3f} +{f['mu_err_up']:.3f} -{f['mu_err_down']:.3f})" if f else ""]
    lines.append("")
    if v3:
        r3 = v3["fit"]["mcsub"]
        lines += [f"v3 reference (tau_h tau_h only, fiducial signal): mu_Z = {r3['mu']:.3f} +{r3['mu_err_up']:.3f} -{r3['mu_err_down']:.3f}, "
                  f"sigma(60-120) = {r3['sigma_60_120_pb']['value']:.0f} +{r3['sigma_60_120_pb']['err_up']:.0f} -{r3['sigma_60_120_pb']['err_down']:.0f} pb", ""]
    if yields:
        lines += ["## Prefit yields per region", "", "| region | " + " | ".join(sorted({k for r in yields["regions"].values() for k in r})) + " |"]
        keys = sorted({k for r in yields["regions"].values() for k in r})
        lines.append("|---|" + "---:|" * len(keys))
        for r, y in yields["regions"].items():
            lines.append(f"| {r} | " + " | ".join(f"{y.get(k, 0):.0f}" for k in keys) + " |")
        lines.append("")
    for ch, fk in fakes.items():
        if not fk:
            continue
        if ch == "emu":
            lines += [f"e mu multijet: OS/SS (SB1) {np.round(fk['SB1']['ratio'], 2).tolist()} vs SB2 {np.round(fk['SB2']['ratio'], 2).tolist()} per dR bin; "
                      f"SS region {fk['ss_region']['data']} data, {fk['ss_region']['mc']:.0f} MC -> {fk['sr_multijet']:.0f} multijet events in the SR", ""]
        else:
            lines += [f"{ch} fakes: C_OS/SS = {fk['osss']['C']:.3f} +- {fk['osss']['stat']:.3f}; r_W = " + ", ".join(f"{g} {v['r']:.2f}" for g, v in fk["w_mt"].items())
                      + f"; same-sign closure obs/pred = {fk['closure_ss']['ratio']:.3f} +- {fk['closure_ss']['stat']:.3f}; SR fakes {fk['sr_fakes']:.0f}", ""]
    (config.OUTPUT_DIR_V4 / "RESULTS.md").write_text("\n".join(lines) + "\n")
    # summary plot: combined and per-channel sigma(60-120), v3 reference
    rows = []
    if f:
        rows.append(("combined (4 channels)", f["sigma_60_120_pb"]))
    for ch, fc in per:
        rows.append((ch, fc["sigma_60_120_pb"]))
    if v3:
        rows.append(("v3 tau_h tau_h (fiducial POI)", v3["fit"]["mcsub"]["sigma_60_120_pb"]))
    if rows:
        plt.style.use("default")
        fig, ax = plt.subplots(figsize=(8, 1.0 + 0.7 * len(rows)))
        pred = rows[0][1]["prediction"]
        ax.axvspan(pred * 0.96, pred * 1.04, color="#ffcc66", alpha=0.5); ax.axvline(pred, color="#c08a00")
        for i, (lab, s) in enumerate(rows):
            yv = len(rows) - i
            ax.errorbar([s["value"]], [yv], xerr=[[s["err_down"]], [s["err_up"]]], fmt="o", color="black", capsize=4, lw=1.5)
            if s.get("stat"):
                ax.errorbar([s["value"]], [yv], xerr=[[s["stat"]], [s["stat"]]], fmt="none", color="tab:red", lw=4)
            ax.text(1230, yv + 0.3, lab, fontsize=10); ax.text(3270, yv + 0.3, f"{s['value']:.0f} +{s['err_up']:.0f} -{s['err_down']:.0f} pb", fontsize=10, ha="right")
        ax.set_xlim(1200, 3300); ax.set_ylim(0.3, len(rows) + 0.9); ax.set_yticks([])
        ax.set_xlabel(r"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau)$, $60<m<120$ GeV [pb]")
        ax.set_title(r"Z$\rightarrow\tau\tau$ v4, CMS Open Data 2016 G+H, 16.4 fb$^{-1}$ (red: stat.)", fontsize=11)
        plotting.fig_tag(fig)
        fig.savefig(config.PLOT_DIR_V4 / "step6_summary.png", bbox_inches="tight", dpi=130); plt.close(fig)
    if f:
        gi = {k: v for k, v in f["grouped_impact"].items() if k != "FullSyst"}
        if f["mu_stat"]:
            gi["Data statistics"] = f["mu_stat"]
        items = sorted(gi.items(), key=lambda kv: kv[1])
        fig, ax = plt.subplots(figsize=(7, 0.35 * len(items) + 1.5))
        ax.barh([k for k, _ in items], [100 * v for _, v in items], color="#4a90d9")
        for i, (_, val) in enumerate(items):
            ax.text(100 * val + 0.1, i, f"{100 * val:.2f}%", va="center", fontsize=9)
        ax.set_xlabel(r"impact on $\mu_Z$ [%]"); ax.set_title("Grouped uncertainties, four-channel fit", fontsize=11)
        plotting.fig_tag(fig)
        fig.savefig(config.PLOT_DIR_V4 / "step6_impacts.png", bbox_inches="tight", dpi=130); plt.close(fig)
    print((config.OUTPUT_DIR_V4 / "RESULTS.md").read_text())


if __name__ == "__main__":
    main()
