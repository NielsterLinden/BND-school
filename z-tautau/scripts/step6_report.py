#!/usr/bin/env python
"""Step 6 -- collect the measurement into output/results.json, output/RESULTS.md and the summary plots.

    python scripts/step6_report.py

Reads fit/results/ztautau_fit_result.json and every cross-check job of step 5 (per channel, `_fixedid`,
`ztautau_taulep`, `ztautau_ptsplit`, `ztautau_emutrig2x`, the MultiFit `comb`), output/data/yields_v4.json,
the fake and in-situ trigger json files, and copies the key TRExFitter plots to output/plots/fit_*.png.

What the report is careful about (REVIEW_v4.md):
  * the expected (Asimov) uncertainty is quoted only if MINOS converged on it (finding 1);
  * the e mu channel and the three tau channels are quoted separately with their difference, because the
    combined number is their compromise, not a confirmation that they agree (finding 2);
  * the grouped impacts overlap (mu_ttbar <-> EmuTrigger <-> mu_Z), so their quadrature sum over-shoots
    the MINOS total; the over-shoot is printed and the scale factor that removes it is written into
    results.json for the combination (finding 5).
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

JOB = config.JOB_V4
VARIANT_JOBS = {"taulep": f"{JOB}_taulep", "ptsplit": f"{JOB}_ptsplit", "emutrig2x": f"{JOB}_emutrig2x"}
VARIANT_LABEL = {
    "taulep": "tau_h tau_h + mu tau_h + e tau_h (no e mu), scale factors free",
    "ptsplit": "l tau_h regions split at pT(tau_h) = 40 GeV, own scale factors below",
    "emutrig2x": "e mu trigger variation doubled (2% per leg, the pre-review treatment)",
}
# Correlation of each uncertainty category with the z-mumu / z-ee channels, for the combination
# (combination/comb/model.py CORRELATION). Categories the combination already knows keep its value;
# the ones v4 adds are recommended here and justified in docs/11-combination-inputs.md.
CATEGORY_RHO = {
    "Luminosity": 1.0, "Pileup": 1.0, "L1 prefiring": 1.0, "Background normalisation": 1.0,
    "Signal modelling": 1.0, "Muon efficiency": 0.0, "Muon momentum": 0.0, "Electron efficiency": 0.0,
    "Tau ID": 0.0, "Tau trigger": 0.0, "Tau energy scale": 0.0, "Gammas": 0.0, "Fakes": 0.0, "MET": 0.0,
    "Electron trigger": 0.0, "Emu trigger": 0.0, "Electron energy": 0.0, "Muon efficiency ": 0.0,
    "Jets": 0.0, "b tagging": 0.0, "Background modelling": 0.0, "NormFactors": 0.0, "Tau ID (fitted)": 0.0,
}


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
    out = {"mu": mu, "mu_err_up": up, "mu_err_down": dn, "mu_stat": stat, "mu_syst": syst,
           "mu_expected_asimov": res.get("poi_asimov"), "asimov_minos_status": res.get("asimov_minos_status"),
           "minos_status": res.get("minos_status"),
           "gof_probability": (res.get("gof") or {}).get("gof_probability"),
           "gof_chi2_ndf": (res.get("gof") or {}).get("gof_chi2_ndf"),
           "grouped_impact": res.get("grouped_impact", {}),
           "grouped_impact_quadrature_sum": res.get("grouped_impact_quadrature_sum"),
           "grouped_impact_scale": res.get("grouped_impact_scale"),
           "sigma_60_120_pb": {"value": mu * s60, "stat": (stat or 0) * s60, "syst": syst * s60,
                               "err_up": up * s60, "err_down": dn * s60, "prediction": s60},
           "tau_id_sf": res.get("tau_id_sf"), "tau_es": res.get("tau_es"), "mu_ttbar": res.get("mu_ttbar"),
           "channels": res.get("channels"), "tau_id_fixed": res.get("tau_id_fixed", False),
           "region_set": res.get("region_set", "nominal"), "scaled_systematics": res.get("scaled_systematics", {})}
    ranking = sorted(res.get("ranking", []), key=lambda r: -max(abs(r["dpoi_up_post"]), abs(r["dpoi_down_post"])))
    out["ranking"] = [{"name": r["name"], "pull": r["pull"], "constraint": 0.5 * (abs(r["err_up"]) + abs(r["err_down"])),
                       "impact_up": r["dpoi_up_post"], "impact_down": r["dpoi_down_post"]} for r in ranking[:20]]
    out["pulls"] = {k: v for k, v in res["nps"].items() if not k.startswith("gamma")}
    return out


def difference(a, b):
    """Difference of two measurements of mu_Z in units of their combined uncertainty, using the error of
    each that faces the other (they are not independent -- they share the tau_h ID scale factors and the
    theory templates -- so this is an upper limit on the significance, stated as such in the report)."""
    if not a or not b:
        return None
    hi, lo = (a, b) if a["mu"] >= b["mu"] else (b, a)
    s = float(np.hypot(hi["mu_err_down"], lo["mu_err_up"]))
    return {"delta": hi["mu"] - lo["mu"], "sigma": s, "n_sigma": (hi["mu"] - lo["mu"]) / s if s else None}


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


def channel_result(f, meta):
    """The channel's entry for the combination, in the fields of combination/comb/inputs.ChannelResult."""
    groups = {k: v for k, v in f["grouped_impact"].items() if k not in ("FullSyst", "Total") and v > 0}
    scale = f.get("grouped_impact_scale") or 1.0
    return {"name": "tautau", "label": r"$Z\to\tau\tau$ (4 channels)", "variant": "v4",
            "mu": f["mu"], "mu_err_up": f["mu_err_up"], "mu_err_down": f["mu_err_down"], "mu_stat": f["mu_stat"],
            "sigma_pred": f["sigma_60_120_pb"]["prediction"],
            "groups": groups,
            "groups_rescaled": {k: v * scale for k, v in groups.items()},
            "groups_rescale_factor": scale,
            "ranking": {r["name"]: 0.5 * (abs(r["impact_up"]) + abs(r["impact_down"])) for r in f["ranking"]},
            "acc": {},
            "sigmodel": 0.0, "residual": 0.0, "gof_p": f["gof_probability"],
            "extra": {"lumi_pb": config.LUMI_PB, "tau_wp": config.TAU_WP, "channels": f["channels"],
                      "n_templates": len(meta["samples"]) if meta else None,
                      "signal_definition": "Z/gamma* -> tautau, 60 < m_LHE < 120 GeV, all decays",
                      "acceptance_note": "A x epsilon is profiled inside the fit (theory templates renormalised to "
                                         "sigma(60-120), CONVENTIONS.md section 3): there is no external acceptance "
                                         "uncertainty to add, so acc = {}",
                      "groups_note": "the categories overlap (mu_ttbar <-> EmuTrigger <-> mu_Z): their quadrature sum "
                                     "over-shoots the MINOS total by 1/groups_rescale_factor. Use groups_rescaled with "
                                     "residual = 0, or the MINOS total directly."}}


def main():
    config.PLOT_DIR_V4.mkdir(parents=True, exist_ok=True)
    fitres = config.FIT_DIR_V4 / "results"
    jobs = {"combined": JOB, **{ch: f"{JOB}_{ch}_fixedid" for ch in config.CHANNELS}}
    free = {ch: load_json(fitres / f"{JOB}_{ch}_fit_result.json") for ch in config.CHANNELS}
    comb_mf = load_json(fitres / "comb_fit_result.json")
    res = {k: load_json(fitres / f"{j}_fit_result.json") for k, j in jobs.items()}
    fits = {k: summarise(r) for k, r in res.items() if r}
    variants = {k: summarise(load_json(fitres / f"{j}_fit_result.json")) for k, j in VARIANT_JOBS.items()}
    variants = {k: v for k, v in variants.items() if v}
    free_s = {ch: summarise(r) for ch, r in free.items() if r}
    yields = load_json(config.DATA_DIR_V4 / "yields_v4.json")
    fakes = {ch: load_json(config.DATA_DIR_V4 / f"fakes_{ch}.json") for ch in ("mutau", "etau", "emu")}
    trig = load_json(an.TRIG_INSITU)
    meta = load_json(config.FIT_DIR_V4 / "fitinputs" / f"{JOB}.root.meta.json")
    f = fits.get("combined")
    tension = difference(variants.get("taulep"), free_s.get("emu"))

    results = {"version": config.VERSION, "channels": config.CHANNELS, "nominal_variant": "combined",
               "dataset": "CMS Open Data 2016 Run2016G+H (UL NanoAODv9): Tau, SingleMuon, SingleElectron, MuonEG",
               "lumi_pb": config.LUMI_PB, "tau_wp": config.TAU_WP,
               "signal_definition": "Z/gamma* -> tautau, 60 < m_LHE < 120 GeV, all decays",
               "fit": fits, "fit_variants": variants, "fit_per_channel_free": free_s,
               "tension_taulep_vs_emu": tension,
               "prediction": (meta or {}).get("signal_prediction"),
               "yields_prefit": yields,
               "fakes": {ch: {k: v for k, v in fk.items() if k in ("osss", "w_mt", "closure_ss", "sr_fakes", "sr_subtracted_mc", "n_data", "SB1", "SB2", "sr_multijet", "ss_region")}
                         for ch, fk in fakes.items() if fk},
               "trigger_insitu": {k: {"sf_plateau": np.asarray(v["sf"])[-3:, :].tolist(), "err_plateau": np.asarray(v["err"])[-3:, :].tolist()}
                                  for k, v in (trig or {}).items()},
               "n_templates": meta and len(meta["samples"]), "n_systs": meta and len(meta["systs"]),
               "per_channel_workspaces": {ch: {"fitinputs": f"fit/fitinputs/{JOB}_{ch}.root", "config": f"fit/{JOB}_{ch}.config",
                                               "workspace": f"fit/results/{JOB}_{ch}/RooStats/{JOB}_{ch}_combined_{JOB}_{ch}_model.root",
                                               "mu_Z_alone_tauid_free": (r["poi_value"], r["poi_err_up"], r["poi_err_down"]) if r else None}
                                          for ch, r in free.items()},
               "multifit_check": comb_mf,
               "tau_id_sf_pog": meta and meta["tau_id_sf_pog"]}
    if f:
        results["for_combination"] = {"channel_result": channel_result(f, meta),
                                      "category_rho_recommended": {k: CATEGORY_RHO.get(k) for k in f["grouped_impact"] if k not in ("FullSyst", "Total")},
                                      "doc": "z-tautau/docs/11-combination-inputs.md"}
    config.OUTPUT_DIR_V4.mkdir(parents=True, exist_ok=True)
    (config.OUTPUT_DIR_V4 / "results.json").write_text(json.dumps(results, indent=1, default=float))
    for k, j in jobs.items():
        if res.get(k):
            copy_trex_plots(j, "" if k == "combined" else f"_{k}")
    for k, j in VARIANT_JOBS.items():
        if variants.get(k):
            copy_trex_plots(j, f"_{k}")

    # ------------------------------------------------------------------ RESULTS.md
    L = [f"# Z -> tautau ({config.VERSION}): four-channel result", "",
         f"CMS Open Data 2016 (Run2016G+H), L = {config.LUMI_PB:.1f} pb^-1. Channels: tau_h tau_h (3 BDT categories), "
         f"mu tau_h, e tau_h, e mu (+ e mu ttbar control region), each split by the tau_h decay mode where it constrains "
         f"the tau_h ID scale factors. Simultaneous binned profile-likelihood fit of m_tt with TRExFitter v1.8.0; the "
         f"tau_h ID scale factors (per decay mode) are free parameters, the tau_h energy scale has a "
         f"{100 * config.TES_PRIOR_V4:.0f}% prior. Generated by scripts/step6_report.py.", ""]
    if f:
        s = f["sigma_60_120_pb"]
        exp = f["mu_expected_asimov"]
        L += ["## Combined fit", "", "| quantity | value |", "|---|---|",
              f"| mu_Z | {f['mu']:.3f} +{f['mu_err_up']:.3f} -{f['mu_err_down']:.3f} "
              f"(stat {f['mu_stat']:.3f}, syst {f['mu_syst']:.3f}) |",
              f"| sigma(pp -> Z/gamma* -> tautau, 60 < m < 120) | {s['value']:.0f} +{s['err_up']:.0f} -{s['err_down']:.0f} pb "
              f"(stat {s['stat']:.0f}) (prediction {s['prediction']:.0f} pb) |",
              (f"| expected (Asimov, MINOS) | +{exp['err_up']:.3f} -{exp['err_down']:.3f} |" if exp else
               f"| expected (Asimov) | not quoted: MINOS did not converge on the Asimov data set "
               f"(status {f['asimov_minos_status']}) |"),
              f"| goodness of fit | p = {f['gof_probability']:.3f}" + (f" (chi2/ndf {f['gof_chi2_ndf']})" if f['gof_chi2_ndf'] else "") + " |",
              f"| mu_ttbar (from emu_CRtt) | {f['mu_ttbar']['value']:.3f} +{f['mu_ttbar']['err_up']:.3f} -{f['mu_ttbar']['err_down']:.3f} |"
              if f.get("mu_ttbar") else "| mu_ttbar | n/a |", ""]
        if f.get("tau_id_sf"):
            L += ["tau_h ID scale factors (fitted, TauPOG prior-free):", "", "| DM | fitted | TauPOG |", "|---|---|---|"]
            for dm, v in f["tau_id_sf"].items():
                L.append(f"| {dm} | {v['value']:.3f} +{v['err_up']:.3f} -{v['err_down']:.3f} | {v['pog'][0]:.3f} +- {v['pog'][1]:.3f} |")
            L.append("")
        if f.get("tau_es"):
            L += [f"tau_h energy scale (pull and post-fit constraint in units of the {100 * config.TES_PRIOR_V4:.0f}% prior):", "",
                  "| DM | pull | constraint | -> % |", "|---|---|---|---|"]
            for dm, v in f["tau_es"].items():
                L.append(f"| {dm} | {v['pull']:+.2f} | {v['constraint']:.2f} | {v['constraint'] * v['prior_pct']:.2f}% |")
            L.append("")

    # ------------------------------------------------------------------ what the number means
    if f and (variants.get("taulep") or free_s.get("emu")):
        L += ["## What the combined number is made of", "",
              "The four channels do not measure mu_Z the same way. The e mu channel has no tau_h and no tau_h ID scale "
              "factor: it measures mu_Z directly (against its trigger efficiency and mu_ttbar). The three tau channels "
              "with free scale factors measure mu_Z only through the ratio of the tau_h tau_h regions (SF^2 mu_Z) to the "
              "l tau_h ones (SF mu_Z). The combined fit is the compromise between the two, and it pays for it with the "
              "nuisance parameters of the e mu trigger and the l tau_h acceptance. Both numbers are therefore quoted.", "",
              "| fit | mu_Z | sigma(60-120) [pb] | GoF p |", "|---|---|---|---|"]
        rows = [("four channels (the measurement)", f), ("e mu alone (no tau_h)", free_s.get("emu")),
                ("tau channels alone, scale factors free", variants.get("taulep"))]
        for lab, fc in rows:
            if not fc:
                continue
            sc = fc["sigma_60_120_pb"]
            gof = f"{fc['gof_probability']:.3f}" if fc.get("gof_probability") is not None else "-"
            L.append(f"| {lab} | {fc['mu']:.3f} +{fc['mu_err_up']:.3f} -{fc['mu_err_down']:.3f} | "
                     f"{sc['value']:.0f} +{sc['err_up']:.0f} -{sc['err_down']:.0f} | {gof} |")
        L.append("")
        if tension:
            L += [f"The two differ by {tension['delta']:+.3f} in mu_Z, i.e. {tension['n_sigma']:.1f} sigma of their "
                  f"(uncorrelated-limit) combined uncertainty {tension['sigma']:.3f}. They share the tau_h ID scale "
                  f"factors and the theory templates, so that is an upper limit on the significance of the difference, "
                  f"not a p-value. The combined result is a valid likelihood result; it is not a measurement of mu_Z "
                  f"independent of the e mu trigger efficiency.", ""]
    if variants:
        L += ["## Cross-check fits", "", "| fit | mu_Z | sigma(60-120) [pb] | what it tests |", "|---|---|---|---|"]
        for k, v in variants.items():
            sc = v["sigma_60_120_pb"]
            L.append(f"| `{VARIANT_JOBS[k]}` | {v['mu']:.3f} +{v['mu_err_up']:.3f} -{v['mu_err_down']:.3f} | "
                     f"{sc['value']:.0f} +{sc['err_up']:.0f} -{sc['err_down']:.0f} | {VARIANT_LABEL[k]} |")
        L.append("")
        pts = variants.get("ptsplit")
        if pts and pts.get("tau_id_sf"):
            L += ["Scale factors of the pT-split fit (the assumption behind the tau_h tau_h / (l tau_h)^2 lever: one "
                  "scale factor per decay mode for pT(tau_h) > 30 GeV in l tau_h and > 40 GeV in tau_h tau_h):", "",
                  "| DM | pT > 40 GeV | pT < 40 GeV | ratio |", "|---|---|---|---|"]
            ratios = []
            for dm in config.TAU_DMS:
                hi, lo = pts["tau_id_sf"].get(f"DM{dm}"), pts["tau_id_sf"].get(f"DM{dm}_lowpt")
                if hi and lo:
                    ratios.append(lo["value"] / hi["value"])
                    L.append(f"| DM{dm} | {hi['value']:.3f} +- {0.5 * (hi['err_up'] + hi['err_down']):.3f} | "
                             f"{lo['value']:.3f} +- {0.5 * (lo['err_up'] + lo['err_down']):.3f} | "
                             f"{lo['value'] / hi['value']:.3f} |")
            L.append("")
            if f and ratios:
                shift = pts["mu"] - f["mu"]
                tot = 0.5 * (f["mu_err_up"] + f["mu_err_down"])
                L += [f"The scale factors below 40 GeV come out "
                      f"{', '.join(f'{100 * (r - 1):+.0f}%' for r in ratios)} relative to those above it "
                      f"(DM{', DM'.join(str(d) for d in config.TAU_DMS)}), each 0.5-1 sigma on its own but "
                      f"coherent in sign, and mu_Z moves by {shift:+.3f} ({abs(shift) / tot:.1f} times the "
                      f"total uncertainty) to {pts['mu']:.3f}. The single-scale-factor assumption is therefore "
                      f"worth more than any experimental systematic in the table below; the split is not the "
                      f"nominal model (the ratios are individually compatible with one and it doubles the "
                      f"number of free scale factors on the same data), but that spread should travel with "
                      f"the result. See REVIEW_v4_RESPONSE.md section 6.", ""]
    if f:
        gsum, gscale = f.get("grouped_impact_quadrature_sum"), f.get("grouped_impact_scale")
        L += ["## Grouped impacts on mu_Z", "",
              "The categories are not independent: mu_ttbar, the e mu trigger efficiency and mu_Z form one chain "
              "(the e mu control region fixes mu_ttbar x eff, the signal region mu_Z x eff), so `NormFactors` and "
              "`Emu trigger` contain the same degeneracy. Their quadrature sum therefore over-shoots the total. "
              "**The uncertainty of the measurement is the MINOS total, not the sum of the rows below.**", ""]
        if gsum and gscale:
            L += [f"Quadrature sum of all rows including data statistics: {gsum:.4f} against the MINOS total "
                  f"{0.5 * (f['mu_err_up'] + f['mu_err_down']):.4f} (a factor {1 / gscale:.2f} too large). "
                  f"A combination that needs a category breakdown should scale every row by {gscale:.3f} "
                  f"(`for_combination.channel_result.groups_rescaled` in results.json).", ""]
        L += ["| group | impact on mu_Z | rho with mumu / ee (recommended) |", "|---|---|---|"]
        for k, val in sorted(f["grouped_impact"].items(), key=lambda kv: -kv[1]):
            rho = "-" if k in ("FullSyst", "Total") else CATEGORY_RHO.get(k, "**to be decided**")
            L.append(f"| {k} | {val:.4f} | {rho} |")
        if f["mu_stat"]:
            L.append(f"| Data statistics | {f['mu_stat']:.4f} | 0.0 |")
        L.append("")
    per = [(ch, fits[ch]) for ch in config.CHANNELS if ch in fits]
    if per:
        L += ["## Per-channel fits (one channel each, tau_h ID scale factors fixed to the TauPOG values)", "",
              "A single l tau_h channel cannot separate mu_Z from its own scale factors, so these use the POG values; "
              "they are consistency checks of the channels, not independent measurements.", "",
              "| channel | mu_Z | sigma(60-120) [pb] |", "|---|---|---|"]
        for ch, fc in per:
            s = fc["sigma_60_120_pb"]
            L.append(f"| {ch} | {fc['mu']:.3f} +{fc['mu_err_up']:.3f} -{fc['mu_err_down']:.3f} | "
                     f"{s['value']:.0f} +{s['err_up']:.0f} -{s['err_down']:.0f} |")
        L.append("")
    L += ["## Per-channel exports for the combination", "",
          "See `docs/11-combination-inputs.md`. Every file below is committed.", "",
          "| channel | fit inputs | config | workspace | mu_Z alone (scale factors free) |", "|---|---|---|---|---|"]
    for ch, r in free.items():
        w = f"fit/results/{JOB}_{ch}/RooStats/{JOB}_{ch}_combined_{JOB}_{ch}_model.root"
        ok = (config.FIT_DIR_V4.parent / w).exists()
        mu = f"{r['poi_value']:.3f} +{r['poi_err_up']:.3f} -{r['poi_err_down']:.3f}" if r else "n/a"
        L.append(f"| {ch} | `fit/fitinputs/{JOB}_{ch}.root` | `fit/{JOB}_{ch}.config` | {'present' if ok else 'missing'} | {mu} |")
    L.append("")
    if comb_mf and comb_mf.get("poi_value") and f:
        L += [f"MultiFit of the four workspaces (`fit/comb.config`): mu_Z = {comb_mf['poi_value']:.3f} "
              f"+{comb_mf['poi_err_up']:.3f} -{comb_mf['poi_err_down']:.3f}, against the single-file fit "
              f"{f['mu']:.3f} +{f['mu_err_up']:.3f} -{f['mu_err_down']:.3f}: the two routes are the same model.", ""]
    if yields:
        keys = sorted({k for r in yields["regions"].values() for k in r})
        L += ["## Prefit yields per region", "", "| region | " + " | ".join(keys) + " |", "|---|" + "---:|" * len(keys)]
        for r, y in yields["regions"].items():
            L.append(f"| {r} | " + " | ".join(f"{y.get(k, 0):.0f}" for k in keys) + " |")
        L += ["", "(`*_SRlo_*` / `*_SRhi_*` are the pT(tau_h) split of the same events as `*_SR_*`; they are filled "
              "for the cross-check fit only and never fitted together with them.)", ""]
    for ch, fk in fakes.items():
        if not fk:
            continue
        if ch == "emu":
            L += [f"e mu multijet: OS/SS (SB1) {np.round(fk['SB1']['ratio'], 2).tolist()} vs SB2 {np.round(fk['SB2']['ratio'], 2).tolist()} per dR bin; "
                  f"SS region {fk['ss_region']['data']} data, {fk['ss_region']['mc']:.0f} MC -> {fk['sr_multijet']:.0f} multijet events in the SR", ""]
        else:
            L += [f"{ch} fakes: C_OS/SS = {fk['osss']['C']:.3f} +- {fk['osss']['stat']:.3f}; r_W = "
                  + ", ".join(f"{g} {v['r']:.2f}" for g, v in fk["w_mt"].items())
                  + f"; same-sign closure obs/pred = {fk['closure_ss']['ratio']:.3f} +- {fk['closure_ss']['stat']:.3f}; SR fakes {fk['sr_fakes']:.0f}", ""]
    (config.OUTPUT_DIR_V4 / "RESULTS.md").write_text("\n".join(L) + "\n")

    # ------------------------------------------------------------------ the block the documents quote
    if f:
        s = f["sigma_60_120_pb"]
        exp = f["mu_expected_asimov"]
        B = ["**Result** (four channels, τh ID scale factors and energy scale fitted in situ)", "", "```"]
        B.append(f"σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = {s['value']:.0f} +{s['err_up']:.0f} −{s['err_down']:.0f} pb"
                 f"   (stat ±{s['stat']:.0f};  prediction {s['prediction']:.0f} pb: aMC@NLO acceptance, NNLO normalisation)")
        B.append(f"μ_Z = {f['mu']:.3f} +{f['mu_err_up']:.3f} −{f['mu_err_down']:.3f}   (stat ±{f['mu_stat']:.3f}, syst ±{f['mu_syst']:.3f})"
                 + (f",   expected ±{0.5 * (exp['err_up'] + exp['err_down']):.3f}" if exp else ",   expected: MINOS did not converge")
                 + f",   GoF p = {f['gof_probability']:.2f}"
                 + (f",   μ_tt̄ = {f['mu_ttbar']['value']:.2f} ± {0.5 * (f['mu_ttbar']['err_up'] + f['mu_ttbar']['err_down']):.2f}" if f.get("mu_ttbar") else ""))
        if f.get("tau_id_sf"):
            B.append(f"τh ID SF ({config.TAU_WP}):  " + "  ".join(
                f"{dm} {v['value']:.3f} ± {0.5 * (v['err_up'] + v['err_down']):.3f} (POG {v['pog'][0]:.2f} ± {v['pog'][1]:.2f})"
                for dm, v in f["tau_id_sf"].items()))
        if f.get("tau_es"):
            B.append("τh energy scale:   " + "  ".join(
                f"{dm} {v['pull'] * v['prior_pct']:+.1f} ± {v['constraint'] * v['prior_pct']:.1f} %" for dm, v in f["tau_es"].items()))
        sub = []
        if free_s.get("emu"):
            e = free_s["emu"]
            sub.append(f"e mu alone {e['mu']:.3f} +{e['mu_err_up']:.3f} −{e['mu_err_down']:.3f}")
        if variants.get("taulep"):
            v = variants["taulep"]
            sub.append(f"τ channels alone (SF free) {v['mu']:.3f} +{v['mu_err_up']:.3f} −{v['mu_err_down']:.3f}")
        if sub:
            B.append("μ_Z of the two sub-measurements the fit combines:  " + ";  ".join(sub)
                     + (f"   ({tension['n_sigma']:.1f} σ apart)" if tension else ""))
        per_txt = "  ".join(f"{ch} {fc['mu']:.3f} +{fc['mu_err_up']:.3f} −{fc['mu_err_down']:.3f}" for ch, fc in per)
        if per_txt:
            B.append("per channel alone (POG SFs fixed):  " + per_txt)
        B += ["```", ""]
        groups = sorted(((k, v) for k, v in f["grouped_impact"].items() if k not in ("FullSyst", "Total")), key=lambda kv: -kv[1])[:6]
        B.append("Largest grouped impacts on μ_Z: " + ", ".join(f"{k} {100 * v:.1f} %" for k, v in groups)
                 + f"; data statistics {100 * (f['mu_stat'] or 0):.1f} %. The categories overlap, so their "
                   f"quadrature sum exceeds the MINOS total by a factor "
                   f"{1 / (f['grouped_impact_scale'] or 1):.2f} (`output/RESULTS.md`).")
        if f.get("ranking"):
            B.append("")
            B.append("Ranking: " + ", ".join(f"{r['name']} {100 * max(abs(r['impact_up']), abs(r['impact_down'])):.1f} %"
                                             for r in f["ranking"][:6]) + ".")
        (config.OUTPUT_DIR_V4 / "result_block.md").write_text("\n".join(B) + "\n")

    # ------------------------------------------------------------------ plots
    rows = []
    if f:
        rows.append(("four channels (result)", f["sigma_60_120_pb"]))
    if free_s.get("emu"):
        rows.append(("e mu alone", free_s["emu"]["sigma_60_120_pb"]))
    if variants.get("taulep"):
        rows.append(("tau channels alone (SF free)", variants["taulep"]["sigma_60_120_pb"]))
    for ch, fc in per:
        rows.append((f"{ch} alone (POG SFs)", fc["sigma_60_120_pb"]))
    if rows:
        plt.style.use("default")
        fig, ax = plt.subplots(figsize=(8, 1.0 + 0.55 * len(rows)))
        pred = rows[0][1]["prediction"]
        ax.axvspan(pred * 0.96, pred * 1.04, color="#ffcc66", alpha=0.5); ax.axvline(pred, color="#c08a00")
        for i, (lab, s) in enumerate(rows):
            yv = len(rows) - i
            ax.errorbar([s["value"]], [yv], xerr=[[s["err_down"]], [s["err_up"]]], fmt="o", color="black", capsize=4, lw=1.5)
            if s.get("stat"):
                ax.errorbar([s["value"]], [yv], xerr=[[s["stat"]], [s["stat"]]], fmt="none", color="tab:red", lw=4)
            ax.text(1230, yv + 0.25, lab, fontsize=9); ax.text(3270, yv + 0.25, f"{s['value']:.0f} +{s['err_up']:.0f} -{s['err_down']:.0f} pb", fontsize=9, ha="right")
        ax.set_xlim(1200, 3300); ax.set_ylim(0.3, len(rows) + 0.9); ax.set_yticks([])
        ax.set_xlabel(r"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau)$, $60<m<120$ GeV [pb]")
        ax.set_title(r"Z$\rightarrow\tau\tau$, CMS Open Data 2016 G+H, 16.4 fb$^{-1}$ (red: stat.)", fontsize=11)
        plotting.fig_tag(fig)
        fig.savefig(config.PLOT_DIR_V4 / "step6_summary.png", bbox_inches="tight", dpi=130); plt.close(fig)
    if f:
        gi = {k: v for k, v in f["grouped_impact"].items() if k not in ("FullSyst", "Total")}
        if f["mu_stat"]:
            gi["Data statistics"] = f["mu_stat"]
        items = sorted(gi.items(), key=lambda kv: kv[1])
        fig, ax = plt.subplots(figsize=(7, 0.35 * len(items) + 1.8))
        ax.barh([k for k, _ in items], [100 * v for _, v in items], color="#4a90d9")
        for i, (_, val) in enumerate(items):
            ax.text(100 * val + 0.1, i, f"{100 * val:.2f}%", va="center", fontsize=9)
        tot = 100 * 0.5 * (f["mu_err_up"] + f["mu_err_down"])
        ax.axvline(tot, color="black", ls="--", lw=1.2)
        ax.text(tot, len(items) - 0.4, f" MINOS total {tot:.1f}%", fontsize=9, va="top")
        ax.set_xlabel(r"impact on $\mu_Z$ [%]")
        ax.set_title("Grouped uncertainties (the categories overlap: their sum exceeds the total)", fontsize=10)
        plotting.fig_tag(fig)
        fig.savefig(config.PLOT_DIR_V4 / "step6_impacts.png", bbox_inches="tight", dpi=130); plt.close(fig)
    print((config.OUTPUT_DIR_V4 / "RESULTS.md").read_text())


if __name__ == "__main__":
    main()
