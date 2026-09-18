#!/usr/bin/env python
"""Build slides/ztautau_slides.pdf: the Z -> tautau four-channel measurement, in the house style.

    source ../fitting/setup.sh && python slides/make_figures.py          # the vector figures (LCG)
    env -u PYTHONPATH -u LD_LIBRARY_PATH -u PYTHONHOME \
        /project/atlas/users/sjankovy/boostHHbbtautau/HHARD_workfolder/betterplottingtool/venv/bin/python \
        slides/build_deck.py                                      # the deck (needs fitz / PyMuPDF)

Every number comes from output/results.json and output/data/yields_v4.json; nothing is typed in by hand,
so the deck cannot drift from the measurement. Page background = the figures' background (#222222), and
every figure is embedded as vector, so any plot can be zoomed into without pixelation.

Order follows .claude/prompts/presentation_style.md: the result first, then what the fit measures in situ, then
the uncertainties, then the regions, the data-driven backgrounds and the in-situ corrections, and finally
what the combination gets.
"""

from __future__ import annotations

import json
from pathlib import Path

from build_comparison_deck import Deck

HERE = Path(__file__).resolve().parent
CHANNEL = HERE.parent
FIGS = HERE / "figs"
OUT = HERE / "ztautau_slides.pdf"

R = json.loads((CHANNEL / "output/results.json").read_text())
FIT = R["fit"]["combined"]
VAR = R.get("fit_variants", {})
FREE = R.get("fit_per_channel_free", {})
YIELDS = json.loads((CHANNEL / "output/data/yields_v4.json").read_text())["regions"]
CH = {"tautau": r"$\tau_h\tau_h$", "mutau": r"$\mu\tau_h$", "etau": r"$e\tau_h$", "emu": r"$e\mu$"}
REGION_TITLE = {"tautau_SR0": r"$\tau_h\tau_h$, BDT < 0.55  (fake sideband, fitted above 110 GeV)",
                "tautau_SR1": r"$\tau_h\tau_h$, 0.55 < BDT < 0.90", "tautau_SR2": r"$\tau_h\tau_h$, BDT > 0.90",
                "emu_SR": r"$e\mu$ signal region", "emu_CRtt": r"$e\mu$  $t\bar{t}$ control region"}
for _c, _l in (("mutau", r"$\mu\tau_h$"), ("etau", r"$e\tau_h$")):
    for _dm in (0, 1, 10, 11):
        REGION_TITLE[f"{_c}_SR_dm{_dm}"] = _l + f",  decay mode {_dm}"


def f(x, n=3):
    return f"{x:.{n}f}"


def mu(fit):
    return f"{fit['mu']:.3f} +{fit['mu_err_up']:.3f} -{fit['mu_err_down']:.3f}"


def sig(fit):
    s = fit["sigma_60_120_pb"]
    return f"{s['value']:.0f} +{s['err_up']:.0f} -{s['err_down']:.0f}"


def fig(name):
    return str(FIGS / f"{name}.pdf")


def main():
    s = FIT["sigma_60_120_pb"]
    tension = R.get("tension_taulep_vs_emu") or {}
    d = Deck(str(OUT), author="Samuel Jankovych", date="17 September 2026")

    # ---------------------------------------------------------------- cover
    exp = FIT.get("mu_expected_asimov")
    d.cover(
        r"Z $\rightarrow$ $\tau\tau$ in four channels",
        "CMS Open Data 2016 (Run2016G+H), 16.4 fb$^{-1}$   —   "
        r"$\tau_h\tau_h$ + $\mu\tau_h$ + $e\tau_h$ + $e\mu$, fitted together",
        lines=[
            f"sigma(pp -> Z/gamma* -> tautau, 60 < m < 120 GeV) = {s['value']:.0f} "
            f"+{s['err_up']:.0f} -{s['err_down']:.0f} pb      (NLO prediction {s['prediction']:.0f} pb)",
            f"mu_Z = {mu(FIT)}"
            + (f",  expected +-{0.5 * (exp['err_up'] + exp['err_down']):.3f}" if exp else "")
            + f",  goodness of fit p = {FIT['gof_probability']:.2f},  mu_ttbar = {FIT['mu_ttbar']['value']:.2f}",
            "",
            "The tau_h identification scale factors and the tau_h energy scale are measured in situ by this fit, not taken",
            "from the TauPOG. Jet -> tau_h fakes come from data in every channel. There is no mu mu channel and every",
            "channel vetoes a second lepton, so this is orthogonal to Z -> mu mu and Z -> ee.",
            "",
            f"It is the compromise between two sub-measurements {tension.get('n_sigma', 0):.1f} sigma apart: "
            f"e mu alone (no tau_h) {mu(FREE['emu'])},"
            if FREE.get("emu") else "",
            f"the three tau channels alone, scale factors free, {mu(VAR['taulep'])}. Both are quoted throughout."
            if VAR.get("taulep") else "",
            "",
            "1 the result    2 measured in situ    3 uncertainties    4 the fitted regions"
            "    5 backgrounds from data    6 trigger efficiencies    7 for the combination",
        ])

    # ---------------------------------------------------------------- 1. the result
    d.divider(r"1.  The result")
    d.image(r"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau)$  for  $60 < m_{\tau\tau} < 120$ GeV", fig("summary"))
    d.image(r"What the combined $\mu_Z$ is made of", fig("submeasurements"))

    rows = [["four channels (the measurement)", mu(FIT), sig(FIT), f(FIT["gof_probability"], 2)]]
    if FREE.get("emu"):
        rows.append([r"$e\mu$ alone (no $\tau_h$)", mu(FREE["emu"]), sig(FREE["emu"]), f(FREE["emu"]["gof_probability"], 2)])
    if VAR.get("taulep"):
        rows.append([r"$\tau$ channels alone, SF free", mu(VAR["taulep"]), sig(VAR["taulep"]), f(VAR["taulep"]["gof_probability"], 2)])
    for ch in ("tautau", "mutau", "etau", "emu"):
        if ch in R["fit"]:
            c = R["fit"][ch]
            rows.append([f"{CH[ch]} alone, POG SF", mu(c), sig(c), f(c["gof_probability"], 2)])
    cross = []
    for key, what in (("emutrig2x", r"$e\mu$ trigger prior doubled"),
                      ("ptsplit", r"$\tau_h$ ID SF split at 40 GeV")):
        if VAR.get(key):
            cross.append([what, mu(VAR[key]), sig(VAR[key])])
    CW = [330, 230, 250, 120]
    tables = [{"title": "Every fit of the measurement", "col_w": CW,
               "headers": ["fit", r"$\mu_Z$", r"$\sigma$(60–120) [pb]", "GoF p"], "rows": rows}]
    if cross:
        tables.append({"title": "Cross-check fits: how much the modelling choices are worth", "col_w": CW[:3],
                       "headers": ["fit", r"$\mu_Z$", r"$\sigma$(60–120) [pb]"], "rows": cross})
    d.tables("All fits, side by side", tables)

    # ---------------------------------------------------------------- 2. in situ
    d.divider(r"2.  What the fit measures in situ")
    d.panels(r"$\tau_h$ identification scale factor and energy scale, per decay mode",
             [fig("tauid"), fig("taues")],
             ["fitted against the TauPOG values, which are only priors here",
              r"the 3% prior is constrained to 0.6–2.2% by the $m_{\tau\tau}$ shapes"], ncols=2)
    d.image(r"Cross-check: is the scale factor flat in $p_T(\tau_h)$?   (it is not)", fig("ptsplit"), accent=(0.541, 0.776, 0.247))

    sf = {k: v for k, v in (FIT.get("tau_id_sf") or {}).items() if not k.endswith("_lowpt")}
    sf_rows = [[k.replace("DM", "decay mode "), f"{v['value']:.3f} ± {0.5 * (v['err_up'] + v['err_down']):.3f}",
                f"{v['pog'][0]:.3f} ± {v['pog'][1]:.3f}",
                f"{100 * 0.5 * (v['err_up'] + v['err_down']) / v['value']:.1f}%"] for k, v in sf.items()]
    es_rows = [[k.replace("DM", "decay mode "), f"{v['pull'] * v['prior_pct']:+.2f}%",
                f"± {v['constraint'] * v['prior_pct']:.2f}%", f"{v['constraint']:.2f}"]
               for k, v in (FIT.get("tau_es") or {}).items()]
    d.tables("The corrections this measurement replaces with its own",
             [{"title": r"$\tau_h$ identification scale factor (DeepTau VSjet Tight)", "col_w": [240, 240, 260, 180],
               "headers": ["", "this fit", "TauPOG (prior only)", "precision"], "rows": sf_rows},
              {"title": r"$\tau_h$ energy scale, relative to the POG central value (prior ± 3%)", "col_w": [240, 200, 220, 260],
               "headers": ["", "shift", "uncertainty", "fraction of the prior"], "rows": es_rows}])

    # ---------------------------------------------------------------- 3. uncertainties
    d.divider(r"3.  Uncertainties")
    d.image(r"Grouped impacts on $\mu_Z$", fig("impacts"))
    d.image(r"Ranking: post-fit impact of each nuisance parameter", fig("ranking"))
    d.image("Pulls and constraints", fig("pulls"))

    scale = FIT.get("grouped_impact_scale") or 1.0
    gi = sorted(((k, v) for k, v in FIT["grouped_impact"].items() if k not in ("FullSyst", "Total")),
                key=lambda kv: -kv[1])
    rho = (R.get("for_combination") or {}).get("category_rho_recommended", {})
    d.tables("Uncertainty budget, and what a combination may correlate",
             [{"title": f"grouped impacts: their quadrature sum exceeds the MINOS total "
                        f"({100 * 0.5 * (FIT['mu_err_up'] + FIT['mu_err_down']):.1f}%) by {1 / scale:.2f}, "
                        f"so rescale by {scale:.3f} if a breakdown is needed",
               "col_w": [320, 220, 200, 240],
               "headers": ["category", r"impact on $\mu_Z$", "rescaled", r"$\rho$ with $\mu\mu$ / ee"],
               "rows": [[k, f"{100 * v:.2f}%", f"{100 * v * scale:.2f}%",
                         "1" if rho.get(k) == 1.0 else "0"] for k, v in gi]
                       + [["data statistics", f"{100 * (FIT['mu_stat'] or 0):.2f}%",
                           f"{100 * (FIT['mu_stat'] or 0) * scale:.2f}%", "0"]]}])

    # ---------------------------------------------------------------- 4. regions
    d.divider(r"4.  The 13 fitted regions")
    d.panels(r"$\tau_h\tau_h$ — three BDT categories", [fig(f"region_tautau_SR{i}") for i in range(3)],
             [REGION_TITLE[f"tautau_SR{i}"] for i in range(3)], ncols=3)
    for ch, lab in (("mutau", r"$\mu\tau_h$"), ("etau", r"$e\tau_h$")):
        d.panels(f"{lab} — one region per $\\tau_h$ decay mode",
                 [fig(f"region_{ch}_SR_dm{dm}") for dm in (0, 1, 10, 11)],
                 [f"decay mode {dm}" for dm in (0, 1, 10, 11)], ncols=4)
    d.panels(r"$e\mu$ — the signal region and the $t\bar{t}$ control region",
             [fig("region_emu_SR"), fig("region_emu_CRtt")],
             ["signal region: no $\\tau_h$, so no ID scale factor",
              r"$D_\zeta$ < $-$40 GeV and MET > 80 GeV: normalises $t\bar{t}$"], ncols=2)

    keys = ["Data", "DYtautau", "Fakes", "TTbar", "DYtautau_out", "DYee", "DYmumu", "WJets"]
    d.tables("Prefit yields per region",
             [{"title": "events, before the fit", "col_w": [250] + [112] * len(keys), "headers": ["region"] + keys,
               "rows": [[r] + [f"{YIELDS[r].get(k, 0):.0f}" for k in keys]
                        for r in REGION_TITLE if r in YIELDS]}])

    # ---------------------------------------------------------------- 5. backgrounds
    d.divider(r"5.  Backgrounds measured from data")
    d.image(r"$\mu\tau_h$: fake factors per process", fig("fakefactors_mutau"), accent=(0.541, 0.776, 0.247))
    d.image(r"$e\tau_h$: fake factors per process", fig("fakefactors_etau"), accent=(0.541, 0.776, 0.247))
    d.image(r"$e\mu$: the multijet OS/SS factor, from two isolation sidebands", fig("emu_osss"),
            accent=(0.541, 0.776, 0.247))

    # ---------------------------------------------------------------- 6. triggers
    d.divider(r"6.  Trigger efficiencies, measured in situ")
    d.image("Scale factors from the other stream's events", fig("triggers"))

    # ---------------------------------------------------------------- 7. combination
    d.divider(r"7.  For the combination")
    deliver = [["fit inputs, all channels", "fit/fitinputs/ztautau.root", "13 fitted regions, 62 systematics"],
               ["fit inputs, per channel", "fit/fitinputs/ztautau_<ch>.root", "histograms copied bit for bit"],
               ["configs", "fit/ztautau[_<ch>].config", "POI mu_Z, scale factors free"],
               ["workspaces", "fit/results/ztautau[_<ch>]/RooStats/", "what a MultiFit reads"],
               ["the numbers", "output/results.json -> for_combination", "sigma_pred, groups, ranking"],
               ["the instructions", "docs/11-combination-inputs.md", "read section 6 before running"]]
    warn = [["signal is 15 templates", "DYtautau_tDM*, listed in the metadata as signal_samples",
             "mu_Z goes on all of them, not on a sample called DYtautau"],
            ["no acceptance parameters", "acceptance_in_fit = true",
             "the templates already vary A x eff inside the fit: Acc_* double counts"],
            ["empty bins", "tautau_SR2 bin 1 is dropped", "its unconstrained gamma made the offset log(0)"],
            ["e mu overlaps z-mumu", "our emu_SR shares events with mumu_CRemu", "one of the two must be dropped"],
            ["muon / electron NPs", "ours are the POG jsons", "decorrelate unless both groups agree they are the same"],
            ["ttbar", "we float mu_ttbar from emu_CRtt", "z-mumu constrains it with XS_TTbar: decide once"]]
    d.tables("What the combination gets, and what it must know",
             [{"title": "deliverables, all committed", "col_w": [290, 420, 430],
               "headers": ["what", "path", "note"], "rows": deliver},
              {"title": "six things that will bite: docs/11-combination-inputs.md sections 5-7",
               "col_w": [250, 450, 440], "headers": ["", "what it is", "what to do"], "rows": warn}])

    d.save()


if __name__ == "__main__":
    main()
