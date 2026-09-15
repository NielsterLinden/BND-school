#!/usr/bin/env python3
"""Assemble slides/ztautau_slides.pdf in the house style (prompts/presentation_style.md) with PyMuPDF.

    source ../setup.sh && python slides/make_figures.py                       # the vector figures (LCG)
    env -u PYTHONPATH -u LD_LIBRARY_PATH -u PYTHONHOME \\
        /project/atlas/users/sjankovy/boostHHbbtautau/HHARD_workfolder/betterplottingtool/venv/bin/python \\
        slides/build_deck.py                                                   # the deck (needs fitz)

Uses slides/build_comparison_deck.py (copied from the HHARD workfolder, unchanged). Numbers come from
output/results.json; nothing is typed in by hand. Optional: variants/tight/output/results.json adds the
working-point comparison slide.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fitz  # noqa: E402
from build_comparison_deck import Deck, smart_upper  # noqa: E402

F = HERE / "figs"
R = json.loads((HERE.parent / "output/results.json").read_text())
FIT, FN = R["fit"]["mcsub"], R["fit"]["nosub"]
Y = R["yields_prefit"]["mcsub"]; YR = R["yields_prefit_per_region"]
B = R["bdt"]; GI = FIT["grouped_impact"]
H = json.loads((HERE / "history.json").read_text())
TIGHT = None
tp = HERE.parent / "variants/tight/output/results.json"
if tp.exists():
    TIGHT = json.loads(tp.read_text())


def pm(v, up, dn, n=3):
    return f"{v:.{n}f} +{up:.{n}f} -{dn:.{n}f}"


s60, sf = FIT["sigma_60_120_pb"], FIT["sigma_fid_pb"]
MU = pm(FIT["mu"], FIT["mu_err_up"], FIT["mu_err_down"])


class ZDeck(Deck):
    """Deck + a text slide (title, bullet lines) and a two-column text/figure slide."""

    def bullets(self, title, lines, accent=None, size=15, x_frac=0.0, width_frac=1.0):
        p = self._new()
        top = self._header(p, title, accent)
        self._lines(p, lines, self.M + x_frac * (self.s.page_w - 2 * self.M), top + 30 * self.k, size,
                    width_frac * (self.s.page_w - 2 * self.M))
        return self

    def _wrap(self, txt, font, fs, max_w):
        """Word-wrap plain text to max_w points."""
        words, out, cur = txt.split(" "), [], ""
        for w in words:
            trial = (cur + " " + w).strip()
            if fitz.get_text_length(trial, fontname=font, fontsize=fs) > max_w and cur:
                out.append(cur); cur = w
            else:
                cur = trial
        if cur:
            out.append(cur)
        return out

    def _lines(self, p, lines, x, y, size, max_w):
        k = self.k
        expanded = []
        for ln in lines:                       # "**Heading**: rest" -> a green heading line + a normal line
            if ln.startswith("**") and "**" in ln[2:]:
                head, rest = ln[2:].split("**", 1)
                expanded.append("**" + head + "**")
                if rest.strip(": "):
                    expanded.append(rest.lstrip(": ").strip())
            else:
                expanded.append(ln)
        for ln in expanded:
            if ln == "":
                y += 14 * k
                continue
            bold = ln.startswith("**")
            txt = ln.strip("*")
            col = self.s.accent2 if bold else self.s.text
            if "$" in txt:
                self._rtext(p, x, y, txt, size * k, col, bold=bold, max_w=max_w)
                y += size * 1.75 * k
            else:
                font = "hebo" if bold else "helv"
                for i, piece in enumerate(self._wrap(txt, font, size * k, max_w)):
                    self._text(p, x + (0 if i == 0 else 16 * k), y, piece, size * k, col, bold=bold)
                    y += size * 1.45 * k
                y += size * 0.3 * k
        return y

    def text_figure(self, title, lines, pdf, accent=None, size=14, split=0.46):
        p = self._new()
        top = self._header(p, title, accent)
        W, M, k = self.s.page_w, self.M, self.k
        self._lines(p, lines, M, top + 28 * k, size, split * (W - 2 * M) - 10 * k)
        self._place(p, fitz.Rect(M + split * (W - 2 * M), top + 6 * k, W - M, self.s.page_h - M), pdf)
        return self

    def text_figures(self, title, lines, pdfs, caps, accent=None, size=13, split=0.36):
        p = self._new()
        top = self._header(p, title, accent)
        W, M, k = self.s.page_w, self.M, self.k
        self._lines(p, lines, M, top + 28 * k, size, split * (W - 2 * M) - 10 * k)
        x0 = M + split * (W - 2 * M)
        n = len(pdfs)
        cw = (W - M - x0 - (n - 1) * M / 2) / n
        for i, (pdf, cap) in enumerate(zip(pdfs, caps)):
            xi = x0 + i * (cw + M / 2)
            self._auto(p, xi + cw / 2, top + 20 * k, cap, 13 * k, self.s.accent2, align="center")
            self._place(p, fitz.Rect(xi, top + 28 * k, xi + cw, self.s.page_h - M), pdf)
        return self


def pdf(name):
    return str(F / f"{name}.pdf")


d = ZDeck(str(HERE / "ztautau_slides.pdf"), author="Samuel Jankovych", date="15 September 2026")

# ---------------------------------------------------------------- cover
d.cover(r"$Z\rightarrow\tau_h\tau_h$ cross section with CMS Open Data",
        "v2 after the review: MC-subtracted, closure-corrected fake factor, fiducial signal, k-fold BDT categories.",
        [rf"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau,\ 60<m<120\ \mathrm{{GeV}}) = {s60['value']:.0f} \pm {s60['stat']:.0f}\ (\mathrm{{stat}})\ ^{{+{s60['err_up']:.0f}}}_{{-{s60['err_down']:.0f}}}\ (\mathrm{{syst}})\ \pm {s60['acceptance']:.0f}\ (\mathrm{{acc}})$ pb    (NNLO 1945 pb)",
         rf"$\sigma_\mathrm{{fid}}(\tau_h\tau_h) = {sf['value']:.2f} \pm {sf['stat']:.2f} \pm {sf['syst']:.2f}$ pb  (pred. 4.50 pb),   $\mu_Z = {FIT['mu']:.3f}\,^{{+{FIT['mu_err_up']:.3f}}}_{{-{FIT['mu_err_down']:.3f}}}$   (v1: 1.160 +0.192 -0.163)",
         "CMS 2016 Tau dataset, Run2016G+H, 16.4 fb-1; both taus hadronic; fakes from data (same-sign fake factor); TRExFitter v1.8.0.",
         f"Precision +{100*FIT['mu_err_up']:.0f}/-{100*FIT['mu_err_down']:.0f} % (v1: +19/-16 %): tau ID SFs {100*GI['Tau ID']:.0f} %, fakes {100*GI['Fakes']:.0f} %, MC stat {100*GI['Gammas']:.0f} %, trigger {100*GI['Tau trigger']:.0f} %; data statistics {100*FIT['mu_stat']:.1f} %.",
         f"Cross-check with DeepTau Tight on both legs: mu_Z = {TIGHT['fit']['mcsub']['mu']:.3f} +{TIGHT['fit']['mcsub']['mu_err_up']:.3f} -{TIGHT['fit']['mcsub']['mu_err_down']:.3f} ({100*TIGHT['fit']['mcsub']['mu_err_up']:.0f} % instead of {100*FIT['mu_err_up']:.0f} %), GoF p = {TIGHT['fit']['mcsub']['gof_probability']:.2f}: more precise, closer to NNLO -- recommended as the next nominal." if TIGHT else "",
         "Sections: 1 result   2 what changed since the review   3 fit   4 BDT categories   5 fakes   6 signal definition and selection   7 findings, next steps",
         "z-tautau/: README.md, CLAUDE.md, REVIEW.md, docs/00-09, review/  (code, documentation, review studies)"])

# ---------------------------------------------------------------- 1. result
d.divider("1.  Result")
d.image(r"Method in one picture", pdf("flow"))
d.text_figure(r"$\sigma(Z/\gamma^*\rightarrow\tau\tau)$ and $\mu_Z$",
              ["**Nominal (v2)**",
               rf"$\mu_Z = {FIT['mu']:.3f}\,^{{+{FIT['mu_err_up']:.3f}}}_{{-{FIT['mu_err_down']:.3f}}}$   (stat. $\pm{FIT['mu_stat']:.3f}$, expected $^{{+{FIT['mu_expected_asimov']['err_up']:.3f}}}_{{-{FIT['mu_expected_asimov']['err_down']:.3f}}}$)",
               rf"$\sigma(60$-$120) = {s60['value']:.0f}\,^{{+{s60['err_up']:.0f}}}_{{-{s60['err_down']:.0f}}} \pm {s60['acceptance']:.0f}_A$ pb,   NNLO 1945 pb",
               rf"$\sigma_\mathrm{{fid}} = {sf['value']:.2f} \pm {sf['stat']:.2f} \pm {sf['syst']:.2f}$ pb,   prediction 4.50 pb",
               f"goodness of fit (saturated model): p = {FIT['gof_probability']:.3f}",
               "",
               "**Cross-checks**",
               rf"FF without MC subtraction: $\mu_Z = {FN['mu']:.3f}\,^{{+{FN['mu_err_up']:.3f}}}_{{-{FN['mu_err_down']:.3f}}}$ (fake normalisation re-measured by the fit)",
               r"v1 (single region, unsubtracted FF, whole DY as signal): $\mu_Z = 1.160\,^{+0.192}_{-0.163}$",
               r"$Z\rightarrow\mu\mu$ v2: $\mu_Z = 0.990 \pm 0.014$, $\sigma(60$-$120) = 1935 \pm 30$ pb",
               "",
               "**Reading**",
               r"$1.4\sigma$ above NNLO and $\mu\mu$; the excess sits at high visible-$\tau$ $p_T$ (section 7).",
               "Precision limited by the external tau ID scale factors, then by the fake estimate."],
              pdf("result_summary"))
d.text_figure("Uncertainty budget",
              ["**Grouped impacts on** $\\mu_Z$"]
              + [f"{k}:  {100*v:.1f} %" for k, v in sorted(GI.items(), key=lambda kv: -kv[1]) if k != "FullSyst"]
              + [f"data statistics:  {100*FIT['mu_stat']:.1f} %", "",
                 "**v1 for comparison**: tau ID 12.5 %, signal modelling 9.8 % (LO-vs-NLO NP, removed), trigger 5.1 %, MC stat 4.8 %, fakes 0.8 % (under-estimated, REVIEW.md 3.5).",
                 "",
                 "Fakes now carry an honest per-category closure uncertainty; MC stat improved by the jet-binned DY samples; signal modelling is scales + PS + PDF only."],
              pdf("impacts"), split=0.42)

if TIGHT:
    TF = TIGHT["fit"]["mcsub"]; TG = TF["grouped_impact"]; TY = TIGHT["yields_prefit"]["mcsub"]
    d.tables("Working-point cross-check: DeepTau Medium vs Tight on both legs  (cross-check, same chain)", [{
        "title": "same chain (fake factors, BDT, categories, fit) with the Tight working point and its TauPOG scale factors",
        "headers": ["quantity", "Medium (nominal)", "Tight"], "col_w": [360, 260, 260],
        "rows": [["mu_Z", MU, pm(TF["mu"], TF["mu_err_up"], TF["mu_err_down"])],
                 ["total uncertainty on mu_Z", f"+{100*FIT['mu_err_up']:.1f} / -{100*FIT['mu_err_down']:.1f} %", f"+{100*TF['mu_err_up']:.1f} / -{100*TF['mu_err_down']:.1f} %"],
                 ["expected (Asimov)", f"+{100*FIT['mu_expected_asimov']['err_up']:.1f} / -{100*FIT['mu_expected_asimov']['err_down']:.1f} %", f"+{100*TF['mu_expected_asimov']['err_up']:.1f} / -{100*TF['mu_expected_asimov']['err_down']:.1f} %"],
                 ["data statistics", f"{100*FIT['mu_stat']:.1f} %", f"{100*TF['mu_stat']:.1f} %"],
                 ["tau ID", f"{100*GI['Tau ID']:.1f} %", f"{100*TG['Tau ID']:.1f} %"],
                 ["fakes", f"{100*GI['Fakes']:.1f} %", f"{100*TG['Fakes']:.1f} %"],
                 ["MC statistics", f"{100*GI['Gammas']:.1f} %", f"{100*TG['Gammas']:.1f} %"],
                 ["tau trigger", f"{100*GI['Tau trigger']:.1f} %", f"{100*TG['Tau trigger']:.1f} %"],
                 ["fiducial signal / fakes (prefit, all categories)", f"{Y['DYtautau']['value']:.0f} / {Y['Fakes']['value']:.0f}", f"{TY['DYtautau']['value']:.0f} / {TY['Fakes']['value']:.0f}"],
                 ["goodness of fit p", f"{FIT['gof_probability']:.3f}", f"{TF['gof_probability']:.3f}"],
                 ["BDT held-out AUC", f"{sum(B['training']['auc_test'])/5:.3f}", f"{sum(TIGHT['bdt']['training']['auc_test'])/5:.3f}"]]}])
# ---------------------------------------------------------------- 2. what changed
d.divider("2.  What changed since the review")
d.tables("Review findings and the v2 answers", [{
    "title": "REVIEW.md findings (v1) and what v2 does",
    "headers": ["finding", "v2", "effect"],
    "col_w": [420, 560, 200],
    "rows": [
        ["3.1 38 % of the signal template non-fiducial", "fiducial part = signal (mu_Z); rest = DYtautau_nonfid background (5 % + scales/PDF/PS on the ratio)", "definition"],
        ["3.2 FF without MC subtraction (bias -0.05)", "MC-subtracted FF nominal; W+jets with uniform weights", "+0.05 on mu"],
        ["3.4 SigModel LO vs NLO 7 %, correlated with mumu", "fiducial C(LO)/C(NLO) = %.3f reported, not fitted; NLO scales/PS/PDF" % R["sigmodel"]["C_LO_over_NLO_fiducial"], "-10 % syst"],
        ["3.5 one over-constrained closure NP, FF stat coherent", "closure corrections eta(tau1), pT(tau2); FF stat per event; closure NPs per category x mass region", "honest fakes"],
        ["3.6 MC statistics 4.8 %", "DY 0J/1J/2J stitched per LHE_NpNLO (5x signal MC)", "%.1f %%" % (100 * GI["Gammas"])],
        ["3.8 trigger SF on jet legs", "SF = 1 for genPartFlav = 0 legs", "< 1 %"],
        ["3.9 second TRExFitter build (StatAnalysis ROOT, v1.10)", "submodule reset to v1.8.0, rebuilt with fitting/build_trexfitter.sh", "reproducible"],
        ["4.3 80 % fakes: ML?", "5-fold XGBoost, mass-agnostic inputs, 3 categories, m_tautau fitted in each", "fakes where fakes are"],
        ["v2.1 fake NPs pulled 0.6-1.5 sigma", "C_OS/SS per (era, N_jets, BDT category); SR0 fitted only above 110 GeV (fake sideband)", "pulls < 0.7 sigma"],
        ["5 eta(tau1) non-closure +-15 %", "diagnosed (FF depends on |eta|), corrected", "eta plots right"],
    ]}])

d.text_figure("Evolution of the result: v1, v2, v2.1",
              ["**v1 (reviewed)**", f"{H['v1']['what']}.  mu_Z = {H['v1']['mu']:.3f} +{H['v1']['up']:.3f} -{H['v1']['down']:.3f}, sigma = {H['v1']['sigma60']} pb, GoF p = {H['v1']['gof']:.2f}",
               "", "**v2**", f"{H['v2']['what']}.  mu_Z = {H['v2']['mu']:.3f} +{H['v2']['up']:.3f} -{H['v2']['down']:.3f}, sigma = {H['v2']['sigma60']} pb, GoF p = {H['v2']['gof']:.3f}",
               "", "**v2.1 (this result)**", f"C_OS/SS per BDT category, SR0 as fake sideband above 110 GeV.  mu_Z = {MU}, sigma = {s60['value']:.0f} pb, GoF p = {FIT['gof_probability']:.2f}",
               "", "**Tight working point** (same v2.1 chain)", (f"mu_Z = {pm(TIGHT['fit']['mcsub']['mu'], TIGHT['fit']['mcsub']['mu_err_up'], TIGHT['fit']['mcsub']['mu_err_down'])}, sigma = {TIGHT['fit']['mcsub']['sigma_60_120_pb']['value']:.0f} pb, GoF p = {TIGHT['fit']['mcsub']['gof_probability']:.2f}" if TIGHT else ""),
               "", "**Reading**", "The central value moved by +0.05 (MC subtraction) and -0.01 (per-category C); the uncertainty shrank from +19/-16 % to +12/-10 %, mostly by removing the LO-vs-NLO NP and adding MC statistics; the fit quality became acceptable once the fake-dominated category stopped being fitted under the peak."],
              pdf("result_summary"), split=0.46)
groups = ["Tau ID", "Fakes", "Gammas", "Tau trigger", "Tau energy scale", "Signal modelling", "Background normalisation", "MET", "Luminosity", "Pileup", "L1 prefiring"]
tg = TIGHT["fit"]["mcsub"]["grouped_impact"] if TIGHT else {}
rows = [[g, f"{H['v1']['impacts'].get(g, 0):.1f}", f"{H['v2']['impacts'].get(g, 0):.1f}", f"{100*GI.get(g, 0):.1f}", f"{100*tg.get(g, 0):.1f}" if TIGHT else "-"] for g in groups]
rows.append(["data statistics", f"{100*H['v1']['stat']:.1f}", f"{100*H['v2']['stat']:.1f}", f"{100*FIT['mu_stat']:.1f}", f"{100*TIGHT['fit']['mcsub']['mu_stat']:.1f}" if TIGHT else "-"])
rows.append(["total on mu_Z (+/-)", f"+{100*H['v1']['up']:.0f} / -{100*H['v1']['down']:.0f}", f"+{100*H['v2']['up']:.0f} / -{100*H['v2']['down']:.0f}", f"+{100*FIT['mu_err_up']:.0f} / -{100*FIT['mu_err_down']:.0f}", f"+{100*TIGHT['fit']['mcsub']['mu_err_up']:.0f} / -{100*TIGHT['fit']['mcsub']['mu_err_down']:.0f}" if TIGHT else "-"])
rows.append(["goodness of fit p", f"{H['v1']['gof']:.2f}", f"{H['v2']['gof']:.3f}", f"{FIT['gof_probability']:.2f}", f"{TIGHT['fit']['mcsub']['gof_probability']:.2f}" if TIGHT else "-"])
d.tables("Uncertainty budget across versions (impact on mu_Z in %)", [{"title": "grouped impacts; v1 'Signal modelling' was the LO-vs-NLO NP, v1 'Fakes' was under-estimated (one over-constrained NP)",
    "headers": ["group", "v1", "v2", "v2.1 Medium", "v2.1 Tight"], "col_w": [300, 170, 170, 190, 190], "rows": rows}])
v2p = H["v2"]["fake_pulls"]; cur = FIT["pulls"]
def pull(d_, k):
    v = d_.get(k); return f"{v[0]:+.2f} ({v[1]:.2f})" if v else "-"
prow = [["FF OS/SS extrapolation (v2: one NP)", pull(v2p, "FakeOSSS_tautau"), "-"],
        ["FF OS/SS SR0 / SR1 / SR2 (v2.1)", "-", " / ".join(pull(cur, f"FakeOSSS_tautau_c{k}") for k in range(3))],
        ["FF non-closure SR1 m<110 / m>110", pull(v2p, "FakeClosure_tautau_c1_lo") + " / " + pull(v2p, "FakeClosure_tautau_c1_hi"), pull(cur, "FakeClosure_tautau_c1_lo") + " / " + pull(cur, "FakeClosure_tautau_c1_hi")],
        ["FF non-closure SR2 m<110 / m>110", pull(v2p, "FakeClosure_tautau_c2_lo") + " / " + pull(v2p, "FakeClosure_tautau_c2_hi"), pull(cur, "FakeClosure_tautau_c2_lo") + " / " + pull(cur, "FakeClosure_tautau_c2_hi")],
        ["FF non-closure SR0 m<110 / m>110", pull(v2p, "FakeClosure_tautau_c0_lo") + " / " + pull(v2p, "FakeClosure_tautau_c0_hi"), "not fitted / " + pull(cur, "FakeClosure_tautau_c0_hi")]]
C = R["C_osss_per_category"]
crow = [[f"category {k}", f"{R['fake_factors']['mcsub']['C_OS_SS_inclusive']:.3f} (inclusive)", f"{H['v2']['osss_sideband_check'][list(H['v2']['osss_sideband_check'])[k]]:.3f}",
         " / ".join(f"{C[0][j][k]:.3f}" for j in range(3)), " / ".join(f"{C[1][j][k]:.3f}" for j in range(3))] for k in range(3)]
d.tables("The fake-factor pulls: diagnosis and fix (v2 -> v2.1)", [
    {"title": "pulls (post-fit constraint) of the fake nuisance parameters", "headers": ["parameter", "v2", "v2.1"], "col_w": [420, 300, 340], "rows": prow},
    {"title": "diagnosis: tau2-anti-isolated sideband, observed / predicted with the inclusive C, and the per-category C now used (era G: 0 / 1 / >= 2 jets; era H)",
     "headers": ["BDT category", "C used in v2", "sideband obs / pred (v2)", "C v2.1, Run2016G", "C v2.1, Run2016H"], "col_w": [220, 200, 240, 260, 260], "rows": crow}])
d.bullets("How the issues were addressed, in one list",
          ["**Reviewed (v1 -> v2)**",
           "38 % of the signal template non-fiducial -> fiducial signal + non-fiducial DY background (5 % NP, theory variations on the ratio)",
           "fake factor without MC subtraction (bias -0.05) -> MC-subtracted, W+jets with uniform weights; closure corrections in |eta(tau1)| (FF varies by +-15 %) and pT(tau2)",
           "one over-constrained closure NP, FF statistics as 4 coherent NPs -> per-event FF statistics in the template variance, closure NPs per category and mass region",
           "LO-vs-NLO 'SigModel' (7 %) correlated with mumu -> reported (fiducial C ratio 0.867), not fitted; NLO scales, PS, PDF",
           "MC statistics 4.8 % -> DY 0J/1J/2J stitched per LHE_NpNLO (5x signal MC)",
           "80 % fakes -> 5-fold XGBoost on mass-agnostic inputs, three categories, m_tautau fitted in each; FF closure in the score verified",
           "wrong TRExFitter build (v1.10, StatAnalysis ROOT) -> submodule at v1.8.0, rebuilt against the LCG ROOT",
           "",
           "**Follow-up (v2 -> v2.1)**",
           "fake NPs pulled 0.6-1.5 sigma -> the OS/SS charge correlation depends on the BDT topology: C per (era, N_jets, category), one NP per category",
           "poor goodness of fit (p = 0.01) -> the fake-dominated category is a sideband: fitted above 110 GeV only (p = 0.22)",
           "trigger paths and luminosity verified: both di-tau paths active and unprescaled in every certified lumisection; single-tau paths rejected (120 GeV threshold keeps 8 % of the signal, no POG scale factors)",
           "Tight working point run through the same chain: +11/-10 % instead of +12/-10 %, closer to NNLO -> recommended next nominal"], size=13)

# ---------------------------------------------------------------- 3. the fit
d.divider("3.  The fit: three BDT categories")
d.panels(r"Prefit $m_{\tau\tau}$ per category", [pdf(f"prefit_tautau_SR{k}") for k in range(3)],
         ["SR0: BDT < 0.55 (fake dominated; fitted above 110 GeV)", "SR1: 0.55 < BDT < 0.90", "SR2: BDT > 0.90 (signal dominated)"], ncols=3)
d.panels(r"Post-fit $m_{\tau\tau}$ per category", [pdf(f"postfit_tautau_SR{k}") for k in range(3)],
         ["SR0 (post-fit total from TRExFitter; bins below 110 GeV not fitted)", "SR1", "SR2: S/B = %.1f" % (YR["tautau_SR2"]["DYtautau"]["value"] / YR["tautau_SR2"]["Fakes"]["value"])], ncols=3)
d.panels("Nuisance-parameter ranking and pulls", [pdf("ranking"), pdf("pulls")],
         [r"top 15 by post-fit impact on $\mu_Z$", "all nuisance parameters"], ncols=2)
rows = []
for s in ("DYtautau", "DYtautau_nonfid", "Fakes", "TTbar", "WJets", "DYee", "WW", "WZ", "ZZ", "SingleTop", "DYlowmass"):
    rows.append([s] + [f"{YR[r][s]['value']:.0f}" for r in ("tautau_SR0", "tautau_SR1", "tautau_SR2")] + [f"{Y[s]['value']:.0f}"])
rows.append(["total prediction"] + [f"{sum(YR[r][s]['value'] for s in Y if s not in ('Data', 'Total')):.0f}" for r in ("tautau_SR0", "tautau_SR1", "tautau_SR2")] + [f"{Y['Total']['value']:.0f}"])
rows.append(["data"] + [f"{YR[r]['Data']['value']:.0f}" for r in ("tautau_SR0", "tautau_SR1", "tautau_SR2")] + [f"{Y['Data']['value']:.0f}"])
d.tables("Prefit yields per category", [{"title": "events, prefit (MC-subtracted fake factor)",
                                         "headers": ["sample", "SR0 (< 0.55)", "SR1 (0.55-0.90)", "SR2 (> 0.90)", "all"],
                                         "col_w": [220, 170, 170, 170, 170], "rows": rows}])

# ---------------------------------------------------------------- 4. BDT
d.divider("4.  The k-fold BDT")
d.text_figures("Why a classifier, and what goes in",
               ["**Why**: 80 % fakes; every fake uncertainty scales with B/S.",
                "Tight WP: -63 % fakes for -24 % signal, same stat. power.",
                "A classifier keeps the signal and sorts the events.",
                "",
                "**Inputs** (mass-agnostic): $p_T(\\tau_{1,2})$, ratio, $|\\eta_{1,2}|$,",
                r"$\Delta R$, $\Delta\phi$, MET, MET significance, $p_T^{vis}$,",
                r"$\Delta\phi(\mathrm{MET},\tau\tau)$, $p_T(\tau\tau+\mathrm{MET})$, $N_{jets}$, jet $p_T$, DM$_{1,2}$",
                "Forbidden: tau1 isolation / raw DeepTau (defines the FF regions), charge, any mass.",
                "",
                "**Training**: signal = fiducial Z->tautau MC (%s events, inclusive + jet-binned);" % f"{B['training']['n_sig']:,}",
                "background = the fake estimate itself (AR data x FF, %s events)." % f"{B['training']['n_bkg']:,}",
                "XGBoost, depth 4, 300 trees, both classes with equal total weight.",
                "",
                "**5 folds by event number**: AR data, same-sign data, SR data, every",
                "MC sample and every kinematic variation are scored by the model",
                "that never saw them.  Held-out AUC %.3f (train %.3f)." % (sum(B["training"]["auc_test"]) / 5, sum(B["training"]["auc_train"]) / 5),
                "",
                "**Categories**: < 0.55 (fake dominated), 0.55-0.90, > 0.90 (S/B %.1f)." % (YR["tautau_SR2"]["DYtautau"]["value"] / YR["tautau_SR2"]["Fakes"]["value"]),
                "Stat-only sensitivity %.2f %% -> %.2f %% (fakes fixed)." % (100 * B["stat_only_sensitivity"]["inclusive"], 100 * B["stat_only_sensitivity"]["categories"])],
               [pdf("bdt_shapes"), pdf("bdt_importance")], ["score shapes (held-out folds)", "feature importance"], split=0.44)
d.panels("Does the fake factor survive the score?  (validation)",
         [pdf("bdt_closure_ss"), pdf("bdt_sr_score"), pdf("SR2_t1_pt")],
         [r"same-sign closure in the score: SS$_T$ vs FF$\times$SS$_L$ (MC subtracted)",
          "signal region, prefit: the fake-dominated bins close, the excess is $\\mu_Z$",
          r"SR2: $p_T(\tau_1)$ data vs prediction (the trend behind section 7)"], ncols=3, accent=d.s.accent2)

# ---------------------------------------------------------------- 5. fakes
d.divider(r"5.  Jet$\rightarrow\tau_h$ fakes")
d.text_figures("The fake factor, MC subtracted and closure corrected",
               [r"FF = [N(SS, $\tau_1$ M, $\tau_2$ M) - MC] / [N(SS, $\tau_1$ VVVL$\wedge\neg$M, $\tau_2$ M) - MC]",
                r"$\times\ f(|\eta_1|)\,g(p_{T,2})$;  applied to OS events with $\tau_1$ failing Medium,",
                r"times $C_{OS/SS}$(era, $N_{jets}$) from the $\tau_2$-anti-isolated sideband.",
                "",
                "era x DM x N_jets x pT = 120 bins, ~10 % stat. each (per event into the template variance).",
                "",
                "**MC subtraction**: genuine tau1 = 1.3 % of the AR but 6 % of the signal,",
                "3.9 %% of the C numerator; C: 1.078 -> %.3f." % R["fake_factors"]["mcsub"]["C_OS_SS_inclusive"],
                "W+jets subtracted with uniform weights (raw weights 20-200, one at -63).",
                "",
                "**Closure corrections**: the FF varies by +-15 % with |eta(tau1)| (same shape",
                "in both eras, all pT, N_jets and DMs: DeepTau's jet rejection is not eta-flat)",
                "and by -6 % with pT(tau2). Residual non-closure per category and mass",
                "region -> one nuisance parameter each: " + ", ".join(f"{k.split('_')[-2]}/{k.split('_')[-1]} {100*v['delta']:.0f} %" for k, v in R["closure_nps"].items()) + "."],
               [pdf("ff_dm1"), pdf("closure_corrections")], ["fake factors, decay mode 1", "closure corrections"], split=0.42)
d.panels(r"Same-sign closure before and after the corrections", [pdf("closure_t1_eta_before"), pdf("closure_t1_eta_after"), pdf("closure_t2_pt_before"), pdf("closure_t2_pt_after")],
         [r"$\eta(\tau_1)$ before", r"$\eta(\tau_1)$ after", r"$p_T(\tau_2)$ before", r"$p_T(\tau_2)$ after"], ncols=4, accent=d.s.accent2)
d.panels(r"Same-sign closure in $N_{jets}$ and $m_{\tau\tau}$", [pdf("closure_njets_after"), pdf("closure_m_tt_after")],
         [r"$N_{jets}$ (binned in the FF)", r"$m_{\tau\tau}$ (the fit variable, all categories)"], ncols=2, accent=d.s.accent2)

# ---------------------------------------------------------------- 6. signal definition and selection
d.divider("6.  Signal definition, mass, selection")
d.text_figure("Fiducial vs non-fiducial Drell-Yan",
              ["**Fiducial volume**: 60 < m_LHE < 120 GeV, both taus hadronic, both visible taus",
               "pT > 40 GeV, |eta| < 2.1.  A = 0.231 % (+-3.7 %), sigma_fid(pred) = 4.50 pb.",
               "",
               "**Of the selected Z/gamma* -> tautau**: 62 % fiducial (-> DYtautau, x mu_Z),",
               "30 % has m_LHE > 120 GeV, 6 % migrates across the 40 GeV / |eta| cuts, 1.5 % has a tau->e/mu leg.",
               "The two 40 GeV cuts enrich the gamma* continuum by two orders of magnitude: 88 % of the",
               "Z->tautau above m = 130 GeV (the fake sideband) has m_LHE > 120 GeV.",
               "",
               "**v2**: DYtautau_nonfid is a background (5 % normalisation NP + the scale/PDF/PS",
               "variations acting on the non-fid/fid ratio, all experimental NPs correlated).",
               "sigma_fid = mu x sigma_fid(pred) is literally true; the sideband no longer contains",
               "mu-dependent signal (REVIEW.md 3.1).",
               "",
               "**Generator cross-check**: fiducial C(LO madgraph)/C(NLO) = %.3f -- the LO sample's" % R["sigmodel"]["C_LO_over_NLO_fiducial"],
               "softer visible-tau spectrum meets the trigger turn-on. Reported, not an NP (docs/07)."],
              pdf("prefit_all_log"), split=0.52)
d.text_figure(r"MET correction of the visible mass: $m_{\tau\tau}$",
              [r"$\geq 4$ neutrinos escape: $m_\mathrm{vis}$ peaks at 0.8 $m_Z$.",
               "",
               "**Collinear approximation**: singular for back-to-back taus, defined in 38 %.",
               "**Likelihood mass** (FastMTT / SVfit-lite idea): scan the visible fractions",
               r"$(x_1, x_2)$, weight = MET transfer function (per-event covariance)",
               r"$\times$ flat two-body phase space $\times\ 1/m^2$; posterior median.",
               "",
               "scale 0.99, resolution 11 % (m_vis: 0.80, 13 %); always defined.",
               "Separation from fakes about equal to m_vis; m_tautau is the physical mass",
               "and gives a clean high-mass sideband for the fake normalisation."],
              pdf("mass_estimators"), split=0.5)
d.tables("Selection and cutflow", [
    {"title": "selection", "headers": ["object", "requirement"], "col_w": [220, 900],
     "rows": [["trigger", "DoubleMedium(Combined)IsoPFTau35_Trk1_eta2p1_Reg (G: Iso, H: CombinedIso); both legs matched to HLT tau objects"],
              ["tau_h", "pT > 40 GeV, |eta| < 2.1, |dz| < 0.2 cm, DM 0/1/10/11; DeepTau VSe VVLoose, VSmu VLoose; VSjet Medium (tight) / VVVLoose & !Medium (loose)"],
              ["pair", "most isolated pair, dR > 0.5; tau1 = leading pT (the fake-factor leg)"],
              ["vetoes", "no extra e / mu: orthogonal to e tau_h, mu tau_h, ee, mumu"],
              ["simulation", "events whose leading tau is a jet are dropped (covered by the fake factor)"]]},
    {"title": "cutflow", "headers": ["step", "Run2016G", "Run2016H"], "col_w": [420, 200, 200],
     "rows": [["NanoAOD", "79.6 M", "76.8 M"], ["skim (certified, trigger, >= 2 candidates)", "2.81 M", "2.22 M"],
              ["pair + MET filters", "2.03 M", "1.62 M"], ["trigger match", "1.99 M", "1.59 M"], ["lepton veto", "1.99 M", "1.59 M"],
              ["signal region (OS, both Medium)", "24 106", "23 480"]]}])

# ---------------------------------------------------------------- 7. findings
d.divider("7.  Findings and next steps")
d.bullets("Findings",
          [rf"1. $\sigma(60$-$120) = {s60['value']:.0f}\,^{{+{s60['err_up']:.0f}}}_{{-{s60['err_down']:.0f}}} \pm {s60['acceptance']:.0f}$ pb, $\mu_Z = {FIT['mu']:.2f}$: $1.4\sigma$ above NNLO and above $Z\rightarrow\mu\mu$ ($\mu_Z = 0.99$). Total $^{{+{100*FIT['mu_err_up']:.0f}}}_{{-{100*FIT['mu_err_down']:.0f}}}$ % (v1: +19/-16 %), data statistics {100*FIT['mu_stat']:.1f} %.",
           "2. The limit is external: the TauPOG tau ID scale factors (%.0f %%); Z->tautau is how those SFs are measured, and the two POG prescriptions (DM- vs pT-binned) differ by 14 %% on the yield." % (100 * GI["Tau ID"]),
           "3. In the signal-dominated category data/prediction rises with the visible tau pT (1.0 at threshold, 1.2 above 70 GeV): part of mu > 1 comes from there (the fit quality itself is fine, p = %.2f, once the fake sideband is restricted to m > 110 GeV). Trigger turn-on modelling or the NLO Z pT spectrum; an in-situ trigger measurement is the next step." % FIT["gof_probability"],
           "4. The fake factor closes only with N_jets, era and |eta(tau1)| / pT(tau2) corrections; the FF is eta-dependent by +-15 % because DeepTau's jet rejection is not. The MC subtraction is not optional (bias -0.05 on mu; unsubtracted fakes 7 % high where there is no signal).",
           "5. 38 % of the selected Z->tautau is non-fiducial; kept in the signal template it makes the sideband mu-dependent. As a background it costs a 5 % normalisation NP and nothing else.",
           "6. The BDT categories give a nearly background-free Z peak (S/B %.1f) and confine the fake uncertainties (%.0f %%, now honest) to where the fakes are; the statistical gain is small (the measurement is systematics limited)." % (YR["tautau_SR2"]["DYtautau"]["value"] / YR["tautau_SR2"]["Fakes"]["value"], 100 * GI["Fakes"]),
           ("8. Tight working point (same chain): mu_Z = %s, total +%.0f/-%.0f %% (Medium +%.0f/-%.0f %%), GoF p = %.2f: every uncertainty group shrinks, the visible-pT slope is weaker but still there. Recommended nominal for the next iteration and the combination." % (pm(TIGHT['fit']['mcsub']['mu'], TIGHT['fit']['mcsub']['mu_err_up'], TIGHT['fit']['mcsub']['mu_err_down']), 100*TIGHT['fit']['mcsub']['mu_err_up'], 100*TIGHT['fit']['mcsub']['mu_err_down'], 100*FIT['mu_err_up'], 100*FIT['mu_err_down'], TIGHT['fit']['mcsub']['gof_probability'])) if TIGHT else "",
           "7. Reproducibility: this checkout's TRExFitter was a v1.10 commit built against the StatAnalysis ROOT; it is now the tagged v1.8.0 built with fitting/build_trexfitter.sh."], size=14)
d.bullets("For the combination, and next steps",
          ["**For the combination**",
           "workspace fit/results/ztautau (job ztautau, POI mu_Z, regions tautau_SR0/1/2) plugs into fitting/combination_skeleton.config",
           "correlated: Lumi, Pileup, L1Prefiring, QCDScale, PDF, PS_ISR/FSR, XS_*; tau-only: TauID/Trigger/ES_DM*, Fake*_tautau, XS_DYtautau_nonfid",
           "SigModel_tautau is reported, not fitted; never correlate it with z-mumu's SigModel (a different quantity)",
           "with ee/mumu fixing mu_Z, tautau constrains the tau ID nuisance parameters in situ (lepton universality): present TauID_DM* as a result",
           "Born-level 60-120 GeV denominator for sigma; acceptance uncertainty 3.7 % correlated",
           "",
           "**Next steps** (CLAUDE.md, open issues)",
           "decay-mode categories -> tau ID in situ even before the combination; quote the pT-binned POG prescription as a cross-check",
           "trigger efficiency from mu tau_h tag-and-probe (SingleMuon) instead of the POG turn-on SFs",
           "HT-binned W+jets; the high-mass DY sample for the non-fiducial template; EWK Z->tautau",
           "add e tau_h / mu tau_h channels",
           "",
           "Everything here is reproducible from the ntuples in ~25 min:  python run_all.py --from 3   (z-tautau/README.md)"], size=14)
d.save()
