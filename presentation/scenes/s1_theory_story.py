"""Section 1 (theory): why we still care about the Z, told through its cross section.

Brief: ``setup_and_reference/briefs/theory.md``. Every printed number is read from
``presentation/data/theory_reference.json`` (``data/extract_theory_reference.py``) and checked at import.
Five scenes, eight clips (manim sections, delivered by ``tools/deliver_chain.py``); a scene with two clips
is a chain, the second clip opens on the first one's last frame.

    TheoryProcess       sm_to_process        the SM table builds; quarks, Z/gamma and leptons turn into the
                                             Drell-Yan diagram (1-01's exact frame); g_V, g_A at both vertices
                        factorisation        protons around the diagram, x1 P and x2 P; sigma = sum int f f
                                             sigma-hat; the alpha_s series with its years; the x1 x2 = m_Z^2/s
                                             line with the |y| < 2.4 part lit
    TheoryPrediction    prediction           sigma axis; the aMC@NLO line; its band grows scale -> alpha_s -> PDF;
                                             four NNLO+NNLL predictions by PDF set, their 49 pb spread
                        lumi_limited         the CMS 2024 point; its error splits into stat / syst / lumi, the
                                             lumi part lights up above the PDF spread
    TheoryEverywhere    z_everywhere         H -> tautau under the Z -> tautau tail (area ratio 1/560, schematic
                                             shapes); Z -> mumu + jet becomes Z -> nunu + jet in a small slice
    TheoryUniversality  lepton_universality  1-01's frame; the lepton legs split into e, mu, tau; the three sit on a
                                             log mass axis (x207, x17), fly into the (g_A, g_V) plane and land
                                             on one point as sin^2 theta runs to 0.2315; equal BR bars; LEP ratios
    TheoryMeasurement   mass_window          the generator-level m_ll spectrum of our aMC@NLO sample; the 60-120 GeV
                                             window shades; its area counts up to the prediction
                        final_frame          the window becomes the prediction line of the final plot's sigma axis
                                             (1650-2400 pb): band, rows Z->ee, Z->mumu, Z->tautau | Z->ll, empty

Schematic on purpose: the Feynman diagram, the protons, the m_tautau shapes of z_everywhere (Gaussians at the
typical resolution; only their area ratio is a number) and the event of the small slice. Layout keeps the title
band (y > 2.7) and the top-left block (x < -5.85, y > 0.22) empty.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    DOWN, LEFT, ORIGIN, RIGHT, UP, ChangeDecimalToValue, Create, DashedLine, DecimalNumber, Dot, Ellipse,
    FadeIn, FadeOut, FadeTransform, GrowFromCenter, GrowFromEdge, LaggedStart, Line, MathTex, Polygon, Rectangle,
    ReplacementTransform, RoundedRectangle, Scene, Square, SurroundingRectangle, Transform, TransformFromCopy,
    ValueTracker, VGroup, always_redraw, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403
from style.bnd_style import _cap_scale  # noqa: E402
from s1_drell_yan import DrellYan  # noqa: E402

# ---------------------------------------------------------------------------
# frozen numbers + anchors (a mismatch stops the render; never adjusted here)
# ---------------------------------------------------------------------------

TH = load_data("theory_reference")
PRED = TH["prediction"]
PDFS = TH["pdf_sets"]["sets"]
CMS24 = TH["cms_2024"]
XR = TH["x_range"]
LEP = TH["leptons"]
SPEC = TH["spectrum"]
HIG = TH["higgs"]

assert abs(PRED["value"] - 1953.93) < 0.01 and round(PRED["err_up"], 1) == 56.1 and round(PRED["err_down"], 1) == 81.7
assert [s["added"] for s in PRED["band_steps"]] == ["scale", "alphas", "pdf"]
assert [p["label"] for p in PDFS] == ["CT18", "MSHT20", "NNPDF3.1", "NNPDF4.0"]
assert [p["value"] for p in PDFS] == [1921, 1935, 1940, 1970] and TH["pdf_sets"]["spread_pb"] == 49
assert [CMS24[k] for k in ("value", "stat", "syst", "lumi")] == [1952, 4, 18, 45]
assert [o["year"] for o in TH["orders"]] == [1970, 1979, 1991, 2022]
assert abs(XR["x_min"] - 6.36e-4) < 1e-5 and abs(XR["x_max"] - 0.0773) < 1e-4
assert round(LEP["ratio_mu_e"]) == 207 and round(LEP["ratio_tau_mu"]) == 17
assert abs(LEP["br_ll"]["value"] - 3.3658) < 1e-9 and abs(LEP["couplings"]["l"]["g_V"] + 0.03694) < 1e-5
assert (LEP["lep_ratios"]["mumu_over_ee"]["value"], LEP["lep_ratios"]["tautau_over_ee"]["value"]) == (1.0009, 1.0019)
assert abs(SPEC["sigma_60_120_pb"] - PRED["value"]) < 1e-6 and round(SPEC["fraction_60_120"], 3) == 0.965
assert SPEC["edges_GeV"][0] == 50 and SPEC["edges_GeV"][-1] == 200
assert round(HIG["z_over_h_tautau"], -1) == 560

EASE = rate_functions.ease_in_out_sine
TH_DARK = darken(THEORY, 0.30)        # predictions drawn on the purple band
CY_DARK = darken(DETECTOR_ACCENT, 0.25)
LEP_COL = {"e": CHANNEL_LINE["ee"], "mu": CHANNEL_LINE["mumu"], "tau": CHANNEL_LINE["tautau"]}


def tex(expr: str, h: float = 0.26, color=INK) -> "MathTex":
    """MathTex whose digits are ``h`` high."""
    m = mathtex(expr, color=color)
    return m.scale(_cap_scale("tex", h))


def txt(s: str, h: float = 0.22, color=INK):
    return text(s, color=color).scale(_cap_scale("text", h))


def right_at(m, x: float, y: float):
    """Right edge at x, vertical centre at y (left-hand label columns)."""
    m.move_to([x - m.width / 2, y, 0.0])
    return m


def pulse(scene, m, width: float = 14.0, run_time: float = 0.7):
    """Emphasis without touching ``m`` (recipes trap 2): a ghost copy widens and fades."""
    ghost = m.copy()
    scene.add(ghost)
    scene.play(ghost.animate.set_stroke(width=width, opacity=0.0), run_time=run_time)
    scene.remove(ghost)


# ---------------------------------------------------------------------------
# 1-02 / 1-03: the SM table -> the Drell-Yan diagram -> factorisation
# ---------------------------------------------------------------------------

SM_COLS = (-2.3, -1.3, -0.3, 1.0, 2.3)
SM_ROWS = (1.55, 0.55, -0.45, -1.45)
SM_CELLS = [   # (tex, column, row, kind)
    ("u", 0, 0, "quark"), ("c", 1, 0, "quark"), ("t", 2, 0, "quark"),
    ("d", 0, 1, "quark"), ("s", 1, 1, "quark"), ("b", 2, 1, "quark"),
    ("e", 0, 2, "lepton"), (r"\mu", 1, 2, "lepton"), (r"\tau", 2, 2, "lepton"),
    (r"\nu_e", 0, 3, "nu"), (r"\nu_\mu", 1, 3, "nu"), (r"\nu_\tau", 2, 3, "nu"),
    ("g", 3, 0, "gluon"), (r"\gamma", 3, 1, "photon"), ("Z", 3, 2, "Z"), ("W", 3, 3, "W"), ("H", 4, 0, "H"),
]


def sm_cell(sym: str, kind: str) -> VGroup:
    boson = kind in ("gluon", "photon", "Z", "W")
    stroke = THEORY if kind in ("gluon", "photon", "Z", "W", "H") else INK
    fill = lighten(THEORY, 0.88) if boson else (lighten(THEORY, 0.72) if kind == "H" else WHITE)
    box = RoundedRectangle(width=0.86, height=0.86, corner_radius=0.1, stroke_color=col(stroke), stroke_width=2.2,
                           fill_color=col(fill), fill_opacity=1.0)
    lab = tex(sym, 0.30, INK).move_to(box)
    g = VGroup(box, lab)
    g.kind = kind
    return g


class TheoryProcess(DrellYan):
    def construct(self):
        white_background(self)
        clip_open(self, "sm_to_process")

        # the table
        cells = []
        for sym, c, r, kind in SM_CELLS:
            cells.append(sm_cell(sym, kind).move_to([SM_COLS[c], SM_ROWS[r], 0.0]))
        table = VGroup(*cells)
        self.play(LaggedStart(*[FadeIn(c, scale=0.75) for c in cells], lag_ratio=0.07), run_time=1.9)
        self.remove(*self.mobjects)
        self.add(table)
        self.wait(0.3)
        keep = [c for c in cells if c.kind in ("quark", "lepton", "photon", "Z")]
        dim = [c for c in cells if c not in keep]
        self.play(*[c.animate.set_opacity(0.18) for c in dim], run_time=0.7)
        self.wait(0.3)

        # the pieces become the diagram
        P = self.build_parts()
        quarks = VGroup(*[c for c in cells if c.kind == "quark"])
        leptons = VGroup(*[c for c in cells if c.kind == "lepton"])
        bosons = VGroup(*[c for c in cells if c.kind in ("photon", "Z")])
        self.play(FadeTransform(quarks.copy(), P["lab_q"]), FadeTransform(quarks.copy(), P["lab_qb"]),
                  FadeTransform(bosons.copy(), P["lab_b"]),
                  FadeTransform(leptons.copy(), P["lab_lm"]), FadeTransform(leptons.copy(), P["lab_lp"]),
                  FadeOut(table), run_time=1.5, rate_func=EASE)
        self.remove(table)
        self.play(Create(P["q"]), Create(P["qb"]), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(P["q_tip"]), FadeIn(P["qb_tip"]), FadeIn(P["v1"]), run_time=0.3)
        self.play(Create(P["boson"]), run_time=0.8, rate_func=EASE)
        self.play(FadeIn(P["v2"]), run_time=0.2)
        self.play(LaggedStart(Create(P["lm"]), Create(P["lp"]), lag_ratio=0.15), run_time=0.8, rate_func=EASE)
        self.play(FadeIn(P["lm_tip"]), FadeIn(P["lp_tip"]), run_time=0.3)
        self.wait(0.3)

        # the couplings at both vertices
        gq = tex(r"g^{\,q}_{V,A}", 0.24, THEORY).next_to(P["v1"], DOWN, buff=0.38)
        gl = tex(r"g^{\,\ell}_{V,A}", 0.24, THEORY).next_to(P["v2"], DOWN, buff=0.38)
        law = tex(r"g_V^f = T_3^f - 2Q_f\sin^2\theta_W \qquad g_A^f = T_3^f", 0.26, THEORY).move_to([0.0, -2.75, 0.0])
        self.play(FadeIn(gq, shift=UP * 0.15), FadeIn(gl, shift=UP * 0.15), run_time=0.6)
        self.play(FadeIn(law, shift=UP * 0.15), run_time=0.7)

        clip_cut(self, "factorisation")

        # ---- factorisation: the diagram moves up, protons around it
        self.play(FadeOut(gq), FadeOut(gl), FadeOut(law), run_time=0.5)
        diagram = VGroup(*P.values())
        self.play(diagram.animate.scale(0.62, about_point=ORIGIN).shift([-0.4, 1.25, 0.0]), run_time=1.1, rate_func=EASE)
        q_in, qb_in = P["q"].get_start(), P["qb"].get_start()
        protons, links, plabs = VGroup(), VGroup(), VGroup()
        for p_in in (q_in, qb_in):
            c = p_in + LEFT * 1.3
            blob = Ellipse(width=1.0, height=0.6, stroke_color=col(DETECTOR_ACCENT), stroke_width=2.5,
                           fill_color=lighten(DETECTOR_ACCENT, 0.85), fill_opacity=1.0).move_to(c)
            valence = VGroup(*[Dot(c + d, radius=0.045, color=col(GREY))
                               for d in (np.array([-0.2, 0.08, 0]), np.array([0.0, -0.1, 0]), np.array([0.2, 0.07, 0]))])
            protons.add(VGroup(blob, valence))
            links.add(fline(c + RIGHT * 0.5, p_in, sw=4.5))
            plabs.add(tex("p", 0.24, INK).next_to(blob, LEFT, buff=0.15))
        x1 = tex(r"x_1 P", 0.2, CY_DARK).next_to(links[0], DOWN, buff=0.1)
        x2 = tex(r"x_2 P", 0.2, CY_DARK).next_to(links[1], UP, buff=0.1)
        new_q = P["lab_q"].copy().move_to([links[0].get_center()[0], q_in[1] + 0.27, 0])
        new_qb = P["lab_qb"].copy().move_to([links[1].get_center()[0], qb_in[1] - 0.27, 0])
        self.play(FadeIn(protons), FadeIn(plabs), Create(links), Transform(P["lab_q"], new_q), Transform(P["lab_qb"], new_qb),
                  run_time=1.0)
        self.play(FadeIn(x1), FadeIn(x2), run_time=0.5)

        # the formula
        F = MathTex(r"\sigma", r"=", r"\sum_{q}\int\! dx_1\, dx_2\;", r"f_{q}(x_1)\,", r"f_{\bar q}(x_2)\;",
                    r"\hat\sigma_{q\bar q}(x_1 x_2 s)", color=col(INK))
        F.scale(_cap_scale("tex", 0.27)).move_to([-0.3, -1.0, 0.0])
        F[3].set_color(col(CY_DARK)); F[4].set_color(col(CY_DARK)); F[5].set_color(col(THEORY))
        F[5].shift(RIGHT * 0.16)
        self.play(LaggedStart(*[FadeIn(part, shift=UP * 0.1) for part in F], lag_ratio=0.25), run_time=1.5)
        self.remove(*F.submobjects)
        self.add(F)
        box_f = SurroundingRectangle(VGroup(F[3], F[4]), color=col(DETECTOR_ACCENT), buff=0.08, corner_radius=0.08, stroke_width=2.5)
        box_s = SurroundingRectangle(F[5], color=col(THEORY), buff=0.08, corner_radius=0.08, stroke_width=2.5)
        self.play(Create(box_s), run_time=0.5)
        self.play(Create(box_f), run_time=0.5)
        self.wait(0.3)

        # the perturbative series, with the year of each order
        S = MathTex(r"\hat\sigma", r"=", r"\hat\sigma^{(0)}", r"+", r"\alpha_s\,\hat\sigma^{(1)}", r"+",
                    r"\alpha_s^2\,\hat\sigma^{(2)}", r"+", r"\alpha_s^3\,\hat\sigma^{(3)}", r"+\cdots", color=col(THEORY))
        S.scale(_cap_scale("tex", 0.25)).move_to([0.0, -2.15, 0.0])
        orders, years = VGroup(), VGroup()
        for term, o in zip((S[2], S[4], S[6], S[8]), TH["orders"]):
            name = {"LO": r"\mathrm{LO}", "NLO": r"\mathrm{NLO}", "NNLO": r"\mathrm{NNLO}", "N3LO": r"\mathrm{N^3LO}"}[o["order"]]
            orders.add(tex(name, 0.19, INK).move_to([term.get_center()[0], -2.72, 0.0]))
            years.add(txt(str(o["year"]), 0.17, GREY).move_to([term.get_center()[0], -3.12, 0.0]))
        self.play(TransformFromCopy(F[5], S[0]), FadeIn(S[1]), run_time=0.8)
        for k, i in enumerate((2, 4, 6, 8)):
            anims = [FadeIn(S[i], shift=RIGHT * 0.1), FadeIn(orders[k]), FadeIn(years[k])]
            if i > 2:
                anims.insert(0, FadeIn(S[i - 1]))
            self.play(*anims, run_time=0.45)
        self.play(FadeIn(S[9]), run_time=0.3)
        self.remove(*S.submobjects)
        self.add(S)

        # inset: the line x1 x2 = m_Z^2 / s, the part at |y| < 2.4 lit
        O, L = np.array([4.35, 0.2, 0.0]), 2.1

        def mp(lx1, lx2):
            return O + np.array([(lx1 + 4) / 4 * L, (lx2 + 4) / 4 * L, 0.0])

        lc = math.log10(XR["tau"])
        ax = VGroup(Line(O, O + RIGHT * L, stroke_color=col(INK), stroke_width=2.2),
                    Line(O, O + UP * L, stroke_color=col(INK), stroke_width=2.2))
        ticks = VGroup()
        for lv, s in ((-4, r"10^{-4}"), (-2, r"10^{-2}"), (0, r"1")):
            p = mp(lv, -4)
            ticks.add(Line(p, p + DOWN * 0.07, stroke_color=col(INK), stroke_width=2.0),
                      tex(s, 0.13, INK).next_to(p + DOWN * 0.07, DOWN, buff=0.06))
            q = mp(-4, lv)
            ticks.add(Line(q, q + LEFT * 0.07, stroke_color=col(INK), stroke_width=2.0))
        t1 = tex(r"x_1", 0.18).next_to(O + RIGHT * L, RIGHT, buff=0.12)
        t2 = tex(r"x_2", 0.18).next_to(O + UP * L, UP, buff=0.08)
        hyper = Line(mp(-4, lc + 4), mp(lc + 4, -4), stroke_color=col(GREY), stroke_width=2.5)
        lo, hi = math.log10(XR["x_min"]), math.log10(XR["x_max"])
        lit = Line(mp(lo, lc - lo), mp(hi, lc - hi), stroke_color=col(THEORY), stroke_width=7)
        centre = Dot(mp(lc / 2, lc / 2), radius=0.05, color=col(INK))
        lab_h = tex(r"x_1 x_2 = m_Z^2/s", 0.17, THEORY).move_to(O + np.array([0.68 * L, 0.84 * L, 0.0]))
        mant, expo = f"{XR['x_min']:.0e}".split("e")
        lab_r = tex(r"x\in[%s\cdot10^{%d};\ %s]" % (mant, int(expo), f"{XR['x_max']:.1g}"), 0.17, THEORY).move_to([5.45, -0.33, 0.0])
        self.play(Create(ax), FadeIn(ticks), FadeIn(t1), FadeIn(t2), run_time=0.6)
        self.play(Create(hyper), FadeIn(lab_h), run_time=0.7)
        self.play(Create(lit), FadeIn(centre), FadeIn(lab_r), run_time=0.8)
        self.wait(0.1)


# ---------------------------------------------------------------------------
# 1-04 / 1-05: the prediction and why luminosity limits every sigma(Z)
# ---------------------------------------------------------------------------

S_LO, S_HI = 1850.0, 2050.0
AX_X0, AX_X1, AX_Y, BAND_TOP = -3.2, 6.6, -2.35, 2.1
LABEL_X = -3.45                  # right edge of the left-hand label column
ROW_Y = (1.45, 0.6, -0.25, -1.1)


def sx(s: float) -> float:
    return AX_X0 + (float(s) - S_LO) / (S_HI - S_LO) * (AX_X1 - AX_X0)


def sigma_axis(lo, hi, x0, x1, y, ticks, title_y):
    axis = Line([x0, y, 0], [x1, y, 0], stroke_color=col(INK), stroke_width=2.5)
    tk = VGroup()
    for t in ticks:
        x = x0 + (t - lo) / (hi - lo) * (x1 - x0)
        tk.add(Line([x, y, 0], [x, y - 0.09, 0], stroke_color=col(INK), stroke_width=2.5),
               txt(f"{t:.0f}", 0.19).move_to([x, y - 0.33, 0]))
    title = tex(r"\sigma(pp\to Z/\gamma^*\to\ell\ell),\ 60<m_{\ell\ell}<120\ \mathrm{GeV}\ \ [\mathrm{pb}]", 0.22)
    title.move_to([(x0 + x1) / 2, title_y, 0])
    return axis, tk, title


def hbar(s, err_down, err_up, y, color, sw=4.0, marker="dot", r=0.08):
    cc = col(color)
    bar = Line([sx(s - err_down), y, 0], [sx(s + err_up), y, 0], stroke_color=cc, stroke_width=sw)
    m = (Dot([sx(s), y, 0], radius=r, color=cc) if marker == "dot"
         else Square(side_length=2 * r, stroke_width=0, fill_color=cc, fill_opacity=1.0).move_to([sx(s), y, 0]))
    g = VGroup(bar, m)
    g.bar, g.marker = bar, m
    return g


class TheoryPrediction(Scene):
    def construct(self):
        white_background(self)
        clip_open(self, "prediction")
        axis, ticks, title = sigma_axis(S_LO, S_HI, AX_X0, AX_X1, AX_Y, [1850, 1900, 1950, 2000, 2050], -3.25)
        self.play(Create(axis), run_time=0.6)
        self.play(FadeIn(ticks), FadeIn(title), run_time=0.5)

        v = PRED["value"]
        line = Line([sx(v), AX_Y, 0], [sx(v), BAND_TOP, 0], stroke_color=col(THEORY), stroke_width=4)
        name = txt("aMC@NLO", 0.18, THEORY)
        val = tex(r"%.1f^{+%.1f}_{-%.1f}\ \mathrm{pb}" % (v, PRED["err_up"], PRED["err_down"]), 0.2, THEORY)
        name.move_to([4.75 + name.width / 2, 2.35, 0])
        val.move_to([4.75 + val.width / 2, 1.85, 0])
        self.play(Create(line), FadeIn(name), FadeIn(val), run_time=0.9)

        def band_rect(lo_pb, hi_pb):
            return Rectangle(width=max(sx(hi_pb) - sx(lo_pb), 1e-3), height=BAND_TOP - AX_Y, stroke_width=0,
                             fill_color=col(THEORY), fill_opacity=0.15).move_to([(sx(lo_pb) + sx(hi_pb)) / 2, (BAND_TOP + AX_Y) / 2, 0])

        band = band_rect(v - 0.01, v + 0.01)
        self.add(band)
        self.bring_to_front(line)
        comp = PRED["components_rel"]
        step_labels = [
            tex(r"\mathrm{scale}\ \ {}^{+%.1f}_{-%.1f}\,\%%" % (100 * comp["scale_up"], -100 * comp["scale_down"]), 0.21, THEORY),
            tex(r"\oplus\ \alpha_s\ \ %.1f\,\%%" % (100 * comp["alphas"]), 0.21, THEORY),
            tex(r"\oplus\ \mathrm{PDF}\ \ %.1f\,\%%" % (100 * comp["pdf"]), 0.21, THEORY),
        ]
        for lab, y in zip(step_labels, (1.9, 1.35, 0.8)):
            right_at(lab, LABEL_X, y)
        for step, lab in zip(PRED["band_steps"], step_labels):
            self.play(Transform(band, band_rect(step["lo_pb"], step["hi_pb"])), FadeIn(lab), run_time=0.8, rate_func=EASE)
            self.wait(0.25)
        self.wait(0.3)

        # four NNLO+NNLL predictions, one per PDF set
        head = right_at(tex(r"\mathrm{NNLO{+}NNLL}", 0.2, INK), LABEL_X, 2.1)
        rows = VGroup()
        for p, y in zip(PDFS, ROW_Y):
            rows.add(VGroup(right_at(txt(p["label"], 0.2), LABEL_X, y), hbar(p["value"], p["down"], p["up"], y, TH_DARK)))
        self.play(FadeOut(VGroup(*step_labels)), FadeIn(head), run_time=0.5)
        self.play(LaggedStart(*[FadeIn(r, shift=DOWN * 0.1) for r in rows], lag_ratio=0.3), run_time=1.6)
        self.remove(*self.mobjects)
        self.add(axis, ticks, title, band, line, name, val, head, rows)
        lo_s, hi_s = min(p["value"] for p in PDFS), max(p["value"] for p in PDFS)
        yb = -1.8
        bracket = VGroup(Line([sx(lo_s), yb, 0], [sx(hi_s), yb, 0], stroke_color=col(INK), stroke_width=2.5),
                         Line([sx(lo_s), yb - 0.1, 0], [sx(lo_s), yb + 0.1, 0], stroke_color=col(INK), stroke_width=2.5),
                         Line([sx(hi_s), yb - 0.1, 0], [sx(hi_s), yb + 0.1, 0], stroke_color=col(INK), stroke_width=2.5))
        spread = tex(r"%d\ \mathrm{pb}" % TH["pdf_sets"]["spread_pb"], 0.2).next_to(bracket, RIGHT, buff=0.2)
        self.play(GrowFromCenter(bracket), FadeIn(spread), run_time=0.7)

        clip_cut(self, "lumi_limited")

        # ---- the best published sigma(Z) at 13 TeV: luminosity dominates
        self.play(FadeOut(rows), FadeOut(head), run_time=0.6)
        c = CMS24
        y0 = ROW_Y[0]
        cms_lab = right_at(txt("CMS 2024", 0.2), LABEL_X, y0)
        cms = hbar(c["value"], c["total"], c["total"], y0, SLATE, sw=4.5, marker="square", r=0.09)
        self.play(FadeIn(cms_lab), FadeIn(cms.marker), GrowFromCenter(cms.bar), run_time=0.9)
        self.wait(0.3)
        parts = []
        for key, y, colour in (("stat", ROW_Y[1], SLATE), ("syst", ROW_Y[2], SLATE), ("lumi", ROW_Y[3], HIGHLIGHT)):
            bar = Line([sx(c["value"] - c[key]), y, 0], [sx(c["value"] + c[key]), y, 0], stroke_color=col(colour), stroke_width=5)
            lab = right_at(tex(r"\mathrm{%s}\ \pm%d" % (key, c[key]), 0.21, colour), LABEL_X, y)
            parts.append((bar, lab))
        self.play(*[TransformFromCopy(cms.bar, b) for b, _ in parts], *[FadeIn(lab) for _, lab in parts], run_time=1.2, rate_func=EASE)
        self.wait(0.3)
        pulse(self, parts[2][0], width=16, run_time=0.8)
        self.wait(0.1)


# ---------------------------------------------------------------------------
# 1-06: the Z under the discoveries
# ---------------------------------------------------------------------------

class TheoryEverywhere(Scene):
    def construct(self):
        white_background(self)
        clip_open(self, "z_everywhere")

        # left: m_tautau, the Higgs under the Z tail (schematic shapes, area ratio from the JSON)
        e_lo, e_hi = -3.5, 0.25
        dax = DataAxes([40, 200, 40], [e_lo, e_hi, 1], 4.7, 3.8, x_ticks=[50, 100, 150, 200], y_ticks=[-3, -2, -1, 0],
                       y_log=True, x_title=r"m_{\tau\tau}\ [\mathrm{GeV}]", tick_label_h=0.19)
        dax.move_frame_to((-2.3, -0.25))
        floor = 10 ** e_lo
        ms = np.arange(40.0, 200.01, 1.0)
        sz, sh = 13.0, 17.0                                   # typical di-tau mass resolution at 91 and 125 GeV
        z = np.exp(-0.5 * ((ms - 91.0) / sz) ** 2)
        h_peak = sz / sh / HIG["z_over_h_tautau"]             # equal-width Gaussians: area ratio = 1 / (Z/H)
        h = h_peak * np.exp(-0.5 * ((ms - 125.0) / sh) ** 2)
        zpoly = Polygon(dax.c2p(40, floor), *[dax.c2p(m, max(v, floor)) for m, v in zip(ms, z)], dax.c2p(200, floor),
                        stroke_color=col(CHANNEL_LINE["tautau"]), stroke_width=2.0, fill_color=col(CHANNEL["tautau"]), fill_opacity=0.75)
        hm = ms[(ms >= 88) & (ms <= 162)]
        hv = h_peak * np.exp(-0.5 * ((hm - 125.0) / sh) ** 2)
        hpoly = Polygon(dax.c2p(hm[0], floor), *[dax.c2p(m, max(v, floor)) for m, v in zip(hm, hv)], dax.c2p(hm[-1], floor),
                        stroke_color=col(SLATE), stroke_width=3.0, fill_color=col(SLATE), fill_opacity=0.25)
        zlab = tex(r"Z\to\tau\tau", 0.24, CHANNEL_LINE["tautau"]).move_to(dax.c2p(140, 10 ** -0.15))
        hlab = tex(r"H\to\tau\tau", 0.24, SLATE)
        hlab.move_to(dax.c2p(172, 10 ** -1.6))
        ratio = tex(r"\sigma\mathcal{B}\,/\,%d" % round(HIG["z_over_h_tautau"], -1), 0.2, SLATE).next_to(hlab, DOWN, buff=0.12)
        hptr = Line(hlab.get_bottom() + DOWN * 0.42 + LEFT * 0.35, dax.c2p(128, h_peak * 1.25), stroke_color=col(SLATE), stroke_width=2.0)
        self.play(Create(dax.ax), FadeIn(VGroup(dax.x_ticks_v, dax.x_labels, dax.x_title, dax.y_ticks_v, dax.y_labels)), run_time=0.7)
        self.play(FadeIn(zpoly, shift=UP * 0.2), FadeIn(zlab), run_time=1.0)
        self.wait(0.3)
        self.play(FadeIn(hpoly, shift=UP * 0.1), run_time=0.8)
        self.play(FadeIn(hlab), FadeIn(ratio), Create(hptr), run_time=0.7)
        self.wait(0.5)

        # right: Z -> mumu + jet, seen as Z -> nunu + jet
        det = mini_slice(0.6, (3.95, -0.1))
        phi_z, phi_jet = math.radians(18), math.radians(200)
        mu1 = signature(det, "mu", math.radians(-12), charge=-1, kappa=0.45, sw=3.0)
        mu2 = signature(det, "mu", math.radians(52), charge=+1, kappa=0.45, sw=3.0)
        jet = signature(det, "jet", phi_jet, sw=3.0, seed=4)
        lab_mm = tex(r"Z\to\mu\mu+\mathrm{jet}", 0.24).move_to([3.95, -2.45, 0])
        self.play(FadeIn(det), run_time=0.6)
        self.play(Create(mu1[0]), Create(mu2[0]), FadeIn(jet), run_time=0.9)
        self.play(FadeIn(VGroup(*mu1[1:])), FadeIn(VGroup(*mu2[1:])), FadeIn(lab_mm), run_time=0.5)
        self.wait(0.5)
        nu1 = signature(det, "nu", math.radians(-12))
        nu2 = signature(det, "nu", math.radians(52))
        met = signature(det, "met", phi_z)
        lab_nn = tex(r"Z\to\nu\nu+\mathrm{jet}", 0.24).move_to(lab_mm)
        self.play(FadeOut(mu1), FadeOut(mu2), FadeIn(nu1), FadeIn(nu2), ReplacementTransform(lab_mm, lab_nn), run_time=1.0)
        self.play(FadeIn(met), run_time=0.6)
        self.wait(0.1)


# ---------------------------------------------------------------------------
# 1-07: three flavours, one Z
# ---------------------------------------------------------------------------

M_X0, M_Y = -4.4, 1.35           # mass axis: 2 units per decade from 0.1 MeV


def mx(m_mev: float) -> float:
    return M_X0 + (math.log10(m_mev) + 1.0) * 2.0


CP, CU = np.array([-3.1, -1.6, 0.0]), 2.0      # coupling plane centre, units per coupling unit


def cpt(ga, gv):
    return CP + np.array([ga * CU, gv * CU, 0.0])


class TheoryUniversality(DrellYan):
    def construct(self):
        white_background(self)
        P = self.build_parts()
        for k in ("q", "qb", "q_tip", "qb_tip", "v1", "boson", "v2", "lm", "lp", "lm_tip", "lp_tip",
                  "lab_q", "lab_qb", "lab_b", "lab_lm", "lab_lp"):
            self.add(P[k])
        clip_open(self, "lepton_universality")

        # the lepton legs split into the three flavours
        v2 = P["v2"].get_center()
        a_m = math.atan2(*(P["lm"].get_end() - v2)[1::-1])
        a_p = math.atan2(*(P["lp"].get_end() - v2)[1::-1])
        L = float(np.linalg.norm(P["lm"].get_end() - v2))
        legs, labels = VGroup(), VGroup()
        for k, (fl, sym, da) in enumerate((("e", "e^-", 0.24), ("mu", r"\mu^-", 0.0), ("tau", r"\tau^-", -0.24))):
            up = fline(v2, v2 + L * unit(a_m + da), color=LEP_COL[fl], sw=4.5)
            dn = fline(v2, v2 + L * unit(a_p - da), color=LEP_COL[fl], sw=4.5)
            legs.add(up, dn)
            labels.add(tex(sym, 0.3, LEP_COL[fl]).next_to(up.get_end(), RIGHT, buff=0.18))
        self.play(FadeOut(VGroup(P["lm"], P["lp"], P["lm_tip"], P["lp_tip"], P["lab_lm"], P["lab_lp"])),
                  LaggedStart(*[Create(l) for l in legs], lag_ratio=0.1), FadeIn(labels), run_time=1.2)
        self.wait(0.4)

        # the log mass axis: 207 and 17 between them
        axis = Line([M_X0, M_Y, 0], [mx(10 ** 4), M_Y, 0], stroke_color=col(INK), stroke_width=2.5)
        ticks = VGroup()
        for k, s in zip(range(-1, 5), ("0.1", "1", "10", "10^2", "10^3", "10^4")):
            x = mx(10 ** k)
            ticks.add(Line([x, M_Y, 0], [x, M_Y - 0.09, 0], stroke_color=col(INK), stroke_width=2.5),
                      tex(s, 0.17).move_to([x, M_Y - 0.33, 0]))
        m_title = tex(r"m\ [\mathrm{MeV}]", 0.2).next_to(axis, RIGHT, buff=0.2)
        dots, targets = VGroup(), []
        for fl, lab in zip(("e", "mu", "tau"), labels):
            x = mx(LEP["mass_MeV"][fl])
            dots.add(Dot([x, M_Y, 0], radius=0.09, color=col(LEP_COL[fl])))
            targets.append(lab.copy().move_to([x, M_Y + 0.5, 0]))
        rest = VGroup(*[P[k] for k in ("q", "qb", "q_tip", "qb_tip", "v1", "boson", "v2", "lab_q", "lab_qb", "lab_b")], legs)
        self.play(FadeOut(rest), *[Transform(lab, t) for lab, t in zip(labels, targets)], Create(axis), FadeIn(ticks),
                  FadeIn(m_title), run_time=1.4, rate_func=EASE)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in dots], lag_ratio=0.2), run_time=0.6)
        r1 = tex(r"\times%d" % round(LEP["ratio_mu_e"]), 0.2, GREY).move_to([(dots[0].get_x() + dots[1].get_x()) / 2, M_Y + 0.5, 0])
        r2 = tex(r"\times%d" % round(LEP["ratio_tau_mu"]), 0.2, GREY).move_to([(dots[1].get_x() + dots[2].get_x()) / 2, M_Y + 0.5, 0])
        self.play(FadeIn(r1), FadeIn(r2), run_time=0.5)
        self.wait(0.4)

        # the (g_A, g_V) plane: sin^2 theta runs, the three leptons stay on one point
        R = 1.3
        plane = VGroup(Line(CP + LEFT * R, CP + RIGHT * R, stroke_color=col(INK), stroke_width=2.2),
                       Line(CP + DOWN * R, CP + UP * R, stroke_color=col(INK), stroke_width=2.2))
        gA = tex(r"g_A", 0.2).next_to(CP + RIGHT * R, RIGHT, buff=0.1)
        gV = tex(r"g_V", 0.2).next_to(CP + UP * R, UP, buff=0.06)
        pticks = VGroup()
        for v in (-0.5, 0.5):
            p = cpt(v, 0)
            pticks.add(Line(p + DOWN * 0.06, p + UP * 0.06, stroke_color=col(INK), stroke_width=2.0))
            q = cpt(0, v)
            pticks.add(Line(q + LEFT * 0.06, q + RIGHT * 0.06, stroke_color=col(INK), stroke_width=2.0),
                       tex(r"%s\tfrac{1}{2}" % ("-" if v < 0 else "+"), 0.13).next_to(q, RIGHT, buff=0.1))
        s2 = ValueTracker(0.0)
        C = LEP["couplings"]

        def gv(f):
            return C[f]["T3"] - 2 * C[f]["Q"] * s2.get_value()

        def mover(f, dx=0.0, colour=GREY, r=0.075):
            return always_redraw(lambda: Dot(cpt(C[f]["g_A"], gv(f)) + RIGHT * dx, radius=r, color=col(colour)))

        others = VGroup(mover("nu"), mover("u"), mover("d"))
        o_labs = VGroup(
            always_redraw(lambda: tex(r"\nu", 0.17, GREY).next_to(cpt(0.5, gv("nu")), RIGHT, buff=0.12)),
            always_redraw(lambda: tex(r"u", 0.17, GREY).next_to(cpt(0.5, gv("u")), RIGHT, buff=0.12)),
            always_redraw(lambda: tex(r"d", 0.17, GREY).next_to(cpt(-0.5, gv("d")), LEFT, buff=0.12)),
        )
        triad = VGroup(mover("l", -0.08, LEP_COL["e"]), mover("l", 0.0, LEP_COL["mu"]), mover("l", 0.08, LEP_COL["tau"]))
        readout_l = tex(r"\sin^2\theta_W=", 0.2)
        num = DecimalNumber(0.0, num_decimal_places=4, color=col(INK))
        num.scale_to_fit_height(tex("0.2315", 0.2).height)
        readout = VGroup(readout_l, num).arrange(RIGHT, buff=0.1).move_to([CP[0], -3.3, 0])
        self.play(Create(plane), FadeIn(gA), FadeIn(gV), FadeIn(pticks), run_time=0.6)
        self.play(*[TransformFromCopy(d, t) for d, t in zip(dots, triad)], FadeIn(others), FadeIn(o_labs), FadeIn(readout),
                  run_time=1.0, rate_func=EASE)
        num.add_updater(lambda m: m.set_value(s2.get_value()))
        self.play(s2.animate.set_value(LEP["sin2_eff"]["value"]), run_time=2.2, rate_func=EASE)
        num.clear_updaters()
        self.remove(s2)
        for grp in (others, o_labs, triad):
            for m in grp:
                m.clear_updaters()
        gvl = tex(r"g_V^{\ell}=%.3f" % C["l"]["g_V"], 0.19, INK).next_to(cpt(-0.5, C["l"]["g_V"]), UP, buff=0.22)
        self.play(FadeIn(gvl), run_time=0.4)
        self.wait(0.3)

        # equal branching fractions; the LEP + SLD ratios
        base, per = -2.75, 0.55
        yax = Line([0.25, base, 0], [0.25, base + 4.2 * per, 0], stroke_color=col(INK), stroke_width=2.2)
        xax = Line([0.25, base, 0], [3.9, base, 0], stroke_color=col(INK), stroke_width=2.2)
        yt = VGroup()
        for k in range(5):
            y = base + k * per
            yt.add(Line([0.25, y, 0], [0.16, y, 0], stroke_color=col(INK), stroke_width=2.0),
                   tex(str(k) + (r"\,\%" if k == 4 else ""), 0.15).next_to([0.16, y, 0], LEFT, buff=0.08))
        br = LEP["br_ll"]["value"]
        bars, blabs = VGroup(), VGroup()
        for x, fl, sym in ((0.9, "e", r"e^+e^-"), (2.1, "mu", r"\mu^+\mu^-"), (3.3, "tau", r"\tau^+\tau^-")):
            fill = CHANNEL[{"e": "ee", "mu": "mumu", "tau": "tautau"}[fl]]
            bars.add(Rectangle(width=0.72, height=br * per, stroke_color=col(LEP_COL[fl]), stroke_width=2.0,
                               fill_color=col(fill), fill_opacity=0.9).move_to([x, base + br * per / 2, 0]))
            blabs.add(tex(sym, 0.18, LEP_COL[fl]).move_to([x, base - 0.3, 0]))
        brline = DashedLine([0.3, base + br * per, 0], [3.85, base + br * per, 0], dash_length=0.1, stroke_color=col(THEORY), stroke_width=2.5)
        brlab = tex(r"\mathcal{B}(Z\to\ell\ell)=%.3f\,\%%" % br, 0.2, THEORY).move_to([2.1, base + br * per + 0.42, 0])
        self.play(Create(yax), Create(xax), FadeIn(yt), FadeIn(blabs), run_time=0.6)
        self.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars], lag_ratio=0.2), run_time=1.1)
        self.play(Create(brline), FadeIn(brlab), run_time=0.6)
        R2 = LEP["lep_ratios"]
        lep_head = txt("LEP + SLD", 0.17, GREY).move_to([5.45, -0.8, 0])
        rat1 = tex(r"\Gamma_{\mu\mu}/\Gamma_{ee}=%.4f\pm%.4f" % (R2["mumu_over_ee"]["value"], R2["mumu_over_ee"]["err"]), 0.17).move_to([5.45, -1.35, 0])
        rat2 = tex(r"\Gamma_{\tau\tau}/\Gamma_{ee}=%.4f\pm%.4f" % (R2["tautau_over_ee"]["value"], R2["tautau_over_ee"]["err"]), 0.17).move_to([5.45, -1.95, 0])
        self.play(FadeIn(lep_head), FadeIn(rat1, shift=UP * 0.1), FadeIn(rat2, shift=UP * 0.1), run_time=0.8)
        self.wait(0.1)


# ---------------------------------------------------------------------------
# 1-08 / 1-09: the window, and the empty frame of the final plot
# ---------------------------------------------------------------------------

F_LO, F_HI, F_X0, F_X1, F_Y = 1650.0, 2400.0, -2.7, 6.4, -2.55
F_ROWS = ((r"Z\to ee", CHANNEL_LINE["ee"], 1.65), (r"Z\to\mu\mu", CHANNEL_LINE["mumu"], 0.8),
          (r"Z\to\tau\tau", CHANNEL_LINE["tautau"], -0.05), (r"Z\to\ell\ell", COMBINED, -1.55))
F_DIVIDER = -0.8


def fx(s: float) -> float:
    return F_X0 + (float(s) - F_LO) / (F_HI - F_LO) * (F_X1 - F_X0)


class TheoryMeasurement(Scene):
    def construct(self):
        white_background(self)
        clip_open(self, "mass_window")
        edges = np.asarray(SPEC["edges_GeV"], dtype=float)
        dens = np.asarray(SPEC["dsigma_dm_pb_per_GeV"], dtype=float)
        dax = DataAxes([50, 200, 25], [-2, 3, 1], 9.8, 4.4, x_ticks=[50, 75, 100, 125, 150, 175, 200], y_log=True,
                       x_title=r"m_{\ell\ell}\ [\mathrm{GeV}]", y_title=r"d\sigma/dm_{\ell\ell}\ [\mathrm{pb/GeV}]", tick_label_h=0.19)
        dax.move_frame_to((1.1, -0.15))
        spec = step_hist(dax, edges, dens, color=INK, stroke_width=3.0)
        self.play(Create(dax.ax), FadeIn(VGroup(dax.x_ticks_v, dax.y_ticks_v, dax.x_labels, dax.y_labels, dax.x_title, dax.y_title)),
                  run_time=0.8)
        self.play(Create(spec), run_time=2.6, rate_func=rate_functions.linear)
        zlab = tex(r"Z/\gamma^*", 0.22).next_to(dax.c2p(91, dens.max()), RIGHT, buff=0.35)
        self.play(FadeIn(zlab), run_time=0.4)
        self.wait(0.3)

        sel = (edges[:-1] >= 60) & (edges[:-1] < 120)
        win_edges = np.append(edges[:-1][sel], 120.0)
        window = step_hist(dax, win_edges, dens[sel], color=THEORY, stroke_width=0, fill_opacity=0.35)
        v60, v120 = vref_line(dax, 60, color=GREY), vref_line(dax, 120, color=GREY)
        self.play(FadeIn(v60), FadeIn(v120), run_time=0.4)
        self.play(FadeIn(window, shift=UP * 0.15), run_time=0.8)
        self.bring_to_front(spec)

        lhs = tex(r"\sigma(60<m_{\ell\ell}<120)=", 0.22, THEORY)
        num = DecimalNumber(0.0, num_decimal_places=1, group_with_commas=False, color=col(THEORY))
        num.scale_to_fit_height(tex("1953.9", 0.22).height)
        unit_pb = tex(r"\mathrm{pb}", 0.22, THEORY)
        lhs.move_to([2.0 + lhs.width / 2, 1.55, 0])
        num.next_to(lhs, RIGHT, buff=0.12)
        unit_pb.next_to(num, RIGHT, buff=0.12)
        num.add_updater(lambda m: m.next_to(lhs, RIGHT, buff=0.12))
        unit_pb.add_updater(lambda m: m.next_to(num, RIGHT, buff=0.12))
        frac = tex(r"=%.3f\ \sigma(m_{\ell\ell}>50\ \mathrm{GeV})" % SPEC["fraction_60_120"], 0.22, THEORY)
        frac.move_to([lhs.get_left()[0] + frac.width / 2, 0.95, 0])
        self.play(FadeIn(lhs), FadeIn(num), FadeIn(unit_pb), run_time=0.4)
        self.play(ChangeDecimalToValue(num, SPEC["sigma_60_120_pb"]), run_time=1.3, rate_func=EASE)
        num.clear_updaters(); unit_pb.clear_updaters()
        self.play(FadeIn(frac, shift=UP * 0.1), run_time=0.5)

        clip_cut(self, "final_frame")

        # ---- the window becomes the prediction line of the final plot
        v = PRED["value"]
        axis, ticks, title = sigma_axis(F_LO, F_HI, F_X0, F_X1, F_Y, list(range(1700, 2401, 100)), -3.45)
        pline = Line([fx(v), F_Y, 0], [fx(v), 2.2, 0], stroke_color=col(THEORY), stroke_width=4)
        plot_parts = VGroup(dax, spec, zlab, v60, v120, lhs, num, unit_pb, frac)
        self.play(FadeOut(plot_parts), ReplacementTransform(window, pline), run_time=1.3, rate_func=EASE)
        self.play(Create(axis), FadeIn(ticks), FadeIn(title), run_time=0.8)
        band0 = Rectangle(width=1e-3, height=2.2 - F_Y, stroke_width=0, fill_color=col(THEORY), fill_opacity=0.15).move_to([fx(v), (2.2 + F_Y) / 2, 0])
        band = Rectangle(width=fx(v + PRED["err_up"]) - fx(v - PRED["err_down"]), height=2.2 - F_Y, stroke_width=0,
                         fill_color=col(THEORY), fill_opacity=0.15).move_to([(fx(v + PRED["err_up"]) + fx(v - PRED["err_down"])) / 2, (2.2 + F_Y) / 2, 0])
        name = txt("aMC@NLO", 0.18, THEORY).move_to([fx(v), 2.45, 0])
        self.add(band0)
        self.bring_to_front(pline)
        self.play(Transform(band0, band), FadeIn(name), run_time=0.8)
        rows = VGroup()
        for sym, colour, y in F_ROWS:
            guide = DashedLine([F_X0, y, 0], [F_X1, y, 0], dash_length=0.12, stroke_color=col(LIGHT_GREY), stroke_width=2.0)
            rows.add(VGroup(guide, right_at(tex(sym, 0.26, colour), F_X0 - 0.25, y)))
        divider = Line([F_X0 - 2.3, F_DIVIDER, 0], [F_X1, F_DIVIDER, 0], stroke_color=col(GREY), stroke_width=1.8)
        self.play(LaggedStart(*[FadeIn(r, shift=RIGHT * 0.15) for r in rows[:3]], lag_ratio=0.3), run_time=1.2)
        self.play(Create(divider), FadeIn(rows[3], shift=RIGHT * 0.15), run_time=0.7)
        self.wait(0.1)
