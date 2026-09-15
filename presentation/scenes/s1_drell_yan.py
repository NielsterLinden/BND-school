"""Section 1 (theory) -- the Drell-Yan process q qbar -> Z/gamma* -> l+ l-.

One tree-level diagram built left to right in reading order: incoming quark
lines converge on the production vertex, the boson propagates as a wavy line,
the lepton pair opens up at the decay vertex. The clip ends on a clean held
frame; the user adds the title and the cross-section formula in PowerPoint.

Physics-symbol labels only (q, qbar, Z/gamma*, l+, l-): exempt from the
animation-only rule. Lepton flavour is deliberately generic (l) so the same
clip serves the ee / mumu / tautau openers; the per-channel subclasses below
(chapter colours, CHANNEL_LINE for strokes and glyphs) are reused by the
section 3-5 process clips through ``build_parts`` + ``SHIFT``.

Source of the topology: any textbook (Drell & Yan 1970); the fitting contract
in fitting/CONVENTIONS.md calls the signal DYee / DYmumu / DYtautau.
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    DOWN, LEFT, RIGHT, UP,
    Create, FadeIn, Flash, LaggedStart, Scene, VGroup,
    rate_functions,
)

from style.bnd_style import (  # noqa: E402
    INK, HIGHLIGHT, PARTICLE, CHANNEL_LINE,
    arrow_tip_on, fline, mathtex, vertex_dot, wavy, white_background,
)

# every vertex and label anchor in one dict (scene units). Tune here only.
_FEY = {
    "q_in":    (-4.6,  1.6),
    "qb_in":   (-4.6, -1.6),
    "v1":      (-2.0,  0.0),      # production vertex
    "v2":      ( 2.0,  0.0),      # decay vertex
    "l_out":   ( 4.6,  1.6),
    "lb_out":  ( 4.6, -1.6),
}
SW = 4.5


class DrellYan(Scene):
    LEPTON_TEX = (r"\ell^{-}", r"\ell^{+}")
    LEPTON_COLOR = INK
    BOSON_TEX = r"Z/\gamma^{*}"

    SHIFT = (0.0, 0.0)          # channel clips move the diagram to make room

    def build_parts(self) -> dict:
        F = {k: np.array([v[0] + self.SHIFT[0], v[1] + self.SHIFT[1], 0.0]) for k, v in _FEY.items()}
        self.F = F
        q = fline(F["q_in"], F["v1"], sw=SW)
        qb = fline(F["qb_in"], F["v1"], sw=SW)
        q_tip = arrow_tip_on(q, at=0.55)
        qb_tip = arrow_tip_on(qb, at=0.55)            # antiquark: arrow against flow
        qb_tip.rotate(np.pi, about_point=qb.point_from_proportion(0.55))
        boson = wavy(F["v1"], F["v2"], sw=SW - 0.5, amplitude=0.13, wavelength=0.6)
        lm = fline(F["v2"], F["l_out"], color=self.LEPTON_COLOR, sw=SW)
        lp = fline(F["v2"], F["lb_out"], color=self.LEPTON_COLOR, sw=SW)
        lm_tip = arrow_tip_on(lm, color=self.LEPTON_COLOR, at=0.55)
        lp_tip = arrow_tip_on(lp, color=self.LEPTON_COLOR, at=0.55)
        lp_tip.rotate(np.pi, about_point=lp.point_from_proportion(0.55))
        v1, v2 = vertex_dot(F["v1"]), vertex_dot(F["v2"])

        lab_q = mathtex("q").scale(1.1).next_to(F["q_in"], LEFT, buff=0.25)
        lab_qb = mathtex(r"\bar{q}").scale(1.1).next_to(F["qb_in"], LEFT, buff=0.25)
        lab_b = mathtex(self.BOSON_TEX).scale(1.05).next_to(boson, UP, buff=0.3)
        lab_lm = mathtex(self.LEPTON_TEX[0], color=self.LEPTON_COLOR).scale(1.1)\
            .next_to(F["l_out"], RIGHT, buff=0.25)
        lab_lp = mathtex(self.LEPTON_TEX[1], color=self.LEPTON_COLOR).scale(1.1)\
            .next_to(F["lb_out"], RIGHT, buff=0.25)
        return dict(q=q, qb=qb, q_tip=q_tip, qb_tip=qb_tip, v1=v1, boson=boson, v2=v2,
                    lm=lm, lp=lp, lm_tip=lm_tip, lp_tip=lp_tip,
                    lab_q=lab_q, lab_qb=lab_qb, lab_b=lab_b, lab_lm=lab_lm, lab_lp=lab_lp)

    def construct(self):
        white_background(self)
        P = self.build_parts()
        ease = rate_functions.ease_in_out_sine

        # 1. the incoming quarks
        self.play(FadeIn(P["lab_q"], shift=RIGHT * 0.2), FadeIn(P["lab_qb"], shift=RIGHT * 0.2),
                  run_time=0.6)
        self.play(Create(P["q"]), Create(P["qb"]), run_time=1.2, rate_func=ease)
        self.play(FadeIn(P["q_tip"]), FadeIn(P["qb_tip"]), FadeIn(P["v1"]), run_time=0.4)
        # 2. annihilation -> boson
        self.play(Flash(P["v1"].get_center(), color=HIGHLIGHT, line_length=0.2,
                        flash_radius=0.4, run_time=0.5))
        self.play(Create(P["boson"]), run_time=1.2, rate_func=ease)
        self.play(FadeIn(P["lab_b"], shift=DOWN * 0.15), FadeIn(P["v2"]), run_time=0.5)
        # 3. the lepton pair
        self.play(LaggedStart(Create(P["lm"]), Create(P["lp"]), lag_ratio=0.15),
                  run_time=1.2, rate_func=ease)
        self.play(FadeIn(P["lm_tip"]), FadeIn(P["lp_tip"]),
                  FadeIn(P["lab_lm"], shift=LEFT * 0.2), FadeIn(P["lab_lp"], shift=LEFT * 0.2),
                  run_time=0.5)
        self.wait(0.1)   # render the true final state (trap 1 in the recipes)


class DrellYanEE(DrellYan):
    LEPTON_TEX = ("e^{-}", "e^{+}")
    LEPTON_COLOR = CHANNEL_LINE["ee"]


class DrellYanMuMu(DrellYan):
    LEPTON_TEX = (r"\mu^{-}", r"\mu^{+}")
    LEPTON_COLOR = CHANNEL_LINE["mumu"]


class DrellYanTauTau(DrellYan):
    LEPTON_TEX = (r"\tau^{-}", r"\tau^{+}")
    LEPTON_COLOR = CHANNEL_LINE["tautau"]
