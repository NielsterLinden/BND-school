"""Section 2: the CMS detector grows out of the CMS logo.

    CMSLogoToSlice   (a) the logo square appears, its four muons draw in;
                     (b) the square expands toward the top right until its
                         corner sits at the detector centre; frame and
                         off-white background dissolve;
                     (c) the quarter layers "round up": every layer sweeps from
                         90 to 360 degrees about the corner; the logo muons fade;
                     (d) the detail fades in (beam pipe, tracker rings, ECAL
                         cells, HCAL towers, muon stations and yoke) so the last
                         frame is exactly ``CMSSlice()`` as used by every other
                         detector clip.

The logo geometry lives in style/bnd_style.py (``LOGO_R``, ``logo_rings``,
``logo_muons``, ``CMSLogo``); this file is only the choreography.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    BOLD, TAU, UP, LEFT, RIGHT, Create, FadeIn, FadeOut, LaggedStart, Scene, ValueTracker, VGroup,
    RoundedRectangle, always_redraw, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403

LOGO_SIDE0 = 2.2                                    # the logo as it first appears
CORNER0 = DET_CENTER + np.array([-LOGO_SIDE0 / 2, -LOGO_SIDE0 / 2, 0.0])   # square centred on the slide
SIDE1 = R_DET / LOGO_R["muon"]                      # side that gives the delivered radius


class CMSLogoToSlice(Scene):

    def construct(self):
        white_background(self)
        ease = rate_functions.ease_in_out_sine
        side = ValueTracker(LOGO_SIDE0)
        cx, cy = ValueTracker(CORNER0[0]), ValueTracker(CORNER0[1])
        ang = ValueTracker(TAU / 4)
        box_op = ValueTracker(1.0)

        def corner():
            return np.array([cx.get_value(), cy.get_value(), 0.0])

        def box(fill):
            s = side.get_value()
            r = RoundedRectangle(width=s, height=s, corner_radius=0.02 * s,
                                 fill_color=col(CMS["bg"]), fill_opacity=box_op.get_value() if fill else 0,
                                 stroke_color=col(INK), stroke_width=2.0,
                                 stroke_opacity=0 if fill else box_op.get_value())
            return r.move_to(corner() + s / 2 * (RIGHT + UP))

        bg = always_redraw(lambda: box(True))
        rings = always_redraw(lambda: logo_rings(corner(), side.get_value(), ang.get_value(), 0.0))
        muons = logo_muons(CORNER0, LOGO_SIDE0)
        frame = always_redraw(lambda: box(False))
        word_op = ValueTracker(0.0)

        def wordmark_at():
            s = side.get_value()
            w = text("CMS", color=CMS["font"], weight=BOLD).scale_to_fit_height(0.135 * s)
            w.move_to(corner() + np.array([0.058 * s, 0.888 * s, 0.0]), aligned_edge=UP + LEFT)
            return w.set_opacity(word_op.get_value())

        wordmark = always_redraw(wordmark_at)

        # (a) the logo
        self.add(bg, rings)
        self.play(FadeIn(bg), FadeIn(rings), run_time=0.6)
        self.add(frame)
        self.play(LaggedStart(*[Create(m, lag_ratio=0.0) for m in muons], lag_ratio=0.2),
                  run_time=1.4, rate_func=ease)
        self.add(wordmark)
        self.play(word_op.animate.set_value(1.0), run_time=0.4)
        self.wait(0.4)

        # (b) expand toward the top right; the square dissolves
        muons_dyn = always_redraw(lambda: logo_muons(corner(), side.get_value()))
        self.remove(*muons); self.add(muons_dyn)
        self.play(side.animate.set_value(SIDE1),
                  cx.animate.set_value(DET_CENTER[0]), cy.animate.set_value(DET_CENTER[1]),
                  box_op.animate.set_value(0.0), word_op.animate.set_value(0.0),
                  run_time=1.6, rate_func=ease)
        self.remove(bg, frame, wordmark)

        # (c) round up
        self.play(ang.animate.set_value(TAU), FadeOut(muons_dyn, run_time=0.7),
                  run_time=2.0, rate_func=ease)
        self.remove(rings)
        det = CMSSlice()
        self.add(det.base)

        # (d) the detail: what makes the logo a detector
        P = dict(det.parts)
        detail_groups = [
            P["tracker"][1],          # no beam pipe: the tracker disc covers it in the composite
            VGroup(*det.ecal_cells), VGroup(*det.hcal_towers),
            VGroup(P["mb1"], P["yoke1"]), VGroup(P["mb2"], P["yoke2"]),
            VGroup(P["mb3"], P["yoke3"]), P["mb4"],
        ]
        self.play(LaggedStart(*[FadeIn(g) for g in detail_groups], lag_ratio=0.18), run_time=2.0)
        self.remove(det.base, *detail_groups)
        self.add(det)
        self.wait(0.1)
