"""Section 6: combination (skeleton).

    ThreeToOne   three schematic measurements (ee, mumu, tautau: a marker with
                 an error bar each, in the chapter colours) sit on one axis
                 with the theory line; they slide together into one combined
                 point with a shorter bar. Schematic positions only: no
                 numbers (the frozen results come in a later pass through
                 presentation/data/).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from manim import (  # noqa: E402
    DOWN, LEFT, RIGHT, UP, Create, DashedLine, Dot, FadeIn, FadeOut, LaggedStart, Line,
    ReplacementTransform, Scene, VGroup, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403

X0, X1 = -4.0, 4.0          # the sigma axis (schematic units)
ROWS = [  # (label, colour, centre offset, half-width, y)
    (r"e^{+}e^{-}",       CHANNEL_LINE["ee"],     0.25, 0.9, 1.4),
    (r"\mu^{+}\mu^{-}",   CHANNEL_LINE["mumu"],  -0.35, 0.7, 0.5),
    (r"\tau^{+}\tau^{-}", CHANNEL_LINE["tautau"], 0.6, 1.9, -0.4),
]
# combined: the inverse-variance mean of the three and its (not much smaller)
# uncertainty, so the merge cannot be read as "beating statistics"
COMB = (-0.04, 0.56, -1.6)   # centre offset, half-width, y
EASE = rate_functions.ease_in_out_sine


def point(x, hw, y, color, r=0.09, sw=4.0) -> VGroup:
    cc = col(color)
    return VGroup(Line([x - hw, y, 0], [x + hw, y, 0], stroke_color=cc, stroke_width=sw),
                  Dot([x, y, 0], radius=r, color=cc))


class ThreeToOne(Scene):

    def construct(self):
        white_background(self)
        axis = Line([X0, -2.4, 0], [X1, -2.4, 0], stroke_color=col(INK), stroke_width=2.5)
        theory = DashedLine([0, -2.4, 0], [0, 2.0, 0], dash_length=0.16,
                            stroke_color=col(THEORY), stroke_width=3.0)
        self.play(Create(axis), run_time=0.6)
        self.play(Create(theory), run_time=0.8)
        pts, labs = VGroup(), VGroup()
        for tex, color, dx, hw, y in ROWS:
            pts.add(point(dx, hw, y, color))
            labs.add(mathtex(tex, color=color).scale(0.95).next_to([X0, y, 0], RIGHT, buff=0.0)
                     .shift(LEFT * 1.3))
        for p, l in zip(pts, labs):
            self.play(FadeIn(l), Create(p[0]), FadeIn(p[1]), run_time=0.7, rate_func=EASE)
        self.wait(0.4)
        comb = point(COMB[0], COMB[1], COMB[2], COMBINED, r=0.11, sw=5.0)
        ghosts = VGroup(*[p.copy().set_opacity(0.25) for p in pts])
        self.add(ghosts)
        self.play(*[ReplacementTransform(p, comb.copy()) for p in pts], run_time=1.4, rate_func=EASE)
        self.add(comb)
        self.play(FadeIn(mathtex(r"Z \to \ell\ell", color=COMBINED).scale(0.95)
                         .next_to([X0, COMB[2], 0], RIGHT, buff=0.0).shift(LEFT * 1.3)), run_time=0.4)
        self.wait(0.1)
