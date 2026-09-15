"""Section 3: Z -> ee.

    ZeeProcess    the Drell-Yan diagram in the ee flavour; then each electron
                  leg radiates (e -> e gamma), the photons convert
                  (gamma -> e+ e-) and the cascade repeats: an electromagnetic
                  shower drawn in Feynman style, inside a faint cone.
    ZeeDetector   starts on the last frame of ZeeProcess, dissolves to the
                  CMS slice: e- and e+ tracks bend opposite ways in the
                  tracker and stop in the ECAL as compact clusters (the shower
                  of the first act); one soft bremsstrahlung photon; nothing in
                  the HCAL or the muon system.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from manim import Polygon  # noqa: E402
from channel_common import *  # noqa: E402,F401,F403

SHOWER = dict(depth=3, length=0.85, spread=0.34, shrink=0.72, sw=4.0)
BEND = 0.30     # showers lean toward the beam axis (the electrons are boosted)


class ZeeProcess(ChannelDiagram, Scene):
    LEPTON_TEX = ("e^{-}", "e^{+}")
    LEPTON_COLOR = CHANNEL_LINE["ee"]

    def end_state(self):
        """Everything on screen at the end of the clip (also the first frame
        of ZeeDetector)."""
        P = self.build_parts()
        a_m, a_p = self.leg_angles(P)
        showers, cones = VGroup(), VGroup()
        for leg, ang, seed in ((P["lm"], a_m, 3), (P["lp"], a_p, 5)):
            root = leg.get_end()
            ang = ang * BEND
            tree = shower_tree(root, ang, seed=seed, **SHOWER)
            L = sum(SHOWER["length"] * SHOWER["shrink"] ** g for g in range(SHOWER["depth"])) * 1.05
            half = 0.45 * SHOWER["spread"] * SHOWER["depth"]
            cone = Polygon(root, root + L * unit(ang + half), root + L * unit(ang - half),
                           fill_color=col(CHANNEL["ee"]), fill_opacity=0.10, stroke_width=0)
            showers.add(tree); cones.add(cone)
        return P, showers, cones

    def construct(self):
        white_background(self)
        P, showers, cones = self.end_state()
        self.play_diagram(P)
        self.wait(0.3)
        # the showers, generation by generation, both legs together
        self.play(FadeIn(cones), run_time=0.5)
        for g in range(SHOWER["depth"]):
            self.play(LaggedStart(*[Create(showers[i][g], lag_ratio=0.0) for i in range(2)],
                                  lag_ratio=0.1), run_time=0.9 - 0.15 * g, rate_func=EASE)
        self.wait(0.1)


class ZeeDetector(ZeeProcess, MovingCameraScene):

    def construct(self):
        white_background(self)
        P, showers, cones = self.end_state()
        old = VGroup(cones, self.diagram_group(P), showers)
        self.add(old)
        self.wait(0.3)
        det = CMSSlice()
        crossfade_to_detector(self, old, det)
        r_lab = det.radii["hcal"][1] + 0.38
        show_signature(self, det, "e", math.radians(32), -1, "e^-", CHANNEL_LINE["ee"],
                       r_label=r_lab, brem=True, sw=4.5)
        show_signature(self, det, "e", math.radians(207), +1, "e^+", CHANNEL_LINE["ee"],
                       r_label=r_lab, sw=4.5)
        zoom_in(self, det, 0.30)
        self.wait(0.1)
