"""Section 4: Z -> mumu (deliberately the simplest of the three).

    ZmumuProcess    the Drell-Yan diagram in the mumu flavour; the muon legs
                    simply continue, straight and clean, to the frame edge:
                    no shower, nothing else happens.
    ZmumuIntoDetector  starts on the last frame of ZmumuProcess: the slice
                    fades in behind the diagram while the diagram shrinks into
                    the interaction point (the collision happens *there*); a
                    flash, then the two muons bend through the tracker, leave
                    only minimum-ionising dots in both calorimeters and end as
                    stubs in the four muon stations. Its last frame is the one
                    s4_zmumu_story.detector_end_state() rebuilds (the chain
                    opens on it): change one, change the other.
                    (Replaces ZmumuDetector, the cross-fade of 15 Sep 2026.)
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from channel_common import *  # noqa: E402,F401,F403
from manim import Transform  # noqa: E402  (not an introducer: keeps the slice *behind* the diagram, FadeIn would re-add it on top)

EXT = 2.4      # how far the legs continue (scene units)


class ZmumuProcess(ChannelDiagram, Scene):
    LEPTON_TEX = (r"\mu^{-}", r"\mu^{+}")
    LEPTON_COLOR = CHANNEL_LINE["mumu"]

    def end_state(self):
        P = self.build_parts()
        a_m, a_p = self.leg_angles(P)
        ext = VGroup()
        for leg, ang in ((P["lm"], a_m), (P["lp"], a_p)):
            e = leg.get_end()
            ln = fline(e, e + EXT * unit(ang), color=self.LEPTON_COLOR, sw=4.5)
            ext.add(VGroup(ln, arrow_tip_on(ln, color=self.LEPTON_COLOR, at=0.8)))
        return P, ext

    def construct(self):
        white_background(self)
        P, ext = self.end_state()
        self.play_diagram(P)
        self.wait(0.3)
        self.play(*[Create(e[0]) for e in ext], run_time=1.0, rate_func=EASE)
        self.play(*[FadeIn(e[1]) for e in ext], run_time=0.3)
        self.wait(0.1)


SHRINK = 0.05  # the diagram at the interaction point: a speck, fully transparent by then


class ZmumuIntoDetector(ZmumuProcess):

    def construct(self):
        white_background(self)
        P, ext = self.end_state()
        old = VGroup(self.diagram_group(P), ext)
        det = CMSSlice()
        self.add(det, old)                       # the slice sits behind the diagram, invisible until the Transform below
        det.set_opacity(0.0)
        self.wait(0.3)
        shown = CMSSlice()
        self.play(Transform(det, shown), run_time=1.0, rate_func=EASE)
        self.play(old.animate.scale(SHRINK).move_to(det.c).set_opacity(0.0), run_time=1.5,
                  rate_func=rate_functions.ease_in_cubic)
        self.remove(old)
        self.play(Flash(det.c, color=HIGHLIGHT, line_length=0.22, flash_radius=0.3, run_time=0.45))
        show_signature(self, det, "mu", math.radians(35), -1, r"\mu^-", CHANNEL_LINE["mumu"], kappa=0.28)
        show_signature(self, det, "mu", math.radians(212), +1, r"\mu^+", CHANNEL_LINE["mumu"], kappa=0.28)
        self.wait(0.1)
