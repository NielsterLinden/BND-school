"""Shared choreography of the section 3-5 channel clips.

Every channel has a *process* clip (Feynman act) and a *detector* clip
(event-display act). The detector clip starts on the last frame of the
process clip, so both build the same final Feynman state through
``process_end_state`` and the seam is checked with tools/framediff.py.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    DOWN, LEFT, RIGHT, UP, Create, FadeIn, FadeOut, Flash, LaggedStart, MovingCameraScene,
    Scene, VGroup, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403
from s1_drell_yan import DrellYan  # noqa: E402

DIAGRAM_SHIFT = (-1.6, -0.3)    # the Drell-Yan diagram sits left and a bit low: room for the legs
EASE = rate_functions.ease_in_out_sine


class ChannelDiagram(DrellYan):
    """DrellYan with the diagram shifted and the lepton labels beside the legs
    (not at their ends, where the process clips continue the legs)."""
    SHIFT = DIAGRAM_SHIFT

    def build_parts(self) -> dict:
        P = super().build_parts()
        P["lab_lm"].next_to(P["lm"].point_from_proportion(0.55), UP + LEFT, buff=0.12)
        P["lab_lp"].next_to(P["lp"].point_from_proportion(0.55), DOWN + LEFT, buff=0.12)
        return P

    def diagram_group(self, P) -> VGroup:
        return VGroup(*P.values())

    def play_diagram(self, P, fast: bool = True):
        """The Drell-Yan diagram, built left to right (~3 s)."""
        f = 0.7 if fast else 1.0
        self.play(FadeIn(P["lab_q"], shift=RIGHT * 0.2), FadeIn(P["lab_qb"], shift=RIGHT * 0.2),
                  run_time=0.4 * f)
        self.play(Create(P["q"]), Create(P["qb"]), run_time=1.0 * f, rate_func=EASE)
        self.play(FadeIn(P["q_tip"]), FadeIn(P["qb_tip"]), FadeIn(P["v1"]), run_time=0.3 * f)
        self.play(Flash(P["v1"].get_center(), color=HIGHLIGHT if self.LEPTON_COLOR != CHANNEL_LINE["tautau"]
                        else SLATE, line_length=0.2, flash_radius=0.4, run_time=0.5 * f))
        self.play(Create(P["boson"]), run_time=1.0 * f, rate_func=EASE)
        self.play(FadeIn(P["lab_b"], shift=DOWN * 0.15), FadeIn(P["v2"]), run_time=0.4 * f)
        self.play(LaggedStart(Create(P["lm"]), Create(P["lp"]), lag_ratio=0.15),
                  run_time=1.0 * f, rate_func=EASE)
        self.play(FadeIn(P["lm_tip"]), FadeIn(P["lp_tip"]),
                  FadeIn(P["lab_lm"]), FadeIn(P["lab_lp"]), run_time=0.4 * f)

    def leg_angles(self, P):
        """Direction (rad) of the l- and l+ legs."""
        d1 = P["lm"].get_end() - P["lm"].get_start()
        d2 = P["lp"].get_end() - P["lp"].get_start()
        return float(np.arctan2(d1[1], d1[0])), float(np.arctan2(d2[1], d2[0]))


def crossfade_to_detector(scene, old: VGroup, det: CMSSlice, run_time: float = 1.2):
    """The Feynman act dissolves while the detector fades in behind it."""
    scene.play(FadeOut(old), FadeIn(det), run_time=run_time)
    scene.remove(old)


def show_signature(scene, det, kind, phi, charge, tex, tex_color, hits=True, r_label=None,
                   label_phi=None, **kw):
    """Standard reveal: track drawn, tracker hits, deposits, then the symbol
    just outside the ring. Returns the VGroup of everything added."""
    sig = signature(det, kind, phi, charge=charge, **kw)
    trk, rest = sig[0], VGroup(*sig[1:])
    scene.play(Create(trk, lag_ratio=0.0), run_time=0.9, rate_func=EASE)
    added = VGroup(trk)
    if hits and kind in ("e", "mu", "tau_h"):
        t0 = trk if kind != "tau_h" else trk[len(trk) // 2]
        th = tracker_hits(det, t0, color=CHANNEL[{"e": "ee", "mu": "mumu", "tau_h": "tautau"}[kind]])
        scene.play(FadeIn(th), run_time=0.3)
        added.add(th)
    if len(rest):
        scene.play(FadeIn(rest), run_time=0.45)
        added.add(rest)
    if tex is not None:
        r_label = det.outer_radius + 0.42 if r_label is None else r_label
        lab = mathtex(tex, color=tex_color).scale(1.1).move_to(
            det.point_at(r_label, phi if label_phi is None else label_phi))
        scene.play(FadeIn(lab), run_time=0.3)
        added.add(lab)
    return added


def zoom_in(scene: MovingCameraScene, det: CMSSlice, scale: float, run_time: float = 1.4, fade=()):
    """End of an event-display clip: the camera closes in on the inner
    detector (the logo geometry keeps the tracker small, so the tracks and
    calorimeter clusters deserve a closer look on the last frame). ``fade``:
    mobjects (outer labels) that would be cut by the new frame."""
    anims = [scene.camera.frame.animate.scale(scale).move_to(det.c)]
    anims += [FadeOut(m) for m in fade]
    scene.play(*anims, run_time=run_time, rate_func=EASE)
