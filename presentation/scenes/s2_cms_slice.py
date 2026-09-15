"""Section 2: the CMS detector in the transverse plane, official logo colours.

    CMSSliceBuild   the detector drawn inside out (beam pipe -> muon system)
    CMSSignatures   the detector already there; e, mu, tau_h (+nu), jet, gamma
                    signatures appear one after another (the legend clip)

The detector itself (``CMSSlice``: logo radii, logo colours, centre (0,-0.25),
outer radius 2.95) and ``signature()`` live in style/bnd_style.py and are
the basis of the per-channel event-display clips in sections 3-5.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from manim import (  # noqa: E402
    Create, FadeIn, GrowFromCenter, LaggedStart, Scene, VGroup,
)
from style.bnd_style import *  # noqa: E402,F401,F403


def build_detector() -> CMSSlice:
    return CMSSlice()


class CMSSliceBuild(Scene):
    """Inside-out build-up, ~7 s. Ends on the complete detector."""

    def construct(self):
        white_background(self)
        det = build_detector()
        P = dict(det.parts)
        self.play(GrowFromCenter(P["beampipe"]), run_time=0.4)
        self.play(GrowFromCenter(P["tracker"][0]), run_time=0.5)
        self.play(LaggedStart(*[Create(r) for r in P["tracker"][1]], lag_ratio=0.25), run_time=1.0)
        self.play(GrowFromCenter(P["tob"]), run_time=0.4)
        self.play(FadeIn(P["ecal"][0]), run_time=0.3)
        self.play(LaggedStart(*[FadeIn(c) for c in det.ecal_cells], lag_ratio=0.5 / N_ECAL), run_time=1.0)
        self.play(FadeIn(P["band"]), run_time=0.2)
        self.play(FadeIn(P["hcal"][0]), run_time=0.3)
        self.play(LaggedStart(*[FadeIn(c) for c in det.hcal_towers], lag_ratio=0.5 / N_HCAL), run_time=1.0)
        self.play(GrowFromCenter(P["solenoid"]), run_time=0.6)
        self.play(GrowFromCenter(P["muon"]), run_time=0.6)
        muon_layers = [P[n] for n in ("mb1", "yoke1", "mb2", "yoke2", "mb3", "yoke3", "mb4")]
        self.play(LaggedStart(*[FadeIn(g) for g in muon_layers], lag_ratio=0.3), run_time=1.4)
        self.remove(*[m for _, g in det.parts for m in [g]])
        self.add(det)
        self.wait(0.1)


# kind, angle (deg), symbol, charge, extra kwargs
SIGNATURES = [
    ("e",     28,  r"e^-",    -1, {}),
    ("mu",    150, r"\mu^+",  +1, {}),
    ("tau_h", 212, r"\tau_h", -1, {"prongs": 3, "pi0": False}),
    ("jet",   275, None,      +1, {}),
    ("gamma", 335, r"\gamma",  0, {}),
]
LABEL_COLOR = {"e": CHANNEL_LINE["ee"], "mu": CHANNEL_LINE["mumu"],
               "tau_h": CHANNEL_LINE["tautau"], "gamma": PARTICLE["photon"]}


class CMSSignatures(Scene):
    """Detector present from frame one; each signature is added with its
    symbol just outside the detector. The tau_h comes with its neutrino."""

    def construct(self):
        white_background(self)
        det = build_detector()
        self.add(det)
        self.wait(0.3)
        r_lab = det.outer_radius + 0.42
        for kind, deg, tex, q, kw in SIGNATURES:
            phi = math.radians(deg)
            sig = signature(det, kind, phi, charge=q, seed=7, **kw)
            trk, deposits = sig[0], VGroup(*sig[1:])
            self.play(Create(trk, lag_ratio=0.0), run_time=0.8)
            if len(deposits):
                self.play(FadeIn(deposits), run_time=0.4)
            if tex is not None:
                lab = mathtex(tex, color=LABEL_COLOR[kind]).scale(1.1).move_to(det.point_at(r_lab, phi))
                self.play(FadeIn(lab), run_time=0.3)
            if kind == "tau_h":
                nu = signature(det, "nu", phi + math.radians(16))
                self.play(Create(nu[0]), FadeIn(nu[1]), run_time=0.6)
                lab = mathtex(r"\nu_\tau", color=PARTICLE["nu"]).scale(1.0).move_to(
                    det.point_at(r_lab, phi + math.radians(16)))
                self.play(FadeIn(lab), run_time=0.3)
        self.wait(0.1)
