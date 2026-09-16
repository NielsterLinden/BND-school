"""Section 5: Z -> tautau, both taus hadronic (tau_h tau_h).

    TauDecay          one tau: tau- -> nu_tau W*-, the virtual W turns into a
                      d ubar pair that hadronises into pions (pi- pi0, decay
                      mode 1, the most common). The ONLY neutrino is the
                      nu_tau from the tau vertex; the W* gives hadrons only.
    TauJet            starts on the last frame of TauDecay: the tau is boosted,
                      so the pions collapse into a narrow cone (the tau jet);
                      the pi0 splits into two photons.
    ZtautauDetector   starts on the last frame of TauJet, dissolves to the
                      CMS slice: a 1-prong tau_h (track + HCAL, gamma gamma in
                      the ECAL) and a 3-prong tau_h, both neutrinos leave
                      unseen, and the missing transverse momentum recoils.

Decay modes 0, 1, 10, 11 are the ones the analysis uses
(z-tautau/docs/02-selection.md); the clip shows DM1 (Feynman) and DM1 + DM10
(detector).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from manim import Polygon, Transform  # noqa: E402
from channel_common import *  # noqa: E402,F401,F403

TAU_COL = CHANNEL_LINE["tautau"]
# vertices of the tau decay diagram (scene units); tune here only
_TD = {
    "tau_in": (-5.6, 0.7),
    "v1":     (-2.8, 0.7),      # tau -> nu W*
    "nu_out": (-0.9, 2.3),
    "v2":     (-0.7, -0.5),     # W* -> d ubar
    "d_out":  ( 1.1,  0.4),
    "u_out":  ( 1.1, -1.4),
    "blob":   ( 1.7, -0.5),
    "pi_ang": (0.28, -0.32),    # pi-, pi0 directions (rad)
    "pi_len": 1.9,
}
SW = 4.5
TD_SHIFT = (0.8, 0.0)       # clears the chapter-identifier corner (x < -5.85, y > 0.22) for the tau label
# 5-03 ends zoomed (scale 0.5) on this centre, a little left of the detector centre so the outer ring
# stays out of the chapter-identifier corner; s5_ztautau_story.A["zoom_503"] opens on it
ZOOM_503 = (-0.2, -0.25)


def _pt(k):
    return np.array([_TD[k][0] + TD_SHIFT[0], _TD[k][1] + TD_SHIFT[1], 0.0])


class TauDecay(Scene):
    """(5-01) The diagram sits TD_SHIFT to the right of its original place so the tau label,
    sliding in from the left, stays out of the chapter-identifier corner."""

    def build_parts(self) -> dict:
        tau = fline(_pt("tau_in"), _pt("v1"), color=TAU_COL, sw=SW)
        tau_tip = arrow_tip_on(tau, color=TAU_COL, at=0.55)
        v1 = vertex_dot(_pt("v1"))
        nu = dashed(_pt("v1"), _pt("nu_out"), color=PARTICLE["nu"], sw=SW - 1)
        nu_tip = arrow_tip_on(fline(_pt("v1"), _pt("nu_out")), color=PARTICLE["nu"], at=0.6)
        w = wavy(_pt("v1"), _pt("v2"), color=PARTICLE["W"], sw=SW - 0.5, amplitude=0.13, wavelength=0.6)
        v2 = vertex_dot(_pt("v2"))
        d = fline(_pt("v2"), _pt("d_out"), color=PARTICLE["quark"], sw=SW - 0.5)
        u = fline(_pt("v2"), _pt("u_out"), color=PARTICLE["quark"], sw=SW - 0.5)
        d_tip = arrow_tip_on(d, at=0.55)
        u_tip = arrow_tip_on(u, at=0.55)
        u_tip.rotate(np.pi, about_point=u.point_from_proportion(0.55))   # antiquark
        blob = hadron_blob(_pt("blob"), width=0.55, height=1.9, angle=0.0)
        p0 = _pt("blob") + RIGHT * 0.2
        pions = VGroup(pion_lines(p0, _TD["pi_ang"][:1], length=_TD["pi_len"], color=TAU_COL, tip=True)[0],
                       pion_lines(p0, _TD["pi_ang"][1:], length=_TD["pi_len"], color=PARTICLE["pion"],
                                  tip=False)[0])                       # neutral pion: grey, no arrow
        lab = dict(
            tau=mathtex(r"\tau^{-}", color=TAU_COL).scale(1.1).next_to(_pt("tau_in"), LEFT, buff=0.25),
            nu=mathtex(r"\nu_\tau", color=PARTICLE["nu"]).scale(1.0).next_to(_pt("nu_out"), UP + RIGHT, buff=0.12),
            w=mathtex(r"W^{*-}").scale(1.0).next_to(w, LEFT, buff=0.22),
            d=mathtex(r"d").scale(0.9).next_to(_pt("d_out"), UP, buff=0.12),
            u=mathtex(r"\bar{u}").scale(0.9).next_to(_pt("u_out"), DOWN, buff=0.12),
            pim=mathtex(r"\pi^{-}", color=TAU_COL).scale(1.0).next_to(pions[0][0].get_end(), RIGHT, buff=0.15),
            pi0=mathtex(r"\pi^{0}", color=PARTICLE["pion"]).scale(1.0).next_to(pions[1][0].get_end(), RIGHT, buff=0.15),
        )
        return dict(tau=tau, tau_tip=tau_tip, v1=v1, nu=nu, nu_tip=nu_tip, w=w, v2=v2, d=d, u=u,
                    d_tip=d_tip, u_tip=u_tip, blob=blob, pions=pions, lab=lab)

    @staticmethod
    def group(P) -> VGroup:
        g = VGroup(*[v for k, v in P.items() if k != "lab"], *P["lab"].values())
        return g

    def construct(self):
        white_background(self)
        P = self.build_parts()
        L = P["lab"]
        self.play(FadeIn(L["tau"], shift=RIGHT * 0.2), run_time=0.4)
        self.play(Create(P["tau"]), run_time=1.0, rate_func=EASE)
        self.play(FadeIn(P["tau_tip"]), FadeIn(P["v1"]), run_time=0.3)
        self.play(Flash(P["v1"].get_center(), color=SLATE, line_length=0.2, flash_radius=0.4, run_time=0.5))
        self.play(Create(P["nu"]), Create(P["w"]), run_time=1.0, rate_func=EASE)
        self.play(FadeIn(P["nu_tip"]), FadeIn(L["nu"]), FadeIn(L["w"]), FadeIn(P["v2"]), run_time=0.4)
        self.play(Create(P["d"]), Create(P["u"]), run_time=0.8, rate_func=EASE)
        self.play(FadeIn(P["d_tip"]), FadeIn(P["u_tip"]), FadeIn(L["d"]), FadeIn(L["u"]), run_time=0.3)
        self.play(FadeIn(P["blob"], scale=0.6), run_time=0.5)
        self.play(LaggedStart(Create(P["pions"][0][0]), Create(P["pions"][1][0]), lag_ratio=0.2),
                  run_time=0.9, rate_func=EASE)
        self.play(FadeIn(P["pions"][0][1]), FadeIn(L["pim"]), FadeIn(L["pi0"]), run_time=0.4)
        self.wait(0.1)


JET_ORIGIN = np.array([-3.2, 0.0, 0.0])
JET_LEN = 5.0
JET_ANG = (0.09, -0.11)     # collimated pi-, pi0 directions
GAMMA_ANG = (0.07, -0.20)   # the two photons of the pi0


class TauJet(TauDecay):
    """The boosted tau: the pions of TauDecay close up into a narrow cone."""

    def end_state(self):
        origin = JET_ORIGIN
        pim = fline(origin, origin + JET_LEN * unit(JET_ANG[0]), color=TAU_COL, sw=SW)
        pim_tip = arrow_tip_on(pim, color=TAU_COL, at=0.6, tip_length=0.18)
        pi0_end = origin + 0.45 * JET_LEN * unit(JET_ANG[1])
        pi0 = fline(origin, pi0_end, color=PARTICLE["pion"], sw=SW - 0.5)
        gam = VGroup(*[wavy(pi0_end, pi0_end + 0.55 * JET_LEN * unit(a), color=PARTICLE["photon"],
                            sw=SW - 1.5, amplitude=0.07, wavelength=0.4) for a in GAMMA_ANG])
        half = 0.30
        cone = Polygon(origin, origin + 1.08 * JET_LEN * unit(half), origin + 1.08 * JET_LEN * unit(-half),
                       fill_color=col(CHANNEL["tautau"]), fill_opacity=0.08, stroke_width=0)
        blob = hadron_blob(origin, width=0.45, height=1.2)
        lab = dict(
            pim=mathtex(r"\pi^{-}", color=TAU_COL).scale(1.0).next_to(pim.get_end(), RIGHT, buff=0.15),
            g1=mathtex(r"\gamma", color=PARTICLE["photon"]).scale(1.0).next_to(gam[0].get_end(), RIGHT, buff=0.15),
            g2=mathtex(r"\gamma", color=PARTICLE["photon"]).scale(1.0).next_to(gam[1].get_end(), RIGHT, buff=0.15),
        )
        return dict(cone=cone, blob=blob, pim=pim, pim_tip=pim_tip, pi0=pi0, gam=gam, lab=lab)

    def construct(self):
        white_background(self)
        P = self.build_parts()
        self.add(self.group(P))
        self.wait(0.3)
        # keep only the hadronisation and its pions, slide them to the left
        keep = VGroup(P["blob"], P["pions"], P["lab"]["pim"], P["lab"]["pi0"])
        rest = VGroup(*[m for m in self.group(P) if m not in keep])
        self.play(FadeOut(rest), run_time=0.6)
        E = self.end_state()
        self.play(keep.animate.shift(JET_ORIGIN - P["blob"].get_center()), run_time=0.8, rate_func=EASE)
        # the boost: pions collimate and stretch, the cone appears
        pim_new = VGroup(E["pim"], E["pim_tip"])
        pi0_new = VGroup(E["pi0"])
        self.add(E["cone"]); E["cone"].set_opacity(0)
        self.play(Transform(P["blob"], E["blob"]),
                  Transform(P["pions"][0], pim_new), Transform(P["pions"][1], pi0_new),
                  Transform(P["lab"]["pim"], E["lab"]["pim"]),
                  P["lab"]["pi0"].animate.next_to(E["pi0"].get_end(), DOWN, buff=0.12),
                  E["cone"].animate.set_opacity(0.08),
                  run_time=1.4, rate_func=EASE)
        # pi0 -> gamma gamma
        self.play(FadeOut(P["lab"]["pi0"]), run_time=0.3)
        self.play(Create(E["gam"][0]), Create(E["gam"][1]), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(E["lab"]["g1"]), FadeIn(E["lab"]["g2"]), run_time=0.3)
        self.wait(0.1)


class ZtautauDetector(TauJet, MovingCameraScene):

    def construct(self):
        white_background(self)
        E = self.end_state()
        old = VGroup(E["cone"], E["blob"], E["pim"], E["pim_tip"], E["pi0"], E["gam"], *E["lab"].values())
        self.add(old)
        self.wait(0.3)
        det = CMSSlice()
        crossfade_to_detector(self, old, det)
        # the taus are boosted (gamma ~ 25 for a Z at rest): each nu_tau leaves
        # almost collinear with its tau_h; p_T^miss is the (small) vector sum
        a1, a2 = math.radians(40), math.radians(218)
        d_nu, d_lab = math.radians(6), math.radians(7)
        r_lab = det.outer_radius + 0.42
        g1 = show_signature(self, det, "tau_h", a1, -1, r"\tau_h", TAU_COL, prongs=1, pi0=True, seed=2,
                            label_phi=a1 - d_lab)
        nu1 = signature(det, "nu", a1 + d_nu)
        self.play(Create(nu1[0]), FadeIn(nu1[1]), run_time=0.5)
        g2 = show_signature(self, det, "tau_h", a2, +1, r"\tau_h", TAU_COL, prongs=3, pi0=False, seed=4,
                            label_phi=a2 + d_lab)
        nu2 = signature(det, "nu", a2 - d_nu)
        self.play(Create(nu2[0]), FadeIn(nu2[1]), run_time=0.5)
        labs = VGroup(
            mathtex(r"\nu_\tau", color=PARTICLE["nu"]).move_to(det.point_at(r_lab, a1 + d_nu + d_lab)),
            mathtex(r"\nu_\tau", color=PARTICLE["nu"]).move_to(det.point_at(r_lab, a2 - d_nu - d_lab)))
        self.play(FadeIn(labs), run_time=0.3)
        # the neutrinos carry momentum away: missing transverse momentum
        phi_met = 0.5 * ((a1 + d_nu) + (a2 - d_nu))          # bisector of the two nu directions
        met = signature(det, "met", phi_met, r_end=det.radii["hcal"][0] + 0.1)
        self.play(Create(met[0]), FadeIn(met[1]), run_time=0.7, rate_func=EASE)
        met_lab = mathtex(r"p_T^{\,\mathrm{miss}}", color=PARTICLE["met"]).scale(0.9).move_to(
            det.point_at(det.radii["solenoid"][1] + 0.42, phi_met))
        self.play(FadeIn(met_lab), run_time=0.3)
        zoom_in(self, det, 0.5, fade=(g1[-1], g2[-1], labs, met_lab), at=ZOOM_503)
        self.wait(0.1)
