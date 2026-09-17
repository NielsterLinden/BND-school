"""Section 5: why tau reconstruction is hard, and the three Z -> tautau final states.

    TauChallenges     (slide 1) the tau flies c*tau = 87 um and decays: its tracks start
                      almost on top of the primary vertex; the nu_tau leaves unseen;
                      the tau splits 65 % hadronic / 35 % leptonic; a q/g jet looks like
                      a (narrow) tau_h jet; the background processes as symbols.
    TautauChannels    (slide 2) the three final states side by side, each a Feynman
                      diagram Z -> tau tau -> nu W* -> (l nu | q qbar'), with its share of
                      Z -> tautau (12 / 46 / 42 %), its neutrino count and its dominant
                      backgrounds. Chained clips a -> b -> c (one diagram each):
                      clip k opens on the final frame of clip k-1 (``end_state(k)``).

Branching fractions: PDG, B(tau -> l nu nu) = 17.82 % + 17.39 % = 35.2 %, hadronic 64.8 %;
pairs: ll = 0.352^2 = 12.4 %, lh = 2 * 0.352 * 0.648 = 45.6 %, hh = 0.648^2 = 42.0 %.
Lifetime tau = 2.903e-13 s, c*tau = 87.03 um. No narrative text (anchors A1): the words
are on the PowerPoint slide.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from manim import AnimationGroup, Arc, DashedLine, GrowFromEdge, Indicate, Line, Polygon, Rectangle, Transform  # noqa: E402
from channel_common import *  # noqa: E402,F401,F403

TAU_COL = CHANNEL_LINE["tautau"]
NU_COL = PARTICLE["nu"]
SW = 4.0

# PDG constants (not analysis results)
BR_LEP, BR_HAD = 35.2, 64.8
FRACTIONS = {"ll": 12, "lh": 46, "hh": 42}
assert abs(BR_LEP ** 2 / 100 - FRACTIONS["ll"]) < 0.6
assert abs(2 * BR_LEP * BR_HAD / 100 - FRACTIONS["lh"]) < 0.6
assert abs(BR_HAD ** 2 / 100 - FRACTIONS["hh"]) < 0.6


def _p(x, y):
    return np.array([x, y, 0.0])


# ---------------------------------------------------------------------------
# Slide 1: TauChallenges
# ---------------------------------------------------------------------------
PV = _p(-4.6, -2.2)          # primary vertex (zoomed view: not to scale)
SV = _p(-2.2, -1.5)          # tau decay vertex, 87 um away


class TauChallenges(Scene):

    def lifetime_parts(self) -> dict:
        beam = DashedLine(_p(-5.6, -2.2), _p(-1.4, -2.2), dash_length=0.12,
                          stroke_color=col(GREY), stroke_width=2.0)
        pv = vertex_dot(PV, radius=0.09)
        tau = fline(PV, SV, color=TAU_COL, sw=SW + 1)
        tau_tip = arrow_tip_on(tau, color=TAU_COL, at=0.55)
        sv = vertex_dot(SV, color=TAU_COL, radius=0.08)
        lab_tau = mathtex(r"\tau^-", color=TAU_COL).scale(1.0).next_to(tau.point_from_proportion(0.5), UP + LEFT, buff=0.1)
        # ruler between the two vertices
        d = SV - PV
        n = np.array([d[1], -d[0], 0.0]) / np.linalg.norm(d)
        off = 0.35 * n
        ruler = VGroup(Line(PV + off, SV + off, stroke_color=col(INK), stroke_width=2.0),
                       Line(PV + off * 0.6, PV + off * 1.4, stroke_color=col(INK), stroke_width=2.0),
                       Line(SV + off * 0.6, SV + off * 1.4, stroke_color=col(INK), stroke_width=2.0))
        lab_ct = mathtex(r"c\tau \approx 87\,\mu\mathrm{m}").scale(0.8).move_to(
            (PV + SV) / 2 + 0.75 * n).rotate(float(np.arctan2(d[1], d[0])))
        lab_t = mathtex(r"\tau_\tau = 2.9\times 10^{-13}\,\mathrm{s}").scale(0.8).move_to(_p(-4.0, -3.35))
        # decay products: three charged tracks fanning from SV, one neutrino
        prongs = VGroup(*[fline(SV, SV + 3.0 * unit(a), color=TAU_COL, sw=SW - 1) for a in (0.42, 0.30, 0.18)])
        nu = dashed(SV, SV + 2.3 * unit(0.95), color=NU_COL, sw=SW - 1)
        nu_lab = mathtex(r"\nu_\tau", color=NU_COL).next_to(nu.get_end(), UP + RIGHT, buff=0.08)
        # other tracks from the primary vertex (pileup / underlying event)
        others = VGroup(*[fline(PV, PV + 2.2 * unit(a), color=PARTICLE["jet"], sw=2.5) for a in (0.75, 1.15, -0.45, 2.5)])
        return dict(beam=beam, pv=pv, tau=tau, tau_tip=tau_tip, sv=sv, lab_tau=lab_tau, ruler=ruler,
                    lab_ct=lab_ct, lab_t=lab_t, prongs=prongs, nu=nu, nu_lab=nu_lab, others=others)

    def br_parts(self) -> dict:
        x0, x1, y = -5.2, 1.2, 1.9
        w = x1 - x0
        xs = x0 + w * BR_HAD / 100
        had = Rectangle(width=xs - x0, height=0.55, fill_color=col(CHANNEL["tautau"]), fill_opacity=0.9,
                        stroke_width=0).move_to(_p((x0 + xs) / 2, y))
        lep = Rectangle(width=x1 - xs, height=0.55, fill_color=col(tint(SLATE, 0.55)), fill_opacity=0.9,
                        stroke_width=0).move_to(_p((xs + x1) / 2, y))
        lab_had = mathtex(r"\tau_h\ \ 65\%", color=WHITE).scale(0.75).move_to(had)
        lab_lep = mathtex(r"\ell\nu\bar\nu\ \ 35\%", color=WHITE).scale(0.75).move_to(lep)
        return dict(had=had, lep=lep, lab_had=lab_had, lab_lep=lab_lep)

    def jet_parts(self) -> dict:
        def cone(o, half, length, color, op):
            return Polygon(o, o + length * unit(half), o + length * unit(-half),
                           fill_color=col(color), fill_opacity=op, stroke_color=col(color), stroke_width=1.5)
        o1, o2 = _p(2.4, 0.9), _p(2.4, -1.6)
        tau_cone = cone(o1, 0.10, 3.4, CHANNEL["tautau"], 0.25)
        jet_cone = cone(o2, 0.28, 3.4, PARTICLE["jet"], 0.25)
        tau_trk = VGroup(*[fline(o1, o1 + 3.2 * unit(a), color=TAU_COL, sw=2.5) for a in (0.04, -0.03)])
        jet_trk = VGroup(*[fline(o2, o2 + 3.2 * unit(a), color=PARTICLE["jet"], sw=2.5)
                           for a in (0.22, 0.12, 0.03, -0.08, -0.19, -0.25)])
        lab_tau = mathtex(r"\tau_h", color=TAU_COL).next_to(tau_cone, RIGHT, buff=0.15)
        lab_jet = mathtex(r"q,\,g").next_to(jet_cone, RIGHT, buff=0.15)
        fake = mathtex(r"\mathrm{jet}\to\tau_h", color=INK).scale(0.8).move_to(_p(4.2, -0.35))
        return dict(tau_cone=tau_cone, jet_cone=jet_cone, tau_trk=tau_trk, jet_trk=jet_trk,
                    lab_tau=lab_tau, lab_jet=lab_jet, fake=fake)

    def bkg_row(self) -> VGroup:
        items = [r"W\!+\!\mathrm{jets}", r"\mathrm{QCD\ multijet}", r"Z/\gamma^*\!\to\ell\ell",
                 r"t\bar t", r"VV"]
        g = VGroup(*[mathtex(s).scale(0.8) for s in items]).arrange(RIGHT, buff=0.6)
        return g.move_to(_p(2.7, -3.35))

    def construct(self):
        white_background(self)
        L = self.lifetime_parts()
        # 1) the tau flies 87 um and decays
        self.play(FadeIn(L["beam"]), FadeIn(L["pv"], scale=0.5), run_time=0.5)
        self.play(Create(L["others"]), run_time=0.6)
        self.play(Create(L["tau"]), FadeIn(L["lab_tau"]), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(L["tau_tip"]), FadeIn(L["sv"], scale=0.5), run_time=0.3)
        self.play(Flash(L["sv"].get_center(), color=SLATE, line_length=0.2, flash_radius=0.4, run_time=0.5))
        self.play(Create(L["ruler"]), FadeIn(L["lab_ct"]), FadeIn(L["lab_t"]), run_time=0.7)
        self.play(LaggedStart(*[Create(p) for p in L["prongs"]], lag_ratio=0.15), Create(L["nu"]),
                  run_time=1.0, rate_func=EASE)
        self.play(FadeIn(L["nu_lab"]), run_time=0.3)
        # the neutrino is invisible: it fades
        self.play(L["nu"].animate.set_stroke(opacity=0.25), run_time=0.6)
        # 2) the branching bar
        B = self.br_parts()
        self.play(GrowFromEdge(B["had"], LEFT), run_time=0.6)
        self.play(GrowFromEdge(B["lep"], LEFT), FadeIn(B["lab_had"]), run_time=0.5)
        self.play(FadeIn(B["lab_lep"]), run_time=0.3)
        # 3) a q/g jet looks like a narrow tau_h jet
        J = self.jet_parts()
        self.play(FadeIn(J["tau_cone"]), Create(J["tau_trk"]), FadeIn(J["lab_tau"]), run_time=0.8)
        self.play(FadeIn(J["jet_cone"]), Create(J["jet_trk"]), FadeIn(J["lab_jet"]), run_time=0.8)
        self.play(J["jet_cone"].animate.stretch(0.4, 1, about_point=J["jet_cone"].get_left()),
                  *[J["jet_trk"][j].animate.set_opacity(0.0) for j in (0, 1, 4, 5)],
                  FadeIn(J["fake"]), run_time=1.0, rate_func=EASE)
        # 4) the backgrounds
        row = self.bkg_row()
        self.play(LaggedStart(*[FadeIn(m, shift=UP * 0.15) for m in row], lag_ratio=0.2), run_time=1.2)
        self.wait(0.1)


# ---------------------------------------------------------------------------
# Slide 2: TautauChannels a/b/c
# ---------------------------------------------------------------------------
# one diagram in local units (tikz-feynman layout of the request), then scaled into a column
_V = dict(z_in=(0.0, 0.0), z=(1.5, 0.0), tm=(3.3, 1.2), tp=(3.3, -1.2),
          nu_m=(4.9, 1.8), w_m=(4.6, 0.7), a_m=(6.1, 1.25), b_m=(6.1, 0.45),
          nu_p=(4.9, -1.8), w_p=(4.6, -0.7), a_p=(6.1, -0.45), b_p=(6.1, -1.25))
DIAG_SCALE = 0.52
COL_X = (-3.5, 0.6, 4.7)       # column centres; left column clears x < -5.85
DIAG_Y = 0.95
KINDS = ("ll", "lh", "hh")
LEGS = {"ll": ("lep", "lep"), "lh": ("lep", "had"), "hh": ("had", "had")}
FINAL = {"ll": r"\tau_\ell\,\tau_\ell", "lh": r"\tau_\ell\,\tau_h", "hh": r"\tau_h\,\tau_h"}
N_NU = {"ll": 4, "lh": 3, "hh": 2}
BKG = {"ll": [r"t\bar t", r"VV", r"Z/\gamma^*\!\to ee,\mu\mu"],
       "lh": [r"Z/\gamma^*\!\to ee,\mu\mu", r"W\!+\!\mathrm{jets}", r"\mathrm{QCD}"],
       "hh": [r"\mathrm{QCD\ multijet}"]}
LS = 0.62   # label scale


def _map(k, i):
    x, y = _V[k]
    c = _p(COL_X[i], DIAG_Y)
    return c + DIAG_SCALE * _p(x - 3.05, y)


def diagram(kind: str, i: int) -> dict:
    P = lambda k: _map(k, i)  # noqa: E731
    sw = 3.2
    parts = {}
    parts["z"] = wavy(P("z_in"), P("z"), color=PARTICLE["Z"], sw=sw, amplitude=0.07, wavelength=0.3)
    parts["lab_z"] = mathtex(r"Z").scale(LS).next_to(P("z_in"), LEFT, buff=0.1)
    parts["v0"] = vertex_dot(P("z"), radius=0.045)
    tm, tp = fline(P("z"), P("tm"), color=TAU_COL, sw=sw), fline(P("z"), P("tp"), color=TAU_COL, sw=sw)
    tm_tip = arrow_tip_on(tm, color=TAU_COL, tip_length=0.15, at=0.55)
    tp_tip = arrow_tip_on(fline(P("tp"), P("z")), color=TAU_COL, tip_length=0.15, at=0.45)  # anti-fermion
    parts["taus"] = VGroup(tm, tp)
    parts["tau_tips"] = VGroup(tm_tip, tp_tip)
    parts["lab_taus"] = VGroup(
        mathtex(r"\tau^-", color=TAU_COL).scale(LS).next_to(tm.point_from_proportion(0.5), UP + LEFT, buff=0.04),
        mathtex(r"\tau^+", color=TAU_COL).scale(LS).next_to(tp.point_from_proportion(0.5), DOWN + LEFT, buff=0.04))
    legs, tips, labs, vtx = VGroup(), VGroup(), VGroup(), VGroup()
    nus = VGroup()
    for side, leg in zip(("m", "p"), LEGS[kind]):
        t, nu_end, w, a, b = P("t" + side), P("nu_" + side), P("w_" + side), P("a_" + side), P("b_" + side)
        top = side == "m"
        nu = dashed(t, nu_end, color=NU_COL, sw=sw - 0.5, dash_length=0.09)
        nus.add(nu)
        labs.add(mathtex(r"\nu_\tau" if top else r"\bar\nu_\tau", color=NU_COL).scale(LS).next_to(nu_end, RIGHT, buff=0.08))
        wl = wavy(t, w, color=PARTICLE["W"], sw=sw - 0.5, amplitude=0.05, wavelength=0.22)
        legs.add(wl)
        labs.add(mathtex(r"W^{-*}" if top else r"W^{+*}").scale(LS * 0.85).next_to(
            wl.get_center(), DOWN + LEFT if top else UP + LEFT, buff=0.02))
        vtx.add(vertex_dot(t, radius=0.04), vertex_dot(w, radius=0.04))
        if leg == "lep":
            # W-* -> l- nubar_l (top), W+* -> nu_l' l'+ (bottom)
            if top:
                la, lb = (r"\ell^-", INK, "line"), (r"\bar\nu_\ell", NU_COL, "nu")
            else:
                la, lb = (r"\nu_{\ell'}", NU_COL, "nu"), (r"\ell'^+", INK, "line")
        else:
            la, lb = ((r"q", INK, "line"), (r"\bar q'", INK, "line")) if top else \
                     ((r"q'", INK, "line"), (r"\bar q", INK, "line"))
        for end, (tex, c, style) in ((a, la), (b, lb)):
            if style == "nu":
                ln = dashed(w, end, color=NU_COL, sw=sw - 0.5, dash_length=0.09)
                nus.add(ln)
            else:
                ln = fline(w, end, color=c, sw=sw - 0.5)
                legs.add(ln)
            labs.add(mathtex(tex, color=c).scale(LS).next_to(end, RIGHT, buff=0.08))
    parts["legs"], parts["nus"], parts["vtx"], parts["labs"] = legs, nus, vtx, labs
    # column annotations
    parts["final"] = mathtex(FINAL[kind], color=TAU_COL).scale(0.95).move_to(_p(COL_X[i], 2.3))
    parts["frac"] = mathtex(rf"{FRACTIONS[kind]}\%").scale(1.1).move_to(_p(COL_X[i], -1.0))
    parts["nnu"] = mathtex(rf"{N_NU[kind]}\,\nu", color=NU_COL).scale(0.85).move_to(_p(COL_X[i], -1.75))
    bk = VGroup(*[mathtex(s).scale(0.62) for s in BKG[kind]]).arrange(DOWN, buff=0.14)
    parts["bkg"] = bk.next_to(_p(COL_X[i], -2.25), DOWN, buff=0.0)
    # share bar: width proportional to the fraction, under the percentage
    wmax = 3.6
    parts["bar"] = Rectangle(width=wmax * FRACTIONS[kind] / 46, height=0.16,
                             fill_color=col(CHANNEL["tautau"]), fill_opacity=0.85,
                             stroke_width=0).move_to(_p(COL_X[i], -0.55))
    return parts


def diagram_group(d: dict) -> VGroup:
    return VGroup(*d.values())


def end_state(k: int) -> VGroup:
    """Everything on screen at the end of clip k (1..3)."""
    return VGroup(*[diagram_group(diagram(KINDS[i], i)) for i in range(k)])


class TautauChannels(Scene):
    K = 1

    def play_diagram(self, d):
        self.play(FadeIn(d["final"], shift=DOWN * 0.15), run_time=0.4)
        self.play(FadeIn(d["lab_z"]), Create(d["z"]), run_time=0.6, rate_func=EASE)
        self.play(FadeIn(d["v0"]), Create(d["taus"]), run_time=0.7, rate_func=EASE)
        self.play(FadeIn(d["tau_tips"]), FadeIn(d["lab_taus"]), run_time=0.3)
        self.play(Create(d["legs"]), Create(d["vtx"]), run_time=0.9, rate_func=EASE)
        self.play(Create(d["nus"]), run_time=0.6, rate_func=EASE)
        self.play(FadeIn(d["labs"]), run_time=0.4)
        self.play(GrowFromEdge(d["bar"], DOWN), FadeIn(d["frac"], scale=0.8), run_time=0.6)
        self.play(Indicate(d["nus"], color=SLATE, scale_factor=1.0), FadeIn(d["nnu"]), run_time=0.7)
        self.play(LaggedStart(*[FadeIn(m, shift=UP * 0.1) for m in d["bkg"]], lag_ratio=0.25), run_time=0.7)

    def construct(self):
        white_background(self)
        if self.K > 1:
            self.add(end_state(self.K - 1))
            self.wait(0.2)
        self.play_diagram(diagram(KINDS[self.K - 1], self.K - 1))
        self.wait(0.1)


class TautauChannelsA(TautauChannels):
    K = 1


class TautauChannelsB(TautauChannels):
    K = 2


class TautauChannelsC(TautauChannels):
    K = 3


# ---------------------------------------------------------------------------
# Section opener: TauChallenges + TautauChannels a/b/c merged into one quicker clip
# ---------------------------------------------------------------------------
class TauIntro(Scene):
    """(5-01) The challenges slide at ~0.55x its run time, cleared, then the three final states
    drawn together (each step on all three columns at once, lagged left to right)."""

    def construct(self):
        white_background(self)
        TC = TauChallenges
        L = TC.lifetime_parts(self)
        self.play(FadeIn(L["beam"]), FadeIn(L["pv"], scale=0.5), Create(L["others"]), run_time=0.4)
        self.play(Create(L["tau"]), FadeIn(L["lab_tau"]), run_time=0.5, rate_func=EASE)
        self.play(FadeIn(L["tau_tip"]), FadeIn(L["sv"], scale=0.5),
                  Flash(L["sv"].get_center(), color=SLATE, line_length=0.2, flash_radius=0.4), run_time=0.35)
        self.play(Create(L["ruler"]), FadeIn(L["lab_ct"]), FadeIn(L["lab_t"]), run_time=0.45)
        self.play(LaggedStart(*[Create(p) for p in L["prongs"]], lag_ratio=0.15), Create(L["nu"]),
                  FadeIn(L["nu_lab"]), run_time=0.6, rate_func=EASE)
        self.play(L["nu"].animate.set_stroke(opacity=0.25), run_time=0.35)
        B = TC.br_parts(self)
        self.play(GrowFromEdge(B["had"], LEFT), GrowFromEdge(B["lep"], LEFT), run_time=0.45)
        self.play(FadeIn(B["lab_had"]), FadeIn(B["lab_lep"]), run_time=0.25)
        J = TC.jet_parts(self)
        self.play(FadeIn(J["tau_cone"]), Create(J["tau_trk"]), FadeIn(J["lab_tau"]),
                  FadeIn(J["jet_cone"]), Create(J["jet_trk"]), FadeIn(J["lab_jet"]), run_time=0.6)
        self.play(J["jet_cone"].animate.stretch(0.4, 1, about_point=J["jet_cone"].get_left()),
                  *[J["jet_trk"][j].animate.set_opacity(0.0) for j in (0, 1, 4, 5)],
                  FadeIn(J["fake"]), run_time=0.6, rate_func=EASE)
        row = TC.bkg_row(self)
        self.play(LaggedStart(*[FadeIn(m, shift=UP * 0.15) for m in row], lag_ratio=0.2), run_time=0.7)
        self.wait(0.3)
        # clear, then the three final states together
        self.play(FadeOut(VGroup(*L.values(), *B.values(), *J.values(), row)), run_time=0.4)
        D = [diagram(k, i) for i, k in enumerate(KINDS)]

        def step(anims, rt, lag=0.12):
            self.play(LaggedStart(*anims, lag_ratio=lag), run_time=rt)

        step([FadeIn(d["final"], shift=DOWN * 0.15) for d in D], 0.35)
        step([AnimationGroup(FadeIn(d["lab_z"]), Create(d["z"])) for d in D], 0.45)
        step([AnimationGroup(FadeIn(d["v0"]), Create(d["taus"]), FadeIn(d["tau_tips"]), FadeIn(d["lab_taus"]))
              for d in D], 0.5)
        step([AnimationGroup(Create(d["legs"]), Create(d["vtx"]), Create(d["nus"])) for d in D], 0.7)
        step([FadeIn(d["labs"]) for d in D], 0.3)
        step([AnimationGroup(GrowFromEdge(d["bar"], DOWN), FadeIn(d["frac"], scale=0.8)) for d in D], 0.5)
        step([AnimationGroup(Indicate(d["nus"], color=SLATE, scale_factor=1.0), FadeIn(d["nnu"])) for d in D], 0.6)
        step([FadeIn(d["bkg"], shift=UP * 0.1) for d in D], 0.5)
        self.wait(0.1)
