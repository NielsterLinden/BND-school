"""Section 2: how an analysis works, as one map with two parallel rows (schematic).

Re-cut of 17 Sep 2026 (user's notes; replaces the one-row pipeline of s2_pipeline.py, which is kept
as a library: its clips were removed, the map imports its icons, camera and masked_play):

    theory row (top, purple)       q qbar -> Z/gamma* -> l l  ->  [MC]  ->  simulated detector  ->  files
    experiment row (bottom)        the CMS slice -> L1 -> HLT -> RAW -> NanoAOD   (s2_trigger.py, 0.6x)
    the two NanoAOD files meet     -> [selection] -> [corrections] -> [fit -> sigma]

The selection node (a funnel, 17 Sep 2026 late, deck owner) is only passed in section 2; the Z -> ee chapter
(scenes/s3_zee_story.py) zooms into it and into the fit.

Chained scenes (every clip opens on the previous clip's final frame; pure builders + ORDER):

    MapBuild        map_build           opens on CMSSlice() = the last frame of 2-01 cms_logo_to_slice;
                                        the slice parks bottom left, the theory row draws, then the
                                        experiment row, then the merge
    MapTheory       theory_prediction   zoom onto the theory row (it becomes a full-size row along the bottom,
                                        a purple panel above, as in the trigger chain): the lineshape
                                        d sigma / d m_ll and its area sigma_theory
                    theory_generator    Monte Carlo: points thrown under the curve, one kept point is an
                                        l+ l- event, the kept points fill a histogram that follows the curve
                    theory_simulation   the event enters the (same) slice: purple tracks, a purple file
    TriggerInMap    trig_a ... trig_i   s2_trigger.TriggerChain with its enter/leave hooks: zoom from the
                                        map onto the experiment row = her spine; at the end the panel folds
                                        into NanoAOD and the camera zooms back out
    MapCorrections  corr_compare        both files on one m_ll axis: data points and the purple prediction
                    corr_factor         a control sample: the same quantity in data and in simulation, k = ratio
                    corr_apply          x k_1, x k_2, x k_3 move the prediction; the data never move
    MapFit          fit_model           data points, the purple line mu x prediction, the mu slider
                    fit_scan            mu scans, the line scales, -2 Delta ln L(mu) draws a parabola, mu-hat +- 1
                    fit_sigma           sigma = mu-hat sigma_theory
    MapThree        map_three           the map shrinks and repeats in the ee / mumu / tautau colours

Physics behind the pictures (symbols only on screen): simulation runs through the same detector
description and reconstruction into the same NanoAOD format (z-mumu/docs/10-skims.md); correction factors
come from data vs simulation in a control sample (tag and probe, docs/12; pile-up, docs/11) or from data
alone (fake factor, docs/13) and are applied to the prediction; the fit scales the prediction by mu,
sigma = mu sigma_pred (docs/CONVENTIONS.md section 2). All shapes are schematic (seeded); no result numbers.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    DOWN, LEFT, RIGHT, UP, AnimationGroup, Circle, Create, CurvedArrow, DashedLine, DashedVMobject, Dot,
    FadeIn, FadeOut, Flash, GrowFromCenter, GrowFromEdge, GrowFromPoint, LaggedStart, Line, Rectangle,
    ReplacementTransform, RoundedRectangle, Scene, Succession, Transform, VGroup, VMobject, ValueTracker,
    always_redraw, rate_functions,
)

from style.bnd_style import *  # noqa: E402,F401,F403
from s2_pipeline import P, camera, dy_icon, icon_fit, icon_funnel, masked_play, tex, win  # noqa: E402
from s2_trigger import (  # noqa: E402
    ORDER_TRIG, SPINE_Y as TRIG_Y, TriggerChain, arrow, file_glyph, node_file, stream, trigger_spine,
)

EASE = rate_functions.ease_in_out_sine
ACC = DETECTOR_ACCENT
TAIL = 0.1

# ---------------------------------------------------------------------------
# geometry (scene units). Title band y > 2.7 and the top-left block x < -5.85, y > 0.22 stay empty.
# Both rows are built at full size along TRIG_Y (= how they look zoomed) and mapped onto the map by
# the camera F(p) = to + Z_MAP (p - about); zooming onto a row is the inverse map.
# ---------------------------------------------------------------------------
Z_MAP = 0.6
B_ABOUT, B_TO = (-5.1, TRIG_Y), (-6.25, -1.95)       # experiment row: her slice centre -> map
T_X = (-3.8, -0.8, 2.0, 4.4)                         # theory row, full size (her L1/HLT/RAW/NanoAOD spacing)
T_ABOUT, T_TO = (T_X[0], TRIG_Y), (-4.57, 0.85)      # theory row -> map (columns over L1 / HLT / RAW / NanoAOD)
M_SEL, M_CORR, M_FIT, M_BOX = (2.85, -0.55), (4.45, -0.55), (6.05, -0.55), 1.0
ZOOM, DETAIL = 5.0, (0.0, -0.35)
T_PANEL, T_PANEL_WH = (0.3, 1.15), (10.2, 2.7)
T_SIM_R = 1.0

TOP = ("t_dy", "t_arr0", "t_mc", "t_arr1", "t_sim", "t_arr2", "t_file")
BOT = dict(zip(("b_det", "b_arr0", "b_l1", "b_arr1", "b_hlt", "b_arr2", "b_raw", "b_raw_lab", "b_arr3",
                "b_nano", "b_nano_lab"), ORDER_TRIG))
MERGE = ("m_link_t", "m_link_b", "m_sel", "m_arr_s", "m_corr", "m_arr", "m_fit")
ORDER_MAP = TOP + tuple(BOT) + MERGE
ALL = ("theory", "trigger", "sel", "corr", "fit")
CHANNELS = ("ee", "mumu", "tautau")
BOX_TINT = {"ee": 0.80, "mumu": 0.55, "tautau": 0.80}   # the pale gold line tint reads grey

# ---------------------------------------------------------------------------
# the schematic spectra (seeded, lepton-generic)
# ---------------------------------------------------------------------------
M_Z = 91.19


def shape(m):
    """Z peak with its radiative left tail plus a small gamma* continuum, height ~1."""
    m = np.asarray(m, dtype=float)
    peak = np.vectorize(lambda t: dscb(t, 0.9, 3.0, 1.6, 5.0))((m - M_Z) / 2.6)
    return peak + 0.015 * (M_Z / m) ** 2


EDGES = np.arange(60.0, 121.0, 2.0)
XC = 0.5 * (EDGES[:-1] + EDGES[1:])
MS = np.linspace(60.0, 120.0, 241)
TRUTH = shape(XC)
ERR = 0.035 * np.sqrt(np.maximum(TRUTH, 0.02))
DATA = TRUTH + ERR * np.random.default_rng(32).standard_normal(len(XC))
TILT = 1.0 + 0.06 * (XC - 90.0) / 30.0
PRED = [1.12 * TRUTH * TILT]                      # the uncorrected simulation: high and tilted
PRED.append(PRED[0] * 0.94)                       # x k_1 (a normalisation-like factor)
PRED.append(PRED[1] / TILT)                       # x k_2 (a shape factor)
PRED.append(PRED[2] * 0.978)                      # x k_3
FIT_NORM, MU_HAT, D_MU = 1.03, 0.97, 0.07         # final prediction = 1.03 x data -> mu-hat ~ 0.97
K_BARS = ((1.03, 1.10), (1.16, 1.10), (1.08, 1.10))   # (data, simulation) bar heights per factor


# ---------------------------------------------------------------------------
# row builders (full size)
# ---------------------------------------------------------------------------

def theory_box(center, done=False, line=THEORY, fill=THEORY, tint=0.85) -> VGroup:
    """The MC generator node: a box around a curve with points kept under it."""
    b = RoundedRectangle(width=1.2, height=1.2, corner_radius=0.12, stroke_color=col(line), stroke_width=2.5,
                         fill_color=lighten(fill, tint) if done else col(WHITE), fill_opacity=1.0).move_to(P(center))
    c = P(center)
    xs = np.linspace(-0.42, 0.42, 61)
    f = shape(60.0 + (xs + 0.42) / 0.84 * 60.0)
    curve = VMobject(stroke_color=col(THEORY), stroke_width=3.0)
    curve.set_points_as_corners([c + P((x, -0.30 + 0.62 * y)) for x, y in zip(xs, f)])
    base = Line(c + P((-0.44, -0.30)), c + P((0.44, -0.30)), stroke_color=col(INK), stroke_width=2.0)
    kept = VGroup(*[Dot(c + P(p), radius=0.035, color=col(THEORY))
                    for p in ((-0.03, -0.05), (0.02, 0.18), (0.08, -0.18), (-0.12, -0.24))])
    thrown = VGroup(*[Dot(c + P(p), radius=0.035, color=col(GREY)) for p in ((-0.30, 0.20), (0.28, 0.05), (0.33, 0.30))])
    g = VGroup(b, VGroup(base, curve, thrown, kept))
    g.box = b
    return g


def sim_slice(center) -> VGroup:
    """The simulated detector: the same slice inside a dashed purple ring."""
    det = mini_slice(T_SIM_R / R_DET, P(center))
    ring = DashedVMobject(Circle(radius=T_SIM_R + 0.12, arc_center=P(center), stroke_color=col(THEORY),
                                 stroke_width=2.5), num_dashes=36)
    g = VGroup(det, ring)
    g.det = det
    return g


def sim_file(center, w=0.42, h=0.56) -> VGroup:
    f = file_glyph(w, h, fill=lighten(THEORY, 0.78), stroke=THEORY)
    lines = VGroup(*[Line(P((-0.3 * w, y * h)), P((0.25 * w, y * h)), stroke_color=col(THEORY), stroke_width=1.4)
                     for y in (0.15, -0.05, -0.25)])
    return VGroup(f, lines).move_to(P(center))


def theory_row(done=False, line=THEORY, fill=THEORY, tint=0.85) -> dict:
    y = TRIG_Y
    dy = dy_icon((T_X[0], y), s=0.23, color=THEORY, sw=4.0)
    mc = theory_box((T_X[1], y), done, line, fill, tint)
    sim = sim_slice((T_X[2], y))
    f = sim_file((T_X[3], y))
    r = T_SIM_R + 0.12
    arrs = (arrow((dy.get_right()[0] + 0.1, y), (T_X[1] - 0.68, y), color=THEORY),
            arrow((T_X[1] + 0.68, y), (T_X[2] - r - 0.1, y), color=THEORY),
            arrow((T_X[2] + r + 0.1, y), (T_X[3] - 0.32, y), color=THEORY))
    return {"t_dy": dy, "t_arr0": arrs[0], "t_mc": mc, "t_arr1": arrs[1], "t_sim": sim, "t_arr2": arrs[2], "t_file": f}


def small_arrow(a, b, sw=1.8, tip=0.12) -> VGroup:
    ln = Line(P(a), P(b), stroke_color=col(INK), stroke_width=sw)
    return VGroup(ln, arrow_tip_on(ln, color=INK, at=1.0, tip_length=tip))


def merge_box(center, icon, done=False, line=ACC, fill=ACC, tint=0.80) -> VGroup:
    b = RoundedRectangle(width=M_BOX, height=M_BOX, corner_radius=0.1, stroke_color=col(line), stroke_width=2.0,
                         fill_color=lighten(fill, tint) if done else col(WHITE), fill_opacity=1.0).move_to(P(center))
    f = min(1.0, 0.78 * M_BOX / icon.width, 0.78 * M_BOX / icon.height)
    icon.scale(f).move_to(P(center))
    g = VGroup(b, icon)
    g.box = b
    return g


def icon_corr(line=ACC) -> VGroup:
    """Two bars (data grey, simulation purple) and the factor k."""
    base = Line(P((-0.40, -0.28)), P((0.10, -0.28)), stroke_color=col(INK), stroke_width=2.0)
    d = Rectangle(width=0.16, height=0.42, stroke_width=1.2, stroke_color=col(INK), fill_color=col(GREY),
                  fill_opacity=1.0).move_to(P((-0.26, -0.07)))
    s = Rectangle(width=0.16, height=0.50, stroke_width=1.2, stroke_color=col(THEORY), fill_color=lighten(THEORY, 0.5),
                  fill_opacity=1.0).move_to(P((-0.04, -0.03)))
    k = tex("k", h=0.30, color=line).move_to(P((0.30, 0.0)))
    return VGroup(base, d, s, k)


def merge_row(done=(), ch=None) -> dict:
    line = ACC if ch is None else CHANNEL_LINE[ch]
    fill = ACC if ch is None else CHANNEL[ch]
    tint = 0.80 if ch is None else BOX_TINT[ch]
    ex, ey = M_SEL[0] - M_BOX / 2 - 0.08, M_SEL[1]
    fx = T_TO[0] + Z_MAP * (T_X[3] - T_X[0]) + Z_MAP * 0.21 + 0.08       # right edge of the two NanoAOD files
    return {
        "m_link_t": small_arrow((fx, T_TO[1]), (ex, ey + 0.2)),
        "m_link_b": small_arrow((fx, B_TO[1]), (ex, ey - 0.2)),
        "m_sel": merge_box(M_SEL, icon_funnel(line, fill), "sel" in done, line, fill, tint),
        "m_arr_s": small_arrow((M_SEL[0] + M_BOX / 2 + 0.08, ey), (M_CORR[0] - M_BOX / 2 - 0.08, ey)),
        "m_corr": merge_box(M_CORR, icon_corr(line), "corr" in done, line, fill, tint),
        "m_arr": small_arrow((M_CORR[0] + M_BOX / 2 + 0.08, ey), (M_FIT[0] - M_BOX / 2 - 0.08, ey)),
        "m_fit": merge_box(M_FIT, icon_fit(line), "fit" in done, line, fill, tint),
    }


def map_state(done=(), ch=None) -> dict:
    """The map with the parts in ``done`` (theory, trigger, sel, corr, fit) tinted; ``ch`` recolours the
    boxes for one channel (the hand-off). Keys in z-order = ``ORDER_MAP``."""
    line = THEORY if ch is None else CHANNEL_LINE[ch]
    fill = THEORY if ch is None else CHANNEL[ch]
    tint = 0.85 if ch is None else BOX_TINT[ch]
    st = {k: camera(v, Z_MAP, T_ABOUT, T_TO)
          for k, v in theory_row("theory" in done or ch is not None, line, fill, tint).items()}
    sp = trigger_spine()
    for k in ("l1", "hlt"):
        b = sp[k].box
        if ch is not None:
            b.set_stroke(col(CHANNEL_LINE[ch])).set_fill(lighten(CHANNEL[ch], BOX_TINT[ch]), opacity=1.0)
        elif "trigger" in done:
            b.set_fill(lighten(ACC, 0.8), opacity=1.0)
    st.update({k: camera(sp[v], Z_MAP, B_ABOUT, B_TO) for k, v in BOT.items()})
    st.update(merge_row(done, ch))
    return st


def to_full(mob, row):
    """Map -> full-size row (the zoom-in camera) for row 'theory' | 'trigger'. In place."""
    about, to = (T_TO, T_ABOUT) if row == "theory" else (B_TO, B_ABOUT)
    return camera(mob, 1.0 / Z_MAP, about, to)


def to_map(mob, row):
    about, to = (T_ABOUT, T_TO) if row == "theory" else (B_ABOUT, B_TO)
    return camera(mob, Z_MAP, about, to)


def reset_state(scene, st, order) -> None:
    """Drop leftovers (finished Successions, faded helpers) and re-add the state in z-order."""
    scene.remove(*scene.mobjects)
    add_state(scene, st, order)
    check_order(scene, st, order)


# ---------------------------------------------------------------------------
# the chained scene base
# ---------------------------------------------------------------------------

class MapScene(Scene):
    DONE: tuple = ()

    def open_map(self, clip: str):
        white_background(self)
        self.st = map_state(self.DONE)
        add_state(self, self.st, ORDER_MAP)
        clip_open(self, clip)

    def zoom_in(self, z, about, to, keep=(), detail=None, run_time=1.4):
        anims = []
        for k in ORDER_MAP:
            tgt = camera(self.st[k].copy(), z, about, to)
            if k not in keep:
                tgt.set_opacity(0.0)
            anims.append(Transform(self.st[k], tgt))
        if detail is not None:
            anims.append(ReplacementTransform(camera(detail.copy(), 1.0 / z, to, about).set_opacity(0.0), detail))
        masked_play(self, *anims, run_time=run_time, rate_func=EASE)

    def zoom_out(self, done, z, about, to, run_time=1.4):
        """The inverse of ``zoom_in(z, about, to)``: the map returns to ``map_state(done)``; everything
        that is not the map is carried back through the inverse camera and fades."""
        end = map_state(done)
        keep = {id(v) for v in self.st.values()}
        det = [m for m in self.mobjects
               if id(m) not in keep and not isinstance(m, ValueTracker) and m.family_members_with_points()]
        for m in det:
            for f in m.get_family():
                f.clear_updaters()
        anims = [Transform(self.st[k], end[k]) for k in ORDER_MAP]
        if det:
            g = VGroup(*det)
            anims.append(Transform(g, camera(g.copy(), 1.0 / z, to, about).set_opacity(0.0), remover=True))
        masked_play(self, *anims, run_time=run_time, rate_func=EASE)
        reset_state(self, self.st, ORDER_MAP)
        self.wait(TAIL)


# ---------------------------------------------------------------------------
# A: the map
# ---------------------------------------------------------------------------

class MapBuild(Scene):
    """map_build (~9 s). Opens on the full slice (last frame of 2-01)."""

    def construct(self):
        white_background(self)
        big = CMSSlice()
        self.add(big)
        clip_open(self, "map_build")
        st = map_state(())
        self.play(ReplacementTransform(big, st["b_det"]), run_time=1.5, rate_func=EASE)

        def grow(parts, t0, step, T, arrow_keys):
            anims = []
            for i, k in enumerate(parts):
                t = t0 + step * i
                m = st[k]
                if k in arrow_keys:
                    anims += [Create(m[0], lag_ratio=0.0, rate_func=win(t, 0.45, T, rate_functions.linear)),
                              FadeIn(m[1], rate_func=win(t + 0.35, 0.2, T))]
                else:
                    anims.append(FadeIn(m, scale=0.6, rate_func=win(t, 0.6, T)))
            return anims

        arrows = {"t_arr0", "t_arr1", "t_arr2", "b_arr0", "b_arr1", "b_arr2", "b_arr3", "m_link_t", "m_link_b", "m_arr_s",
                  "m_arr"}
        T = 3.0
        self.play(*grow(TOP, 0.05, 0.40, T, arrows), run_time=T)
        bot = [k for k in BOT if k != "b_det"]
        T = 3.0
        self.play(*grow(bot, 0.05, 0.26, T, arrows), run_time=T)
        T = 2.4
        self.play(*grow(MERGE, 0.05, 0.28, T, arrows), run_time=T)
        reset_state(self, st, ORDER_MAP)
        self.wait(TAIL)


# ---------------------------------------------------------------------------
# B: the theory row
# ---------------------------------------------------------------------------

def t_panel() -> RoundedRectangle:
    w, h = T_PANEL_WH
    return RoundedRectangle(width=w, height=h, corner_radius=0.18, stroke_color=col(THEORY), stroke_width=2.5,
                            fill_color=lighten(THEORY, 0.94), fill_opacity=1.0).move_to(P(T_PANEL))


def t_callout(center, top, half=0.6) -> VGroup:
    px, py = T_PANEL
    pw, ph = T_PANEL_WH
    x = center
    xl, xr = max(px - pw / 2 + 0.2, x - 1.6), min(px + pw / 2 - 0.2, x + 1.6)
    return VGroup(*[Line(P((x + s * half, top)), P((xe, py - ph / 2)), stroke_color=col(THEORY), stroke_width=1.8)
                    for s, xe in ((-1, xl), (1, xr))])


PLOT_T = (-1.3, 1.42)


def with_y_title(dax: DataAxes, expr: str, h: float, buff: float = 0.18) -> DataAxes:
    """Rotated y title beside the plot frame (DataAxes anchors its own on the tick labels, which a
    schematic axis does not have, so it would land inside the plot)."""
    dax.add(tex(expr, h=h).rotate(np.pi / 2).next_to(dax.frame, LEFT, buff=buff))
    return dax


def theory_axes() -> DataAxes:
    dax = DataAxes([60, 120, 10], [0, 1.15, 0.5], 4.8, 1.6, x_ticks=[70, 80, 90, 100, 110], y_ticks=[],
                   show_y_labels=False, tick_label_h=0.16, x_title=r"m_{\ell\ell}\ [\mathrm{GeV}]",
                   title_h=0.2, title_buff=0.1).move_frame_to(P(PLOT_T))
    return with_y_title(dax, r"d\sigma / dm_{\ell\ell}", h=0.2)


def event_glyph(center) -> VGroup:
    c = P(center)
    a1, a2 = 0.35, 0.35 + np.pi + 0.3
    v = Dot(c, radius=0.06, color=col(INK))
    arrs = VGroup()
    labs = VGroup()
    for a, lab in ((a1, r"\ell^{-}"), (a2, r"\ell^{+}")):
        end = c + 0.95 * np.array([np.cos(a), np.sin(a), 0.0])
        ln = Line(c, end, stroke_color=col(THEORY), stroke_width=4.0)
        arrs.add(VGroup(ln, arrow_tip_on(ln, color=THEORY, at=1.0, tip_length=0.18)))
        labs.add(tex(lab, h=0.26).move_to(c + 1.25 * np.array([np.cos(a), np.sin(a), 0.0])))
    return VGroup(arrs, v, labs)


class MapTheory(MapScene):
    """theory_prediction, theory_generator, theory_simulation (~15 s), then the zoom back out
    (opens trig_a_collisions)."""
    DONE = ()

    def construct(self):
        self.open_map("theory_prediction")
        self.zoom_in(1.0 / Z_MAP, T_TO, T_ABOUT, keep=TOP)
        s = self.st
        dy, mc, sim, fnode = s["t_dy"], s["t_mc"], s["t_sim"], s["t_file"]

        # -- the prediction: lineshape and its area ---------------------------------------------
        pan = t_panel()
        call = t_callout(T_X[0], dy.get_top()[1] + 0.05, half=0.8)
        self.play(FadeIn(pan, scale=0.92), Create(call), run_time=0.7, rate_func=EASE)
        dax = theory_axes()
        curve = data_trace(dax, MS, shape(MS), color=THEORY, stroke_width=4.0)
        area = data_band(dax, MS, np.zeros_like(MS), shape(MS), color=THEORY, opacity=0.22)
        self.play(FadeIn(dax), run_time=0.5)
        self.play(Create(curve), run_time=1.3, rate_func=rate_functions.linear)
        sig = tex(r"\sigma_{\mathrm{theory}}", h=0.34, color=THEORY).move_to(dax.c2p(109.0, 0.62))
        self.add(area)
        self.bring_to_front(curve)
        area.set_opacity(0.0)
        self.play(area.animate.set_fill(col(THEORY), opacity=0.22), FadeIn(sig, shift=LEFT * 0.15), run_time=0.8)

        # -- the generator: hit or miss ---------------------------------------------------------
        clip_cut(self, "theory_generator")
        call2 = t_callout(T_X[1], mc.get_top()[1], half=0.6)
        self.play(ReplacementTransform(call, call2), mc.box.animate.set_fill(lighten(THEORY, 0.85)),
                  FadeOut(area), run_time=0.6, rate_func=EASE)
        rng = np.random.default_rng(41)
        pts = np.column_stack([rng.uniform(60, 120, 90), rng.uniform(0.0, 1.08, 90)])
        pts[0] = (91.5, 0.45)                                        # the first throw is kept: the event
        kept = [bool(y < shape(x)) for x, y in pts]
        dots = [Dot(dax.c2p(x, y), radius=0.04, color=col(THEORY if k else GREY)) for (x, y), k in zip(pts, kept)]
        first = dots[0]
        self.play(FadeIn(first, scale=2.0), run_time=0.3)
        glyph = event_glyph((3.9, 1.25))
        self.play(Flash(first.get_center(), color=col(THEORY), line_length=0.12, flash_radius=0.16, num_lines=8),
                  GrowFromPoint(glyph, first.get_center()), run_time=0.9, rate_func=EASE)
        self.play(LaggedStart(*[FadeIn(d, scale=1.6) for d in dots[1:]], lag_ratio=0.06), run_time=1.6)
        hist0 = step_hist(dax, EDGES, np.zeros(len(XC)), color=THEORY, stroke_width=3.0, fill_opacity=0.3)
        hist = step_hist(dax, EDGES, shape(XC), color=THEORY, stroke_width=3.0, fill_opacity=0.3)
        self.add(hist0)
        drops = []
        for d, (x, _), k in zip(dots, pts, kept):
            if k:
                j = min(int((x - 60.0) // 2.0), len(XC) - 1)
                drops.append(Succession(d.animate(run_time=0.5, rate_func=rate_functions.ease_in_quad)
                                        .move_to(dax.c2p(x, shape(XC)[j])), FadeOut(d, run_time=0.15)))
            else:
                drops.append(FadeOut(d, run_time=0.4))
        self.play(LaggedStart(*drops, lag_ratio=0.004), Transform(hist0, hist, rate_func=win(0.3, 0.9, 1.3)),
                  run_time=1.3)
        self.bring_to_front(curve, sig)

        # -- the simulated detector ---------------------------------------------------------------
        clip_cut(self, "theory_simulation")
        panel_stuff = VGroup(pan, call2, dax, curve, sig, hist0)
        c0 = sim.get_center()
        self.play(FadeOut(panel_stuff), glyph.animate.scale(0.25).move_to(c0), run_time=0.9, rate_func=EASE)
        geo = mini_slice(T_SIM_R / R_DET, c0)                          # geometry only
        r_trk = geo.radii["tob"][1]
        trks = VGroup(*[track(geo, phi, k, r_end=r_trk, color=THEORY, sw=4.0, outline=WHITE, outline_sw=2.5)
                        for phi, k in ((0.35, 1.2), (0.35 + np.pi + 0.3, -1.2))])
        self.play(FadeOut(glyph, scale=0.3), Flash(c0, color=col(THEORY), line_length=0.12, flash_radius=0.2,
                                                   num_lines=10), run_time=0.35)
        self.play(Create(trks, lag_ratio=0.0), run_time=0.45, rate_func=rate_functions.linear)
        fly = sim_file(c0, 0.42, 0.56).scale(0.3)
        self.play(fly.animate.scale(1 / 0.3).move_to(fnode.get_center()).set_opacity(0.0), FadeOut(trks),
                  run_time=0.8, rate_func=EASE)
        a1, a2 = s["t_arr1"][0], s["t_arr2"][0]
        self.play(stream(a1.get_start(), a1.get_end(), 8, 1.6, seed=12, color=THEORY),
                  stream(a2.get_start(), a2.get_end(), 8, 1.6, seed=13, color=THEORY))
        clip_cut(self, "trig_a_collisions")
        self.zoom_out(("theory",), 1.0 / Z_MAP, T_TO, T_ABOUT)
        check_order(self, self.st, ORDER_MAP)


# ---------------------------------------------------------------------------
# C: the experiment row = the trigger chain (s2_trigger.py), zoomed in from the map
# ---------------------------------------------------------------------------

class TriggerInMap(TriggerChain):
    """trig_a_collisions ... trig_i_nanoaod: Ella's chain, entered from and left to the map."""

    def enter(self, st: dict) -> None:
        white_background(self)
        self.spine = st
        mp = map_state(("theory",))
        add_state(self, mp, ORDER_MAP)
        clip_open(self, "trig_a_collisions")
        anims = []
        for k in ORDER_MAP:
            if k in BOT:
                anims.append(ReplacementTransform(mp[k], st[BOT[k]]))
            else:
                anims.append(Transform(mp[k], to_full(mp[k].copy(), "trigger").set_opacity(0.0), remover=True))
        masked_play(self, *anims, run_time=1.5, rate_func=EASE)

    def leave(self) -> None:
        clip_cut(self, "corr_compare")
        st = self.spine
        keep = {id(v) for v in st.values()}
        left = [m for m in self.mobjects
                if id(m) not in keep and not isinstance(m, ValueTracker) and m.family_members_with_points()]
        for m in left:
            for f in m.get_family():
                f.clear_updaters()
        up = [m for m in left if m.get_center()[1] > -0.45]           # the panel and what is in it
        up_ids = {id(m) for m in up}
        rest = [m for m in left if id(m) not in up_ids]                # rates, clocks, callout
        nano_c = st["nano_icon"].get_center()
        fold = VGroup(*up)
        anims = [Transform(fold, fold.copy().scale(0.04).move_to(nano_c).set_opacity(0.0), remover=True)]
        if rest:
            anims.append(FadeOut(VGroup(*rest)))
        self.play(*anims, run_time=0.9, rate_func=EASE)
        self.remove(*left)
        end = map_state(("theory", "trigger", "sel"))
        anims = []
        for k in ORDER_MAP:
            if k in BOT:
                anims.append(ReplacementTransform(st[BOT[k]], end[k]))
            else:
                anims.append(ReplacementTransform(to_full(end[k].copy(), "trigger").set_opacity(0.0), end[k]))
        masked_play(self, *anims, run_time=1.5, rate_func=EASE)
        reset_state(self, end, ORDER_MAP)
        self.wait(TAIL)


# ---------------------------------------------------------------------------
# D: corrections
# ---------------------------------------------------------------------------
PLOT_C, PLOT_C_WH = (1.75, -0.05), (6.7, 3.5)
CR_BOX, CR_BOX_WH = (-4.2, -0.2), (2.6, 3.3)
BAR_X, BAR_BASE, BAR_W, BAR_SCALE = (-4.72, -3.68), -0.05, 0.62, 1.0


def plot_axes(center, wh) -> DataAxes:
    dax = DataAxes([60, 120, 10], [0, 1.28, 0.5], *wh, x_ticks=[70, 80, 90, 100, 110], y_ticks=[],
                   show_y_labels=False, tick_label_h=0.18, x_title=r"m_{\ell\ell}\ [\mathrm{GeV}]",
                   title_h=0.22, title_buff=0.14).move_frame_to(P(center))
    return with_y_title(dax, r"\mathrm{events}", h=0.2)


def data_points(dax, radius=0.045) -> VGroup:
    dots = VGroup(*[data_dot(dax, x, y, color=SAMPLE["Data"], radius=radius) for x, y in zip(XC, DATA)])
    bars = VGroup(*[data_errorbar(dax, x, y, e, color=SAMPLE["Data"], stroke_width=2.0)
                    for x, y, e in zip(XC, DATA, ERR)])
    g = VGroup(bars, dots)
    g.dots, g.bars = dots, bars
    return g


def pred_hist(dax, i) -> VMobject:
    return step_hist(dax, EDGES, PRED[i], color=THEORY, stroke_width=3.0, fill_opacity=0.22)


def k_bars(i) -> VGroup:
    d, s = K_BARS[i]
    out = VGroup()
    for x, hgt, fill, line in ((BAR_X[0], d, GREY, INK), (BAR_X[1], s, lighten(THEORY, 0.5), THEORY)):
        out.add(Rectangle(width=BAR_W, height=hgt * BAR_SCALE, stroke_color=col(line), stroke_width=2.0,
                          fill_color=col(fill), fill_opacity=1.0).move_to(P((x, BAR_BASE + hgt * BAR_SCALE / 2))))
    return out


def k_fraction(i) -> VGroup:
    """k_i = (grey bar) / (purple bar): the two bars shrunk into a fraction."""
    y = -1.22
    lab = tex(r"k_{%d} =" % (i + 1), h=0.34, color=ACC).move_to(P((-4.85, y)))
    bar = Line(P((-4.25, y)), P((-3.35, y)), stroke_color=col(INK), stroke_width=2.5)
    d, s = K_BARS[i]
    num = Rectangle(width=0.5, height=0.22 * d, stroke_color=col(INK), stroke_width=1.5, fill_color=col(GREY),
                    fill_opacity=1.0).next_to(bar, UP, buff=0.07)
    den = Rectangle(width=0.5, height=0.22 * s, stroke_color=col(THEORY), stroke_width=1.5,
                    fill_color=lighten(THEORY, 0.5), fill_opacity=1.0).next_to(bar, DOWN, buff=0.07)
    g = VGroup(lab, bar, num, den)
    g.lab, g.bar, g.num, g.den = lab, bar, num, den
    return g


class MapCorrections(MapScene):
    """corr_compare, corr_factor, corr_apply (~14 s), then the zoom out (opens fit_model)."""
    DONE = ("theory", "trigger", "sel")

    def construct(self):
        self.open_map("corr_compare")
        dax = plot_axes(PLOT_C, PLOT_C_WH)
        self.zoom_in(ZOOM, M_CORR, DETAIL, detail=VGroup(dax))
        f_sim = sim_file((-4.3, 1.0), 0.75, 1.0)
        f_dat = node_file((-4.3, -1.3), 0.75, 1.0)
        self.play(FadeIn(f_sim, shift=RIGHT * 0.5 + DOWN * 0.3), FadeIn(f_dat, shift=RIGHT * 0.5 + UP * 0.3),
                  run_time=0.7, rate_func=EASE)
        pts = data_points(dax)
        src = f_dat.get_center()
        flights = [ReplacementTransform(Dot(src, radius=0.045, color=col(SAMPLE["Data"])), d) for d in pts.dots]
        self.play(LaggedStart(*flights, lag_ratio=0.05), run_time=1.3)
        self.play(FadeIn(pts.bars), run_time=0.4)
        hist = pred_hist(dax, 0)
        self.play(GrowFromEdge(hist, DOWN), f_sim.animate.shift(RIGHT * 0.12), run_time=0.9, rate_func=EASE)
        self.bring_to_front(pts)
        self.play(f_sim.animate.shift(LEFT * 0.12), run_time=0.2)

        # -- a control sample: the same quantity in data and in simulation ------------------------
        clip_cut(self, "corr_factor")
        box = DashedVMobject(RoundedRectangle(width=CR_BOX_WH[0], height=CR_BOX_WH[1], corner_radius=0.2,
                                              stroke_color=col(ACC), stroke_width=3.0).move_to(P(CR_BOX)),
                             num_dashes=48)
        base = Line(P((BAR_X[0] - 0.5, BAR_BASE)), P((BAR_X[1] + 0.5, BAR_BASE)), stroke_color=col(INK), stroke_width=2.5)
        self.play(FadeOut(f_sim, shift=LEFT * 0.3), FadeOut(f_dat, shift=LEFT * 0.3), Create(box), FadeIn(base),
                  run_time=0.8, rate_func=EASE)
        bars = k_bars(0)
        self.play(*[GrowFromEdge(b, DOWN) for b in bars], run_time=0.8, rate_func=EASE)
        frac = k_fraction(0)
        self.play(ReplacementTransform(bars[0].copy(), frac.num), ReplacementTransform(bars[1].copy(), frac.den),
                  Create(frac.bar), FadeIn(frac.lab, shift=RIGHT * 0.15), run_time=1.0, rate_func=EASE)

        # -- apply: the prediction moves, the data do not ------------------------------------------
        clip_cut(self, "corr_apply")
        peak = dax.c2p(M_Z + 7.0, 0.95)
        for i in range(3):
            if i > 0:
                nb, nf = k_bars(i), k_fraction(i)
                self.play(Transform(bars, nb), Transform(frac, nf), run_time=0.5, rate_func=EASE)
            chip = tex(r"\times k_{%d}" % (i + 1), h=0.34, color=ACC).move_to(P((CR_BOX[0], -2.25)))
            self.play(FadeIn(chip, shift=UP * 0.15), run_time=0.3)
            arc = CurvedArrow(chip.get_right() + RIGHT * 0.1, peak, angle=0.5, color=col(ACC), stroke_width=3.0)
            self.play(Create(arc), chip.animate.move_to(peak + UP * 0.28 + RIGHT * 0.35), run_time=0.6, rate_func=EASE)
            self.play(Transform(hist, pred_hist(dax, i + 1)), FadeOut(arc), FadeOut(chip, scale=0.6),
                      run_time=0.7, rate_func=EASE)
            self.bring_to_front(pts)
        clip_cut(self, "fit_model")
        self.zoom_out(("theory", "trigger", "sel", "corr"), ZOOM, M_CORR, DETAIL)
        check_order(self, self.st, ORDER_MAP)


# ---------------------------------------------------------------------------
# E: the fit (one parameter)
# ---------------------------------------------------------------------------
PLOT_F, PLOT_F_WH = (-1.95, 0.05), (6.1, 3.6)
SL_X, SL_Y, SL_LEN = 4.1, 1.85, 3.4
NLL_C, NLL_WH, NLL_TOP = (4.1, -0.45), (3.4, 2.1), 6.5


def nll(mu):
    return ((np.asarray(mu, dtype=float) - MU_HAT) / D_MU) ** 2


def fit_line(dax, mu) -> VMobject:
    return data_trace(dax, MS, mu * FIT_NORM * shape(MS), color=THEORY, stroke_width=4.0)


class MapFit(MapScene):
    """fit_model, fit_scan, fit_sigma (~11 s), then the zoom out (opens map_three)."""
    DONE = ("theory", "trigger", "sel", "corr")

    def construct(self):
        self.open_map("fit_model")
        dax = plot_axes(PLOT_F, PLOT_F_WH)
        self.zoom_in(ZOOM, M_FIT, DETAIL, detail=VGroup(dax))
        pts = data_points(dax)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in pts.dots], lag_ratio=0.05), FadeIn(pts.bars), run_time=1.0)
        mu = ValueTracker(1.0)
        line = fit_line(dax, 1.0)
        lab = tex(r"\mu\,\sigma_{\mathrm{theory}}", h=0.32, color=THEORY).move_to(dax.c2p(106.0, 0.78))
        self.play(Create(line), FadeIn(lab, shift=LEFT * 0.15), run_time=1.0, rate_func=rate_functions.linear)
        self.bring_to_front(pts)
        sl = slider(0.8, 1.2, 1.0, ref=1.0, length=SL_LEN, ticks=[0.8, 1.0, 1.2], fmt="{:.1f}", marker_r=0.09)
        sl.shift(P((SL_X, SL_Y)) - sl.axis.get_center())
        mu_lab = tex(r"\mu", h=0.34).next_to(sl.axis, LEFT, buff=0.3)
        self.play(FadeIn(sl), FadeIn(mu_lab), run_time=0.6)

        # -- the scan ---------------------------------------------------------------------------------
        clip_cut(self, "fit_scan")
        nax = DataAxes([0.8, 1.2, 0.1], [0, NLL_TOP, 1], *NLL_WH, x_ticks=[], y_ticks=[], show_x_labels=False,
                       show_y_labels=False)
        with_y_title(nax.move_frame_to(P(NLL_C)), r"-2\,\Delta\ln L", h=0.22, buff=0.3)
        self.play(FadeIn(nax), run_time=0.4)
        self.remove(line)
        live = always_redraw(lambda: fit_line(dax, mu.get_value()))
        self.add(live)
        self.bring_to_front(pts)
        sl.marker.add_updater(lambda m: m.become(sl.marker_at(mu.get_value())))
        mu_lo, mu_hi = 0.80, MU_HAT + D_MU * np.sqrt(NLL_TOP - 0.6)
        self.play(mu.animate.set_value(mu_lo), run_time=0.7, rate_func=EASE)
        grid = np.linspace(mu_lo, mu_hi, 121)
        par = data_trace(nax, grid, nll(grid), color=INK, stroke_width=3.5)
        ball = always_redraw(lambda: Dot(nax.c2p(mu.get_value(), min(float(nll(mu.get_value())), NLL_TOP)),
                                         radius=0.07, color=col(ACC)))
        self.add(ball)
        self.play(mu.animate.set_value(mu_hi), Create(par, lag_ratio=0.0), run_time=1.8, rate_func=rate_functions.linear)
        self.play(mu.animate.set_value(MU_HAT), run_time=1.0, rate_func=EASE)
        sl.marker.clear_updaters()
        live.clear_updaters()
        ball.clear_updaters()
        self.remove(mu)
        lo, hi = MU_HAT - D_MU, MU_HAT + D_MU
        one = DashedLine(nax.c2p(0.8, 1.0), nax.c2p(1.2, 1.0), dash_length=0.08, stroke_color=col(GREY), stroke_width=2.0)
        drops = VGroup(*[DashedLine(nax.c2p(v, 1.0), P((sl.x_of(v), sl.axis.get_center()[1])), dash_length=0.08,
                                    stroke_color=col(ACC), stroke_width=2.0) for v in (lo, hi)])
        one_lab = tex("1", h=0.22, color=GREY).next_to(nax.c2p(0.8, 1.0), LEFT, buff=0.12)
        self.play(Create(one), FadeIn(one_lab), run_time=0.5)
        self.play(Create(drops), Transform(sl.marker, sl.marker_at(MU_HAT, D_MU)), run_time=0.8, rate_func=EASE)

        # -- sigma ----------------------------------------------------------------------------------
        clip_cut(self, "fit_sigma")
        s = tex(r"\sigma = \hat{\mu}\;\sigma_{\mathrm{theory}}", h=0.42,
                substrings_to_isolate=[r"\sigma_{\mathrm{theory}}"]).move_to(P((SL_X, -2.45)))
        s.get_part_by_tex(r"\sigma_{\mathrm{theory}}").set_color(col(THEORY))
        self.play(FadeIn(s, shift=UP * 0.2), FadeOut(drops), run_time=0.8, rate_func=EASE)
        clip_cut(self, "map_three")
        self.zoom_out(ALL, ZOOM, M_FIT, DETAIL)
        check_order(self, self.st, ORDER_MAP)


# ---------------------------------------------------------------------------
# F: three channels
# ---------------------------------------------------------------------------
H_Z, H_ABOUT, H_X = 0.42, (-0.2, -0.55), 1.05
H_YS = {"ee": 1.35, "mumu": -0.55, "tautau": -2.45}
H_LAB_X = -1.95
CHANNEL_TEX = {"ee": r"e^{+}e^{-}", "mumu": r"\mu^{+}\mu^{-}", "tautau": r"\tau^{+}\tau^{-}"}


def thumb(st: dict, at: str) -> dict:
    return {k: camera(v, H_Z, H_ABOUT, (H_X, H_YS[at])) for k, v in st.items()}


def state_three() -> dict:
    out = {}
    for ch in CHANNELS:
        out.update({f"{ch}:{k}": v for k, v in thumb(map_state(ALL, ch), ch).items()})
    for ch in CHANNELS:
        out[f"lab_{ch}"] = tex(CHANNEL_TEX[ch], h=0.34, color=CHANNEL_LINE[ch]).next_to(
            P((H_LAB_X, H_YS[ch])), LEFT, buff=0.0)
    return out


ORDER_THREE = tuple(f"{ch}:{k}" for ch in CHANNELS for k in ORDER_MAP) + tuple(f"lab_{ch}" for ch in CHANNELS)


class MapThree(Scene):
    """map_three (~6 s): the finished map shrinks, repeats three times and takes the channel colours."""

    def construct(self):
        white_background(self)
        st = map_state(ALL)
        add_state(self, st, ORDER_MAP)
        clip_open(self, "map_three")
        plain = thumb(map_state(ALL), "ee")
        masked_play(self, *[Transform(st[k], plain[k]) for k in ORDER_MAP], run_time=1.5, rate_func=EASE)
        live = {"ee": st}
        for ch in ("mumu", "tautau"):
            live[ch] = {k: v.copy() for k, v in st.items()}
        for ch in ("mumu", "tautau"):
            add_state(self, live[ch], ORDER_MAP)
        tgt = {ch: thumb(map_state(ALL), ch) for ch in ("mumu", "tautau")}
        self.play(*[Transform(live[ch][k], tgt[ch][k]) for ch in ("mumu", "tautau") for k in ORDER_MAP],
                  run_time=1.3, rate_func=EASE)
        end = state_three()
        T = 2.6
        xs = {k: st[k].get_center()[0] for k in ORDER_MAP}
        x0, x1 = min(xs.values()), max(xs.values())
        anims = []
        for r, ch in enumerate(CHANNELS):
            for k in ORDER_MAP:
                t0 = 0.3 * r + 1.4 * (xs[k] - x0) / (x1 - x0)
                anims.append(Transform(live[ch][k], end[f"{ch}:{k}"], rate_func=win(t0, 0.5, T)))
        labs = [FadeIn(end[f"lab_{ch}"], shift=RIGHT * 0.2, rate_func=win(0.3 * r + 0.1, 0.6, T))
                for r, ch in enumerate(CHANNELS)]
        self.play(*anims, *labs, run_time=T)
        final = {f"{ch}:{k}": live[ch][k] for ch in CHANNELS for k in ORDER_MAP}
        final.update({f"lab_{ch}": end[f"lab_{ch}"] for ch in CHANNELS})
        reset_state(self, final, ORDER_THREE)
        self.wait(TAIL)
