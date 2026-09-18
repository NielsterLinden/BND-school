"""Section 2 tail: how an analysis works, as one pipeline spine (schematic; the stacks are real).

Eight chained plain Scenes. ``PipeMap`` opens on an empty white frame; every
later clip opens on the previous clip's final frame (same pure builder, same
ORDER; joints checked with tools/framediff.py):

    PipeMap          pipe_a_map          the seven nodes appear left to right, the
                                         simulated river drops into the files node, sigma last
    PipeData         pipe_b_data         zoom on nodes 1-2: collisions flash in the slice and
                                         stream into a file pile that is squeezed (skim); the
                                         simulated river (Drell-Yan -> the same slice -> files) merges in
    PipeSelection    pipe_c_selection    the funnel: one bar shrinks through trigger, lepton
                                         ID/iso, exactly two leptons, opposite sign, mass window
    PipeTnP          pipe_d_tnp          tag and probe: the probe is tested, pass/fail spectra
                                         grow, eps = pass/(pass+fail), data vs simulation -> SF slider
    PipeBackgrounds  pipe_e_backgrounds  simulated processes stack up; the fakes come from a
                                         same-sign control region via a fake factor (x f): a sliver at the bottom
    PipeCompare      pipe_f_compare      data over the stack, the ratio panel, the systematic
                                         band breathes
    PipeFit          pipe_g_fit          mu slides and tightens, nuisance parameters pull, the
                                         ratio flattens; sigma = mu sigma_pred and
                                         sigma = (N - B)/(A eps L) lit term by term
    PipeThree        pipe_h_three        the spine shrinks and triplicates in the ee / mumu /
                                         tautau channel colours (hand-off to sections 3-5)

What each node is (drawn lepton-generic, grounded in this project's Z->mumu analysis):
  0 collisions / CMS detector            CMSSlice (logo geometry)
  1 recorded + simulated events, skims   z-mumu/docs/10-skims.md
  2 selection                            z-mumu/docs/01-selection.md, docs/14 (mumu_SR definition)
  3 corrections: tag and probe, SF       z-mumu/docs/12-tag-and-probe-fits.md
  4 backgrounds: simulation + Fakes      z-mumu/docs/13-fake-factor.md (same-sign control region)
  5 data vs prediction, ratio, syst.     z-mumu/docs/14-fit-and-systematics.md
  6 fit: mu_Z, nuisance parameters,      docs/CONVENTIONS.md section 2 (sigma_fid = mu x sigma_pred),
    sigma = (N - B)/(A eps L)            z-mumu/docs/06-cross-section.md
Symbols only (user decision, 15 Sep 2026), except the stacked histograms: on the user's
request (16 Sep 2026) the stacks of pipe_e / pipe_f are the real Z->mumu signal region
(data/zmumu_sr_stack.json, nominal prediction + Fakes, rebinned to 2 GeV) on a log axis
with decade and mass tick labels, and the fit panel of pipe_g shows the real 5 GeV
pre-fit -> post-fit data/pred (data/zmumu_fit.json). Everything else is schematic
(seeded RNG for texture); the legend stays lepton-generic.

Delivered clips (manim sections, tools/deliver_chain.py; a zoom-out and the next zoom-in
share one clip name, so they are joined):
    pipe_a_map | pipe_b1_collisions  pipe_b2_skim  pipe_b3_simulation |
    pipe_c1_trigger  pipe_c2_lepton_id  pipe_c3_two_leptons  pipe_c4_mass_window |
    pipe_d1_tag_probe  pipe_d2_efficiency  pipe_d3_scale_factor |
    pipe_e1_simulation_stack  pipe_e2_control_region  pipe_e3_fake_factor |
    pipe_f1_data  pipe_f2_ratio  pipe_f3_uncertainty |
    pipe_g1_fit_model  pipe_g2_fit  pipe_g3_cross_section | pipe_h_three

Zoom mechanics: a zoom is a camera map F(p) = to + z (p - about) applied to the
whole spine (``Transform`` to a mapped copy of the same builder output, strokes
scaled per family member) while the clip's detail group arrives through the
inverse map; the zoom-out maps the detail back onto the node and ``Transform``s
the spine to ``spine_state(done)``, so every clip ends exactly on the builder. The
camera moves run under ``keepout_mask()`` (the zoom sweeps the neighbours through the
title band and the top-left block); every held state keeps both clear by its anchors.
ManimCE 0.20.1 trap: ``VMobject.scale(..., scale_stroke=True)`` on a group sets
*every* member to the group's own width x factor, so ``scale_strokes`` below
rescales each member instead.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    DOWN, LEFT, ORIGIN, RIGHT, UP, AnimationGroup, ArcBetweenPoints, Create, CurvedArrow,
    DashedLine, DashedVMobject, Dot, FadeIn, FadeOut, Flash, GrowFromCenter, LaggedStart, Line,
    Polygon, Rectangle, ReplacementTransform, RoundedRectangle, Scene, Sector, Succession, Transform,
    VGroup, VMobject, ValueTracker, always_redraw, config, rate_functions,
)
from manim.utils.rate_functions import squish_rate_func  # noqa: E402

from style.bnd_style import *  # noqa: E402,F401,F403
from s1_drell_yan import _FEY  # noqa: E402  (the section-1 diagram vertices: the mini icon is that diagram)

EASE = rate_functions.ease_in_out_sine

# ---------------------------------------------------------------------------
# every anchor in one place (scene units; keep the title band y > 2.7 and the top-left
# block x < -5.85, y > 0.22 empty: TITLE_BAND_Y, CORNER_X_MAX, CORNER_Y_MIN)
# ---------------------------------------------------------------------------
NODE_X = (-6.2, -4.13, -2.07, 0.0, 2.07, 4.13, 6.2)
SPINE_Y = -0.55                                  # centres the map in the usable area (y < 2.7)
BOX = 1.24
ZOOM = 5.0
ALL = tuple(range(7))
TAIL = 0.1                     # end on the last change (tools/render.py grabs the true last frame)
DONE = {"a": (), "b": (0, 1), "c": (0, 1, 2), "d": (0, 1, 2, 3),
        "e": (0, 1, 2, 3, 4), "f": (0, 1, 2, 3, 4, 5), "g": ALL}

A = {
    # map ----------------------------------------------------------------
    "river": (-4.13, SPINE_Y + 1.78),            # mini Drell-Yan icon (simulation) above the files node
    "river_arrow": ((-4.13, SPINE_Y + 1.40), (-4.13, SPINE_Y + 0.72)),
    "detail": (0.0, -0.35),                      # where a zoomed node's detail is centred
    # pipe_b: nodes 0 and 1 zoomed 3.2x about their midpoint ------------
    "b_zoom": 3.2, "b_about": (-5.165, SPINE_Y), "b_to": (0.0, -0.62),
    "b_dy": (-3.31, 2.03),                       # the Drell-Yan icon above the zoomed slice
    "b_kappa": 0.85, "b_pile_step": (0.06, 0.12), "b_dy_scale": 0.34,
    "b_events": ((0.55, 3.55, 1), (2.35, 5.45, -1), (4.2, 1.05, 1)),   # (phi_1, phi_2, charge of 1)
    "b_sim_event": (1.35, 4.35, -1),
    # pipe_c: the funnel -------------------------------------------------
    "c_x": 1.0, "c_ys": (1.72, 0.84, -0.04, -0.92, -1.80), "c_bar_h": 0.50, "c_w": 5.4,
    "c_frac": (1.0, 0.55, 0.28, 0.275, 0.25),      # the big drop is "exactly two"; opposite sign removes ~nothing
    "c_rim": (3.55, 2.25), "c_neck": (0.95, -0.55), "c_spout_y": -2.75,   # (half-width, y) of the funnel
    "c_lab_x": -2.55, "c_lab_h": 0.30,
    # pipe_d: tag and probe ----------------------------------------------
    "d_slice": (-4.10, -1.20), "d_slice_s": 0.75, "d_tag_phi": 2.45, "d_probe_phi": -0.40,
    "d_kappa": 0.75, "d_lab_r": 0.62, "d_pass": (0.55, 1.05), "d_fail": (4.25, 1.05), "d_hist_wh": (3.0, 1.8),
    "d_eps": (2.4, -1.05), "d_slider": (3.05, -2.80), "d_slider_len": 4.0,
    # pipe_e: backgrounds ------------------------------------------------
    "e_box": (-4.25, -0.50), "e_box_wh": (2.95, 3.05), "e_slice_s": 0.43, "e_tag_phi": 2.55,
    "e_fake_phi": 0.25, "e_kappa": 0.8, "e_jet_half": 0.32, "e_plot": (2.35, -0.25), "e_plot_wh": (5.5, 3.3),
    "e_key": (4.20, 1.10), "e_arrow": ((-2.65, -1.62), 65.0),    # start (scene), end mass: into the Fakes sliver
    # pipe_f: comparison -------------------------------------------------
    "f_plot": (1.35, 0.95), "f_plot_wh": (7.0, 2.8), "f_ratio": (1.35, -1.40), "f_ratio_h": 1.15,
    "f_key": (-4.80, 1.30), "f_band": 0.02,
    # pipe_g: fit ---------------------------------------------------------
    "g_slider": (-2.60, 1.75), "g_slider_len": 4.4, "g_pulls": (-3.0, -0.75), "g_pull_len": 3.4,
    "g_ratio": (3.45, 1.50), "g_ratio_wh": (5.0, 1.5), "g_sig1": (3.45, -0.50), "g_sig2": (3.45, -2.20),
    # pipe_h: three strips ----------------------------------------------
    "h_scale": 0.55, "h_x": 0.75, "h_ys": {"ee": 1.25, "mumu": -0.25, "tautau": -1.75}, "h_lab_x": -3.2,
    "h_box_tint": {"ee": 0.80, "mumu": 0.55, "tautau": 0.80},
}

# the stacks: the real Z->mumu signal region (nominal prediction + Fakes), 1 GeV -> 2 GeV bins
SR = load_data("zmumu_sr_stack")
FIT = load_data("zmumu_fit")
assert SR["data"]["total"] == 10378567 and len(SR["edges"]) == 61
assert len(FIT["sr_12bin"]["ratio_prefit"]) == 12


def _r2(a) -> np.ndarray:
    return np.asarray(a, dtype=float).reshape(30, 2).sum(axis=1)


def _mc(*names) -> np.ndarray:
    return sum(_r2(SR["mc"][n]["nominal"]) for n in names)


EDGES = np.asarray(SR["edges"], dtype=float)[::2]    # 30 bins of 2 GeV, 60-120
assert len(EDGES) == 31 and EDGES[0] == 60.0 and EDGES[-1] == 120.0
XC = 0.5 * (EDGES[:-1] + EDGES[1:])
SIG_FILL = lighten(DETECTOR_ACCENT, 0.45)             # generic Z -> ll (method chapter tint)
# bottom-up; the legend is lepton-generic (the method), the numbers are the mumu channel's
LAYERS = (("Fakes", _r2(SR["fakes"]["counts"]), SAMPLE["Fakes"]),
          ("VV", _mc("WW", "WZ", "ZZ"), SAMPLE["WZ"]),
          ("TauTau", _mc("DYtautau"), SAMPLE["DYtautau"]),
          ("Top", _mc("TTbar", "SingleTop"), SAMPLE["TTbar"]),
          ("Z", _mc("DYmumu"), SIG_FILL))
TOTAL = sum(c for _, c, _ in LAYERS)
DATA = _r2(SR["data"]["counts"])
assert abs(TOTAL.sum() - SR["totals"]["nominal"]["mc_plus_fakes"]) < 1e-3
Y_EXP = (1.8, 6.8)                                    # log axis 10^1.8 .. 10^6.8 events / 2 GeV
EDGES12 = np.asarray(FIT["sr_12bin"]["edges"], dtype=float)
RATIO12_PRE = np.asarray(FIT["sr_12bin"]["ratio_prefit"], dtype=float)
RATIO12_POST = np.asarray(FIT["sr_12bin"]["ratio_postfit"], dtype=float)

EDGES_T = np.arange(60.0, 121.0, 3.0)                # tag-and-probe spectra (20 bins)
_PK = dict(sigma=4.2, tail=(1.0, 3.0), floor=0.0)
PASS = schematic_zpeak(EDGES_T, height=0.78, seed=5, jitter=0.03, **_PK) + 0.04
FAIL = 0.24 * schematic_zpeak(EDGES_T, height=0.78, seed=6, jitter=0.06, **_PK) + 0.09
PASS_MC = 1.12 * schematic_zpeak(EDGES_T, height=0.78, seed=5, jitter=0.0, **_PK) + 0.035
FAIL_MC = 0.12 * schematic_zpeak(EDGES_T, height=0.78, seed=6, jitter=0.0, **_PK) + 0.06

PULLS = (0.45, -0.35, 0.15, 0.70, -0.20)
CONSTR = (0.55, 0.80, 0.95, 0.45, 0.85)


# ---------------------------------------------------------------------------
# small scene-local helpers
# ---------------------------------------------------------------------------

def P(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    return p if p.shape == (3,) else np.array([p[0], p[1], 0.0])


_TEX_N: list = []


def tex(expr: str, h: float = 0.30, color=INK, **kw):
    """MathTex scaled so that a capital N is ``h`` high (uniform symbol size)."""
    if not _TEX_N:
        _TEX_N.append(mathtex("N").height)
    return mathtex(expr, color=color, **kw).scale(h / _TEX_N[0])


def scale_strokes(mob, f: float):
    """Multiply every family member's own stroke width by ``f`` (what
    ``scale_stroke=True`` should do; ManimCE 0.20.1 overwrites members with the
    group's width instead)."""
    for m in mob.get_family():
        if isinstance(m, VMobject):
            m.set_stroke(width=m.get_stroke_width() * f, family=False)
            m.set_stroke(width=m.get_stroke_width(background=True) * f, background=True, family=False)
    return mob


def camera(mob, z: float, about, to):
    """The zoom map F(p) = to + z (p - about), strokes scaled by z. In place."""
    mob.scale(z, about_point=P(about))
    mob.shift(P(to) - P(about))
    return scale_strokes(mob, z)


def keepout_mask(pad: float = 0.02) -> VGroup:
    """Background-coloured cover of the user's PowerPoint overlays (title band y > TITLE_BAND_Y,
    top-left block x < CORNER_X_MAX, y > CORNER_Y_MIN), held in the foreground during the camera
    moves only (``masked_play``). A zoom about a node sweeps its neighbours and the river
    up and left through the overlays; they now pass under the overlay edge instead of being drawn
    there. Every held frame is clear on its own (the anchors in ``A``); the cover is white on white."""
    out = 0.3
    x0, x1 = -config.frame_width / 2 - out, config.frame_width / 2 + out
    y1 = config.frame_height / 2 + out
    band = Polygon(P((x0, TITLE_BAND_Y - pad)), P((x1, TITLE_BAND_Y - pad)), P((x1, y1)), P((x0, y1)))
    corner = Polygon(P((x0, CORNER_Y_MIN - pad)), P((CORNER_X_MAX + pad, CORNER_Y_MIN - pad)),
                     P((CORNER_X_MAX + pad, y1)), P((x0, y1)))
    return VGroup(band, corner).set_fill(col(BG), opacity=1.0).set_stroke(width=0)


def masked_play(scene, *anims, **kw):
    """``scene.play`` with ``keepout_mask()`` on top of everything for that play only."""
    mask = keepout_mask()
    scene.add_foreground_mobject(mask)
    scene.play(*anims, **kw)
    scene.remove(mask)


def win(t0: float, dur: float, total: float, func=EASE):
    """Rate function that runs ``func`` in the window [t0, t0+dur] of a play of ``total`` s."""
    return squish_rate_func(func, max(0.0, t0 / total), min(1.0, (t0 + dur) / total))


def axis_titles(dax, x_expr=None, y_expr=None, h=0.26, buff=0.28):
    """Axis titles placed from the plot frame (DataAxes anchors its titles on the
    tick groups, which are empty here: no ticks on a schematic axis)."""
    g = VGroup()
    if x_expr is not None:
        g.add(tex(x_expr, h=h).next_to(dax.frame, DOWN, buff=buff))
    if y_expr is not None:
        g.add(tex(y_expr, h=h).next_to(dax.frame, LEFT, buff=buff))
    return g


def schematic_axes(x_length, y_length, y_top=1.15, center=ORIGIN):
    dax = DataAxes([60, 120, 10], [0, y_top, 0.5], x_length, y_length, x_ticks=[], y_ticks=[],
                   show_x_labels=False, show_y_labels=False)
    return dax.move_frame_to(P(center))


def log_axes(x_length, y_length, center, x_labels=True) -> DataAxes:
    """The real-stack axes: log y (decade labels 10^3..10^6; none at the floor, where the
    Fakes sliver and the control-region arrow are), mass ticks without the corner labels."""
    dax = DataAxes([60, 120, 10], [*Y_EXP, 1], x_length, y_length, x_ticks=[70, 80, 90, 100, 110],
                   y_ticks=[3, 4, 5, 6], y_log=True, show_x_labels=x_labels, tick_label_h=0.18,
                   x_title=r"m_{\ell\ell}\ [\mathrm{GeV}]" if x_labels else None,
                   y_title=r"\mathrm{events}\,/\,2\,\mathrm{GeV}", title_h=0.22, title_buff=0.16)
    return dax.move_frame_to(P(center))


def stack_stage(dax, present):
    return stack_hist(dax, EDGES, [(n, c if n in present else np.zeros_like(c), colour)
                                   for n, c, colour in LAYERS])


def step_xy(edges, vals):
    xs, ys = [], []
    for i, v in enumerate(vals):
        xs += [edges[i], edges[i + 1]]
        ys += [v, v]
    return xs, ys


def dy_icon(center, s: float = 0.25, color=THEORY, sw: float = 2.4) -> VGroup:
    """The section-1 Drell-Yan diagram (s1_drell_yan._FEY) shrunk by ``s``, no labels."""
    c = P(center)
    F = {k: c + s * np.array([x, y, 0.0]) for k, (x, y) in _FEY.items()}
    lines = VGroup(fline(F["q_in"], F["v1"], color=color, sw=sw),
                   fline(F["qb_in"], F["v1"], color=color, sw=sw),
                   wavy(F["v1"], F["v2"], color=color, amplitude=0.18 * s, wavelength=0.78 * s, sw=0.9 * sw),
                   fline(F["v2"], F["l_out"], color=color, sw=sw),
                   fline(F["v2"], F["lb_out"], color=color, sw=sw))
    dots = VGroup(vertex_dot(F["v1"], color=color, radius=0.035), vertex_dot(F["v2"], color=color, radius=0.035))
    return VGroup(lines, dots)


# ---------------------------------------------------------------------------
# node icons (built at spine size around ORIGIN; spine() centres them)
# ---------------------------------------------------------------------------
FILE_W, FILE_H, PILE_STEP = 0.30, 0.40, (0.05, 0.06)


def file_glyph(fill, stroke, w=FILE_W, h=FILE_H, sw=1.6) -> VGroup:
    e = 0.34 * w
    body = Polygon(P((-w / 2, -h / 2)), P((w / 2, -h / 2)), P((w / 2, h / 2 - e)), P((w / 2 - e, h / 2)),
                   P((-w / 2, h / 2)), fill_color=col(fill), fill_opacity=1.0,
                   stroke_color=col(stroke), stroke_width=sw)
    fold = Polygon(P((w / 2 - e, h / 2)), P((w / 2 - e, h / 2 - e)), P((w / 2, h / 2 - e)),
                   fill_color=lighten(stroke, 0.45), fill_opacity=1.0, stroke_color=col(stroke), stroke_width=sw)
    lines = VGroup(*[Line(P((-w / 2 + 0.18 * w, y)), P((xr, y)), stroke_color=col(stroke), stroke_width=0.8 * sw)
                     for y, xr in ((0.14 * h, w / 2 - e - 0.10 * w), (-0.06 * h, w / 2 - 0.18 * w),
                                   (-0.26 * h, w / 2 - 0.18 * w))])
    return VGroup(body, fold, lines)


def pile(n, fill, stroke) -> VGroup:
    return VGroup(*[file_glyph(fill, stroke).shift(P((i * PILE_STEP[0], i * PILE_STEP[1])))
                    for i in range(n)]).move_to(ORIGIN)


DATA_FILE = (LIGHT_GREY, INK)
SIM_FILE = (lighten(THEORY, 0.78), THEORY)


def icon_detector() -> CMSSlice:
    return mini_slice(0.16, ORIGIN)


def icon_files() -> VGroup:
    data = pile(2, *DATA_FILE).move_to(P((-0.21, 0.0)))
    sim = pile(2, *SIM_FILE).move_to(P((0.21, 0.0)))
    g = VGroup(data, sim)
    g.data, g.sim = data, sim
    return g


def icon_funnel(line, fill) -> VGroup:
    pts = ((-0.46, 0.40), (0.46, 0.40), (0.10, -0.08), (0.10, -0.40), (-0.10, -0.40), (-0.10, -0.08))
    body = Polygon(*[P(p) for p in pts], stroke_color=col(line), stroke_width=2.5,
                   fill_color=lighten(fill, 0.85), fill_opacity=1.0)
    bars = VGroup(*[Rectangle(width=w, height=0.06, stroke_width=0, fill_color=col(c), fill_opacity=1.0)
                    .move_to(P((0.0, y))) for w, y, c in ((0.56, 0.29, GREY), (0.32, 0.15, GREY), (0.13, 0.01, line))])
    return VGroup(body, bars)


def icon_tnp(line) -> VGroup:
    v = P((-0.06, -0.10))
    tag = ArcBetweenPoints(v, P((-0.42, 0.36)), angle=-0.7, stroke_color=col(INK), stroke_width=4.5)
    probe = DashedVMobject(ArcBetweenPoints(v, P((0.44, -0.34)), angle=-0.7, stroke_color=col(INK),
                                            stroke_width=3.0), num_dashes=5)
    vtx = Dot(v, radius=0.045, color=col(INK))
    eps = tex(r"\varepsilon", h=0.34, color=line).move_to(P((0.24, 0.22)))
    return VGroup(tag, probe, vtx, eps)


def icon_control(line) -> VGroup:
    box = DashedVMobject(RoundedRectangle(width=0.94, height=0.70, corner_radius=0.12, stroke_color=col(line),
                                          stroke_width=2.5), num_dashes=20)
    lab = tex(r"\ell^{\pm}\ell^{\pm}", h=0.22).move_to(ORIGIN)
    return VGroup(box, lab)


def icon_compare(line, fill, sig_fill=None) -> VGroup:
    x0, x1, yb, yt = -0.42, 0.42, -0.06, 0.40
    ink = col(INK)
    frame = VGroup(Line(P((x0, yt)), P((x0, yb)), stroke_color=ink, stroke_width=1.8),
                   Line(P((x0, yb)), P((x1, yb)), stroke_color=ink, stroke_width=1.8))
    bh = 0.07
    peak = lambda x: 0.26 * np.exp(-0.5 * (x / 0.11) ** 2)   # noqa: E731
    xs = np.linspace(x0, x1, 41)
    bkg = Polygon(P((x0, yb)), P((x0, yb + bh)), P((x1, yb + bh)), P((x1, yb)), stroke_width=0,
                  fill_color=col(SAMPLE["TTbar"]), fill_opacity=1.0)
    sig = Polygon(*[P((x, yb + bh + peak(x))) for x in xs], *[P((x, yb + bh)) for x in xs[::-1]],
                  stroke_width=0, fill_color=lighten(fill, 0.45) if sig_fill is None else col(sig_fill), fill_opacity=1.0)
    jit = (0.01, -0.012, 0.006, -0.004, 0.012, -0.008, 0.004)
    dots = VGroup(*[Dot(P((x, yb + bh + peak(x) + j)), radius=0.028, color=ink)
                    for x, j in zip(np.linspace(-0.33, 0.33, 7), jit)])
    ref = DashedLine(P((x0, -0.30)), P((x1, -0.30)), dash_length=0.05, stroke_color=col(GREY), stroke_width=1.5)
    rdots = VGroup(*[Dot(P((x, -0.30 + j)), radius=0.024, color=ink)
                     for x, j in zip(np.linspace(-0.32, 0.32, 5), (-0.03, 0.02, -0.02, 0.03, -0.01))])
    return VGroup(frame, bkg, sig, dots, ref, rdots)


def icon_fit(line) -> VGroup:
    y, ink = 0.20, col(INK)
    axis = Line(P((-0.40, y)), P((0.40, y)), stroke_color=ink, stroke_width=2.0)
    ref = DashedLine(P((0.08, y - 0.05)), P((0.08, y + 0.24)), dash_length=0.05, stroke_color=col(THEORY),
                     stroke_width=2.0)
    bar = Line(P((-0.16, y)), P((0.06, y)), stroke_color=ink, stroke_width=4.0)
    dot = Dot(P((-0.05, y)), radius=0.05, color=ink)
    mu = tex(r"\mu", h=0.24).move_to(P((-0.27, -0.20)))
    ln = Line(P((-0.13, -0.19)), P((0.10, -0.19)), stroke_color=ink, stroke_width=2.0)
    arr = VGroup(ln, arrow_tip_on(ln, color=INK, at=1.0, tip_length=0.08))
    sigma = tex(r"\sigma", h=0.30, color=line).move_to(P((0.27, -0.19)))
    main, sig = VGroup(axis, ref, bar, dot, mu), VGroup(arr, sigma)
    g = VGroup(main, sig)
    g.main, g.sig = main, sig
    return g


def node_icons(line=DETECTOR_ACCENT, fill=DETECTOR_ACCENT, sig_fill=None) -> list:
    return [icon_detector(), icon_files(), icon_funnel(line, fill), icon_tnp(line), icon_control(line),
            icon_compare(line, fill, sig_fill), icon_fit(line)]


def river_group() -> VGroup:
    icon = dy_icon(A["river"])
    a0, a1 = A["river_arrow"]
    ln = Line(P(a0), P(a1), stroke_color=col(THEORY), stroke_width=3.0)
    arrow = VGroup(ln, arrow_tip_on(ln, color=THEORY, at=1.0, tip_length=0.2))
    g = VGroup(icon, arrow)
    g.icon, g.arrow = icon, arrow
    return g


# ---------------------------------------------------------------------------
# builders: the chain states
# ---------------------------------------------------------------------------

def spine_state(done=()) -> dict:
    """End state of every clip a..g: the spine with nodes ``done`` tinted, and the river."""
    return {"spine": spine(node_icons(), NODE_X, y=SPINE_Y, box=(BOX, BOX), done=done),
            "river": river_group()}


ORDER_SPINE = ("spine", "river")


def strip(ch=None, at="ee") -> VGroup:
    """The spine shrunk to a strip at the slot of channel ``at``; coloured for
    channel ``ch`` (None = the method cyan of section 2)."""
    line = DETECTOR_ACCENT if ch is None else CHANNEL_LINE[ch]
    fill = DETECTOR_ACCENT if ch is None else CHANNEL[ch]
    sig = darken(fill, 0.18) if ch == "mumu" else None      # the pale gold peak vanishes on the gold box tint
    sp = spine(node_icons(line, fill, sig), NODE_X, y=SPINE_Y, box=(BOX, BOX), done=ALL, accent=line)
    if ch is not None:        # box tint from the channel fill colour (the pale gold line tint reads grey)
        for b in sp.boxes:
            b.set_fill(lighten(fill, A["h_box_tint"][ch]), opacity=1.0)
    return camera(sp, A["h_scale"], (0.0, SPINE_Y), (A["h_x"], A["h_ys"][at]))


CHANNEL_TEX = {"ee": r"e^{+}e^{-}", "mumu": r"\mu^{+}\mu^{-}", "tautau": r"\tau^{+}\tau^{-}"}


def strip_label(ch) -> VGroup:
    lab = tex(CHANNEL_TEX[ch], h=0.34, color=CHANNEL_LINE[ch])
    return lab.next_to(P((A["h_lab_x"], A["h_ys"][ch])), LEFT, buff=0.0)


def state_h() -> dict:
    """End state of PipeThree: three channel strips and their symbols."""
    st = {ch: strip(ch, ch) for ch in ("ee", "mumu", "tautau")}
    st.update({f"lab_{ch}": strip_label(ch) for ch in ("ee", "mumu", "tautau")})
    return st


ORDER_H = ("ee", "mumu", "tautau", "lab_ee", "lab_mumu", "lab_tautau")


def zoomed_spine(sp, z, about, to, boxes=(), arrows=(), icons=()) -> VGroup:
    """Camera-mapped copy of the spine; parts not listed go to opacity 0."""
    g = camera(sp.copy(), z, about, to)
    for part, keep in ((g[0], boxes), (g[1], arrows), (g[2], icons)):
        for i, m in enumerate(part):
            if i not in keep:
                m.set_opacity(0.0)
    return g


def node_point(k):
    return (NODE_X[k], SPINE_Y)


# ---------------------------------------------------------------------------
# the chained scene base
# ---------------------------------------------------------------------------

class PipeScene(Scene):
    DONE_BEFORE: tuple = ()

    def open_chain(self, clip: str):
        white_background(self)
        self.st = spine_state(self.DONE_BEFORE)
        add_state(self, self.st, ORDER_SPINE)
        clip_open(self, clip)

    def zoom_in(self, detail, z, about, to, spine_target, run_time=1.3, extra=()):
        anims = [Transform(self.st["spine"], spine_target),
                 Transform(self.st["river"], camera(self.st["river"].copy(), z, about, to).set_opacity(0.0))]
        if detail is not None:
            start = camera(detail.copy(), 1.0 / z, to, about).set_opacity(0.0)
            anims.append(ReplacementTransform(start, detail))
        masked_play(self, *anims, *extra, run_time=run_time, rate_func=EASE)

    def detail_mobjects(self) -> list:
        keep = {id(self.st["spine"]), id(self.st["river"])}
        out = []

        def visit(m):
            if id(m) in keep or isinstance(m, ValueTracker):
                return
            if isinstance(m, VMobject):
                if m.has_points() or len(m.get_family()) > 1:
                    out.append(m)
                return
            for s in m.submobjects:
                visit(s)

        for m in list(self.mobjects):
            visit(m)
        return out

    def zoom_out(self, done, z, about, to, run_time=1.3):
        end = spine_state(done)
        det = VGroup(*self.detail_mobjects())
        anims = [Transform(self.st["spine"], end["spine"]), Transform(self.st["river"], end["river"])]
        if len(det):
            anims.append(Transform(det, camera(det.copy(), 1.0 / z, to, about).set_opacity(0.0), remover=True))
        masked_play(self, *anims, run_time=run_time, rate_func=EASE)
        # empty Group wrappers left behind by LaggedStart / Succession (no points in their family)
        keep = {id(self.st["spine"]), id(self.st["river"])}
        stray = [m for m in self.mobjects if id(m) not in keep and not m.family_members_with_points()]
        if stray:
            self.remove(*stray)

    def close_chain(self):
        check_order(self, self.st, ORDER_SPINE)
        self.wait(TAIL)


# ---------------------------------------------------------------------------
# a: the map
# ---------------------------------------------------------------------------

class PipeMap(Scene):
    """pipe_a_map (~8 s). Empty frame -> the spine: each node box grows with its
    icon, the arrow to the next node draws; the simulated river (mini
    Drell-Yan) drops into the files node; the sigma of the fit node comes last.
    Ends on ``spine_state(())``."""

    def construct(self):
        white_background(self)
        st = spine_state(DONE["a"])
        sp, rv = st["spine"], st["river"]
        clip_open(self, "pipe_a_map")
        T = 7.4
        anims = []
        for i in range(7):
            t0 = 0.05 + 0.92 * i
            icon = sp.icons[i]
            anims += [GrowFromCenter(sp.boxes[i], rate_func=win(t0, 0.65, T)),
                      FadeIn(icon.main if i == 6 else icon, scale=0.6, rate_func=win(t0 + 0.15, 0.65, T))]
            if i < 6:
                # squished rate functions need lag_ratio=0 (the lag re-maps alpha per submobject)
                anims += [Create(sp.arrows[i][0], lag_ratio=0.0, rate_func=win(t0 + 0.55, 0.4, T, rate_functions.linear)),
                          FadeIn(sp.arrows[i][1], rate_func=win(t0 + 0.85, 0.2, T))]
        anims.append(FadeIn(sp.icons[6].sig, shift=LEFT * 0.12, rate_func=win(6.75, 0.6, T)))
        river = AnimationGroup(FadeIn(rv.icon, shift=DOWN * 0.15, rate_func=win(1.75, 0.7, T)),
                               Create(rv.arrow[0], lag_ratio=0.0, rate_func=win(2.25, 0.45, T, rate_functions.linear)),
                               FadeIn(rv.arrow[1], rate_func=win(2.6, 0.2, T)), group=rv)
        self.play(AnimationGroup(*anims, group=sp), river, run_time=T)
        check_order(self, st, ORDER_SPINE)
        self.wait(TAIL)


# ---------------------------------------------------------------------------
# b: recorded and simulated events
# ---------------------------------------------------------------------------

class PipeData(PipeScene):
    """pipe_b_data (~10 s). Zoom 3.2x on nodes 1-2 (the slice and the files
    stay as frames). Three collisions flash in the slice, two tracks each, and
    a file flies to the pile; five more stream in; the pile of 8 is squeezed
    to 2 (the skim, docs/10). The simulated river: the Drell-Yan icon drops
    into the same slice, the simulated event gives purple files that land next
    to the data pile. Zoom out to ``spine_state((0, 1))``."""
    DONE_BEFORE = DONE["a"]

    def construct(self):
        self.open_chain("pipe_b1_collisions")
        z, about, to = A["b_zoom"], A["b_about"], A["b_to"]
        fresh = spine_state(())["spine"]
        c0 = P(to) + z * (P(node_point(0)) - P(about))               # zoomed beam spot
        geo = mini_slice(0.16 * z, c0)                                 # geometry only (never added)
        zicon = camera(fresh.icons[1].copy(), z, about, to)            # the files icon as zoomed
        d0, d1 = zicon.data[0], zicon.data[1]
        pile_c = 0.5 * (d0.get_center() + d1.get_center())
        slots = [d0.copy().move_to(pile_c + P(A["b_pile_step"]) * (i - 3.5)) for i in range(8)]

        self.zoom_in(None, z, about, to, zoomed_spine(self.st["spine"], z, about, to,
                                                       boxes=(0, 1), arrows=(0,), icons=(0,)))

        def event_tracks(ev, color):
            p1, p2, q = ev
            k = A["b_kappa"]
            r_trk = geo.radii["tob"][1]                  # end at the tracker: says nothing about the flavour
            return VGroup(*[track(geo, phi, sgn * k, r_end=r_trk, color=color, sw=4.4, outline=WHITE, outline_sw=3.2)
                            for phi, sgn in ((p1, q), (p2, -q))])

        def flight(target):
            return ReplacementTransform(target.copy().scale(0.3).move_to(c0), target, path_arc=-0.7)

        for j, ev in enumerate(A["b_events"]):
            f = 1.0 if j == 0 else 0.65
            self.play(Flash(c0, color=col(INK), line_length=0.14, flash_radius=0.26, num_lines=10,
                            run_time=0.32 * f))
            trks = event_tracks(ev, INK)
            self.play(Create(trks, lag_ratio=0.0), run_time=0.45 * f, rate_func=rate_functions.linear)
            self.play(flight(slots[j]), FadeOut(trks), run_time=0.55 * f, rate_func=EASE)
        self.play(LaggedStart(*[flight(s) for s in slots[3:]], lag_ratio=0.35), run_time=0.9)
        clip_cut(self, "pipe_b2_skim")
        # the skim: 8 -> 2
        self.play(ReplacementTransform(VGroup(*slots[:4]), VGroup(d0.copy())),
                  ReplacementTransform(VGroup(*slots[4:]), VGroup(d1.copy())), run_time=0.9, rate_func=EASE)
        clip_cut(self, "pipe_b3_simulation")
        # the simulated river
        dy = dy_icon(A["b_dy"], s=A["b_dy_scale"], sw=3.0)
        self.play(FadeIn(dy, shift=DOWN * 0.15), run_time=0.4)
        self.play(FadeOut(dy.copy(), scale=0.1, shift=c0 - dy.get_center()), run_time=0.5, rate_func=EASE)
        self.play(Flash(c0, color=col(THEORY), line_length=0.14, flash_radius=0.26, num_lines=10, run_time=0.25))
        strk = event_tracks(A["b_sim_event"], THEORY)
        self.play(Create(strk, lag_ratio=0.0), run_time=0.35, rate_func=rate_functions.linear)
        s0, s1 = zicon.sim[0].copy(), zicon.sim[1].copy()
        self.play(LaggedStart(flight(s0), flight(s1), lag_ratio=0.4), FadeOut(strk), run_time=0.75)
        clip_cut(self, "pipe_c1_trigger")
        self.zoom_out(DONE["b"], z, about, to)
        self.close_chain()


# ---------------------------------------------------------------------------
# c: selection funnel
# ---------------------------------------------------------------------------

def funnel_outline() -> VGroup:
    """A straight cone from the rim to the neck with a soft corner, then a straight neck."""
    x = A["c_x"]
    (r, y_rim), (n, y_neck), y_sp = A["c_rim"], A["c_neck"], A["c_spout_y"]
    k = 14.0
    soft = lambda t: np.log1p(np.exp(k * t)) / k                      # noqa: E731  smooth max(t, 0)
    ys = np.linspace(y_rim, y_sp, 120)
    hw = n + (r - n) * soft((ys - y_neck) / (y_rim - y_neck)) / soft(1.0)
    sides = VGroup()
    for sgn in (-1, 1):
        vm = VMobject(stroke_color=col(DETECTOR_ACCENT), stroke_width=3.5)
        vm.set_points_as_corners([P((x + sgn * h, y)) for h, y in zip(hw, ys)])
        sides.add(vm)
    return sides


def funnel_bar(i, color) -> Rectangle:
    return Rectangle(width=A["c_w"] * A["c_frac"][i], height=A["c_bar_h"], fill_color=col(color),
                     fill_opacity=1.0, stroke_color=darken(color, 0.25), stroke_width=1.5).move_to(
        P((A["c_x"], A["c_ys"][i])))


CUT_TEX = (r"\mathrm{trigger}", r"\ell\ p_{T},\,\mathrm{ID/iso}", r"N_{\ell}=2", r"\ell^{+}\ell^{-}",
           r"m_{\ell\ell}\in[\,m_{\mathrm{min}},\,m_{\mathrm{max}}\,]")


class PipeSelection(PipeScene):
    """pipe_c_selection (~8 s). Zoom on the funnel: one bar enters at the rim and
    shrinks through the cuts of docs/01 (trigger, lepton ID and isolation,
    exactly two leptons, opposite sign, the mass window), each stage leaving a
    pale ghost; the surviving bar is cyan. Proportions only."""
    DONE_BEFORE = DONE["b"]
    CUTS = {1: "pipe_c2_lepton_id", 2: "pipe_c3_two_leptons", 4: "pipe_c4_mass_window"}

    def construct(self):
        self.open_chain("pipe_c1_trigger")
        k, z, to = 2, ZOOM, A["detail"]
        about = node_point(k)
        outline = funnel_outline()
        self.zoom_in(VGroup(outline), z, about, to, zoomed_spine(self.st["spine"], z, about, to))
        labs = [tex(t, h=A["c_lab_h"]).next_to(P((A["c_lab_x"], y)), LEFT, buff=0.0)
                for t, y in zip(CUT_TEX, A["c_ys"])]
        bar = funnel_bar(0, GREY)
        self.play(GrowFromCenter(bar), FadeIn(labs[0], shift=RIGHT * 0.2), run_time=0.65, rate_func=EASE)
        for i in range(1, 5):
            if i in self.CUTS:
                clip_cut(self, self.CUTS[i])
            else:
                self.wait(0.15)
            mover = bar.copy()
            self.add(mover)
            self.play(Transform(mover, funnel_bar(i, DETECTOR_ACCENT if i == 4 else GREY)),
                      bar.animate.set_fill(col(LIGHT_GREY)).set_stroke(darken(LIGHT_GREY, 0.15)),
                      FadeIn(labs[i], shift=RIGHT * 0.2), run_time=1.0, rate_func=EASE)
            bar = mover
        clip_cut(self, "pipe_d1_tag_probe")
        self.zoom_out(DONE["c"], z, about, to)
        self.close_chain()


# ---------------------------------------------------------------------------
# d: tag and probe
# ---------------------------------------------------------------------------

def probe_rain(src, dax, counts, n, seed) -> list:
    rng = np.random.default_rng(seed)
    p = counts / counts.sum()
    out = []
    for _ in range(n):
        j = int(rng.choice(len(counts), p=p))
        x = rng.uniform(EDGES_T[j], EDGES_T[j + 1])
        dot = Dot(P(src), radius=0.05, color=col(SAMPLE["Data"]))
        out.append(Succession(dot.animate(run_time=0.55, rate_func=rate_functions.ease_in_quad)
                              .move_to(dax.c2p(x, counts[j])), FadeOut(dot, run_time=0.15)))
    return out


class PipeTnP(PipeScene):
    """pipe_d_tnp (~9 s). Zoom on the corrections node. A Z -> ll event: the
    tag drawn bold, the probe dashed (docs/12: tag tight + trigger-matched,
    probe tested). Probe copies fall into the pass and fail spectra, which
    grow; eps = N_pass / (N_pass + N_fail); the same measurement in simulation
    (purple dashed) gives eps_MC; the scale-factor slider eps_data / eps_MC
    leaves the reference 1 (symbols only)."""
    DONE_BEFORE = DONE["c"]

    def construct(self):
        self.open_chain("pipe_d1_tag_probe")
        k, z, to = 3, ZOOM, A["detail"]
        about = node_point(k)
        det = mini_slice(A["d_slice_s"], A["d_slice"])
        w, h = A["d_hist_wh"]
        ax_p = DataAxes([60, 120, 10], [0, 1.0, 0.5], w, h, x_ticks=[], y_ticks=[], show_x_labels=False,
                        show_y_labels=False).move_frame_to(P(A["d_pass"]))
        ax_f = DataAxes([60, 120, 10], [0, 1.0, 0.5], w, h, x_ticks=[], y_ticks=[], show_x_labels=False,
                        show_y_labels=False).move_frame_to(P(A["d_fail"]))
        t_p = tex(r"N_{\mathrm{pass}}", h=0.30).next_to(ax_p.frame, UP, buff=0.16)
        t_f = tex(r"N_{\mathrm{fail}}", h=0.30).next_to(ax_f.frame, UP, buff=0.16)
        xt = VGroup(axis_titles(ax_p, r"m_{\ell\ell}", h=0.22, buff=0.18),
                    axis_titles(ax_f, r"m_{\ell\ell}", h=0.22, buff=0.18))
        zero = np.zeros(len(PASS))
        h_p = step_hist(ax_p, EDGES_T, zero, color=SAMPLE["Data"], stroke_width=2.5, fill_opacity=0.35)
        h_f = step_hist(ax_f, EDGES_T, zero, color=SAMPLE["Data"], stroke_width=2.5, fill_opacity=0.35)
        detail = VGroup(det, ax_p, ax_f, t_p, t_f, xt, h_p, h_f)
        self.zoom_in(detail, z, about, to, zoomed_spine(self.st["spine"], z, about, to))

        kap, r_trk = A["d_kappa"], det.radii["tob"][1]       # tracks end at the tracker (flavour-neutral)
        tag = track(det, A["d_tag_phi"], kap, r_end=r_trk, color=INK, sw=6.5, outline=WHITE, outline_sw=3.0)
        probe_line = track(det, A["d_probe_phi"], -kap, r_end=r_trk, color=INK, sw=4.5)
        probe = VGroup(track(det, A["d_probe_phi"], -kap, r_end=r_trk, color=WHITE, sw=7.5),
                       DashedVMobject(probe_line, num_dashes=4))
        end = probe_line.pts[-1]

        def near_end(expr, p):           # symbol just beyond the track end, white halo over the detector
            phi = float(np.arctan2(p[1] - det.c[1], p[0] - det.c[0]))
            lab = tex(expr, h=0.30).move_to(det.point_at(A["d_lab_r"], phi))
            return lab.set_stroke(col(WHITE), width=6, background=True)

        l_tag = near_end(r"\ell_{\mathrm{tag}}", tag.pts[-1])
        l_probe = near_end(r"\ell_{\mathrm{probe}}", end)
        self.play(LaggedStart(Create(tag), FadeIn(l_tag), Create(probe), FadeIn(l_probe), lag_ratio=0.45),
                  run_time=1.4)
        clip_cut(self, "pipe_d2_efficiency")
        rain_p = probe_rain(end, ax_p, PASS, 8, seed=21)
        rain_f = probe_rain(end, ax_f, FAIL, 3, seed=22)
        drops = [rain_p[0], rain_p[1], rain_f[0], rain_p[2], rain_p[3], rain_p[4], rain_f[1], rain_p[5],
                 rain_p[6], rain_f[2], rain_p[7]]
        self.play(LaggedStart(*drops, lag_ratio=0.18),
                  Transform(h_p, step_hist(ax_p, EDGES_T, PASS, color=SAMPLE["Data"], stroke_width=2.5,
                                           fill_opacity=0.35),
                            rate_func=squish_rate_func(EASE, 0.25, 1.0)),
                  Transform(h_f, step_hist(ax_f, EDGES_T, FAIL, color=SAMPLE["Data"], stroke_width=2.5,
                                           fill_opacity=0.35),
                            rate_func=squish_rate_func(EASE, 0.25, 1.0)),
                  run_time=1.7)
        eps = tex(r"\varepsilon = \frac{N_{\mathrm{pass}}}{N_{\mathrm{pass}} + N_{\mathrm{fail}}}", h=0.30)
        eps.move_to(P(A["d_eps"]))
        self.play(FadeIn(eps, shift=UP * 0.2), run_time=0.6, rate_func=EASE)
        clip_cut(self, "pipe_d3_scale_factor")
        mc_p = data_trace(ax_p, *step_xy(EDGES_T, PASS_MC), color=THEORY, stroke_width=3.5, dashed_=True)
        mc_f = data_trace(ax_f, *step_xy(EDGES_T, FAIL_MC), color=THEORY, stroke_width=3.5, dashed_=True)
        self.play(Create(mc_p), Create(mc_f), run_time=0.7, rate_func=rate_functions.linear)
        sl = slider(0.9, 1.1, 1.0, ref=1.0, length=A["d_slider_len"], ticks=[], marker_r=0.09)
        sl.move_to(P(A["d_slider"]) - (sl.axis.get_center() - sl.get_center()))
        sf = tex(r"\varepsilon_{\mathrm{data}} / \varepsilon_{\mathrm{MC}}", h=0.30,
                 substrings_to_isolate=[r"\varepsilon_{\mathrm{MC}}"])
        sf.get_part_by_tex(r"\varepsilon_{\mathrm{MC}}").set_color(col(THEORY))
        sf.next_to(sl.axis, LEFT, buff=0.45)
        self.play(FadeIn(sf), FadeIn(sl), run_time=0.5)
        self.play(Transform(sl.marker, sl.marker_at(0.965, 0.012)), run_time=0.9, rate_func=EASE)
        clip_cut(self, "pipe_e1_simulation_stack")
        self.zoom_out(DONE["d"], z, about, to)
        self.close_chain()


# ---------------------------------------------------------------------------
# e: backgrounds
# ---------------------------------------------------------------------------

def jet_cone(det, phi, half) -> VGroup:
    """A jet drawn the event-display way: a pale cone from the beam spot into the
    coil around ``phi`` (scene-local; not in the style module)."""
    r = 0.5 * sum(det.radii["solenoid"])            # drawn into the coil so it reads at icon size
    return Sector(radius=r, angle=2 * half, start_angle=phi - half, arc_center=det.c,
                  fill_color=col(PARTICLE["hadron"]), fill_opacity=0.55, stroke_color=darken(PARTICLE["hadron"], 0.3),
                  stroke_width=1.5)


KEY_ROWS = ((r"Z/\gamma^{*}\to\ell\ell", SIG_FILL), (r"t\bar{t},\ tW", SAMPLE["TTbar"]),
            (r"VV,\ \tau\tau", SAMPLE["WZ"]), (r"\mathrm{Fakes}", SAMPLE["Fakes"]))
KEY_ROW = {"Z": 0, "Top": 1, "VVTT": 2, "Fakes": 3}


def stack_key() -> VGroup:
    key = colour_key(list(KEY_ROWS), label_h=0.21, swatch=(0.32, 0.21))
    sw = key.rows[KEY_ROW["VVTT"]][1]                 # two layers, one row: a two-stripe swatch
    halves = VGroup(*[Rectangle(width=sw.width / 2, height=sw.height, stroke_width=0, fill_color=col(c),
                                fill_opacity=0.95) for c in (SAMPLE["WZ"], SAMPLE["DYtautau"])])
    halves.arrange(RIGHT, buff=0).move_to(sw)
    frame = Rectangle(width=sw.width, height=sw.height, fill_opacity=0.0, stroke_color=darken(SAMPLE["WZ"], 0.25),
                      stroke_width=1.0).move_to(sw)
    key.rows[KEY_ROW["VVTT"]].submobjects[1] = VGroup(halves, frame)
    return key


class PipeBackgrounds(PipeScene):
    """pipe_e_backgrounds (~9.5 s). Zoom on the backgrounds node. The simulated
    processes stack up layer by layer (VV, ttbar, Z -> ll; small, ordered
    ttbar > VV > Fakes as in the mumu SR); a dashed control-region box opens
    with a same-sign pair on a small slice, one lepton inside a jet. docs/13:
    the same-sign control region measures a fake *factor* f (tight over
    anti-isolated), which is then applied to opposite-sign events with an
    anti-isolated lepton; the arrow carries "x f" (not the control-region
    events) into the plot, where the Fakes sliver grows at the bottom of the
    stack and lifts everything above it."""
    DONE_BEFORE = DONE["d"]

    def construct(self):
        self.open_chain("pipe_e1_simulation_stack")
        k, z, to = 4, ZOOM, A["detail"]
        about = node_point(k)
        dax = log_axes(*A["e_plot_wh"], center=A["e_plot"])
        stack = stack_stage(dax, ())
        detail = VGroup(dax, stack)
        self.zoom_in(detail, z, about, to, zoomed_spine(self.st["spine"], z, about, to))
        key = stack_key().move_to(P(A["e_key"]))
        present = []
        for layers, row in ((("VV", "TauTau"), "VVTT"), (("Top",), "Top"), (("Z",), "Z")):
            present += layers
            self.play(Transform(stack, stack_stage(dax, present)), FadeIn(key.rows[KEY_ROW[row]], shift=LEFT * 0.15),
                      run_time=0.85 if row != "Z" else 1.05, rate_func=EASE)
        clip_cut(self, "pipe_e2_control_region")
        bw, bh = A["e_box_wh"]
        box = DashedVMobject(RoundedRectangle(width=bw, height=bh, corner_radius=0.2, stroke_color=col(DETECTOR_ACCENT),
                                              stroke_width=3.0).move_to(P(A["e_box"])), num_dashes=44)
        geo = mini_slice(A["e_slice_s"], A["e_box"])
        lab = tex(r"\ell^{\pm}\ell^{\pm}", h=0.34).next_to(box, UP, buff=0.18)
        self.play(Create(box), FadeIn(geo), FadeIn(lab), run_time=0.8, rate_func=EASE)
        kap = A["e_kappa"]
        r_trk = geo.radii["tob"][1]                                          # tracks end at the tracker
        tag = track(geo, A["e_tag_phi"], kap, r_end=r_trk, color=INK, sw=4.5, outline=WHITE, outline_sw=3.0)
        fake = track(geo, A["e_fake_phi"], kap, r_end=r_trk, color=INK, sw=4.5)  # same sign; inside the jet cone
        p_end = fake.pts[-1]
        phi_j = float(np.arctan2(p_end[1] - geo.c[1], p_end[0] - geo.c[0]))
        jet = VGroup(jet_cone(geo, phi_j, A["e_jet_half"]), signature(geo, "jet", phi_j, seed=4))
        self.play(Create(tag), run_time=0.5, rate_func=rate_functions.linear)
        self.play(FadeIn(jet), Create(fake), run_time=0.7, rate_func=EASE)
        clip_cut(self, "pipe_e3_fake_factor")
        a0, m_end = A["e_arrow"]
        j = int(np.searchsorted(EDGES, m_end) - 1)
        y_end = np.sqrt(10 ** Y_EXP[0] * LAYERS[0][1][j])            # log-middle of the Fakes sliver
        arrow = CurvedArrow(P(a0), dax.c2p(m_end, y_end), angle=0.45, color=col(DETECTOR_ACCENT), stroke_width=3.5)
        ff = tex(r"\times f", h=0.34, color=DETECTOR_ACCENT).next_to(arrow.point_from_proportion(0.5), DOWN, buff=0.14)
        self.play(Create(arrow), run_time=0.55, rate_func=EASE)
        self.play(FadeIn(ff, shift=UP * 0.15), run_time=0.35, rate_func=EASE)
        self.play(Transform(stack, stack_stage(dax, ("VV", "TauTau", "Top", "Z", "Fakes"))),
                  FadeIn(key.rows[KEY_ROW["Fakes"]], shift=LEFT * 0.15), run_time=0.9, rate_func=EASE)
        clip_cut(self, "pipe_f1_data")
        self.zoom_out(DONE["e"], z, about, to)
        self.close_chain()


# ---------------------------------------------------------------------------
# f: comparison
# ---------------------------------------------------------------------------

class PipeCompare(PipeScene):
    """pipe_f_compare (~7.5 s). Zoom on the comparison node: the stack (with its
    key) arrives; data points appear over it; the ratio panel opens slightly
    below 1; the systematic band breathes (tracker-driven redraw, 1 -> wider ->
    1) as each source wiggles it."""
    DONE_BEFORE = DONE["e"]

    def construct(self):
        self.open_chain("pipe_f1_data")
        k, z, to = 5, ZOOM, A["detail"]
        about = node_point(k)
        dax = log_axes(*A["f_plot_wh"], center=A["f_plot"], x_labels=False)
        stack = stack_stage(dax, ("VV", "TauTau", "Top", "Z", "Fakes"))
        key = stack_key().move_to(P(A["f_key"]))
        self.zoom_in(VGroup(dax, stack, key), z, about, to, zoomed_spine(self.st["spine"], z, about, to))
        dots = VGroup(*[data_dot(dax, x, y, color=SAMPLE["Data"], radius=0.05) for x, y in zip(XC, DATA)])
        self.play(LaggedStart(*[GrowFromCenter(d) for d in dots], lag_ratio=0.12, group=dots), run_time=1.1)
        clip_cut(self, "pipe_f2_ratio")
        rp = ratio_panel(dax, EDGES, DATA, TOTAL, center=A["f_ratio"], y_range=(0.9, 1.1, 0.1),
                         y_length=A["f_ratio_h"], dot_radius=0.045, y_ticks=[0.9, 1.0, 1.1],
                         x_title=r"m_{\ell\ell}\ [\mathrm{GeV}]", y_title=r"\mathrm{data}/\mathrm{pred}",
                         title_h=0.22, title_buff=0.16)
        self.play(FadeIn(VGroup(rp.dax, rp.ref)), run_time=0.5, rate_func=EASE)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in rp.dots], lag_ratio=0.12, group=rp.dots), run_time=0.9)
        clip_cut(self, "pipe_f3_uncertainty")
        wdt, op = ValueTracker(1.0), ValueTracker(0.0)
        xs, _ = step_xy(EDGES, TOTAL)
        _, tot = step_xy(EDGES, TOTAL)
        tot = np.array(tot)
        ones = np.ones_like(tot)
        d = A["f_band"]
        band_m = always_redraw(lambda: data_band(dax, xs, tot * (1 - d * wdt.get_value()), tot * (1 + d * wdt.get_value()),
                                                 color=INK, opacity=0.22 * op.get_value()))
        band_r = always_redraw(lambda: data_band(rp.dax, xs, ones * (1 - d * wdt.get_value()),
                                                 ones * (1 + d * wdt.get_value()), color=INK,
                                                 opacity=0.22 * op.get_value()))
        self.add(band_m, band_r)
        self.bring_to_front(dots, rp.dots)
        self.play(op.animate.set_value(1.0), run_time=0.45)
        for v in (1.9, 1.25, 1.7, 1.0):
            self.play(wdt.animate.set_value(v), run_time=0.45, rate_func=EASE)
        band_m.clear_updaters()
        band_r.clear_updaters()
        self.remove(wdt, op)
        clip_cut(self, "pipe_g1_fit_model")
        self.zoom_out(DONE["f"], z, about, to)
        self.close_chain()


# ---------------------------------------------------------------------------
# g: fit
# ---------------------------------------------------------------------------

def sigma_formula() -> VGroup:
    return tex(r"\sigma = \frac{N - B}{A \, \varepsilon \, L}", h=0.40,
               substrings_to_isolate=["N", "B", "A", r"\varepsilon", "L"])


class PipeFit(PipeScene):
    """pipe_g_fit (~8.5 s). Zoom on the fit node: the mu slider (purple
    reference = the prediction), five schematic nuisance parameters theta_i,
    the pre-fit ratio. One play is the fit: mu slides and its bar shrinks,
    the pulls move and tighten, the ratio flattens. Then sigma = mu sigma_pred
    (CONVENTIONS section 2) and sigma = (N - B)/(A eps L) (docs/06), its terms
    lit one by one in the method cyan."""
    DONE_BEFORE = DONE["f"]

    def construct(self):
        self.open_chain("pipe_g1_fit_model")
        k, z, to = 6, ZOOM, A["detail"]
        about = node_point(k)
        sl = slider(0.8, 1.2, 1.0, err=0.15, ref=1.0, length=A["g_slider_len"], ticks=[], marker_r=0.09)
        sl.shift(P(A["g_slider"]) - sl.axis.get_center())
        mu = tex(r"\mu_{Z}", h=0.36).next_to(sl.axis, LEFT, buff=0.35)
        names = [r"\theta_{%d}" % (i + 1) for i in range(5)]
        pp = pull_plot(names, [0.0] * 5, [1.0] * 5, x_length=A["g_pull_len"], row_h=0.44, label_h=0.26,
                       tex=True).move_to(P(A["g_pulls"]))
        post = pull_plot(names, PULLS, CONSTR, x_length=A["g_pull_len"], row_h=0.44, label_h=0.26,
                         tex=True).move_to(P(A["g_pulls"]))
        rw, rh = A["g_ratio_wh"]
        dax_h = schematic_axes(rw, 1.0, center=(A["g_ratio"][0], 3.5))      # geometry only (never added)
        ones = np.ones(12)          # the real 5 GeV fit bins: data/pred pre-fit -> post-fit
        rp = ratio_panel(dax_h, EDGES12, RATIO12_PRE, ones, center=A["g_ratio"], y_range=(0.9, 1.1, 0.1),
                         y_length=rh, dot_radius=0.05, show_x_labels=False, show_y_labels=False)
        rp_post = ratio_panel(dax_h, EDGES12, RATIO12_POST, ones, center=A["g_ratio"], y_range=(0.9, 1.1, 0.1),
                              y_length=rh, dot_radius=0.05, show_x_labels=False, show_y_labels=False)
        r_t = VGroup(axis_titles(rp.dax, r"m_{\ell\ell}", None, h=0.24, buff=0.22),
                     tex(r"\mathrm{data}/\mathrm{pred}", h=0.22).rotate(np.pi / 2).next_to(rp.dax.frame, LEFT, buff=0.3))
        frame = VGroup(pp.band, pp.zero, pp.baseline, pp.ticks)
        detail = VGroup(sl, mu, frame, rp.dax, rp.ref, rp.dots, r_t)
        self.zoom_in(detail, z, about, to, zoomed_spine(self.st["spine"], z, about, to))
        self.play(LaggedStart(*[FadeIn(r) for r in pp.rows], lag_ratio=0.3, group=pp.rows), run_time=0.8)
        clip_cut(self, "pipe_g2_fit")
        self.play(Transform(sl.marker, sl.marker_at(0.97, 0.03)), Transform(pp.rows, post.rows),
                  Transform(rp.dots, rp_post.dots), run_time=1.8, rate_func=EASE)
        clip_cut(self, "pipe_g3_cross_section")
        s1 = tex(r"\sigma = \mu_{Z} \, \sigma_{\mathrm{pred}}", h=0.40).move_to(P(A["g_sig1"]))
        self.play(FadeIn(s1, shift=UP * 0.2), run_time=0.6, rate_func=EASE)
        s2 = sigma_formula().move_to(P(A["g_sig2"]))
        self.play(FadeIn(s2, shift=UP * 0.2), run_time=0.6, rate_func=EASE)
        for t in ("N", "B", "A", r"\varepsilon", "L"):
            part = s2.get_part_by_tex(t)
            self.play(*[g.animate.set_color(col(DETECTOR_ACCENT)) for g in part.family_members_with_points()],
                      run_time=0.3)
        clip_cut(self, "pipe_h_three")
        self.zoom_out(DONE["g"], z, about, to)
        self.close_chain()


# ---------------------------------------------------------------------------
# h: three channels
# ---------------------------------------------------------------------------

class PipeThree(Scene):
    """pipe_h_three (~6 s). The finished spine shrinks into a strip at the top
    (the river fades); two copies slide down; a colour wave runs left to right
    through the three strips (ee green, mumu gold, tautau red via
    CHANNEL_LINE / CHANNEL) and the channel symbols appear. Ends on
    ``state_h()``: the hand-off to sections 3-5."""

    def construct(self):
        white_background(self)
        st = spine_state(DONE["g"])
        add_state(self, st, ORDER_SPINE)
        clip_open(self, "pipe_h_three")
        sp, rv = st["spine"], st["river"]
        shrink = (A["h_scale"], (0.0, SPINE_Y), (A["h_x"], A["h_ys"]["ee"]))
        masked_play(self, Transform(sp, strip(None, "ee")),       # node 0 grazes the top-left block
                    Transform(rv, camera(rv.copy(), *shrink).set_opacity(0.0), remover=True),
                    run_time=1.6, rate_func=EASE)
        mm, tt = sp.copy(), sp.copy()
        self.add(mm, tt)
        self.play(Transform(mm, strip(None, "mumu")), Transform(tt, strip(None, "tautau")), run_time=1.4,
                  rate_func=EASE)
        end = state_h()
        live = {"ee": sp, "mumu": mm, "tautau": tt}
        T = 2.4
        anims = []
        for r, ch in enumerate(("ee", "mumu", "tautau")):
            s, tgt = live[ch], end[ch]
            for i in range(7):
                t0 = 0.12 * i + 0.3 * r
                anims += [Transform(s[0][i], tgt[0][i], rate_func=win(t0, 0.5, T)),
                          Transform(s[2][i], tgt[2][i], rate_func=win(t0, 0.5, T))]
                if i < 6:
                    anims.append(Transform(s[1][i], tgt[1][i], rate_func=win(t0 + 0.06, 0.5, T)))
        labs = [FadeIn(end[f"lab_{ch}"], shift=RIGHT * 0.2, rate_func=win(0.3 * r + 0.1, 0.6, T))
                for r, ch in enumerate(("ee", "mumu", "tautau"))]
        self.play(*anims, *labs, run_time=T)
        final = {**live, **{f"lab_{ch}": end[f"lab_{ch}"] for ch in ("ee", "mumu", "tautau")}}
        check_order(self, final, ORDER_H)
        self.wait(TAIL)
