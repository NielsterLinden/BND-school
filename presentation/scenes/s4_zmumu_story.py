"""Section 4, real data: the Z -> mumu story, chained on 4-02 (``mumu_into_detector``).

Every number printed here is read from ``presentation/data/zmumu_*.json`` (frozen result, tag
``zmumu-freeze-2026-09-17``) and formatted in code; module-level anchor asserts stop the render
on schema drift. Plain ``Scene``s, each opening on the previous clip's final frame: pure builders
``state_*()`` + ``ORDER_*`` tuples, ``finish()`` (strays -> error, z-order set, ``check_order``)
at the end of every clip, seams checked with ``tools/framediff.py`` at delivery.

Re-cut of 17 Sep 2026 (final pass): play order = clip number; tag-and-probe comes first and is
told in full, the other corrections follow in plain words, then rebin -> fit -> cross sections.

    MumuEvent         mumu_event_to_mass     the schematic tracks of 4-02 fade; a real SR event (run 278969)
                                             is drawn from its muons' phi, charge, p_T; the p_T values
                                             converge to m_mumu
                      mumu_first_entry       the slice parks, *linear* count axes appear, m_mumu becomes
                                             the first entry: a block, not a dot
    MumuFill          mumu_fill              four more real events drop blocks (1,0,1,2,0,0,0,0,1: two share
                                             a bin); without a pause the axis becomes logarithmic, the clock
                                             runs and blocks rain over the whole plot: N = 10,378,567
    MumuStack         mumu_prediction        bars -> points; the *uncorrected* prediction slides in under them
                                             (simulation + the data-driven fakes wedge); key
                      mumu_ratio             the ratio panel opens at data/pred = 0.944
    MumuTagProbe      mumu_tnp_tag_probe     the plot parks top right; a real Z pair: the tag (solid gold,
                                             trigger ring) and the probe (dashed grey); m counts up to
                                             90.45 GeV -> the probe is a real muon and turns gold
                      mumu_tnp_pass_fail     the question, tight ID: the probe reached 3 muon stations ->
                                             pass; a second real pair, tracker-only probe, 1 station -> fail;
                                             panels "pass" (twice as high) and "fail", own linear scales; all
                                             probes of the 40-45 GeV barrel cell rain in; fitted peak on the
                                             grey background; eps_data = N_pass / (N_pass + N_fail)
                      mumu_tnp_grid          the two panels move up = row "data"; the identical procedure on
                                             "simulation" below: pass | fail x data / simulation; eps per row
    MumuTagProbeScan  mumu_tnp_scan          the grid parks left; eps_ID vs p_T opens; the p_T window steps
                                             through the ten barrel cells: the probe's curvature, the four
                                             real spectra and both eps change, two points land per cell
                      mumu_tnp_scale_factor  below the eps plot: each data / simulation pair merges into
                                             SF = eps_data / eps_sim
                      mumu_tnp_table         the SF points become the p_T column of a table; a side view of
                                             the detector steps through the |eta| bins, the points follow,
                                             the table fills to the real 10 x 4 map
                      mumu_tnp_apply         table beside the m_mumu plot: three real simulated events pick
                                             their two cells, w = SF x SF, a reweighted block lands; then the
                                             whole prediction steps 0.944 -> 0.971; same method for isolation,
                                             trigger, reconstruction -> 0.968
    MumuCorrections   mumu_corr_pileup       plain-word row, the prediction steps -> 0.974
                      mumu_corr_prefiring    plain-word row -> 0.994
                      mumu_corrected         the rows collapse to "prediction x 0.950": the final data vs
                                             prediction (= the fit input)
    MumuFit           mumu_rebin             60 x 1 GeV -> 12 x 5 GeV, five adjacent bins summed (v2_5_fit.py)
                      mumu_fit               pulls + mu slider, post-fit ratio flat; mu_Z = 0.988 +- 0.014
                      mumu_sigma_fid         sigma_fid = (N - B) / (C L) with our numbers = 794.5 pb (counting);
                                             the fit: 790.2 +- 0.2 +- 6.1 +- 9.6 pb
                      mumu_sigma_total       side view: muons beyond |eta| = 2.4 or too soft are lost, A = 0.409;
                                             sigma_fid / A = 1931 +- 30 pb beside the prediction 1953.9 pb

Physics honesty: the data points never move, only the prediction and the ratio do. The correction
stages are the frozen per-event weights in the order told here (``recut`` block of zmumu_sr_stack.json,
``final`` == the frozen fit input). The first entries use a linear count axis (a count of 1 is below the
log floor 10^1.5); it becomes the log axis when the rain starts. Track curvature follows kappa_from_pt
(exaggerated scale, honest relative curvature). Pass / fail spectra are the real 0.5 GeV fit inputs
summed to 1 GeV by the extractor, each panel on its own linear scale (the fail peak is ~20x smaller);
the pass panel is drawn twice as high. The side view (``eta_fan``) and the ten pairs of the acceptance
clip are schematic; A itself is the frozen number. Words on screen ("pass", "fail", "data", "simulation",
the correction rows) are an exception the deck owner asked for (briefs/zmumu.md). Layout keeps the
title band (y > 2.7) and the top-left block (x < -5.85, y > 0.22) empty.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    DOWN, LEFT, PI, RIGHT, UP, ChangeDecimalToValue, Circle, Create, DashedLine, DashedVMobject, DecimalNumber,
    Dot, FadeIn, FadeOut, GrowFromCenter, LaggedStart, Line, Rectangle, ReplacementTransform, Scene, Succession,
    SurroundingRectangle, Transform, ValueTracker, VGroup, VMobject, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403
from style.bnd_style import _cap_scale, _is_placeholder  # noqa: E402

# ---------------------------------------------------------------------------
# frozen data + anchors (a mismatch stops the render; never adjusted here)
# ---------------------------------------------------------------------------

SR = load_data("zmumu_sr_stack")
EV = load_data("zmumu_events")
FK = load_data("zmumu_fakes")
CO = load_data("zmumu_corrections")
FIT = load_data("zmumu_fit")
TNP = load_data("zmumu_tnp")
RC = SR["recut"]                                   # the correction stages in the order this chapter tells them

assert SR["data"]["total"] == 10378567
assert len(SR["data"]["counts"]) == 60
assert sum(SR["data"]["counts"]) == SR["data"]["total"]
assert abs(FK["yields"]["sr_fakes"] - 3869.955) < 0.01
assert tuple(SR["stack_order"]) == ("Fakes", "WW", "WZ", "ZZ", "SingleTop", "TTbar", "DYtautau", "DYmumu")
assert abs(SR["fakes"]["total"] - FK["yields"]["sr_fakes"]) < 1e-6
assert len(FIT["sr_12bin"]["data"]) == 12 and len(FIT["sr_12bin"]["edges"]) == 13
# ---- the freeze of 17 Sep 2026 (FREEZE.md)
assert abs(FIT["poi"]["value"] - 0.988286) < 1e-6 and round(FIT["poi"]["err_up"], 3) == round(FIT["poi"]["err_down"], 3) == 0.014
assert abs(FIT["sigma_fid"]["value"] - 790.200) < 0.01 and abs(FIT["sigma_fid"]["syst"] - 6.118) < 0.01
assert abs(FIT["sigma_60_120"]["value"] - 1931.043) < 0.01 and round(FIT["sigma_60_120"]["total"]) == 30
assert abs(FIT["sigma_60_120"]["pred"] - 1953.9) < 0.05
assert abs(FIT["sigma_60_120"]["A"] - 0.409209) < 1e-6
assert abs(FIT["C"] - 0.791567) < 1e-6 and abs(FIT["counting"]["n_bkg"] - 68799.02) < 0.01
assert abs(FIT["counting"]["sigma_fid_pb"] - 794.497) < 0.01
STAGES = tuple(RC["stages"])
assert STAGES == ("raw", "sf_id", "sf_muon", "sf_pileup", "final")
LADDER = {k: float(RC["totals"][k]["data_over_pred"]) for k in STAGES}
assert [round(LADDER[k], 3) for k in STAGES] == [0.944, 0.971, 0.968, 0.974, 0.994]
STEP = RC["step_factors"]
MEANS = RC["raw_weighted_means"]
assert [round(STEP[k], 3) for k in ("muon_efficiency", "pileup", "prefiring", "all")] == [0.976, 0.994, 0.980, 0.950]
assert [round(MEANS[k], 3) for k in ("sf_id_pair", "sf_iso_pair", "sf_trigger_event")] == [0.972, 1.007, 0.996]
# ---- tag and probe (unchanged by the freeze)
assert TNP["n_pairs_data"] == 20050129
assert TNP["cell"]["pt_range"] == [40.0, 45.0] and TNP["cell"]["eta_range"] == [0.0, 0.9]
assert abs(TNP["cell"]["data"]["eps"] - 0.958508) < 1e-6 and abs(TNP["cell"]["mc"]["eps"] - 0.971125) < 1e-6
assert TNP["events"]["pass"]["probe"]["passes_id"] and not TNP["events"]["fail"]["probe"]["passes_id"]
assert round(TNP["id"]["map_mean"], 3) == 0.980
CELLS = TNP["cells_barrel"]
assert len(CELLS) == 10 and CELLS[4]["pt_range"] == [40.0, 45.0] and abs(CELLS[4]["data"]["eps"] - 0.958508) < 1e-6
APPLY = TNP["apply_examples"]["events"]
assert [e["rule"] for e in APPLY] == ["barrel", "forward", "overlap_endcap"]

EDGES = np.asarray(SR["edges"], dtype=float)
XC = 0.5 * (EDGES[:-1] + EDGES[1:])
DATA = np.asarray(SR["data"]["counts"], dtype=float)
TOTAL = int(SR["data"]["total"])
FAKES = np.asarray(SR["fakes"]["counts"], dtype=float)
B12 = FIT["sr_12bin"]
EDGES12 = np.asarray(B12["edges"], dtype=float)
XC12 = 0.5 * (EDGES12[:-1] + EDGES12[1:])
DATA12 = np.asarray(B12["data"], dtype=float)
LAYERS = tuple(SR["stack_order"])                  # bottom-up
# the first five entries: real events of zmumu_events.json (sr.events), chosen so that two share a 1 GeV bin
FIRST_EVENTS = [EV["sr"]["events"][i] for i in (0, 14, 4, 8, 12)]
assert FIRST_EVENTS[0]["event"] == EV["sr"]["chosen"][0]["event"]
FIRST_BINS = [int(math.floor(float(r["mass"]))) for r in FIRST_EVENTS]
assert FIRST_BINS == [89, 86, 88, 89, 94]          # -> 1,0,1,2,0,0,0,0,1 from 86 GeV
FIRST_LEVELS = [FIRST_BINS[:k].count(b) for k, b in enumerate(FIRST_BINS)]
CELL = CELLS[4]                                    # 40-45 GeV, |eta| < 0.9: the most populated cell
TP_EVENTS = TNP["events"]                          # "pass": probe passes tight ID, "fail": it does not
PT_EDGES = [float(v) for v in TNP["pt_edges"]]
ETA_EDGES = [float(v) for v in TNP["eta_edges"]]

# ---------------------------------------------------------------------------
# one coordinate dictionary (title band y > 2.7 and top-left block x < -5.85, y > 0.22 stay empty)
# ---------------------------------------------------------------------------

EASE = rate_functions.ease_in_out_sine
PHI_M, PHI_P, KAPPA_402 = math.radians(35), math.radians(212), 0.28   # 4-02 (s4_zmumu.py)
R_LABEL = R_DET + 0.42                             # channel_common.show_signature
Y_SAFE = 2.42                                      # a label centre above this would enter the title band
DET_PARK = (0.42, (-4.6, -0.7))                    # slice parked left of the plot
PLOT_MAIN_C = (1.85, 0.85)                         # main plot area centre (natural)
RATIO_C = (1.85, -1.55)
PLOT_C = np.array([1.85, -0.25, 0.0])              # natural centre every plot placement scales about
PLOT_MAIN = (1.0, PLOT_C)
PLOT_GHOST = (0.30, (5.55, 1.72))                  # the plot parked top right while tag and probe are introduced
PLOT_LEFT = (0.58, (-3.0, -0.45))
Y_EXP = (1.5, 7)                                   # log axis 10^1.5 .. 10^7
LIN_Y = (0, 4, 1)                                  # the linear count axis of the first five entries
RATIO_RANGE = (0.88, 1.12, 0.1)                    # one symmetric range for every stage (checked by the extractor)
RATIO_TICKS = (0.9, 1.0, 1.1)
KEY_LEFT = (-1.25, -2.92)                          # left end of the one-row colour key
COUNTER_AT = (-0.2, 2.12)                          # bottom-left of "N" (natural)
CLOCK_AT = (-0.72, 2.24)
STAMP_AT = (-4.6, -2.3)
LEFT_COL_X = -5.7                                  # left edge of the correction rows (right of the top-left block)
ROW_Y0, ROW_DY = 2.05, 0.95
ROW_RIGHT = -2.55                                  # the rows end left of the plot's y title
DET_TP = (0.80, (-3.25, -0.5))                     # the slice of the tag-and-probe events
DET_TP_PARK = (0.34, (-4.85, -2.62))               # small, bottom left
MASS_AT = (2.3, -0.4)                              # the m(tag, probe) readout of clip tnp_tag_probe
VERDICT_AT = (-3.25, 2.24)                         # "3 muon stations" + mark, above the slice
# -- the 2 x 2 grid (natural = the grid clip); the pass panels are twice as high as the fail panels
GRID = dict(w=2.6, h_pass=1.9, h_fail=0.95, x_pass=1.55, x_fail=4.75, y_data=0.15, y_mc=-2.55)   # y = panel bottoms
GRID_C = np.array([3.15, -0.25, 0.0])              # natural centre of every grid placement
ROW_D_C = np.array([3.15, 1.10, 0.0])              # centre of the data row
GRID_NAT = (1.0, GRID_C)
TP_BIG = (1.2, (3.3, -1.0))                        # the data row alone, larger (clip tnp_pass_fail); about ROW_D_C
GRID_SCAN = (0.6, (-2.9, 0.62))                    # the grid parked left while p_T is scanned
EPS_LONG_C = (-3.2, 1.0)                           # the eps formula of clip tnp_pass_fail
CELL_LAB = dict(big=(5.3, 0.36), grid=(-3.6, 1.55), scan=(3.35, 2.42))
EPS_SCAN = dict(data=(-2.0, -1.72), mc=(-2.0, -2.28))
EFF_C, EFF_WH = (3.35, 1.0), (5.0, 2.2)
EFF_Y = (0.88, 1.00, 0.04)                         # one eps range for all four |eta| columns (0.893 .. 0.991)
SFP_C, SFP_WH = (3.35, -1.45), (5.0, 1.0)
SF_Y = (0.93, 1.00)                                # one SF range for all four columns (0.938 .. 0.990)
SF_TEX_LEFT = (-3.6, -2.75)                        # left end of the SF formula: right of the parked slice (x < -3.85),
                                                   # below the SF panel's y title, left of its first tick label
FAN_T = (0.62, (-3.2, 1.2))                        # the side view of clip tnp_table
ETA_LAB_C = (-3.2, 2.32)
TAB_T = (0.85, (-3.05, -1.55))                     # the table while it fills
TAB_APPLY = (0.85, (-4.37, -1.25))                 # the table beside the plot: the size of TAB_T, it ends left of the
                                                   # plot's y title (ROW_RIGHT); its top stays under the top-left block
DOT_R, RDOT_R = 0.036, 0.034


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def tex_h(expr: str, h: float, color=INK):
    """MathTex whose digits are ``h`` high (uniform size with or without sub/superscripts)."""
    return mathtex(expr, color=color).scale(_cap_scale("tex", h))


def text_h(s: str, h: float, color=INK):
    return text(s, color=color).scale(_cap_scale("text", h))


def _q(rec) -> int:
    return int(np.sign(float(rec["charge"])))


def _p3(p):
    p = np.asarray(p, dtype=float)
    return p if p.shape == (3,) else np.array([p[0], p[1], 0.0])


def tex_int(n) -> str:
    return f"{int(round(float(n))):,}".replace(",", "{,}")


def phi_end(trk) -> float:
    """phi of the outer end of a track built on a natural slice (``.pts`` stay natural)."""
    d = trk.pts[-1] - DET_CENTER
    return math.atan2(d[1], d[0])


def label_spot(phi: float, placement=(1.0, DET_CENTER)) -> np.ndarray:
    """Where a charge label goes: just outside the ring at ``phi`` (as in 4-02); if
    that would put its centre above ``Y_SAFE`` (the title band), slide it along the
    same circle toward the horizontal until it fits."""
    s, c = placement
    c = _p3(c)
    r = s * R_LABEL

    def at(a):
        return c + r * np.array([math.cos(a), math.sin(a), 0.0])

    p = at(phi)
    if p[1] <= Y_SAFE:
        return p
    a0 = math.asin(max(-1.0, min(1.0, (Y_SAFE - c[1]) / r)))
    ph = phi % (2 * math.pi)
    return at(min((a0, math.pi - a0), key=lambda a: abs(a - ph)))


# -- the verified swap, regrouping, and the end-of-clip guard ----------------------------

def _drawn(s) -> bool:
    """Does this single family member put ink on the frame?"""
    if not s.has_points():
        return False
    if not isinstance(s, VMobject):
        return True
    fo = float(np.max(s.get_fill_opacity()))
    sw = float(np.max(s.get_stroke_width()))
    so = float(np.max(s.get_stroke_opacity()))
    pts = s.points[:, :2]
    area = float(np.prod(pts.max(axis=0) - pts.min(axis=0)))
    return (fo > 0.01 and area > 1e-7) or (sw > 0.01 and so > 0.01)


def _look(m):
    """What a mobject draws: per family member with points, (style, bbox)."""
    out = []
    for s in m.get_family():
        if not _drawn(s):
            continue
        if isinstance(s, VMobject):
            fo = float(np.max(s.get_fill_opacity()))
            sw = float(np.max(s.get_stroke_width()))
            so = float(np.max(s.get_stroke_opacity()))
        else:
            fo, sw, so = 1.0, 0.0, 0.0
        pts = s.points[:, :2]
        lo, hi = pts.min(axis=0), pts.max(axis=0)
        fill_vis = fo > 0.01 and float(np.prod(hi - lo)) > 1e-7
        stroke_vis = sw > 0.01 and so > 0.01
        style = (s.get_fill_color().to_hex() if fill_vis else "-", round(fo, 2) if fill_vis else 0,
                 s.get_stroke_color().to_hex() if stroke_vis else "-",
                 round(sw, 2) if stroke_vis else 0, round(so, 2) if stroke_vis else 0)
        out.append((style, np.concatenate([lo, hi])))
    return out


def _match(a, b, tol):
    if len(a) == 0 or len(b) == 0:
        return len(a) == len(b)
    d = np.abs(a[:, None, :] - b[None, :, :]).max(axis=2)
    return bool(np.all(d.min(axis=1) <= tol) and np.all(d.min(axis=0) <= tol))


def looks_same(live, target, tol=6e-3) -> tuple[bool, str]:
    la, lb = _look(live), _look(target)
    sa, sb = {s for s, _ in la}, {s for s, _ in lb}
    if sa != sb:
        return False, f"styles differ: {sorted(sa ^ sb)[:4]}"
    for st in sa:
        A = np.array([bb for s, bb in la if s == st])
        B = np.array([bb for s, bb in lb if s == st])
        if not _match(A, B, tol):
            return False, f"geometry differs for style {st}"
    return True, ""


def settle(scene, live, target, what=""):
    """Swap a live top-level mobject for the builder's (identity for ``finish`` / the next
    clip) after proving they draw the same thing (same styles, same member boxes within
    6e-3 units; coincident duplicates left by a 5 -> 1 Transform count once)."""
    ok, why = looks_same(live, target)
    if not ok:
        raise AssertionError(f"settle({what}): live mobject does not match the builder: {why}")
    scene.replace(live, target)


def settle_same(scene, prev, end, keys):
    for k in keys:
        settle(scene, prev[k], end[k], k)


def _belongs(m, fam) -> bool:
    if id(m) in fam:
        return True
    pts = [x for x in m.get_family() if x.has_points()]
    return bool(pts) and all(id(x) in fam for x in pts)


def regroup(scene, parent, what=""):
    """Children of ``parent`` revealed one by one (or left over when a sibling was faded
    out) sit top-level in ``scene.mobjects``: replace them by ``parent`` itself, at the
    z-position of the first piece. Every *drawn* member of ``parent`` must be on stage."""
    fam = {id(m) for m in parent.get_family()}
    idx = [i for i, m in enumerate(scene.mobjects) if _belongs(m, fam)]
    if not idx:
        raise AssertionError(f"regroup({what}): nothing of it is in the scene")
    shown = {id(x) for i in idx for x in scene.mobjects[i].get_family()}
    missing = [x for x in parent.get_family() if _drawn(x) and id(x) not in shown]
    if missing:
        raise AssertionError(f"regroup({what}): {len(missing)} drawn members were never revealed")
    drop = set(idx)
    keep = [m for i, m in enumerate(scene.mobjects) if i not in drop]
    keep.insert(idx[0], parent)
    scene.mobjects = keep


def swap_family(scene, old, new, what=""):
    """``old`` was taken apart on stage (some children faded out) while children of ``new`` were
    faded in: drop every top-level piece of either and put ``new`` where ``old`` began."""
    fo, fn = {id(m) for m in old.get_family()}, {id(m) for m in new.get_family()}
    idx = [i for i, m in enumerate(scene.mobjects) if _belongs(m, fo)]
    if not idx:
        raise AssertionError(f"swap_family({what}): nothing of the old object is in the scene")
    keep = []
    for i, m in enumerate(scene.mobjects):
        if i == idx[0]:
            keep.append(new)
        if _belongs(m, fo) or _belongs(m, fn):
            continue
        keep.append(m)
    scene.mobjects = keep


def finish(scene, state, order):
    """End of a clip: every top-level mobject that draws something must be one of the state's
    objects (after ``settle`` / ``regroup``); invisible leftovers are dropped; the z-order is
    set to ``order`` and verified with ``check_order``."""
    want = [state[k] for k in order]
    ids = {id(m) for m in want}
    live = [m for m in scene.mobjects if not _is_placeholder(m) and not isinstance(m, ValueTracker)]
    stray = [m for m in live if id(m) not in ids and _look(m)]
    if stray:
        by_kind = sorted({type(m).__name__ for m in stray})
        raise AssertionError(f"finish: {len(stray)} stray drawn mobject(s) not in the end state: {by_kind}")
    on_stage = {id(m) for m in live}
    missing = [k for k in order if id(state[k]) not in on_stage]
    if missing:
        raise AssertionError(f"finish: end-state objects not on stage: {missing}")
    scene.mobjects = list(want)
    check_order(scene, state, order)


def hidden_rain(*args, **kw) -> list:
    """``rain()`` with each mark invisible until its own flight starts: a short FadeIn
    is prepended to every Succession (LaggedStart begins all children up front, so
    the FadeIn's start state hides the mark; without it all marks sit in the tracker
    as one dark blob from the first frame)."""
    out = []
    for a in rain(*args, **kw):
        move, fade = a.animations
        out.append(Succession(FadeIn(move.mobject, run_time=0.06), move, fade))
    return out


def pulse(scene, mobs, color=DETECTOR_ACCENT, scale=1.7, run_time=0.55):
    """Look here: coloured copies grow and fade (a copy, so nothing is left at peak scale)."""
    ghosts = VGroup(*[m.copy().set_color(col(color)) for m in mobs])
    scene.add(ghosts)
    scene.play(*[g.animate.scale(scale).set_opacity(0.0) for g in ghosts], run_time=run_time)
    scene.remove(ghosts)


# ---------------------------------------------------------------------------
# plot builders (always at the natural position, then place())
# ---------------------------------------------------------------------------

def sr_layers(stage):
    """The eight stack layers bottom-up at a ``recut`` stage; ``None`` = everything at the floor.
    The data-driven fakes are part of the prediction from its first appearance."""
    rows = []
    for n in LAYERS:
        if stage is None:
            c = np.zeros(60)
        elif n == "Fakes":
            c = FAKES
        else:
            c = np.asarray(RC["mc"][n][stage], dtype=float)
        rows.append((n, c, SAMPLE[n]))
    return rows


def fit_layers(kind):
    return [(n, np.asarray(B12[kind][n], dtype=float), SAMPLE[n]) for n in LAYERS]


def main_axes(per5=False, x_labels=False, linear=False) -> DataAxes:
    """The m_mumu axes. ``linear``: counts 0 .. 4 (the first entries). Otherwise log from
    10^1.5; ``per5``: counts per 5 GeV, the range shifted by log10(5) so the rebinned stack
    keeps its height on screen and the relabelled ticks carry the factor 5. ``x_labels``:
    tick labels + title on the main plot itself (before the ratio panel takes them over)."""
    common = dict(show_x_labels=x_labels, title_h=0.2, title_buff=0.16,
                  x_title=r"m_{\mu\mu}\ [\mathrm{GeV}]" if x_labels else None)
    if linear:
        dax = DataAxes([60, 120, 10], list(LIN_Y), 6.4, 2.9, y_ticks=[1, 2, 3, 4],
                       y_title=r"\mathrm{events}\,/\,\mathrm{GeV}", **common)
    else:
        off = math.log10(5.0) if per5 else 0.0
        e0, e1 = Y_EXP[0] + off, Y_EXP[1] + off
        dax = DataAxes([60, 120, 10], [e0, e1, 1], 6.4, 2.9, y_ticks=[k for k in range(1, 9) if e0 < k <= e1 + 1e-9],
                       y_log=True, y_title=r"\mathrm{events}\,/\,5\,\mathrm{GeV}" if per5 else r"\mathrm{events}\,/\,\mathrm{GeV}",
                       **common)
    dax.move_frame_to(PLOT_MAIN_C)
    return dax


def entry_block(dax, m, level=0) -> Rectangle:
    """One event as one block in its 1 GeV bin, ``level`` blocks already below it (linear axis)."""
    return data_bar(dax, math.floor(float(m)) + 0.5, level + 1.0, 0.5, color=SAMPLE["Data"], base=float(level),
                    fill_opacity=0.55, stroke_width=1.0)


def data_bars(dax, visible=True) -> VGroup:
    if visible:
        return VGroup(*[data_bar(dax, XC[i], DATA[i], 0.5, color=SAMPLE["Data"], fill_opacity=0.55,
                                 stroke_width=1.0) for i in range(60)])
    # the same bars at the floor, fully transparent (Transform source)
    return VGroup(*[data_bar(dax, XC[i], 0.0, 0.5, color=SAMPLE["Data"], fill_opacity=0.55,
                             stroke_width=1.0).set_opacity(0.0) for i in range(60)])


def data_dots(dax, bins=60) -> VGroup:
    xc, n = (XC, DATA) if bins == 60 else (XC12, DATA12)
    return VGroup(*[data_dot(dax, xc[i], n[i], color=SAMPLE["Data"], radius=DOT_R) for i in range(len(n))])


def layout_counter(pre, num):
    num.next_to(pre, RIGHT, buff=0.14)
    num.shift(UP * (pre.get_bottom()[1] - num[0].get_bottom()[1]))
    return num


def make_counter(value) -> VGroup:
    pre = tex_h(r"N =", 0.24)
    pre.shift(_p3(COUNTER_AT) - pre.get_corner(DOWN + LEFT))
    num = DecimalNumber(float(value), num_decimal_places=0, group_with_commas=True, color=col(INK))
    num.scale(_cap_scale("tex", 0.24))
    layout_counter(pre, num)
    g = VGroup(pre, num)
    g.pre, g.num = pre, num
    return g


def set_counter(counter, value):
    counter.num.set_value(float(round(value)))
    layout_counter(counter.pre, counter.num)


def make_key() -> VGroup:
    entries = [(r"\mathrm{Data}", SAMPLE["Data"], "dot"), (r"Z\to\mu\mu", SAMPLE["DYmumu"]),
               (r"\tau\tau,\ t\bar{t},\ tW,\ VV", SAMPLE["TTbar"]), (r"\mathrm{Fakes}", SAMPLE["Fakes"])]
    key = colour_key(entries, label_h=0.2)
    # the background row: one stripe per simulated background colour
    row = key.rows[2]
    sw = row[1]
    names = ("DYtautau", "TTbar", "SingleTop", "WW", "WZ", "ZZ")
    w, h = sw.width, sw.height
    stripes = VGroup(*[Rectangle(width=w / len(names), height=h, stroke_width=0, fill_color=col(SAMPLE[n]),
                                 fill_opacity=0.95) for n in names]).arrange(RIGHT, buff=0).move_to(sw)
    frame = Rectangle(width=w, height=h, fill_opacity=0.0, stroke_color=darken(SAMPLE["TTbar"], 0.25),
                      stroke_width=1.0).move_to(sw)
    row.submobjects[1] = VGroup(stripes, frame)
    key.rows.arrange(RIGHT, buff=0.42)
    key.rows.shift(_p3(KEY_LEFT) - key.rows.get_left())
    return key


def plot_parts(placement=PLOT_MAIN, *, stage=None, data="dots", bins=60, fit="prefit", counter=True, ratio=True,
               rlabel=True, key=True, entries=(), x_labels=False, linear=False) -> dict:
    """The plot family of the chain, built at the natural position and then placed about
    PLOT_C: dax, stack (8 layers bottom-up), data (blocks / bars / dots), counter, ratio
    (panel with the data stat error bars), rlabel (total data/pred of the stage), key."""
    P = {}
    dax = main_axes(per5=bins == 12, x_labels=x_labels, linear=linear)
    P["dax"] = dax
    if bins == 60:
        P["stack"] = stack_hist(dax, EDGES, sr_layers(stage))
    else:
        P["stack"] = stack_hist(dax, EDGES12, fit_layers(fit))
    if data == "blocks":
        P["data"] = VGroup(*[entry_block(dax, m, lv) for m, lv in entries])
    elif data == "bars":
        P["data"] = data_bars(dax)
    else:
        P["data"] = data_dots(dax, bins)
    if counter:
        P["counter"] = make_counter(TOTAL)
    if ratio:
        kw = dict(y_range=RATIO_RANGE, y_length=0.9, y_ticks=list(RATIO_TICKS), dot_radius=RDOT_R,
                  x_title=r"m_{\mu\mu}\ [\mathrm{GeV}]", y_title=r"\mathrm{data/pred.}",
                  title_h=0.2, title_buff=0.16)
        if bins == 60:
            P["ratio"] = ratio_panel(dax, EDGES, DATA, P["stack"].total, RATIO_C, num_err=np.sqrt(DATA), **kw)
        else:
            r12 = np.asarray(B12[f"ratio_{fit}"], dtype=float)
            P["ratio"] = ratio_panel(dax, EDGES12, r12, np.ones(12), RATIO_C, num_err=r12 / np.sqrt(DATA12), **kw)
        if rlabel and bins == 60:
            v = LADDER[stage]
            lab = tex_h(f"{v:.3f}", 0.2)
            lab.next_to(P["ratio"].dax.c2p(120, v), RIGHT, buff=0.22)
            P["rlabel"] = lab
    if key:
        P["key"] = make_key()
    for m in P.values():
        place(m, placement, PLOT_C)
    return P


PLOT_KEYS = ("dax", "stack", "data", "counter", "ratio", "rlabel", "key")


def move_ratio(scene, live, stage, run_time=1.6):
    """The prediction steps to ``stage``: stack, ratio dots + error bars and the label. The data never move."""
    S = plot_parts(stage=stage, key=False, counter=False)
    scene.play(Transform(live["stack"], S["stack"]), Transform(live["ratio"].dots, S["ratio"].dots),
               Transform(live["ratio"].errs, S["ratio"].errs), Transform(live["rlabel"], S["rlabel"]),
               run_time=run_time, rate_func=EASE)


# ---------------------------------------------------------------------------
# detector builders
# ---------------------------------------------------------------------------

def parked_slice(placement=DET_PARK) -> CMSSlice:
    return place(CMSSlice(), placement, DET_CENTER)


def event_group(rec) -> VGroup:
    """A real SR event on a natural slice: the muons (event_from_json) and any FSR
    photon attached to them (straight, to the ECAL). Place it afterwards."""
    det0 = CMSSlice()
    ev = event_from_json(det0, rec)
    photons = VGroup(*[signature(det0, "gamma", float(ph["phi"])) for ph in rec.get("fsr_photons", [])])
    g = VGroup(ev, photons)
    g.ev, g.photons = ev, photons
    return g


def run_stamp(rec) -> VGroup:
    s = text_h(f"Run {rec['run']}   Event {rec['event']}", 0.14, color=GREY)
    return s.move_to(_p3(STAMP_AT))


def mass_tex(rec):
    return tex_h(rf"m_{{\mu\mu}} = {float(rec['mass']):.1f}\ \mathrm{{GeV}}", 0.28)


# ---------------------------------------------------------------------------
# the correction rows: plain words first, the factor below, the technical name small
# ---------------------------------------------------------------------------

def mini_table(n_rows=5, n_cols=4, w=0.72, h=0.5, color=SLATE, strength=0.35) -> VGroup:
    """A table as a glyph (no values): the icon of "a scale-factor table"."""
    cw, ch = w / n_cols, h / n_rows
    g = VGroup()
    for i in range(n_rows):
        for j in range(n_cols):
            a = strength * (0.55 + 0.45 * ((i * 7 + j * 3) % 5) / 4.0)
            g.add(Rectangle(width=cw, height=ch, stroke_color=col(WHITE), stroke_width=0.8,
                            fill_color=col(mix(WHITE, col(color).to_hex(), a)), fill_opacity=1.0)
                  .move_to([(j - (n_cols - 1) / 2) * cw, ((n_rows - 1) / 2 - i) * ch, 0.0]))
    return g


def icon_pileup(w=0.72) -> VGroup:
    """Several collisions in one bunch crossing: one is ours (gold), the others grey."""
    line = Line(LEFT * w / 2, RIGHT * w / 2, stroke_color=col(GREY), stroke_width=2.0)
    xs = np.linspace(-0.4, 0.4, 7) * w
    dots = VGroup(*[Dot([x, 0.0, 0.0], radius=0.045, color=col(GREY)) for x in xs])
    dots[3].set_color(col(CHANNEL_LINE["mumu"])).scale(1.5)
    return VGroup(line, dots)


def icon_prefiring(w=0.72) -> VGroup:
    """Bunch crossings as ticks on a time line; the trigger fires one tick too early."""
    line = Line(LEFT * w / 2, RIGHT * w / 2, stroke_color=col(GREY), stroke_width=2.0)
    xs = np.array([-0.3, 0.0, 0.3]) * w
    ticks = VGroup(*[Line([x, -0.09, 0], [x, 0.09, 0], stroke_color=col(GREY), stroke_width=2.0) for x in xs])
    ours = Dot([xs[1], 0.0, 0.0], radius=0.055, color=col(CHANNEL_LINE["mumu"]))
    fired = Circle(radius=0.085, stroke_color=col(DETECTOR_ACCENT), stroke_width=3.0).move_to([xs[0], 0.0, 0.0])
    return VGroup(line, ticks, ours, fired)


def word_row(y, words, factor, note=None, icon=None) -> VGroup:
    """[icon] words / (technical name, small) / x factor: left-aligned at LEFT_COL_X, never wider than
    ROW_RIGHT (the plot's y title begins there)."""
    g = VGroup()
    x = LEFT_COL_X
    y_f = y - (0.62 if note else 0.42)
    if icon is not None:
        icon.move_to([x + 0.36, 0.5 * (y + y_f), 0.0])
        g.add(icon)
        x += 0.72 + 0.18
    w = text_h(words, 0.19)
    w.shift(np.array([x, y, 0.0]) - w.get_left())
    f = tex_h(rf"\times\,{factor}", 0.24)
    f.shift(np.array([x, y_f, 0.0]) - f.get_left())
    g.add(w, f)
    if note:
        nt = text_h(note, 0.13, color=GREY)
        nt.shift(np.array([x, y - 0.28, 0.0]) - nt.get_left())
        g.add(nt)
    assert g.get_right()[0] <= ROW_RIGHT + 1e-6, f"row {words!r} reaches x = {g.get_right()[0]:.2f}, past {ROW_RIGHT}"
    g.words, g.factor = w, f
    return g


def row_y(slot: int) -> float:
    return ROW_Y0 - ROW_DY * slot


def row_mu():
    return word_row(row_y(0), "muon efficiency", f"{STEP['muon_efficiency']:.3f}", "tag and probe", mini_table())


def row_pu():
    return word_row(row_y(1), "extra collisions", f"{STEP['pileup']:.3f}", "pile-up", icon_pileup())


def row_l1():
    return word_row(row_y(2), "trigger fired early", f"{STEP['prefiring']:.3f}", "prefiring", icon_prefiring())


def row_all():
    return word_row(row_y(0), "prediction", f"{STEP['all']:.3f}")


SUB_DY = 0.8


def sub_rows() -> list:
    """The four tag-and-probe measurements behind "muon efficiency", same method each."""
    reco = TNP["reco"]["sf_reco_per_event"]
    spec = [("identification", f"{MEANS['sf_id_pair']:.3f}"), ("isolation", f"{MEANS['sf_iso_pair']:.3f}"),
            ("trigger", f"{MEANS['sf_trigger_event']:.3f}"), ("reconstruction", f"{reco:.4f}")]
    return [word_row(ROW_Y0 - SUB_DY * k, w, f, icon=mini_table()) for k, (w, f) in enumerate(spec)]


# ---------------------------------------------------------------------------
# tag and probe
# ---------------------------------------------------------------------------

def tp_slice(placement=DET_TP) -> CMSSlice:
    return place(CMSSlice(), placement, DET_CENTER)


def tp_pair(rec, placement=DET_TP) -> VGroup:
    """A real opposite-sign tag + probe pair (zmumu_tnp.json events) on a slice: two
    ``muon_pieces``, each keeping only the muon stations it was matched in (``nStations``:
    the tracker-only probe has one). ``.tag``, ``.probe``, ``.tag_phi``, ``.probe_phi``
    (natural), ``.record``."""
    det0 = CMSSlice()
    out = []
    for m in (rec["tag"], rec["probe"]):
        p = muon_pieces(det0, float(m["phi"]), _q(m), kappa_from_pt(m["pt"]))
        stubs = p.deposits[1]
        stubs.submobjects = stubs.submobjects[:max(0, min(4, int(m["nStations"])))]
        out.append(p)
    t, p = out
    g = VGroup(t, p)
    g.tag, g.probe, g.record = t, p, rec
    g.tag_phi, g.probe_phi = phi_end(t.trk), phi_end(p.trk)
    place(g, placement, DET_CENTER)
    return g


def tp_label(pair, which, color, placement=DET_TP):
    phi = pair.tag_phi if which == "tag" else pair.probe_phi
    lab = tex_h(rf"\mu_{{\mathrm{{{which}}}}}", 0.24, color=color)
    spot = label_spot(phi, placement)
    c = _p3(placement[1])
    u = (spot - c) / max(np.linalg.norm(spot - c), 1e-9)
    return lab.move_to(spot + u * 0.5 * max(lab.width - 0.3, 0.0))


def probe_asked(pair) -> DashedVMobject:
    """The probe before anything is known about it: a dashed grey line along its track."""
    core = pair.probe.trk[-1].copy().set_stroke(color=col(GREY), width=3.0)
    return DashedVMobject(core, num_dashes=26, dashed_ratio=0.55)


def trigger_ring(tag) -> Circle:
    """The tag fired the trigger: a cyan ring on its outermost muon-station stub."""
    stub = tag.deposits[1][-1]
    return Circle(radius=0.16, stroke_color=col(DETECTOR_ACCENT), stroke_width=4.0).move_to(stub.get_center())


def mass_readout(value) -> VGroup:
    """m(tag, probe) as a running number (VGroup(pre, num, unit), ``.num`` a DecimalNumber)."""
    pre = tex_h(r"m_{\mu\mu} =", 0.28)
    num = DecimalNumber(float(value), num_decimal_places=2, color=col(INK)).scale(_cap_scale("tex", 0.28))
    unit = tex_h(r"\mathrm{GeV}", 0.28)
    num.next_to(pre, RIGHT, buff=0.16)
    num.shift(UP * (pre[0][-1].get_center()[1] - num.get_center()[1]))
    unit.next_to(num, RIGHT, buff=0.16).align_to(num, DOWN)
    g = VGroup(pre, num, unit)
    g.num = num
    return g.move_to(_p3(MASS_AT))


def verdict(n_stations: int, passes: bool) -> VGroup:
    """Why the probe passes or fails tight ID, as seen in the slice: the muon stations it reached."""
    words = text_h(f"{n_stations} muon station" + ("s" if n_stations != 1 else ""), 0.2)
    mark = tex_h(r"\checkmark" if passes else r"\times", 0.26, color=CHANNEL_LINE["mumu"] if passes else INK)
    g = VGroup(words, mark).arrange(RIGHT, buff=0.22)
    g.words, g.mark = words, mark
    return g.move_to(_p3(VERDICT_AT))


def tp_panel(sample: str, which: str, cell, placement=GRID_NAT, nc=GRID_C) -> VGroup:
    """One m(tag, probe) panel of the grid: the real spectrum of ``cell`` (data: bars, simulation: gold
    outline), the fitted background as a grey band and the fit model (method cyan), on its own linear
    scale; pass panels are twice as high as fail panels. ``.ax``, ``.body``, ``.fit`` (= band, line)."""
    d = cell["data" if sample == "data" else "mc"]
    h = GRID["h_pass" if which == "pass" else "h_fail"]
    cx = GRID["x_pass" if which == "pass" else "x_fail"]
    by = GRID["y_data" if sample == "data" else "y_mc"]
    spec, model, bkg = (np.asarray(d[k], dtype=float) for k in (which, f"model_{which}", f"bkg_{which}"))
    ymax = 1.12 * max(spec.max(), model.max())
    dax = DataAxes([60, 120, 30], [0.0, ymax, ymax], GRID["w"], h, y_ticks=[], show_y_labels=False,
                   tick_label_h=0.15, title_h=0.16, title_buff=0.1, x_title=r"m_{\mu\mu}\ [\mathrm{GeV}]")
    dax.move_frame_to((cx, by + h / 2))
    if sample == "data":
        body = VGroup(*[data_bar(dax, XC[i], spec[i], 0.5, color=SAMPLE["Data"], fill_opacity=0.55, stroke_width=1.0)
                        for i in range(60)])
    else:
        body = step_hist(dax, EDGES, spec, color=CHANNEL_LINE["mumu"], stroke_width=3.0)
    band = data_band(dax, XC, np.zeros(60), bkg, color=GREY, opacity=0.6)
    line = data_trace(dax, XC, model, color=DETECTOR_ACCENT, stroke_width=3.0)
    fit = VGroup(band, line)
    # simulation: the fit describes it so well that the model line would hide the gold outline -> outline on top
    g = VGroup(dax, body, fit) if sample == "data" else VGroup(dax, fit, body)
    g.ax, g.body, g.fit, g.band, g.line = dax, body, fit, band, line
    return place(g, placement, nc)


def flat_bars(panel) -> VGroup:
    """The data bars of ``panel`` at the floor, transparent (Transform source)."""
    out = VGroup()
    for b in panel.body:
        f = b.copy().stretch_to_fit_height(1e-3).align_to(b, DOWN).set_opacity(0.0)
        out.add(f)
    return out


PANELS = (("data", "pass"), ("data", "fail"), ("mc", "pass"), ("mc", "fail"))


def pkey(sample, which) -> str:
    return f"p_{sample}_{which}"


def grid_panels(cell, placement=GRID_NAT, samples=("data", "mc"), nc=GRID_C) -> dict:
    return {pkey(s, w): tp_panel(s, w, cell, placement, nc) for s, w in PANELS if s in samples}


def col_words(placement=GRID_NAT, nc=GRID_C, row="mc") -> dict:
    """"pass" / "fail" under the panels of the lowest row on screen."""
    y = GRID["y_data" if row == "data" else "y_mc"] - 1.0
    out = {}
    for which in ("pass", "fail"):
        w = text_h(which, 0.24).move_to([GRID["x_pass" if which == "pass" else "x_fail"], y, 0.0])
        out[f"w_{which}"] = place(w, placement, nc)
    return out


def row_words(placement=GRID_NAT) -> dict:
    out = {}
    for s, name in (("data", "data"), ("mc", "simulation")):
        h = GRID["h_pass"]
        w = text_h(name, 0.2, color=INK if s == "data" else CHANNEL_LINE["mumu"]).rotate(PI / 2)
        w.move_to([GRID["x_pass"] - GRID["w"] / 2 - 0.42, GRID["y_data" if s == "data" else "y_mc"] + h / 2, 0.0])
        out[f"rw_{s}"] = place(w, placement, GRID_C)
    return out


def tp_count(panel, which: str):
    d = CELL["data"]
    n = d["n_pass_signal"] if which == "pass" else d["n_fail_signal"]
    sub = r"\mathrm{pass}" if which == "pass" else r"\mathrm{fail}"
    return tex_h(rf"N_{{{sub}}} = {tex_int(n)}", 0.2).next_to(panel.ax.frame, UP, buff=0.14)


def cell_label(cell, where: str):
    """The cell being measured, two lines (p_T window, |eta| window)."""
    pt0, pt1 = cell["pt_range"]
    a = tex_h(rf"{pt0:g} < p_T < {pt1:g}\ \mathrm{{GeV}}", 0.19)
    b = tex_h(rf"|\eta| < {cell['eta_range'][1]:g}", 0.19)
    if where == "scan":
        g = VGroup(a, b).arrange(RIGHT, buff=0.5)
    else:
        g = VGroup(a, b).arrange(DOWN, buff=0.14)
    return g.move_to(_p3(CELL_LAB[where]))


def eps_long():
    d = CELL["data"]
    expr = (r"\varepsilon_{\mathrm{data}} = \frac{N_{\mathrm{pass}}}{N_{\mathrm{pass}} + N_{\mathrm{fail}}}"
            rf" = {d['eps']:.4f}")
    return tex_h(expr, 0.22).move_to(_p3(EPS_LONG_C))


def eps_readout(sample: str, value: float, h: float, center) -> VGroup:
    """eps_data / eps_sim = 0.xxxx as VGroup(pre, num): the number can tick (ChangeDecimalToValue)."""
    sub, colour = (r"\mathrm{data}", INK) if sample == "data" else (r"\mathrm{sim}", CHANNEL_LINE["mumu"])
    pre = tex_h(rf"\varepsilon_{{{sub}}} =", h, color=colour)
    num = DecimalNumber(float(value), num_decimal_places=4, color=col(colour)).scale(_cap_scale("tex", h))
    num.next_to(pre, RIGHT, buff=0.55 * h)
    num.shift(UP * (pre[0][-1].get_center()[1] - num.get_center()[1]))
    g = VGroup(pre, num)
    g.num = num
    return g.move_to(_p3(center))


def eps_in_grid(sample, placement=GRID_NAT):
    """The efficiency of a row, in the free space above its (half-height) fail panel."""
    y = GRID["y_data" if sample == "data" else "y_mc"] + GRID["h_fail"] + 0.5 * (GRID["h_pass"] - GRID["h_fail"])
    g = eps_readout(sample, CELL["data" if sample == "data" else "mc"]["eps"], 0.22, (GRID["x_fail"], y))
    return place(g, placement, GRID_C)


def scan_track(i: int) -> VGroup:
    """The probe of p_T cell ``i`` in the parked slice: only its curvature changes."""
    pt = 0.5 * (PT_EDGES[i] + PT_EDGES[i + 1])
    trk = muon_pieces(CMSSlice(), math.radians(68), -1, kappa_from_pt(pt)).trk
    return place(trk, DET_TP_PARK, DET_CENTER)


def eff_axes(x_title=True) -> DataAxes:
    dax = DataAxes([20, 200, 30], list(EFF_Y), EFF_WH[0], EFF_WH[1], x_ticks=[20, 50, 100, 150, 200],
                   y_ticks=[0.88, 0.92, 0.96, 1.00], y_fmt="{:.2f}", tick_label_h=0.16, title_h=0.2, title_buff=0.14,
                   x_title=r"p_T\ [\mathrm{GeV}]", y_title=r"\varepsilon_{\mathrm{ID}}")
    dax.move_frame_to(EFF_C)
    if not x_title:                      # the SF panel below owns the p_T title from then on
        dax.x_title.set_opacity(0.0)
        dax.x_labels.set_opacity(0.0)
    return dax


def _bin_point(dax, i, value, colour) -> VGroup:
    return VGroup(Line(dax.c2p(PT_EDGES[i], value), dax.c2p(PT_EDGES[i + 1], value), stroke_color=col(colour), stroke_width=2.5),
                  Dot(dax.c2p(0.5 * (PT_EDGES[i] + PT_EDGES[i + 1]), value), radius=0.045, color=col(colour)))


def eff_points(dax, sample: str, j: int) -> VGroup:
    E = TNP["id"]["eff_data" if sample == "data" else "eff_mc"]
    colour = SAMPLE["Data"] if sample == "data" else CHANNEL_LINE["mumu"]
    return VGroup(*[_bin_point(dax, i, float(E[i][j]), colour) for i in range(10)])


def eff_key(dax) -> VGroup:
    key = colour_key([(r"\mathrm{data}", SAMPLE["Data"], "dot"), (r"\mathrm{simulation}", CHANNEL_LINE["mumu"], "dot")],
                     label_h=0.16)
    return key.next_to(dax.frame.get_corner(DOWN + LEFT), UP + RIGHT, buff=0.12)


def pt_band(dax, i: int) -> Rectangle:
    """The p_T bin being measured, lit over the full height of the eps plot."""
    lo, hi = dax.c2p(PT_EDGES[i], EFF_Y[0]), dax.c2p(PT_EDGES[i + 1], EFF_Y[1])
    return Rectangle(width=max(hi[0] - lo[0], 0.06), height=hi[1] - lo[1], stroke_width=0, fill_color=col(DETECTOR_ACCENT),
                     fill_opacity=0.22).move_to([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, 0.0])


def sf_axes() -> DataAxes:
    dax = DataAxes([20, 200, 30], [SF_Y[0], SF_Y[1], 0.035], SFP_WH[0], SFP_WH[1], x_ticks=[20, 50, 100, 150, 200],
                   y_ticks=[0.95, 1.00], y_fmt="{:.2f}", tick_label_h=0.16, title_h=0.2, title_buff=0.14,
                   x_title=r"p_T\ [\mathrm{GeV}]", y_title=r"\mathrm{SF}")
    dax.move_frame_to(SFP_C)
    return dax


def sf_points(dax, j: int) -> VGroup:
    S = TNP["id"]["sf"]
    return VGroup(*[_bin_point(dax, i, float(S[i][j]), DETECTOR_ACCENT) for i in range(10)])


def sf_tex():
    d, m = CELL["data"]["eps"], CELL["mc"]["eps"]
    expr = (r"\mathrm{SF} = \frac{\varepsilon_{\mathrm{data}}}{\varepsilon_{\mathrm{sim}}}"
            rf" = \frac{{{d:.4f}}}{{{m:.4f}}} = {d / m:.3f}")
    t = tex_h(expr, 0.2)
    return t.shift(_p3(SF_TEX_LEFT) - t.get_left())


def sf_map(placement=TAB_T) -> VGroup:
    """The real 10 x 4 tight-ID scale-factor map (rows p_T, columns |eta|)."""
    rws = [rf"{PT_EDGES[i]:g}\!-\!{PT_EDGES[i + 1]:g}" for i in range(10)]
    cls = [rf"{ETA_EDGES[j]:g}\!-\!{ETA_EDGES[j + 1]:g}" for j in range(4)]
    vg = value_grid(TNP["id"]["sf"], rws, cls, fmt="{:.3f}", cell=(0.82, 0.29), label_h=0.145, vmin=0.93, vmax=1.0)
    head = tex_h(r"p_T \,\backslash\, |\eta|", 0.145)
    head.next_to(vg.row_labels, UP, buff=0.10).align_to(vg.col_labels, DOWN)
    g = VGroup(vg, head)
    g.move_to(np.zeros(3))
    g.vg, g.head, g.n_col = vg, head, 4
    return place(g, placement, np.zeros(3))


def map_cell(m, i, j) -> VGroup:
    """Cell (i, j) of a ``sf_map``: its rectangle and its number (texts are row-major, none blank)."""
    k = i * m.n_col + j
    return VGroup(m.vg.cells[k], m.vg.texts[k])


def fan_state(j, placement=FAN_T) -> VGroup:
    """The side view with |eta| bin ``j`` lit (``None``: none)."""
    fan = eta_fan(np.zeros(3), eta_edges=tuple(ETA_EDGES))
    if j is not None:
        fan.wedges[j].set_fill(col(DETECTOR_ACCENT), opacity=0.75)
    return place(fan, placement, np.zeros(3))


def eta_label(j: int):
    lo, hi = ETA_EDGES[j], ETA_EDGES[j + 1]
    expr = rf"|\eta| < {hi:g}" if lo == 0 else rf"{lo:g} < |\eta| < {hi:g}"
    return tex_h(expr, 0.19).move_to(_p3(ETA_LAB_C))


def apply_lines(ev) -> VGroup:
    """One simulated event: its two muons (p_T, |eta|) and the weight the table gives it."""
    mu1, mu2 = ev["muons"]
    rows = [tex_h(rf"p_T = {m['pt']:.0f}\ \mathrm{{GeV}},\ \ |\eta| = {abs(m['eta']):.1f}", 0.17) for m in (mu1, mu2)]
    w = tex_h(rf"w = {mu1['sf_id']:.3f} \times {mu2['sf_id']:.3f} = {ev['w_id']:.3f}", 0.18)
    g = VGroup(*rows, w)
    for k, r in enumerate(g):
        r.shift(np.array([LEFT_COL_X, 2.32 - 0.36 * k - (0.08 if k == 2 else 0.0), 0.0]) - r.get_left())
    assert g.get_right()[0] <= ROW_RIGHT + 1e-6, f"apply lines reach x = {g.get_right()[0]:.2f}, past {ROW_RIGHT}"
    g.mu, g.w = rows, w
    return g


# ---------------------------------------------------------------------------
# the fit and the cross sections
# ---------------------------------------------------------------------------

NP_SHOWN = ("SigModel", "MuonScale", "MuonRes", "PDF", "QCDScale", "Lumi")
NP = {n["name"]: n for n in FIT["nps"]}
PULLS_AT = (-4.35, -1.3)
SLIDER_AT = (-4.1, 1.35)


def pulls(post: bool) -> VGroup:
    names = list(NP_SHOWN)
    p = [NP[n]["pull"] for n in names] if post else [0.0] * len(names)
    c = [NP[n]["constraint"] for n in names] if post else [1.0] * len(names)
    pp = pull_plot(names, p, c, x_length=2.4, row_h=0.4, label_h=0.17)
    pp.move_to(_p3(PULLS_AT))
    return pp


def mu_strings():
    poi = FIT["poi"]
    up, dn = f"{poi['err_up']:.3f}", f"{poi['err_down']:.3f}"
    err = rf"\pm {up}" if up == dn else rf"^{{+{up}}}_{{-{dn}}}"
    return f"{poi['value']:.3f}", err


def mu_slider(post: bool):
    poi = FIT["poi"]
    sl = slider(0.95, 1.03, poi["value"] if post else 1.0, err=max(poi["err_up"], poi["err_down"]) if post else None,
                ref=1.0, length=2.8, ticks=[0.96, 1.00, 1.02], fmt="{:.2f}", tick_label_h=0.16)
    sl.move_to(_p3(SLIDER_AT))
    return sl


def mu_line(h=0.30, at=(3.1, 2.15)):
    v, e = mu_strings()
    return tex_h(rf"\mu_{{Z}} = {v} {e}", h).move_to(_p3(at))


def sigma_formula():
    """How a cross section is measured: events = cross section x luminosity x the fraction we keep."""
    return tex_h(r"\sigma_{\mathrm{fid}} = \frac{N - B}{C \cdot L}", 0.26).move_to(np.array([3.1, 1.2, 0.0]))


def sigma_counting():
    c = FIT["counting"]
    expr = (rf"\frac{{{tex_int(c['n_obs'])} - {tex_int(c['n_bkg'])}}}{{{FIT['C']:.3f} \times {FIT['lumi_pb'] / 1000:.2f}\ \mathrm{{fb}}^{{-1}}}}"
            rf" = {c['sigma_fid_pb']:.1f}\ \mathrm{{pb}}")
    line = tex_h(expr, 0.19).move_to(np.array([3.1, 0.2, 0.0]))
    tag = text_h("counting", 0.14, color=GREY).next_to(line, UP, buff=0.12).align_to(line, LEFT)
    g = VGroup(line, tag)
    g.line, g.tag = line, tag
    return g


def sigma_fid_line(y=-0.85):
    s = FIT["sigma_fid"]
    expr = (rf"\sigma_{{\mathrm{{fid}}}} = {s['value']:.1f} \pm {s['stat']:.1f}_{{\mathrm{{stat}}}}"
            rf" \pm {s['syst']:.1f}_{{\mathrm{{syst}}}} \pm {s['lumi']:.1f}_{{\mathrm{{lumi}}}}\ \mathrm{{pb}}")
    line = tex_h(expr, 0.25).move_to(np.array([3.1, y, 0.0]))
    tag = text_h("fit", 0.14, color=GREY).next_to(line, UP, buff=0.12).align_to(line, LEFT)
    g = VGroup(line, tag)
    g.line, g.tag = line, tag
    return g


# ten schematic pairs for the acceptance picture (eta of the two muons, soft?): four of ten are kept, A = 0.409
ACC_PAIRS = ((0.3, -0.6, False), (1.1, -0.2, False), (-1.6, 0.7, False), (1.9, 0.4, False),
             (2.9, 1.3, False), (-3.3, -0.8, False), (3.6, 2.7, False), (-2.8, 0.2, False), (2.6, -1.0, False),
             (0.6, -0.9, True))
ACC_FAN = (0.5, (1.35, 0.2))                       # x 0.17 .. 2.53, y -0.48 .. 0.88: left of A and of the result line


def acceptance_picture() -> VGroup:
    """Side view: muons inside |eta| < 2.4 with enough p_T are measured (gold); a pair with a muon down the
    beam pipe, or too soft, is lost (grey). Schematic directions; A is the frozen number."""
    fan = eta_fan(np.zeros(3), eta_edges=tuple(ETA_EDGES))
    kept, lost = VGroup(), VGroup()
    up = +1
    for e1, e2, soft in ACC_PAIRS:
        ok = abs(e1) < 2.4 and abs(e2) < 2.4 and not soft
        for k, e in enumerate((e1, e2)):
            side = +1 if e >= 0 else -1
            u = up if k == 0 else -up
            if abs(e) < 2.4:
                end = fan.edge_point(e, side=side, up=u)
                if soft:
                    end = fan.c + 0.35 * (end - fan.c)
            else:                                        # inside the empty cone, out through the endcap hole and beyond
                th = eta_theta(abs(e))
                end = fan.c + np.array([side * (fan.half_z + 0.3), u * (fan.half_z + 0.3) * math.tan(th), 0.0])
            ln = Line(fan.c, end, stroke_color=col(CHANNEL_LINE["mumu"] if ok else GREY), stroke_width=3.0 if ok else 2.0)
            (kept if ok else lost).add(ln)
        up = -up
    a = FIT["sigma_60_120"]["A"]
    g = VGroup(fan, lost, kept)
    place(g, ACC_FAN, np.zeros(3))
    a_tex = tex_h(rf"A = {a:.3f}", 0.26).move_to(np.array([4.75, 0.62, 0.0]))
    g.add(a_tex)
    g.fan, g.lost, g.kept, g.a_tex = fan, lost, kept, a_tex
    assert len(kept) == 8 and len(lost) == 12            # 4 of 10 pairs
    return g


SIG_AX = dict(lo=1850.0, hi=2050.0, x0=0.25, length=5.4, y=-2.3)


def sigma_panel() -> VGroup:
    """sigma(60-120) = sigma_fid / A on its axis: the measurement (dot + total error bar, labelled "this
    analysis") beside the dashed prediction it is normalised to (aMC@NLO lineshape, NNLO
    sigma(m > 50) = 6077.22 pb; z-mumu/handoff.md:100), labelled with what it is."""
    s = FIT["sigma_60_120"]
    lo, hi, x0, L, y0 = SIG_AX["lo"], SIG_AX["hi"], SIG_AX["x0"], SIG_AX["length"], SIG_AX["y"]
    ink = col(INK)

    def xs(v):
        return x0 + (float(v) - lo) / (hi - lo) * L

    axis = Line([xs(lo), y0, 0], [xs(hi), y0, 0], stroke_color=ink, stroke_width=2.5)
    ticks, tick_labels = VGroup(), VGroup()
    for v in np.arange(lo, hi + 1, 50.0):
        p = np.array([xs(v), y0, 0.0])
        ticks.add(Line(p, p + DOWN * 0.08, stroke_color=ink, stroke_width=2.5))
        tick_labels.add(text_h(f"{v:.0f}", 0.15).next_to(p + DOWN * 0.08, DOWN, buff=0.08))
    unit_ = tex_h(r"\mathrm{pb}", 0.17).next_to(tick_labels[-1], RIGHT, buff=0.18)
    pred = float(s["pred"])
    theory = DashedLine([xs(pred), y0, 0], [xs(pred), y0 + 0.8, 0], dash_length=0.1,
                        stroke_color=col(THEORY), stroke_width=3.0)
    th_lab = tex_h(rf"{pred:.1f}\ \mathrm{{pb}}", 0.17, color=THEORY).next_to(theory.get_end(), UP, buff=0.08)
    th_name = VGroup(text_h("prediction", 0.15, color=THEORY),
                     text_h("aMC@NLO, NNLO norm.", 0.12, color=THEORY)).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
    th_name.next_to(th_lab, RIGHT, buff=0.22).align_to(th_lab, UP)
    ypt = y0 + 0.4
    cc = col(CHANNEL_LINE["mumu"])
    bar = Line([xs(s["value"] - s["total"]), ypt, 0], [xs(s["value"] + s["total"]), ypt, 0],
               stroke_color=cc, stroke_width=5.0)
    dot = Dot([xs(s["value"]), ypt, 0], radius=0.1, color=cc)
    me_name = text_h("this analysis", 0.15).next_to(bar, LEFT, buff=0.22)
    label = tex_h(rf"\sigma_{{60\text{{--}}120}} = \frac{{\sigma_{{\mathrm{{fid}}}}}}{{A}} = {s['value']:.0f} \pm {s['total']:.0f}\ \mathrm{{pb}}", 0.2)
    label.move_to(np.array([4.85, -0.12, 0.0]))        # right of the side view, under A
    g = VGroup(axis, ticks, tick_labels, unit_, theory, th_lab, th_name, bar, dot, me_name, label)
    g.axis, g.theory, g.th_lab, g.th_name, g.bar, g.dot, g.me_name, g.label = (
        axis, theory, th_lab, th_name, bar, dot, me_name, label)
    g.frame = VGroup(axis, ticks, tick_labels, unit_)
    return g


# ---------------------------------------------------------------------------
# end states (pure) and their z-orders
# ---------------------------------------------------------------------------

def detector_end_state() -> dict:
    """Exact replica of 4-02's last frame (s4_zmumu.ZmumuIntoDetector: CMSSlice +
    channel_common.show_signature for mu- at 35 deg and mu+ at 212 deg)."""
    det = CMSSlice()
    st = {"det": det}
    for sfx, phi, q, tex in (("m", PHI_M, -1, r"\mu^-"), ("p", PHI_P, +1, r"\mu^+")):
        sig = signature(det, "mu", phi, charge=q, kappa=KAPPA_402)
        trk, rest = sig[0], VGroup(*sig[1:])
        th = tracker_hits(det, trk, color=CHANNEL["mumu"])
        lab = mathtex(tex, color=CHANNEL_LINE["mumu"]).scale(1.1).move_to(
            det.point_at(det.outer_radius + 0.42, phi))
        st.update({f"trk_{sfx}": trk, f"hits_{sfx}": th, f"dep_{sfx}": rest, f"lab_{sfx}": lab})
    return st


ORDER_402 = ("det", "trk_m", "hits_m", "dep_m", "lab_m", "trk_p", "hits_p", "dep_p", "lab_p")


def state_a() -> dict:
    """First entry: the parked slice with the real event, the linear count axes, one block."""
    rec = FIRST_EVENTS[0]
    P = plot_parts(stage=None, data="blocks", entries=((rec["mass"], 0),), counter=False, ratio=False, key=False,
                   x_labels=True, linear=True)
    return {"det": parked_slice(), "ev": place(event_group(rec), DET_PARK, DET_CENTER),
            "stamp": run_stamp(rec), "dax": P["dax"], "stack": P["stack"], "data": P["data"]}


ORDER_A = ("det", "ev", "stamp", "dax", "stack", "data")


def state_b() -> dict:
    P = plot_parts(stage=None, data="bars", ratio=False, key=False, x_labels=True)
    return {"det": parked_slice(), **P}


ORDER_B = ("det", "dax", "stack", "data", "counter")


def state_c() -> dict:
    return {"det": parked_slice(), **plot_parts(stage="raw")}


ORDER_C = ("det", *PLOT_KEYS)


def state_t1() -> dict:
    pair = tp_pair(TP_EVENTS["pass"])
    return {**plot_parts(PLOT_GHOST, stage="raw"), "det": tp_slice(),
            "tp_tag": pair.tag, "tp_tag_lab": tp_label(pair, "tag", CHANNEL_LINE["mumu"]),
            "tp_probe": pair.probe, "tp_probe_lab": tp_label(pair, "probe", CHANNEL_LINE["mumu"]),
            "mass": mass_readout(TP_EVENTS["pass"]["mass_bare"])}


ORDER_T1 = (*PLOT_KEYS, "det", "tp_tag", "tp_tag_lab", "tp_probe", "tp_probe_lab", "mass")


def state_t2() -> dict:
    P = grid_panels(CELL, TP_BIG, samples=("data",), nc=ROW_D_C)
    return {**plot_parts(PLOT_GHOST, stage="raw"), "det": tp_slice(DET_TP_PARK), **P,
            **col_words(TP_BIG, ROW_D_C, row="data"), "cell_lab": cell_label(CELL, "big"),
            "pass_n": tp_count(P[pkey("data", "pass")], "pass"), "fail_n": tp_count(P[pkey("data", "fail")], "fail"),
            "eps_long": eps_long()}


ORDER_T2 = (*PLOT_KEYS, "det", pkey("data", "pass"), pkey("data", "fail"), "w_pass", "w_fail", "cell_lab",
            "pass_n", "fail_n", "eps_long")

GRID_KEYS = tuple(pkey(s, w) for s, w in PANELS)


def state_t3() -> dict:
    return {"det": tp_slice(DET_TP_PARK), **grid_panels(CELL), **col_words(), **row_words(),
            "cell_lab": cell_label(CELL, "grid"), "eps_data": eps_in_grid("data"), "eps_mc": eps_in_grid("mc")}


ORDER_T3 = ("det", *GRID_KEYS, "w_pass", "w_fail", "rw_data", "rw_mc", "cell_lab", "eps_data", "eps_mc")


def scan_common(j=0, x_title=True) -> dict:
    ax = eff_axes(x_title=x_title)
    return {"eff_ax": ax, "eff_key": eff_key(ax), "eff_data": eff_points(ax, "data", j), "eff_mc": eff_points(ax, "mc", j)}


def state_t4() -> dict:
    return {"det": tp_slice(DET_TP_PARK), "scan_trk": scan_track(4), **grid_panels(CELL, GRID_SCAN),
            **col_words(GRID_SCAN), **row_words(GRID_SCAN), "cell_lab": cell_label(CELL, "scan"),
            "eps_data": eps_readout("data", CELL["data"]["eps"], 0.24, EPS_SCAN["data"]),
            "eps_mc": eps_readout("mc", CELL["mc"]["eps"], 0.24, EPS_SCAN["mc"]), **scan_common()}


ORDER_T4 = ("det", "scan_trk", *GRID_KEYS, "w_pass", "w_fail", "rw_data", "rw_mc", "cell_lab", "eps_data", "eps_mc",
            "eff_ax", "eff_key", "eff_data", "eff_mc")


def state_t5() -> dict:
    sax = sf_axes()
    return {"det": tp_slice(DET_TP_PARK), "scan_trk": scan_track(4), **grid_panels(CELL, GRID_SCAN),
            **col_words(GRID_SCAN), **row_words(GRID_SCAN), "cell_lab": cell_label(CELL, "scan"),
            "sf_tex": sf_tex(), **scan_common(x_title=False), "sf_ax": sax, "sf_pts": sf_points(sax, 0)}


ORDER_T5 = ("det", "scan_trk", *GRID_KEYS, "w_pass", "w_fail", "rw_data", "rw_mc", "cell_lab", "sf_tex",
            "eff_ax", "eff_key", "eff_data", "eff_mc", "sf_ax", "sf_pts")


def state_t6() -> dict:
    sax = sf_axes()
    return {"fan": fan_state(3), "eta_lab": eta_label(3), "sfmap": sf_map(), **scan_common(j=3, x_title=False),
            "sf_ax": sax, "sf_pts": sf_points(sax, 3)}


ORDER_T6 = ("fan", "eta_lab", "sfmap", "eff_ax", "eff_key", "eff_data", "eff_mc", "sf_ax", "sf_pts")


def state_t7() -> dict:
    return {**plot_parts(stage="sf_muon"), "row_mu": row_mu()}


ORDER_T7 = (*PLOT_KEYS, "row_mu")


def state_g1() -> dict:
    return {**plot_parts(stage="sf_pileup"), "row_mu": row_mu(), "row_pu": row_pu()}


ORDER_G1 = (*PLOT_KEYS, "row_mu", "row_pu")


def state_g2() -> dict:
    return {**plot_parts(stage="final"), "row_mu": row_mu(), "row_pu": row_pu(), "row_l1": row_l1()}


ORDER_G2 = (*PLOT_KEYS, "row_mu", "row_pu", "row_l1")


def state_g3() -> dict:
    return {**plot_parts(stage="final"), "row_all": row_all()}


ORDER_G3 = (*PLOT_KEYS, "row_all")

FIT_KEYS = ("dax", "stack", "data", "ratio", "key")


def left_plot() -> dict:
    P = plot_parts(PLOT_LEFT, bins=12, fit="postfit", counter=False, rlabel=False)
    return {k: P[k] for k in FIT_KEYS}


def state_h3() -> dict:
    return {**left_plot(), "mu_line": mu_line(), "f_sym": sigma_formula(), "f_num": sigma_counting(),
            "sig_fid": sigma_fid_line()}


ORDER_H3 = (*FIT_KEYS, "mu_line", "f_sym", "f_num", "sig_fid")


def state_h4() -> dict:
    return {**left_plot(), "mu_line": mu_line(), "sig_fid": sigma_fid_line(y=1.2), "acc": acceptance_picture(),
            "sig_tot": sigma_panel()}


ORDER_H4 = (*FIT_KEYS, "mu_line", "sig_fid", "acc", "sig_tot")

ALL_STATES = {"a": state_a, "b": state_b, "c": state_c, "t1": state_t1, "t2": state_t2, "t3": state_t3, "t4": state_t4,
              "t5": state_t5, "t6": state_t6, "t7": state_t7, "g1": state_g1, "g2": state_g2, "g3": state_g3,
              "h3": state_h3, "h4": state_h4}


# ---------------------------------------------------------------------------
# the clips
# ---------------------------------------------------------------------------

class MumuEvent(Scene):
    """4-02's schematic tracks fade, the real event is drawn, its p_T values converge to
    m_mumu (one clip); the slice parks, linear count axes appear, first entry: a block."""

    def construct(self):
        white_background(self)
        prev = detector_end_state()
        add_state(self, prev, ORDER_402)
        clip_open(self, "mumu_event_to_mass")
        end = state_a()
        rec = FIRST_EVENTS[0]

        self.play(*[FadeOut(prev[k]) for k in ("trk_m", "hits_m", "dep_m", "trk_p", "hits_p", "dep_p")],
                  run_time=0.5)
        grp = event_group(rec)                      # on a natural slice == the live one
        labs = {-1: prev["lab_m"], +1: prev["lab_p"]}
        pts = VGroup()
        for mu, m in zip(grp.ev, grp.ev.muons):
            lab = labs[_q(m)]
            spot = label_spot(phi_end(mu.trk))
            self.play(Create(mu.trk, lag_ratio=0.0), lab.animate.move_to(spot), run_time=0.9, rate_func=EASE)
            self.play(FadeIn(mu.hits), run_time=0.3)
            self.play(FadeIn(mu.deposits), run_time=0.45)
            pt = tex_h(rf"p_{{T}} = {float(m['pt']):.1f}\ \mathrm{{GeV}}", 0.22)
            ghost = lab.copy().move_to(spot)
            pt.next_to(ghost, RIGHT if spot[0] >= 0 else LEFT, buff=0.2)
            pts.add(pt)
        regroup(self, grp, "event 1")
        self.play(*[FadeIn(p, shift=0.25 * (RIGHT if p.get_center()[0] > 0 else LEFT)) for p in pts],
                  run_time=0.5)
        self.wait(0.4)
        m_tex = mass_tex(rec).move_to(np.array([4.95, 0.2, 0.0]))
        self.remove(*pts)
        self.add(pts)
        self.play(ReplacementTransform(pts, m_tex), run_time=0.9, rate_func=EASE)
        clip_cut(self, "mumu_first_entry")

        m_spot = end["dax"].c2p(float(rec["mass"]), 3.45)           # MAIN placement = natural
        self.play(ReplacementTransform(prev["det"], end["det"]), ReplacementTransform(grp, end["ev"]),
                  FadeOut(prev["lab_m"]), FadeOut(prev["lab_p"]), m_tex.animate.move_to(m_spot),
                  FadeIn(end["stamp"]), FadeIn(end["dax"]), run_time=1.8, rate_func=EASE)
        self.add(end["stack"])                      # 8 degenerate layers at the floor (invisible)
        self.bring_to_front(m_tex)
        self.wait(0.2)
        block = end["data"][0]
        self.play(ReplacementTransform(m_tex, block), run_time=0.8, rate_func=EASE)
        regroup(self, end["data"], "first block")
        finish(self, end, ORDER_A)
        self.wait(0.2)


class MumuFill(Scene):
    """Four more real events drop blocks (two share a bin); without a pause the count axis
    becomes logarithmic (five events are nothing at that scale), the clock runs and blocks
    rain over the whole plot while the 60 bins grow to the real counts (N = 10,378,567)."""

    def construct(self):
        white_background(self)
        prev = state_a()
        add_state(self, prev, ORDER_A)
        clip_open(self, "mumu_fill")
        end = state_b()
        det, dax, stamp, entries = prev["det"], prev["dax"], prev["stamp"], prev["data"]

        self.play(FadeOut(prev["ev"]), run_time=0.35)
        for k, rec in enumerate(FIRST_EVENTS[1:], start=1):
            f = 1.0 if k == 1 else 0.7              # the first one at reading pace, then quicker
            g = place(event_group(rec), DET_PARK, DET_CENTER)
            lines = VGroup(*[mu.trk for mu in g.ev], *[ph[0] for ph in g.photons])
            rest = VGroup(*[VGroup(mu.hits, mu.deposits) for mu in g.ev], *[ph[1] for ph in g.photons])
            self.play(Create(lines, lag_ratio=0.0), Transform(stamp, run_stamp(rec)), run_time=0.7 * f,
                      rate_func=EASE)
            self.play(FadeIn(rest), run_time=0.3 * f)
            blk = entry_block(dax, rec["mass"], FIRST_LEVELS[k])
            fly = blk.copy().scale(0.6).move_to(det.get_center())
            self.play(Transform(fly, blk), FadeOut(lines), FadeOut(rest), run_time=0.65 * f, rate_func=EASE)
            self.remove(fly)
            entries.add(blk)
        self.wait(0.15)

        # -- no pause: the axis becomes logarithmic; at that scale five events are nothing
        dax_log = end["dax"]
        clk = clock(CLOCK_AT)
        counter = make_counter(len(entries))
        gone = [b.copy().stretch_to_fit_height(1e-3).align_to(dax.frame, DOWN).set_opacity(0.0) for b in entries]
        self.play(*[Transform(b, t) for b, t in zip(entries, gone)],
                  FadeOut(dax.y_labels), FadeOut(dax.y_ticks_v), FadeOut(dax.y_title),
                  FadeIn(dax_log.y_labels), FadeIn(dax_log.y_ticks_v), FadeIn(dax_log.y_title),
                  FadeIn(clk), FadeIn(counter), FadeOut(stamp), run_time=0.7, rate_func=EASE)
        self.remove(entries)
        swap_family(self, dax, dax_log, "count axis: linear -> log")

        bars = end["data"].copy()
        bars0 = data_bars(dax_log, visible=False)
        self.add(bars0)
        n = ValueTracker(float(len(entries)))
        counter.num.add_updater(lambda m: (m.set_value(float(round(n.get_value()))),
                                           layout_counter(counter.pre, m)))
        order = np.random.default_rng(20260915).permutation(60)
        height = np.clip(np.log10(np.maximum(DATA, 1.0)) - Y_EXP[0], 0.05, None)    # what the eye sees on the log axis
        drops = hidden_rain(det, dax_log, EDGES, DATA, n=120, seed=3, color=SAMPLE["Data"], p=height, square=True,
                            y_frac=(0.25, 1.0))
        self.play(LaggedStart(*[Transform(bars0[i], bars[i]) for i in order], lag_ratio=0.02, group=bars0),
                  LaggedStart(*drops, lag_ratio=0.03),
                  clk.turns.animate(rate_func=rate_functions.ease_in_quad).set_value(3.0),
                  n.animate(rate_func=rate_functions.ease_in_quad).set_value(float(TOTAL)), run_time=4.2)
        counter.num.clear_updaters()
        set_counter(counter, TOTAL)
        self.remove(n, clk.turns)
        self.play(FadeOut(clk), run_time=0.35)
        settle(self, bars0, end["data"], "bars")
        settle(self, counter, end["counter"], "counter")
        settle_same(self, prev, end, ("det", "stack"))
        finish(self, end, ORDER_B)
        self.wait(0.2)


class MumuStack(Scene):
    """Bars become points; the uncorrected prediction (simulation + the data-driven fakes
    wedge) slides in under them bottom-up; the ratio panel opens at data/pred = 0.944."""

    def construct(self):
        white_background(self)
        prev = state_b()
        add_state(self, prev, ORDER_B)
        clip_open(self, "mumu_prediction")
        end = state_c()

        self.play(ReplacementTransform(prev["data"], end["data"]), run_time=1.0, rate_func=EASE)
        self.wait(0.2)
        stack = prev["stack"]
        dax = prev["dax"]
        base = sr_layers("raw")
        for i, poly in enumerate(stack):            # each (invisible) layer waits at its own base
            rows_ = [(n, c if j < i else np.zeros(60), cc) for j, (n, c, cc) in enumerate(base)]
            poly.become(stack_hist(dax, EDGES, rows_)[i])
        self.play(LaggedStart(*[Transform(stack[i], end["stack"][i]) for i in range(len(stack))],
                              lag_ratio=0.15, group=stack), run_time=2.6, rate_func=EASE)
        settle(self, stack, end["stack"], "stack")
        self.play(FadeIn(end["key"]), run_time=0.6)
        clip_cut(self, "mumu_ratio")
        ratio = end["ratio"]
        self.play(FadeIn(ratio.dax), FadeIn(ratio.ref), dax.x_labels.animate.set_opacity(0.0),
                  dax.x_title.animate.set_opacity(0.0), run_time=0.6)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in ratio.dots], lag_ratio=0.03, group=ratio.dots),
                  FadeIn(ratio.errs), run_time=1.3)
        regroup(self, ratio, "ratio")
        self.wait(0.3)
        self.play(FadeIn(end["rlabel"], shift=LEFT * 0.2), run_time=0.5)
        settle_same(self, prev, end, ("det", "dax", "counter"))
        finish(self, end, ORDER_C)
        self.wait(0.2)


class MumuTagProbe(Scene):
    """Tag and probe with real data, the idea: the tag, the probe, the Z mass as the label; why a probe
    passes or fails tight ID; the pass / fail spectra of one cell and eps_data; the same on simulation."""

    def construct(self):
        white_background(self)
        prev = state_c()
        add_state(self, prev, ORDER_C)
        clip_open(self, "mumu_tnp_tag_probe")

        # -- the plot parks; a real pair: the tag, then the probe, then the mass ---------------
        e1 = state_t1()
        self.play(*[ReplacementTransform(prev[k], e1[k]) for k in PLOT_KEYS], ReplacementTransform(prev["det"], e1["det"]),
                  run_time=1.6, rate_func=EASE)
        self.wait(0.2)
        pair_rec = TP_EVENTS["pass"]
        tag, probe = e1["tp_tag"], e1["tp_probe"]
        self.play(Create(tag.trk, lag_ratio=0.0), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(tag.hits), run_time=0.3)
        self.play(FadeIn(tag.deposits), run_time=0.45)
        regroup(self, tag, "tag")
        ring = trigger_ring(tag)
        self.play(GrowFromCenter(ring), run_time=0.35)
        self.play(ring.animate.scale(2.2).set_stroke(opacity=0.0), run_time=0.5)
        self.remove(ring)
        self.play(FadeIn(e1["tp_tag_lab"]), run_time=0.4)
        self.wait(0.3)
        pair = VGroup(tag, probe)
        pair.probe = probe
        asked = probe_asked(pair)
        grey_lab = e1["tp_probe_lab"].copy().set_color(col(GREY))
        self.play(Create(asked), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(grey_lab), run_time=0.4)
        mass = mass_readout(60.0)
        self.play(FadeIn(mass), run_time=0.4)
        self.play(ChangeDecimalToValue(mass.num, float(pair_rec["mass_bare"])), run_time=1.6, rate_func=EASE)
        self.wait(0.2)
        self.play(Create(probe.trk, lag_ratio=0.0), FadeOut(asked), Transform(grey_lab, e1["tp_probe_lab"]),
                  run_time=0.9, rate_func=EASE)
        self.play(FadeIn(probe.hits), FadeIn(probe.deposits), run_time=0.5)
        regroup(self, probe, "probe")
        settle(self, grey_lab, e1["tp_probe_lab"], "probe label")
        settle(self, mass, e1["mass"], "mass readout")
        finish(self, e1, ORDER_T1)
        clip_cut(self, "mumu_tnp_pass_fail")

        # -- the question: tight ID. Why this probe passes, why another one fails ----------------
        e2 = state_t2()
        settle_same(self, e1, e2, PLOT_KEYS)
        det = e1["det"]
        p_pass, p_fail = e2[pkey("data", "pass")], e2[pkey("data", "fail")]
        self.play(FadeOut(e1["mass"]), run_time=0.4)
        pulse(self, probe.deposits[1])
        v_pass = verdict(int(pair_rec["probe"]["nStations"]), True)
        self.play(FadeIn(v_pass.words, shift=UP * 0.1), run_time=0.4)
        self.play(GrowFromCenter(v_pass.mark), run_time=0.35)
        self.play(FadeIn(p_pass.ax), FadeIn(p_fail.ax), FadeIn(e2["w_pass"]), FadeIn(e2["w_fail"]),
                  FadeIn(e2["cell_lab"]), run_time=0.6)

        def fly(trk, panel, mass_):
            src = trk[-1].get_end()                  # the track core (VGroup(rim, core))
            d = Rectangle(width=0.13, height=0.13, stroke_width=0, fill_color=col(SAMPLE["Data"]), fill_opacity=1.0).move_to(src)
            self.play(GrowFromCenter(d), run_time=0.2)
            self.play(d.animate.move_to(panel.ax.c2p(float(mass_), 0.0)), run_time=0.7, rate_func=rate_functions.ease_in_quad)
            self.play(FadeOut(d), run_time=0.2)

        fly(probe.trk, p_pass, pair_rec["mass_bare"])
        self.play(*[FadeOut(e1[k]) for k in ("tp_tag", "tp_tag_lab", "tp_probe", "tp_probe_lab")], FadeOut(v_pass),
                  run_time=0.4)
        fail_rec = TP_EVENTS["fail"]
        fail = tp_pair(fail_rec)
        flabs = VGroup(tp_label(fail, "tag", CHANNEL_LINE["mumu"]), tp_label(fail, "probe", CHANNEL_LINE["mumu"]))
        lines = VGroup(fail.tag.trk, fail.probe.trk)
        rest = VGroup(fail.tag.hits, fail.tag.deposits, fail.probe.hits, fail.probe.deposits)
        self.play(Create(lines, lag_ratio=0.0), run_time=0.8, rate_func=EASE)
        self.play(FadeIn(rest), FadeIn(flabs), run_time=0.5)
        ring2 = trigger_ring(fail.tag)
        self.play(GrowFromCenter(ring2), run_time=0.3)
        self.play(ring2.animate.scale(2.2).set_stroke(opacity=0.0), run_time=0.4)
        self.remove(ring2)
        pulse(self, fail.probe.deposits[1])
        v_fail = verdict(int(fail_rec["probe"]["nStations"]), False)
        self.play(FadeIn(v_fail.words, shift=UP * 0.1), run_time=0.4)
        self.play(GrowFromCenter(v_fail.mark), run_time=0.35)
        self.wait(0.2)
        fly(fail.probe.trk, p_fail, fail_rec["mass_bare"])
        self.play(FadeOut(lines), FadeOut(rest), FadeOut(flabs), FadeOut(v_fail), run_time=0.4)

        # -- every probe of the cell
        pb0, fb0 = flat_bars(p_pass), flat_bars(p_fail)
        pb, fb = p_pass.body.copy(), p_fail.body.copy()
        self.add(pb0, fb0)
        order = np.random.default_rng(20260916).permutation(60)
        drops = (hidden_rain(det, p_pass.ax, EDGES, CELL["data"]["pass"], n=45, seed=5, color=SAMPLE["Data"], square=True)
                 + hidden_rain(det, p_fail.ax, EDGES, CELL["data"]["fail"], n=25, seed=6, color=SAMPLE["Data"], square=True))
        np.random.default_rng(7).shuffle(drops)
        self.play(LaggedStart(*[Transform(pb0[i], pb[i]) for i in order], lag_ratio=0.02, group=pb0),
                  LaggedStart(*[Transform(fb0[i], fb[i]) for i in order], lag_ratio=0.02, group=fb0),
                  LaggedStart(*drops, lag_ratio=0.05), run_time=3.2)
        settle(self, pb0, p_pass.body, "pass bars")
        settle(self, fb0, p_fail.body, "fail bars")
        self.play(FadeIn(p_pass.band), FadeIn(p_fail.band), run_time=0.5)
        self.play(Create(p_pass.line), Create(p_fail.line), run_time=1.3, rate_func=rate_functions.linear)
        regroup(self, p_pass, "pass panel")
        regroup(self, p_fail, "fail panel")
        self.play(FadeIn(e2["pass_n"], shift=UP * 0.1), FadeIn(e2["fail_n"], shift=UP * 0.1), run_time=0.5)
        self.play(ReplacementTransform(det, e2["det"]), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(e2["eps_long"], shift=UP * 0.15), run_time=0.8)
        finish(self, e2, ORDER_T2)
        clip_cut(self, "mumu_tnp_grid")

        # -- the two panels move up = "data"; the identical procedure on simulation below --------------
        e3 = state_t3()
        settle(self, e2["det"], e3["det"], "det")
        self.play(*[FadeOut(e2[k]) for k in PLOT_KEYS], FadeOut(e2["pass_n"]), FadeOut(e2["fail_n"]), run_time=0.6)
        self.play(*[ReplacementTransform(e2[pkey("data", w)], e3[pkey("data", w)]) for w in ("pass", "fail")],
                  ReplacementTransform(e2["w_pass"], e3["w_pass"]), ReplacementTransform(e2["w_fail"], e3["w_fail"]),
                  ReplacementTransform(e2["cell_lab"], e3["cell_lab"]), FadeOut(e2["eps_long"], shift=UP * 0.2),
                  run_time=1.5, rate_func=EASE)
        self.play(FadeIn(e3["rw_data"]), FadeIn(e3["eps_data"], shift=UP * 0.1), run_time=0.6)
        self.wait(0.3)
        mp, mf = e3[pkey("mc", "pass")], e3[pkey("mc", "fail")]
        self.play(FadeIn(mp.ax), FadeIn(mf.ax), FadeIn(e3["rw_mc"]), run_time=0.6)
        self.play(Create(mp.body), Create(mf.body), run_time=1.3, rate_func=rate_functions.linear)
        self.play(FadeIn(mp.band), FadeIn(mf.band), Create(mp.line), Create(mf.line), run_time=1.0,
                  rate_func=rate_functions.linear)
        regroup(self, mp, "simulation pass panel")
        regroup(self, mf, "simulation fail panel")
        self.play(FadeIn(e3["eps_mc"], shift=UP * 0.1), run_time=0.6)
        finish(self, e3, ORDER_T3)
        self.wait(0.2)


class MumuTagProbeScan(Scene):
    """Tag and probe, the measurement: the p_T scan fills eps_ID(p_T); SF = eps_data / eps_sim below it;
    the SF points become a table, the |eta| scan completes it; the table is applied to the prediction."""

    def construct(self):
        white_background(self)
        prev = state_t3()
        add_state(self, prev, ORDER_T3)
        clip_open(self, "mumu_tnp_scan")

        # -- the grid parks left, the eps plot opens, the p_T window steps through the ten barrel cells ----
        e4 = state_t4()
        moved = [*GRID_KEYS, "w_pass", "w_fail", "rw_data", "rw_mc", "cell_lab", "eps_data", "eps_mc"]
        self.play(*[ReplacementTransform(prev[k], e4[k]) for k in moved], run_time=1.5, rate_func=EASE)
        settle(self, prev["det"], e4["det"], "det")
        ax = e4["eff_ax"]
        trk = scan_track(4)
        self.play(FadeIn(ax), FadeIn(e4["eff_key"]), Create(trk, lag_ratio=0.0), run_time=0.8)
        band = pt_band(ax, 4)
        lab = e4["cell_lab"]
        live = {k: e4[k] for k in GRID_KEYS}
        eps_d, eps_m = e4["eps_data"], e4["eps_mc"]
        self.add(band)
        self.bring_to_back(band)
        landed_d, landed_m = VGroup(), VGroup()

        def show_cell(i, run_time):
            cell = CELLS[i]
            tgt = grid_panels(cell, GRID_SCAN)
            anims = []
            for k in GRID_KEYS:
                anims += [Transform(live[k].body, tgt[k].body), Transform(live[k].fit, tgt[k].fit)]
            new_lab = cell_label(cell, "scan")
            anims += [Transform(band, pt_band(ax, i)), Transform(trk, scan_track(i)),
                      Transform(lab[0], new_lab[0]),
                      ChangeDecimalToValue(eps_d.num, float(cell["data"]["eps"])),
                      ChangeDecimalToValue(eps_m.num, float(cell["mc"]["eps"]))]
            self.play(*anims, run_time=run_time, rate_func=EASE)

        band.set_opacity(0.0)
        for n_step, i in enumerate(range(10)):
            slow = n_step < 2
            if n_step == 0:
                band.become(pt_band(ax, 4).set_opacity(0.0))
                self.play(band.animate.set_fill(opacity=0.22), run_time=0.3)
            show_cell(i, 0.8 if slow else 0.5)
            pd, pm = e4["eff_data"][i], e4["eff_mc"][i]
            self.play(GrowFromCenter(pd), GrowFromCenter(pm), run_time=0.4 if slow else 0.25)
            landed_d.add(pd)
            landed_m.add(pm)
        show_cell(4, 0.7)                                # back to the cell everybody knows
        self.play(FadeOut(band), run_time=0.3)
        regroup(self, e4["eff_data"], "eps data points")
        regroup(self, e4["eff_mc"], "eps simulation points")
        ref = grid_panels(CELL, GRID_SCAN)
        for k in GRID_KEYS:
            ok, why = looks_same(live[k], ref[k])
            assert ok, f"scan: panel {k} did not return to the 40-45 GeV cell: {why}"
        settle(self, trk, e4["scan_trk"], "scan track")
        finish(self, e4, ORDER_T4)
        clip_cut(self, "mumu_tnp_scale_factor")

        # -- SF = eps_data / eps_sim, below the eps plot ---------------------------------------------------
        e5 = state_t5()
        keep = ["det", "scan_trk", *GRID_KEYS, "w_pass", "w_fail", "rw_data", "rw_mc", "cell_lab", "eff_key", "eff_data", "eff_mc"]
        self.play(FadeOut(e4["eps_data"], shift=DOWN * 0.1), FadeOut(e4["eps_mc"], shift=DOWN * 0.1),
                  FadeIn(e5["sf_tex"], shift=UP * 0.1), run_time=0.8)
        sax = e5["sf_ax"]
        self.play(FadeIn(sax), ax.x_title.animate.set_opacity(0.0), ax.x_labels.animate.set_opacity(0.0), run_time=0.7)
        merges, copies = [], []
        for i in range(10):
            a, b = e4["eff_data"][i].copy(), e4["eff_mc"][i].copy()
            self.add(a, b)
            copies += [a, b]
            tgt = e5["sf_pts"][i]
            # the pair lands as one SF point, and that point stays: the real one is revealed on landing (the copies,
            # identical to it by then, lie under it and leave after the play). Fading the copies out on landing made
            # every point vanish until all ten popped back at the end.
            merges.append(Succession(Transform(a, tgt.copy(), run_time=0.55), FadeIn(tgt, run_time=0.02)))
            merges.append(Transform(b, tgt.copy(), run_time=0.55))
        self.play(LaggedStart(*merges, lag_ratio=0.06), run_time=2.0, rate_func=EASE)
        self.remove(*copies)
        regroup(self, e5["sf_pts"], "SF points")
        for k in keep:
            settle(self, e4[k], e5[k], k)
        settle(self, ax, e5["eff_ax"], "eps axes")
        finish(self, e5, ORDER_T5)
        clip_cut(self, "mumu_tnp_table")

        # -- the SF points become the p_T column; the side view steps through |eta|; the table fills ------
        e6 = state_t6()
        gone = ["det", "scan_trk", *GRID_KEYS, "w_pass", "w_fail", "rw_data", "rw_mc", "cell_lab", "sf_tex"]
        fan = fan_state(0)
        eta_lab = eta_label(0)
        m = e6["sfmap"]
        self.play(*[FadeOut(e5[k]) for k in gone], run_time=0.6)
        self.play(FadeIn(fan), FadeIn(eta_lab), FadeIn(m.head), FadeIn(m.vg.row_labels), FadeIn(m.vg.col_labels[0]),
                  run_time=0.8)
        eff_d, eff_m, sfp = e5["eff_data"], e5["eff_mc"], e5["sf_pts"]
        ax5, sax5 = e5["eff_ax"], e5["sf_ax"]

        def fill_column(j, run_time):
            flights = []
            for i in range(10):
                c = map_cell(m, i, j)
                d = sfp[i][1].copy()
                self.add(d)
                flights.append(Succession(d.animate(run_time=0.5).move_to(c.get_center()).set_opacity(0.0),
                                          FadeOut(d, run_time=0.01)))
                flights.append(Succession(FadeIn(c, run_time=0.5)))
            self.play(LaggedStart(*flights, lag_ratio=0.04), run_time=run_time)

        fill_column(0, 1.8)
        for j in (1, 2, 3):
            self.play(Transform(fan, fan_state(j)), Transform(eta_lab, eta_label(j)),
                      Transform(eff_d, eff_points(ax5, "data", j)), Transform(eff_m, eff_points(ax5, "mc", j)),
                      Transform(sfp, sf_points(sax5, j)), FadeIn(m.vg.col_labels[j]), run_time=0.9, rate_func=EASE)
            fill_column(j, 1.2)
        regroup(self, m, "scale-factor table")
        settle(self, fan, e6["fan"], "side view")
        settle(self, eta_lab, e6["eta_lab"], "eta label")
        settle(self, eff_d, e6["eff_data"], "eps data points")
        settle(self, eff_m, e6["eff_mc"], "eps simulation points")
        settle(self, sfp, e6["sf_pts"], "SF points")
        settle_same(self, e5, e6, ("eff_ax", "eff_key", "sf_ax"))
        finish(self, e6, ORDER_T6)
        clip_cut(self, "mumu_tnp_apply")

        # -- the table beside the plot: event by event, then all of the prediction -------------------------
        e7 = state_t7()
        P = plot_parts(stage="raw")
        tab = sf_map(TAB_APPLY)
        assert tab.get_right()[0] <= ROW_RIGHT + 1e-6 and tab.get_top()[1] < 0.2, "the table must clear the plot's y title"
        self.play(*[FadeOut(e6[k]) for k in ("fan", "eta_lab", "eff_ax", "eff_key", "eff_data", "eff_mc", "sf_ax", "sf_pts")],
                  ReplacementTransform(e6["sfmap"], tab), run_time=1.0, rate_func=EASE)
        self.play(*[FadeIn(P[k]) for k in PLOT_KEYS], run_time=0.8)
        top = np.asarray(P["stack"].total, dtype=float)
        for n_ev, ev in enumerate(APPLY):
            f = 1.0 if n_ev == 0 else 0.65
            lines = apply_lines(ev)
            rects = VGroup(*[SurroundingRectangle(tab.vg.cells[mu["ipt"] * 4 + mu["ieta"]], color=col(DETECTOR_ACCENT),
                                                  buff=0.0, stroke_width=3.5) for mu in ev["muons"]])
            for k in (0, 1):                     # left-column items enter from the right: never across the top-left block
                self.play(FadeIn(lines.mu[k], shift=LEFT * 0.15), run_time=0.35 * f)
                self.play(Create(rects[k]), run_time=0.35 * f)
            self.play(FadeIn(lines.w, shift=LEFT * 0.15), run_time=0.45 * f)
            blk = Rectangle(width=0.13, height=0.13, stroke_color=darken(SAMPLE["DYmumu"], 0.4), stroke_width=1.5,
                            fill_color=col(SAMPLE["DYmumu"]), fill_opacity=1.0).move_to(lines.w.get_right() + RIGHT * 0.2)
            k_bin = int(np.clip(math.floor(ev["mass"]) - 60, 0, 59))
            self.play(GrowFromCenter(blk), run_time=0.2 * f)
            self.play(blk.animate.move_to(P["dax"].c2p(float(ev["mass"]), top[k_bin])).scale(ev["w_id"]),
                      run_time=0.75 * f, rate_func=rate_functions.ease_in_quad)
            self.play(FadeOut(blk), FadeOut(lines), FadeOut(rects), run_time=0.3 * f)
        move_ratio(self, P, "sf_id", run_time=1.8)       # every simulated event, both muons: the identification table
        self.wait(0.2)
        subs = sub_rows()
        self.play(ReplacementTransform(tab, subs[0][0]), FadeIn(VGroup(*subs[0][1:]), shift=LEFT * 0.15), run_time=1.0,
                  rate_func=EASE)
        for r in subs[1:]:
            self.play(FadeIn(r, shift=LEFT * 0.2), run_time=0.4, rate_func=EASE)
        move_ratio(self, P, "sf_muon", run_time=1.4)     # isolation, trigger, reconstruction: the same method
        self.wait(0.2)
        self.play(*[FadeOut(VGroup(*r[1:])) for r in subs], *[FadeOut(r[0]) for r in subs[1:]],
                  ReplacementTransform(subs[0][0], e7["row_mu"][0]), FadeIn(VGroup(*e7["row_mu"][1:])), run_time=0.9,
                  rate_func=EASE)
        regroup(self, e7["row_mu"], "muon efficiency row")
        settle_same(self, P, e7, PLOT_KEYS)
        finish(self, e7, ORDER_T7)
        self.wait(0.2)


class MumuCorrections(Scene):
    """The other corrections, in plain words: extra collisions (pile-up) and the trigger that fired one
    crossing early (prefiring) step the prediction to the final one; the data never move."""

    def construct(self):
        white_background(self)
        prev = state_t7()
        add_state(self, prev, ORDER_T7)
        clip_open(self, "mumu_corr_pileup")
        e1 = state_g1()
        # the rows start at LEFT_COL_X, 0.15 right of the top-left block: they enter from the right and leave
        # downward, never from / to the left (a row sliding in from x - 0.25 crosses the block: keepout WARN)
        self.play(FadeIn(e1["row_pu"], shift=LEFT * 0.25), run_time=0.6, rate_func=EASE)
        move_ratio(self, prev, "sf_pileup")
        settle_same(self, prev, e1, (*PLOT_KEYS, "row_mu"))
        finish(self, e1, ORDER_G1)
        clip_cut(self, "mumu_corr_prefiring")

        e2 = state_g2()
        self.play(FadeIn(e2["row_l1"], shift=LEFT * 0.25), run_time=0.6, rate_func=EASE)
        move_ratio(self, e1, "final")
        settle_same(self, e1, e2, (*PLOT_KEYS, "row_mu", "row_pu"))
        finish(self, e2, ORDER_G2)
        clip_cut(self, "mumu_corrected")

        e3 = state_g3()
        self.play(*[FadeOut(e2[k], shift=DOWN * 0.15) for k in ("row_mu", "row_pu", "row_l1")], run_time=0.5)
        self.play(FadeIn(e3["row_all"], shift=LEFT * 0.25), run_time=0.6, rate_func=EASE)
        pulse(self, [e2["rlabel"]], scale=1.5, run_time=0.7)
        settle_same(self, e2, e3, PLOT_KEYS)
        finish(self, e3, ORDER_G3)
        self.wait(0.2)


class MumuFit(Scene):
    """Rebin 60 -> 12 (5 GeV); the pulls and mu_Z move as the fit goes post-fit and the ratio flattens;
    sigma_fid: counting with our numbers, then the fit; sigma(60-120) = sigma_fid / A beside the prediction."""

    def construct(self):
        white_background(self)
        prev = state_g3()
        add_state(self, prev, ORDER_G3)
        clip_open(self, "mumu_rebin")

        self.play(FadeOut(prev["row_all"]), FadeOut(prev["counter"]), FadeOut(prev["rlabel"]), run_time=0.5)
        pre12 = plot_parts(bins=12, fit="prefit", counter=False, rlabel=False)
        anims = [Transform(prev["dax"], pre12["dax"]), Transform(prev["stack"], pre12["stack"])]
        for j in range(12):
            for k in range(5):
                anims.append(Transform(prev["data"][5 * j + k], pre12["data"][j].copy()))
                anims.append(Transform(prev["ratio"].dots[5 * j + k], pre12["ratio"].dots[j].copy()))
                anims.append(Transform(prev["ratio"].errs[5 * j + k], pre12["ratio"].errs[j].copy()))
        self.play(*anims, run_time=1.7, rate_func=EASE)
        for k in FIT_KEYS:
            settle(self, prev[k], pre12[k], k)
        clip_cut(self, "mumu_fit")

        pp0, pp1 = pulls(False), pulls(True)
        sl = mu_slider(False)
        self.play(FadeIn(pp0), FadeIn(sl), run_time=0.7)
        post12 = plot_parts(bins=12, fit="postfit", counter=False, rlabel=False)
        poi = FIT["poi"]
        self.play(Transform(pre12["stack"], post12["stack"]), Transform(pre12["ratio"].dots, post12["ratio"].dots),
                  Transform(pre12["ratio"].errs, post12["ratio"].errs),
                  Transform(pp0.rows, pp1.rows),
                  Transform(sl.marker, sl.marker_at(poi["value"], max(poi["err_up"], poi["err_down"]))),
                  run_time=2.6, rate_func=EASE)
        for k in FIT_KEYS:
            settle(self, pre12[k], post12[k], k)
        v, e = mu_strings()
        seed = tex_h(rf"\mu_{{Z}} = {v} {e}", 0.18).next_to(sl.marker, UP, buff=0.62)
        self.play(FadeIn(seed, shift=UP * 0.15), run_time=0.5)
        clip_cut(self, "mumu_sigma_fid")

        # -- how: events = cross section x luminosity x the fraction we keep; counting, then the fit --------
        e3 = state_h3()
        self.play(*[ReplacementTransform(post12[k], e3[k]) for k in FIT_KEYS], FadeOut(pp0), FadeOut(sl),
                  ReplacementTransform(seed, e3["mu_line"]), run_time=1.8, rate_func=EASE)
        self.play(FadeIn(e3["f_sym"], shift=UP * 0.1), run_time=0.7)
        self.wait(0.3)
        self.play(FadeIn(e3["f_num"].line, shift=UP * 0.1), run_time=0.8)
        self.play(FadeIn(e3["f_num"].tag), run_time=0.3)
        regroup(self, e3["f_num"], "counting line")
        self.wait(0.4)
        self.play(FadeIn(e3["sig_fid"].line, shift=UP * 0.1), run_time=0.8)
        self.play(FadeIn(e3["sig_fid"].tag), run_time=0.3)
        regroup(self, e3["sig_fid"], "fit line")
        finish(self, e3, ORDER_H3)
        clip_cut(self, "mumu_sigma_total")

        # -- the acceptance: where the factor 2.44 comes from; then the result beside the prediction ---------
        e4 = state_h4()
        settle_same(self, e3, e4, (*FIT_KEYS, "mu_line"))
        self.play(FadeOut(e3["f_sym"]), FadeOut(e3["f_num"]), ReplacementTransform(e3["sig_fid"], e4["sig_fid"]),
                  run_time=1.0, rate_func=EASE)
        acc = e4["acc"]
        self.play(FadeIn(acc.fan), run_time=0.6)
        self.play(LaggedStart(*[Create(l) for l in acc.lost], lag_ratio=0.08), run_time=1.1)
        self.play(LaggedStart(*[Create(l) for l in acc.kept], lag_ratio=0.08), run_time=1.0)
        self.play(FadeIn(acc.a_tex, shift=LEFT * 0.15), run_time=0.6)
        regroup(self, acc, "acceptance picture")
        sig = e4["sig_tot"]
        self.play(FadeIn(sig.label, shift=UP * 0.1), run_time=0.8)
        self.wait(0.3)
        self.play(FadeIn(sig.frame), run_time=0.5)
        self.play(Create(sig.theory), FadeIn(sig.th_lab), run_time=0.6)
        self.play(FadeIn(sig.th_name, shift=LEFT * 0.15), run_time=0.5)
        self.play(GrowFromCenter(sig.bar), GrowFromCenter(sig.dot), run_time=0.6)
        self.play(FadeIn(sig.me_name, shift=RIGHT * 0.15), run_time=0.5)
        regroup(self, sig, "sigma panel")
        finish(self, e4, ORDER_H4)
        self.wait(0.2)
