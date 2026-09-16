"""Section 5, real data: the Z -> tautau story, chained on 5-03 (``tautau_detector``).

Every number printed here is read from ``presentation/data/ztautau_*.json`` (v3, DeepTau
Tight, MC-subtracted fake factor) and formatted in code; module-level asserts stop the
render on schema drift. Nine scenes, each opening on the previous scene's final frame
(pure builders ``state_*()`` + ``ORDER_*`` tuples, ``check_order`` at the end).

    TautauEvent (MovingCameraScene)
        tautau_a1_event         camera zooms out of 5-03; the schematic tracks fade; the real SR
                                event (run 281976) is drawn from its tau_h phi / charge / DM / p_T
                                and the real p_T^miss; p_T values; run/event stamp
        tautau_a2_visible_mass  the two p_T values converge to m_vis = 69.3 GeV; a mass scale
                                with the dashed m_Z shows it sits low
    TautauMass
        tautau_b1_neutrinos     two dashed nu arrows grow along the tau_h axes (head-to-tail),
                                their sum lands inside the p_T^miss resolution ellipse
        tautau_b2_likelihood    the event's (x1, x2) posterior lights up cell by cell; the formula
                                m_tautau = m_vis / sqrt(x1 x2); the median cell flashes; 86.6 GeV
        tautau_b3_shapes        the grid parks; simulated m_vis (grey ghost) vs m_tautau (red)
                                shapes with dashed m_Z and median ticks; 0.80 -> 0.99, 13% -> 11%
    TautauRain
        tautau_c1_first_entry   shapes fade, slice parks left, linear m_tautau axes (14 fit bins);
                                the event's m_tautau is the first entry
        tautau_c2_rain          clock, rain, the 14 data bins grow to the real counts; N = 21,160
    TautauStack
        tautau_d1_simulation    bars -> points; the simulation stack (rest / non-fid / fiducial)
                                slides in and reaches a third of the data; key
        tautau_d2_gap           the ratio panel (0-4) opens at data / simulation = 2.97
    TautauFakes
        tautau_e1_same_sign     plot parks; the slice in the dashed same-sign box with a real SS
                                pair: tau2 Tight (red), tau1 a jet with a grey tau_h-like core
        tautau_e2_tally         ten real SS pairs, a ten-box tally: f = N_T / N_L = 2/8
        tautau_e3_ff_map        the real 4 x 5 fake-factor map (period G, 0 jets; decay modes named);
                                the full dependence written as f = f(period, DM, N_jets, p_T)
        tautau_e4_eta_closure   first the disagreement: same-sign data vs FF prediction in |eta(tau1)|
                                with a ratio strip (0.93 ... 1.14 ... 0.84); the ratio becomes the
                                correction f(|eta|) and slides to 1; then the same for pT(tau2), g(pT)
        tautau_e5_osss          C_OS/SS = C(period, N_jets, D_BDT) in [0.96; 1.305]; box turns OS
    TautauTransfer
        tautau_f1_apply         N_fake = C_OS/SS f x N_AR; the plot returns to the main position
        tautau_f2_template      the fake template enters the bottom of the stack; the ratio drops
                                from 2.97 to 1.02; N_fake = 13592 (64%)
    TautauCorrections
        tautau_g1_tau_sf        the fiducial layer inflates to its raw yield, then the tau_h ID /
                                TES / trigger scale factors bring it down
        tautau_g2_event_weights pileup, L1 prefiring, the MC subtraction in the AR (2.5%), the
                                non-fiducial fraction (38%)
    TautauBDT (its own scene: everything of g2 fades to a blank frame first)
        tautau_h1_inputs        the 16 mass-agnostic input symbols, by importance, two columns
        tautau_h2_shapes        the list parks left; for the 8 leading inputs a small panel each:
                                unit-normalised fakes (pale slate fill) vs Z -> tautau (red line)
        tautau_h3_sketch        the list collapses to a bracket; three schematic trees (one
                                highlighted path each) summed into D_BDT
        tautau_h4_score         D_BDT (log y from 10^1.5), 4-layer stack, data, data/pred ratio panel,
                                cuts at 0.55 / 0.90, AUC, key
        tautau_h5_categories    three category panels; SR0 below 110 GeV greyed; yields
    TautauFit
        tautau_i1_fit           the mu_Z slider (top) and the pulls of 8 named NPs (below) go post-fit;
                                the panels too
        tautau_i2_impacts       grouped impacts as horizontal bars ("Gammas" printed as
                                "Template stat. (gamma)"; stat in red)
        tautau_i3_sigma_fid     sigma_fid = 4.82 +- 0.09_stat +- 0.47_syst pb
        tautau_i4_sigma_total   a large 1700-2400 pb axis: this work (red, +222/-194), the published
                                CMS 1952 +- 49 pb (60-120 GeV) and ATLAS 1981 +- 57 pb (66-116 GeV,
                                narrower window, noted here only) as slate points, the 1944.9 pb
                                prediction as a purple dashed line with its +15/-21 pb band (NNLO+NNLL
                                NNPDF3.1 uncertainty of CMS-SMP-20-004 Table 5, ztautau_reference.json)

Reserved areas (06 A1): nothing above y = 2.7 (title band) and nothing at x < -5.85, y > 0.22 (chapter
identifier), except the first frames of a1, which open on 5-03's zoomed event display.

Physics honesty (schematic parts on top of real numbers):
  * Track curvature follows kappa_from_pt (exaggerated scale, honest relative curvature); eta is
    not drawn (r-phi view); the p_T^miss arrow length is proportional to p_T with the event's
    57.4 GeV reaching the HCAL inner radius.
  * b1: the neutrino arrows use the event's posterior mode (x1, x2) of the MET likelihood; this
    event's collinear solution is undefined (m_col null: the p_T^miss points outside the tau_h
    wedge), so the head-to-tail sum lands inside the p_T^miss resolution ellipse (drawn from the
    event's MET covariance), not on the arrow tip. The nu2 arrow (0.4 GeV) is floored to a
    visible length.
  * b3: both shapes are the inclusive aMC@NLO Z -> tautau signal sample alone with unit weights
    (``signal_weighting == "inclusive_unit"``: one clean peak each; the stitched jet-binned mix
    had a shoulder); the medians are marked and the docs/04 scale / resolution numbers printed.
  * e4: the same-sign closure plots draw closure.eta / closure.pt2 obs and pred; the corrected
    prediction is pred x before (= obs by construction), the ratio after is the frozen "after".
  * h2: the input shapes are unit-normalised (fakes = AR data x FF, signal = unit-weight Z -> tautau).
  * h3: the trees are a schematic of a gradient-boosted classifier, not the real 300 trees.
  * The stack draws the clamped fit-input histograms (region-wise negative bins at 0); the
    printed data / simulation = 2.969 and data / pred. = 1.021 are the RESULTS.md totals
    (unclamped, 1 % below the drawn stack's total).
  * g1/g2: the fiducial layer is scaled uniformly by the frozen raw / corrected signal-yield
    ratio and the mean scale factors (shape assumed unchanged); the ratio dots follow.
  * h4: two negative "rest" score bins are drawn at 0; the ratio uses the frozen sr_score.pred.
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
    DOWN, LEFT, ORIGIN, RIGHT, UP, Circle, Create, DashedLine, DashedVMobject, DecimalNumber, Dot, Ellipse,
    FadeIn, FadeOut, Flash, GrowFromCenter, GrowFromEdge, LaggedStart, Line, MovingCameraScene,
    Rectangle, ReplacementTransform, RoundedRectangle, Scene, Succession, Transform, Triangle,
    ValueTracker, VGroup, VMobject, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403
from style.bnd_style import _cap_scale, _is_placeholder  # noqa: E402

# ---------------------------------------------------------------------------
# frozen data + anchors (a mismatch stops the render; never adjusted here)
# ---------------------------------------------------------------------------

EV = load_data("ztautau_events")
MA = load_data("ztautau_mass")
SR = load_data("ztautau_sr_stack")
FK = load_data("ztautau_fakes")
BD = load_data("ztautau_bdt")
CO = load_data("ztautau_corrections")
FIT = load_data("ztautau_fit")
BI = load_data("ztautau_bdt_inputs")
RF = load_data("ztautau_reference")

EVENT = EV["sr"]["chosen"][0]
SS_PAIRS = EV["ss_ff"]["chosen"]
REGIONS = ("tautau_SR0", "tautau_SR1", "tautau_SR2")

assert SR["data"]["total"] == 21160
assert sum(SR["data"]["inclusive"]) == 21160 and len(SR["edges"]) == 15
assert abs(SR["fakes"]["total"] - 13592) < 1 and abs(FK["yields"]["n_fake"]["total"] - 13592) < 1
assert tuple(SR["stack_order"]) == ("Fakes", "rest", "DYtautau_nonfid", "DYtautau")
assert abs(FIT["poi"]["value"] - 1.071) < 5e-4
assert round(FIT["poi"]["err_up"], 3) == 0.114 and round(FIT["poi"]["err_down"], 3) == 0.100
assert abs(FIT["sigma_60_120"]["value"] - 2082) < 1
assert round(FIT["sigma_60_120"]["err_up"]) == 222 and round(FIT["sigma_60_120"]["err_down"]) == 194
assert round(FIT["sigma_60_120"]["pred"], 1) == 1944.9
assert (round(FIT["sigma_fid"]["value"], 2), round(FIT["sigma_fid"]["stat"], 2),
        round(FIT["sigma_fid"]["syst"], 2)) == (4.82, 0.09, 0.47)
assert round(EVENT["m_tt"], 1) == 86.6 and round(EVENT["m_vis"], 1) == 69.3
assert EVENT["run"] == 281976 and EVENT["event"] == 2983852771
assert len(SS_PAIRS) == 10 and sum(bool(p["passes"]) for p in SS_PAIRS) == 2
assert [BD["categories"][r]["data"] for r in REGIONS] == [15742, 2704, 2714]
assert len(MA["signal"]["m_vis"]) == 40 and len(MA["signal"]["m_tt"]) == 40 and len(MA["edges"]) == 41
assert len(FK["closure"]["eta"]["before"]) == 6 and len(FK["closure"]["pt2"]["before"]) == 5
assert abs(FIT["grouped_impact"]["Tau ID"] - 0.084) < 0.001
assert len(FIT["nps_shown"]) == 8 and len(FIT["impact_order"]) == 11
assert np.asarray(EVENT["posterior"]["weights"]).shape == (12, 12)
assert MA["signal_weighting"] == "inclusive_unit"
assert len(BI["features"]) == 16 and len(BI["dm_labels"]) == 4
assert all(f["name"] == n for f, n in zip(BI["features"], BI["feature_order"]))
assert abs(RF["this_work"]["value"] - FIT["sigma_60_120"]["value"]) < 1e-6
assert RF["cms"]["value"] == 1952 and RF["atlas"]["value"] == 1981 and round(RF["theory"]["value"], 1) == 1944.9
assert (round(RF["theory"]["unc_up"]), round(RF["theory"]["unc_down"])) == (15, 21)
assert [round(x, 3) for x in FK["osss"]["range_per_category"]] == [0.960, 1.305]
for r in REGIONS:
    assert FIT["regions"]["prefit"][r]["data"] == SR["data"]["per_region"][r]

EDGES = np.asarray(SR["edges"], dtype=float)
NB = len(EDGES) - 1
XC = 0.5 * (EDGES[:-1] + EDGES[1:])
HW = 0.5 * (EDGES[1:] - EDGES[:-1])
DATA = np.asarray(SR["data"]["inclusive"], dtype=float)
TOTAL = int(SR["data"]["total"])
FAKES = np.asarray(SR["fakes"]["inclusive"], dtype=float)
LAYERS = tuple(SR["stack_order"])                  # bottom-up
GROUP_INCL = {n: np.clip(np.asarray(SR["groups_counts"][n]["inclusive_clamped"], dtype=float), 0, None)
              for n in ("rest", "DYtautau_nonfid", "DYtautau")}
REST_NAMES = tuple(SR["groups"]["rest"])
LAYER_COL = {"Fakes": SAMPLE["Fakes"], "rest": SAMPLE["TTbar"],
             "DYtautau_nonfid": lighten(CHANNEL["tautau"], 0.55), "DYtautau": CHANNEL["tautau"]}
M_Z = float(EV["m_Z"])
MEANS = CO["means_on_fiducial_signal_in_sr"]
K_RAW = float(MEANS["signal_yield"]["raw"]) / float(MEANS["signal_yield"]["corrected"])
K_SF = K_RAW * float(MEANS["sf_id_both_legs"]) * float(MEANS["sf_trig_both_legs"])
DY_TOTAL = float(SR["groups_counts"]["DYtautau"]["total"])
DM_TEX = {str(k): v for k, v in BI["dm_labels"].items()}      # decay modes named, never coded
DMS = ("0", "1", "10", "11")
FEATURES = {f["name"]: f for f in BI["features"]}
SHAPE_FEATURES = ("dr_tt", "t1_pt", "pt_vis", "dphi_tt", "pt_tt", "t2_pt", "t1_abseta", "t2_abseta")
IMPACT_NAME = {"Gammas": "Template stat. (γ)"}
_CAT = {"c0": r"D_{\mathrm{BDT}} < 0.55", "c1": r"0.55 < D_{\mathrm{BDT}} < 0.90", "c2": r"D_{\mathrm{BDT}} > 0.90"}
NP_TEX = {                                                     # pull labels: named, never coded (06 B5)
    **{f"TauID_DM{d}": rf"\tau_h\ \mathrm{{ID}}\ ({DM_TEX[d]})" for d in DMS},
    **{f"TauTrigger_DM{d}": rf"\tau_h\ \mathrm{{trigger}}\ ({DM_TEX[d]})" for d in DMS},
    **{f"FakeOSSS_tautau_{c}": rf"C_{{\mathrm{{OS/SS}}}}\ ({t})" for c, t in _CAT.items()},
    **{f"FakeClosure_tautau_{c}_{r}": rf"f\ \mathrm{{closure}}\ ({t},\ m_{{\tau\tau}} {o} 110)"
       for c, t in _CAT.items() for r, o in (("lo", "<"), ("hi", ">"))},
    "Lumi": r"\mathrm{Luminosity}",
}
assert all(n in NP_TEX for n in FIT["nps_shown"]), [n for n in FIT["nps_shown"] if n not in NP_TEX]
PRED_TOTAL = float(SR["totals"]["pred"])

# ---------------------------------------------------------------------------
# one coordinate dictionary (06 section B3 defaults unless the brief says otherwise)
# ---------------------------------------------------------------------------

EASE = rate_functions.ease_in_out_sine
PLOT_C = np.array([1.85, -0.25, 0.0])              # every plot placement scales about this
A = dict(
    phi_503=(math.radians(40), math.radians(218)), d_nu_503=math.radians(6),
    zoom_503=(0.5, (-0.2, -0.25)),                 # s5_ztautau.ZOOM_503: 5-03 ends zoomed here (corner clear)
    r_label=R_DET + 0.42, y_safe=2.42,
    det_park=(0.42, (-4.6, -0.7)), det_cr=(0.84, (-1.2, -0.45)),
    cr_box=dict(center=(-1.2, -0.45), width=6.2, height=6.2),
    plot_main_c=(1.85, 0.85), ratio_c=(1.85, -1.55),
    plot_main=(1.0, PLOT_C), plot_ghost=(0.44, (4.75, 1.24)), plot_left=(0.58, (-4.2, -1.5)),
    key_left=(-1.75, -2.92), key_buff=0.2, key_label_h=0.17, key_i=(-6.6, -3.3), key_i_scale=0.72,
    counter_at=(-0.2, 2.12), clock_at=(-0.72, 2.24), stamp_at=(-4.6, -2.3),
    pt_at=((4.7, 1.35), (4.7, 0.65), (4.7, -0.05)), mass_at=(4.7, 2.2),
    mscale=dict(x0=3.2, x1=6.2, y=-2.3, lo=40.0, hi=120.0),
    grid_c=(4.7, 0.65), grid_cell=0.2, grid_park=(0.36, (6.5, -2.0)), formula_at=(4.7, -1.35),
    shapes_c=(4.85, 0.6), shapes_len=(3.4, 1.9), perf_at=((4.35, -1.5), (4.35, -2.0)),
    y_top=3500.0, ratio_wide=(0, 4, 2), ratio_zoom=(0.5, 1.5, 0.5),
    tally_c=(4.7, -0.5), ff_tex_c=(4.7, -1.05), ff_tex_end=(4.7, -0.5), ffmap_c=(4.9, -2.45),
    ffdep_at=(4.85, -3.5),
    clos_main=(4.3, -1.6), clos_len=(2.9, 1.0), clos_ratio=(4.3, -2.8), clos_rlen=0.55,   # f(|eta|) label fits right
    osss_at=((4.6, -1.35), (4.6, -1.9)),
    formula_f=(-4.75, -0.35), n_fake_at=(-4.75, -0.35), pct_at=(-4.75, -1.05),
    left_col_x=-6.95, corr_y0=-0.02, corr_dy=0.44, corr_h=0.175,
    inputs_cols=(-4.2, 0.9), inputs_y0=2.25, inputs_dy=0.55, inputs_h=0.24,
    inputs_park=dict(x=-5.75, y0=2.3, dy=0.31, h=0.15),
    fpanel=dict(xs=(-2.85, -0.35, 2.15, 4.65), ys=(1.25, -1.3), len=(2.0, 1.3), key_at=(-0.2, -2.95)),
    tree=dict(xs=(-2.6, 0.3, 3.2), y=0.3, dy=0.95, dx=(0.62, 0.32), r=0.15, bracket_x=-5.5,
              plus_dx=1.45, out_at=(5.9, 0.3)),
    score_y=(1.5, 4.2),
    panel_x=(-0.45, 1.85, 4.15), panel_w=1.75, panel_h=2.5, panel_cy=0.85, panel_lab_y=(2.57, 2.3),
    pulls_at=(-3.0, -1.6), slider_at=(-3.9, 1.35), mu_seed_dy=0.5, mu_seed_i2=(-4.1, -2.62),
    imp_name_r=-4.25, imp_bar_x0=-4.15, imp_y0=1.3, imp_dy=0.3, imp_per_pct=0.18, imp_head=(-3.3, 1.68),
    mu_line_at=(3.0, 2.2), sig_fid_at=(3.0, 1.45), sig_label_at=(3.0, 0.75),
    sig_ax=dict(lo=1700.0, hi=2400.0, x0=-0.85, length=7.1, y=-2.3, rows=(2.0, 1.3, 0.6), height=2.6),
)
P_YTOP = {"tautau_SR0": (3500.0, [0, 1000, 2000, 3000]), "tautau_SR1": (800.0, [0, 400, 800]),
          "tautau_SR2": (1000.0, [0, 500, 1000])}
DOT_R, RDOT_R = 0.036, 0.034


# ---------------------------------------------------------------------------
# small helpers (structure copied from s4_zmumu_story)
# ---------------------------------------------------------------------------

def tex_h(expr: str, h: float, color=INK):
    return mathtex(expr, color=color).scale(_cap_scale("tex", h))


def text_h(s: str, h: float, color=INK):
    return text(s, color=color).scale(_cap_scale("text", h))


def _p3(p):
    p = np.asarray(p, dtype=float)
    return p if p.shape == (3,) else np.array([p[0], p[1], 0.0])


def _q(rec) -> int:
    return int(np.sign(float(rec["charge"])))


def _prongs(rec) -> int:
    return 1 if int(rec["dm"]) in (0, 1) else 3


def _pi0(rec) -> bool:
    return int(rec["dm"]) in (1, 11)


def phi_end(trk) -> float:
    d = trk.pts[-1] - DET_CENTER
    return math.atan2(d[1], d[0])


def label_spot(phi: float, placement=(1.0, DET_CENTER)) -> np.ndarray:
    s, c = placement
    c = _p3(c)
    r = s * A["r_label"]

    def at(a):
        return c + r * np.array([math.cos(a), math.sin(a), 0.0])

    p = at(phi)
    if p[1] <= A["y_safe"]:
        return p
    a0 = math.asin(max(-1.0, min(1.0, (A["y_safe"] - c[1]) / r)))
    ph = phi % (2 * math.pi)
    return at(min((a0, math.pi - a0), key=lambda a: abs(a - ph)))


def scale_strokes(mob, f: float):
    for m in mob.get_family():
        if isinstance(m, VMobject):
            m.set_stroke(width=m.get_stroke_width() * f, family=False)
            m.set_stroke(width=m.get_stroke_width(background=True) * f, background=True, family=False)
    return mob


def _look(m):
    out = []
    for s in m.get_family():
        if not s.has_points():
            continue
        if isinstance(s, VMobject):
            fo = float(np.max(s.get_fill_opacity()))
            sw = float(np.max(s.get_stroke_width()))
            so = float(np.max(s.get_stroke_opacity()))
        else:
            fo, sw, so = 1.0, 0.0, 0.0
        pts = s.points[:, :2]
        lo, hi = pts.min(axis=0), pts.max(axis=0)
        area = float(np.prod(hi - lo))
        fill_vis = fo > 0.01 and area > 1e-7
        stroke_vis = sw > 0.01 and so > 0.01
        if not (fill_vis or stroke_vis):
            continue
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


def looks_same(live, target, tol=6e-3):
    la, lb = _look(live), _look(target)
    sa, sb = {s for s, _ in la}, {s for s, _ in lb}
    if sa != sb:
        return False, f"styles differ: {sorted(sa ^ sb)[:4]}"
    for st in sa:
        Aa = np.array([bb for s, bb in la if s == st])
        Bb = np.array([bb for s, bb in lb if s == st])
        if not _match(Aa, Bb, tol):
            return False, f"geometry differs for style {st}"
    return True, ""


def settle(scene, live, target, what=""):
    ok, why = looks_same(live, target)
    if not ok:
        raise AssertionError(f"settle({what}): live mobject does not match the builder: {why}")
    scene.replace(live, target)


def adopt(scene, parent, what=""):
    fam = {id(m) for m in parent.get_family()}

    def belongs(m):
        if id(m) in fam:
            return True
        pts = [x for x in m.get_family() if x.has_points()]
        return bool(pts) and all(id(x) in fam for x in pts)

    idx = [i for i, m in enumerate(scene.mobjects) if belongs(m)]
    if not idx:
        raise AssertionError(f"adopt({what}): nothing of it is in the scene")
    for m in scene.mobjects[idx[0]:idx[-1] + 1]:
        if not belongs(m) and not _is_placeholder(m) and not isinstance(m, ValueTracker):
            raise AssertionError(f"adopt({what}): {type(m).__name__} interleaved")
    shown = {id(x) for i in idx for x in scene.mobjects[i].get_family()}
    missing = [x for x in parent.get_family() if x.has_points() and id(x) not in shown]
    if missing:
        raise AssertionError(f"adopt({what}): {len(missing)} drawn members were never revealed")
    keep = [m for i, m in enumerate(scene.mobjects) if i not in set(idx)]
    keep.insert(idx[0], parent)
    scene.mobjects = keep


def hidden_rain(*args, **kw) -> list:
    out = []
    for a in rain(*args, **kw):
        move, fade = a.animations
        out.append(Succession(FadeIn(move.mobject, run_time=0.06), move, fade))
    return out


def settle_same(scene, prev, end, keys):
    for k in keys:
        settle(scene, prev[k], end[k], k)


# ---------------------------------------------------------------------------
# detector builders
# ---------------------------------------------------------------------------

def parked_slice(placement=None) -> CMSSlice:
    return place(CMSSlice(), placement or A["det_park"], DET_CENTER)


def tau_pieces(det: CMSSlice, rec, seed: int = 0) -> VGroup:
    """One reconstructed tau_h as 5-03 draws it: VGroup(prong tracks, tracker hits,
    deposits); prongs from the decay mode, bend sign from the charge, curvature from p_T."""
    assert _prongs(rec) == int(rec["prongs"]) and _pi0(rec) == bool(rec["pi0"])
    sig = signature(det, "tau_h", float(rec["phi"]), charge=_q(rec), kappa=kappa_from_pt(rec["pt"]),
                    prongs=_prongs(rec), pi0=_pi0(rec), seed=seed)
    trk, deposits = sig[0], VGroup(*sig[1:])
    hits = tracker_hits(det, trk[len(trk) // 2], color=CHANNEL["tautau"])
    g = VGroup(trk, hits, deposits)
    g.trk, g.hits, g.deposits, g.record = trk, hits, deposits, rec
    g.lines, g.rest = trk, VGroup(hits, deposits)
    return g


def fake_tau_pieces(det: CMSSlice, rec, seed: int = 0) -> VGroup:
    """A tau_h candidate failing Tight: a jet spray with a grey tau_h-like core."""
    rng = np.random.default_rng(seed + 100)
    phi, q, k = float(rec["phi"]), _q(rec), kappa_from_pt(rec["pt"])
    jet = signature(det, "jet", phi, seed=seed)
    cc = col(PARTICLE["hadron"])
    dphis = (0.0,) if _prongs(rec) == 1 else (-0.07, 0.0, 0.07)
    trk = VGroup()
    for i, dphi in enumerate(dphis):
        qq = q if len(dphis) == 1 else (q, -q, q)[i]
        trk.add(track(det, phi + dphi, qq * k * (1.0 + rng.uniform(-0.3, 0.6)), r_end=det.radii["hcal"][0],
                      color=cc, sw=2.8, outline=darken(cc, 0.4), outline_sw=1.8))
    dep = calo_hit(det, "hcal", phi_end(trk[len(dphis) // 2]), width=1 if len(dphis) == 1 else 3,
                   extent=0.85, side_opacity=0.5)
    g = VGroup(jet, trk, dep)
    g.jet, g.trk, g.dep, g.record = jet, trk, dep, rec
    g.lines, g.rest = VGroup(jet[0], trk), VGroup(jet[1], jet[2], dep)
    return g


def met_arrow(det: CMSSlice, pt: float, phi: float) -> VGroup:
    """p_T^miss arrow: length proportional to p_T, the event's value reaching the HCAL inner radius."""
    r_end = det.radii["hcal"][0] * float(pt) / float(EVENT["met"]["pt"])
    return signature(det, "met", float(phi), r_end=r_end)


def r_per_gev(det: CMSSlice) -> float:
    return det.radii["hcal"][0] / float(EVENT["met"]["pt"])


def event_taus() -> tuple:
    det0 = CMSSlice()
    return tau_pieces(det0, EVENT["t1"], seed=11), tau_pieces(det0, EVENT["t2"], seed=12)


def event_met() -> VGroup:
    return met_arrow(CMSSlice(), EVENT["met"]["pt"], EVENT["met"]["phi"])


def nu_arrow(p0, p1, color=PARTICLE["nu"], sw: float = 2.8) -> VGroup:
    line = DashedLine(_p3(p0), _p3(p1), dash_length=0.1, stroke_color=col(color), stroke_width=sw)
    return VGroup(line, arrow_tip_on(Line(_p3(p0), _p3(p1)), color=color, at=1.0, tip_length=0.14))


def nu_geometry(det: CMSSlice) -> dict:
    """Collinear nu momenta at the posterior mode: p_nu = p_T (1 - x) / x along each tau_h;
    the second one is floored to a visible length (0.4 GeV in this event)."""
    x1, x2 = EVENT["posterior"]["mode_x"]
    p1 = float(EVENT["t1"]["pt"]) * (1.0 - x1) / x1
    p2 = float(EVENT["t2"]["pt"]) * (1.0 - x2) / x2
    k = r_per_gev(det)
    L1, L2 = max(p1 * k, 0.12), max(p2 * k, 0.12)
    e1, e2 = unit(float(EVENT["t1"]["phi"])), unit(float(EVENT["t2"]["phi"]))
    P1 = det.c + L1 * e1
    P2 = P1 + L2 * e2
    return dict(P1=P1, P2=P2)


def event_nus() -> VGroup:
    det0 = CMSSlice()
    g = nu_geometry(det0)
    return VGroup(nu_arrow(det0.c, g["P1"]), nu_arrow(g["P1"], g["P2"]))


def met_ellipse() -> Ellipse:
    det0 = CMSSlice()
    m = EVENT["met"]
    cov = np.array([[m["covxx"], m["covxy"]], [m["covxy"], m["covyy"]]], dtype=float)
    lam, vec = np.linalg.eigh(cov)
    k = r_per_gev(det0)
    ell = Ellipse(width=2 * math.sqrt(lam[1]) * k, height=2 * math.sqrt(lam[0]) * k,
                  fill_color=col(PARTICLE["met"]), fill_opacity=0.15, stroke_width=0)
    ell.rotate(math.atan2(vec[1, 1], vec[0, 1]))
    tip = det0.c + det0.radii["hcal"][0] * unit(float(m["phi"]))
    return ell.move_to(tip)


def run_stamp() -> VGroup:
    s = text_h(f"Run {EVENT['run']}   Event {EVENT['event']}", 0.14, color=GREY)
    return s.move_to(_p3(A["stamp_at"]))


def pt_labels() -> dict:
    t1, t2, m = EVENT["t1"], EVENT["t2"], EVENT["met"]
    labs = {}
    for key, expr, at in (("pt1", rf"p_T(\tau_h^{{{'+' if _q(t1) > 0 else '-'}}}) = {t1['pt']:.1f}\ \mathrm{{GeV}}",
                           A["pt_at"][0]),
                          ("pt2", rf"p_T(\tau_h^{{{'+' if _q(t2) > 0 else '-'}}}) = {t2['pt']:.1f}\ \mathrm{{GeV}}",
                           A["pt_at"][1]),
                          ("pt_met", rf"p_T^{{\,\mathrm{{miss}}}} = {m['pt']:.1f}\ \mathrm{{GeV}}", A["pt_at"][2])):
        labs[key] = tex_h(expr, 0.22).move_to(_p3(at))
    return labs


def mass_label(kind: str):
    if kind == "vis":
        return tex_h(rf"m_{{\mathrm{{vis}}}} = {EVENT['m_vis']:.1f}\ \mathrm{{GeV}}", 0.28).move_to(_p3(A["mass_at"]))
    return tex_h(rf"m_{{\tau\tau}} = {EVENT['m_tt']:.1f}\ \mathrm{{GeV}}", 0.28).move_to(_p3(A["mass_at"]))


def mass_scale(value: float) -> VGroup:
    """A short mass axis (40-120 GeV) with the dashed m_Z and a red dot at ``value``."""
    S = A["mscale"]
    ink = col(INK)

    def xs(v):
        return S["x0"] + (float(v) - S["lo"]) / (S["hi"] - S["lo"]) * (S["x1"] - S["x0"])

    y = S["y"]
    axis = Line([S["x0"], y, 0], [S["x1"], y, 0], stroke_color=ink, stroke_width=2.5)
    ticks, labels = VGroup(), VGroup()
    for v in (40, 80, 120):
        p = np.array([xs(v), y, 0.0])
        ticks.add(Line(p, p + DOWN * 0.08, stroke_color=ink, stroke_width=2.5))
        labels.add(text_h(f"{v}", 0.15).next_to(p + DOWN * 0.08, DOWN, buff=0.08))
    mz = DashedLine([xs(M_Z), y, 0], [xs(M_Z), y + 0.45, 0], dash_length=0.08, stroke_color=col(GREY),
                    stroke_width=2.5)
    mz_lab = tex_h(r"m_Z", 0.17, color=GREY).next_to(mz.get_end(), UP, buff=0.06)
    dot = Dot([xs(value), y, 0], radius=0.07, color=col(CHANNEL_LINE["tautau"]))
    g = VGroup(axis, ticks, labels, mz, mz_lab, dot)
    g.frame, g.dot = VGroup(axis, ticks, labels, mz, mz_lab), dot
    return g


# -- the likelihood grid and the mass shapes ----------------------------------------

def heat_grid(parked: bool = False) -> VGroup:
    """The event's 12 x 12 (x1, x2) posterior as a heat map without numbers: column = x1
    bin, row = x2 bin (upward)."""
    W = np.asarray(EVENT["posterior"]["weights"], dtype=float)
    a = W / W.max()
    cw = A["grid_cell"]
    c0 = _p3(A["grid_c"])
    cells = VGroup()
    for i2 in range(12):
        for i1 in range(12):
            fill = lighten(CHANNEL["tautau"], 1.0 - (0.06 + 0.86 * float(a[i1, i2])))
            cells.add(Rectangle(width=cw, height=cw, fill_color=fill, fill_opacity=1.0, stroke_color=col(WHITE),
                                stroke_width=1.0).move_to(c0 + np.array([(i1 - 5.5) * cw, (i2 - 5.5) * cw, 0.0])))
    if parked:
        s, c = A["grid_park"]
        cells.scale(s).move_to(_p3(c))
    return cells


def grid_labels(cells) -> VGroup:
    x1 = tex_h(r"x_1", 0.22).next_to(cells, DOWN, buff=0.12)
    x2 = tex_h(r"x_2", 0.22).next_to(cells, LEFT, buff=0.12)
    return VGroup(x1, x2)


def median_cell_index() -> int:
    """The high-weight cell whose block mass is nearest the posterior median (the flash)."""
    W = np.asarray(EVENT["posterior"]["weights"], dtype=float)
    Mb = np.asarray(EVENT["posterior"]["m_block_median"], dtype=float)
    med = float(EVENT["posterior"]["m_median"])
    best, best_d = None, 1e9
    for i1 in range(12):
        for i2 in range(12):
            if W[i1, i2] >= 0.5 * W.max() and abs(Mb[i1, i2] - med) < best_d:
                best, best_d = (i1, i2), abs(Mb[i1, i2] - med)
    i1, i2 = best
    return i2 * 12 + i1


def mass_formula():
    return tex_h(r"m_{\tau\tau} = m_{\mathrm{vis}} / \sqrt{x_1 x_2}", 0.24).move_to(_p3(A["formula_at"]))


def shapes_parts() -> dict:
    e = np.asarray(MA["edges"], dtype=float)
    vis = np.clip(np.asarray(MA["signal"]["m_vis"], dtype=float), 0, None)
    tt = np.clip(np.asarray(MA["signal"]["m_tt"], dtype=float), 0, None)
    top = float(max(vis.max(), tt.max())) * 1.15
    dax = DataAxes([0, 200, 50], [0, top, top], A["shapes_len"][0], A["shapes_len"][1], y_ticks=[],
                   show_y_labels=False, tick_label_h=0.18, title_h=0.18, title_buff=0.14,
                   x_title=r"m\ [\mathrm{GeV}]")
    dax.move_frame_to(A["shapes_c"])
    ghost = step_hist(dax, e, vis, color=GREY, stroke_width=2.0, fill_opacity=0.25)
    red = step_hist(dax, e, tt, color=CHANNEL_LINE["tautau"], stroke_width=3.0)
    mz = vref_line(dax, M_Z, color=GREY, stroke_width=2.0, num_dashes=18)
    med = VGroup()
    for key, cc in (("m_vis", GREY), ("m_tt", CHANNEL_LINE["tautau"])):
        x = float(MA["signal"]["median_over_mZ"][key]) * M_Z
        med.add(DashedVMobject(Line(dax.c2p(x, 0.0), dax.c2p(x, 0.55 * top), stroke_color=col(cc),
                                    stroke_width=2.5), num_dashes=8))
    lab = VGroup(tex_h(r"m_{\mathrm{vis}}", 0.2, color=GREY).move_to(dax.c2p(38.0, 0.8 * top)),
                 tex_h(r"m_{\tau\tau}", 0.2, color=CHANNEL_LINE["tautau"]).move_to(dax.c2p(150.0, 0.8 * top)))
    pf = MA["performance"]
    perf = VGroup(
        tex_h(rf"\langle m\rangle / m_Z:\ {pf['m_vis']['scale']:.2f} \to {pf['m_tt']['scale']:.2f}", 0.22)
        .move_to(_p3(A["perf_at"][0])),
        tex_h(rf"\sigma_m / m:\ {100 * pf['m_vis']['core_res']:.0f}\% \to {100 * pf['m_tt']['core_res']:.0f}\%", 0.22)
        .move_to(_p3(A["perf_at"][1])))
    return {"shp_dax": dax, "shp_vis": ghost, "shp_tt": red, "shp_mz": mz, "shp_med": med, "shp_lab": lab,
            "perf": perf}


SHAPE_KEYS = ("shp_dax", "shp_vis", "shp_tt", "shp_mz", "shp_med", "shp_lab", "perf")


# ---------------------------------------------------------------------------
# the main m_tautau plot (always at the natural position, then place())
# ---------------------------------------------------------------------------

def sr_layers(fakes: bool, dy_scale: float = 1.0, mc: bool = True):
    """The four layers bottom-up; ``mc=False`` = every layer at zero (the invisible,
    Transform-able stack of the clips before the simulation arrives)."""
    rows = []
    for n in LAYERS:
        if n == "Fakes":
            c = FAKES if fakes else np.zeros(NB)
        elif not mc:
            c = np.zeros(NB)
        elif n == "DYtautau":
            c = GROUP_INCL[n] * dy_scale
        else:
            c = GROUP_INCL[n]
        rows.append((n, c, LAYER_COL[n]))
    return rows


def main_axes(x_labels: bool = False) -> DataAxes:
    dax = DataAxes([0, 350, 50], [0, A["y_top"], 500], 6.4, 2.9, y_ticks=[0, 1000, 2000, 3000],
                   show_x_labels=x_labels, title_h=0.2, title_buff=0.16,
                   x_title=r"m_{\tau\tau}\ [\mathrm{GeV}]" if x_labels else None,
                   y_title=r"\mathrm{events\,/\,bin}")
    dax.move_frame_to(A["plot_main_c"])
    return dax


def entry_dot(dax, m) -> Dot:
    return Dot(dax.c2p(float(m), 0.03 * A["y_top"]), radius=0.05, color=col(SAMPLE["Data"]))


def data_bars(dax, visible=True) -> VGroup:
    g = VGroup()
    for i in range(NB):
        b = data_bar(dax, XC[i], DATA[i] if visible else 0.0, HW[i], color=SAMPLE["Data"], fill_opacity=0.55,
                     stroke_width=1.0)
        if not visible:
            b.set_opacity(0.0)
        g.add(b)
    return g


def data_dots(dax) -> VGroup:
    return VGroup(*[data_dot(dax, XC[i], DATA[i], color=SAMPLE["Data"], radius=DOT_R) for i in range(NB)])


def layout_counter(pre, num):
    num.next_to(pre, RIGHT, buff=0.14)
    num.shift(UP * (pre.get_bottom()[1] - num[0].get_bottom()[1]))
    return num


def make_counter(value) -> VGroup:
    pre = tex_h(r"N =", 0.24)
    pre.shift(_p3(A["counter_at"]) - pre.get_corner(DOWN + LEFT))
    num = DecimalNumber(float(value), num_decimal_places=0, group_with_commas=True, color=col(INK))
    num.scale(_cap_scale("tex", 0.24))
    layout_counter(pre, num)
    g = VGroup(pre, num)
    g.pre, g.num = pre, num
    return g


def set_counter(counter, value):
    counter.num.set_value(float(round(value)))
    layout_counter(counter.pre, counter.num)


def make_key(rows: int) -> VGroup:
    """One-row key: Data, fiducial Z->tautau, non-fiducial Z/gamma*->tautau, the striped
    'rest' swatch and (rows == 5) the fakes. colour_key takes 4 rows, so two are merged."""
    entries = [(r"\mathrm{Data}", SAMPLE["Data"], "dot"), (r"Z\to\tau\tau", LAYER_COL["DYtautau"]),
               (r"Z/\gamma^{*}\to\tau\tau", LAYER_COL["DYtautau_nonfid"]),
               (r"t\bar{t},\ tW,\ VV,\ W", LAYER_COL["rest"])]
    if rows == 5:
        entries.append((r"\mathrm{Fakes}", LAYER_COL["Fakes"]))
    keys = [colour_key(entries[:4], label_h=A["key_label_h"])]
    if rows == 5:
        keys.append(colour_key(entries[4:], label_h=A["key_label_h"]))
    row_list = [r for k in keys for r in list(k.rows)]
    row = row_list[3]
    sw = row[1]
    names = ("TTbar", "SingleTop", "WW", "WZ", "ZZ")
    w, h = sw.width, sw.height
    stripes = VGroup(*[Rectangle(width=w / len(names), height=h, stroke_width=0, fill_color=col(SAMPLE[n]),
                                 fill_opacity=0.95) for n in names]).arrange(RIGHT, buff=0).move_to(sw)
    frame = Rectangle(width=w, height=h, fill_opacity=0.0, stroke_color=darken(SAMPLE["TTbar"], 0.25),
                      stroke_width=1.0).move_to(sw)
    row.submobjects[1] = VGroup(stripes, frame)
    rows_g = VGroup(*row_list).arrange(RIGHT, buff=A["key_buff"])
    rows_g.shift(_p3(A["key_left"]) - rows_g.get_left())
    key = VGroup(rows_g)
    key.rows = rows_g
    return key


def key_i() -> VGroup:
    key = make_key(5)
    key.scale(A["key_i_scale"], about_point=key.get_left())
    key.shift(_p3(A["key_i"]) - key.get_left())
    return key


def ratio_parts(dax, total, wide: bool) -> VGroup:
    """data / prediction panel; on the wide (0-4) panel bins above the top are drawn as
    small triangles at the top edge (off scale)."""
    yr = A["ratio_wide"] if wide else A["ratio_zoom"]
    kw = dict(y_range=(yr[0], yr[1], yr[2]), y_length=0.9, y_ticks=list(np.arange(yr[0], yr[1] + 1e-9, yr[2])),
              y_fmt="{:.0f}" if wide else "{:.1f}", dot_radius=RDOT_R,
              x_title=r"m_{\tau\tau}\ [\mathrm{GeV}]", y_title=r"\mathrm{data/pred.}", title_h=0.2, title_buff=0.16)
    rp = ratio_panel(dax, EDGES, DATA, total, A["ratio_c"], **kw)
    for j, i in enumerate(rp.bins):
        if rp.ratio[i] > yr[1]:
            tri = Triangle(fill_color=col(SAMPLE["Data"]), fill_opacity=1.0, stroke_width=0)
            tri.scale_to_fit_height(0.11).move_to(rp.dax.c2p(XC[i], yr[1]) + DOWN * 0.07)
            rp.dots.submobjects[j] = tri
    return rp


def ratio_value(fakes: bool, dy_scale: float = 1.0) -> float:
    t = SR["totals"]
    if fakes:
        if dy_scale != 1.0:
            return TOTAL / (PRED_TOTAL + DY_TOTAL * (dy_scale - 1.0))
        return float(t["data_over_pred"])
    return float(t["data_over_mc"])


def n_fakes_big() -> VGroup:
    y = FK["yields"]
    return tex_h(rf"N_{{\mathrm{{fake}}}} = {y['n_fake']['total']:.0f}", 0.32).move_to(_p3(A["n_fake_at"]))


def pct_big() -> VGroup:
    return tex_h(rf"{100 * FK['yields']['fraction_of_sr']:.0f}\%", 0.32).move_to(_p3(A["pct_at"]))


def n_fakes_small() -> VGroup:
    y = FK["yields"]
    lab = tex_h(rf"{y['n_fake']['total']:.0f}\ ({100 * y['fraction_of_sr']:.0f}\%)", 0.18)
    dax = main_axes()
    lab.next_to(dax.c2p(350.0, 0.5 * FAKES[-1]), RIGHT, buff=0.18)
    return lab


def plot_parts(placement=None, *, fakes=False, dy_scale=1.0, mc=True, data="dots", counter=True, ratio=True,
               wide=None, rlabel=True, key_rows=4, entries=(), x_labels=False) -> dict:
    placement = placement or A["plot_main"]
    P = {}
    dax = main_axes(x_labels=x_labels)
    P["dax"] = dax
    P["stack"] = stack_hist(dax, EDGES, sr_layers(fakes, dy_scale, mc))
    if data == "entries":
        P["data"] = VGroup(*[entry_dot(dax, m) for m in entries])
    elif data == "bars":
        P["data"] = data_bars(dax)
    elif data == "bars0":
        P["data"] = data_bars(dax, visible=False)
    else:
        P["data"] = data_dots(dax)
    if counter:
        P["counter"] = make_counter(TOTAL)
    if ratio:
        wide = (not fakes) if wide is None else wide
        P["ratio"] = ratio_parts(dax, P["stack"].total, wide)
        if rlabel:
            v = ratio_value(fakes, dy_scale)
            lab = tex_h(f"{v:.3f}", 0.2)
            lab.next_to(P["ratio"].dax.c2p(350.0, min(v, P["ratio"].dax.y_top)), RIGHT, buff=0.22)
            P["rlabel"] = lab
    if key_rows:
        P["key"] = make_key(key_rows)
    for m in P.values():
        place(m, placement, PLOT_C)
    return P


PLOT_KEYS = ("dax", "stack", "data", "counter", "ratio", "rlabel", "key")


# ---------------------------------------------------------------------------
# control region: box, pairs, tally, fake-factor map, closure strips, C_OS/SS
# ---------------------------------------------------------------------------

def cr_box(sym_tex=r"\tau_h^{\pm}\tau_h^{\pm}") -> VGroup:
    B = A["cr_box"]
    rect = RoundedRectangle(width=B["width"], height=B["height"], corner_radius=0.25,
                            stroke_color=col(DETECTOR_ACCENT), stroke_width=3.0, fill_opacity=0.0)
    rect.move_to(_p3(B["center"]))
    dashed_rect = DashedVMobject(rect, num_dashes=70)
    sym = tex_h(sym_tex, 0.24, color=DETECTOR_ACCENT)
    sym.move_to(rect.get_corner(UP + LEFT) + np.array([0.7, -0.36, 0.0]))
    g = VGroup(dashed_rect, sym)
    g.rect, g.sym = dashed_rect, sym
    return g


def ss_pair(i: int) -> VGroup:
    """Real same-sign pair i on a natural slice, placed into the box: tau2 (Tight, red) and
    tau1 (red if it passes Tight, else a jet with a grey tau_h-like core)."""
    rec = SS_PAIRS[i]
    det0 = CMSSlice()
    c1 = tau_pieces(det0, rec["t1"], seed=20 + i) if rec["passes"] else fake_tau_pieces(det0, rec["t1"], seed=20 + i)
    c2 = tau_pieces(det0, rec["t2"], seed=40 + i)
    phi1 = phi_end(c1.trk[len(c1.trk) // 2])
    phi2 = phi_end(c2.trk[len(c2.trk) // 2])
    g = VGroup(c1, c2)
    place(g, A["det_cr"], DET_CENTER)
    g.c1, g.c2, g.record, g.phi1, g.phi2 = c1, c2, rec, phi1, phi2
    return g


def pair_labels(pair) -> VGroup:
    rec = pair.record
    labs = VGroup()
    for who, phi, tight in (("t1", pair.phi1, bool(rec["passes"])), ("t2", pair.phi2, True)):
        sgn = "+" if _q(rec[who]) > 0 else "-"
        lab = mathtex(rf"\tau_h^{{{sgn}}}", color=CHANNEL_LINE["tautau"] if tight else GREY)
        labs.add(lab.move_to(label_spot(phi, A["det_cr"])))
    return labs


def tally_colour(p):
    return CHANNEL["tautau"] if p["passes"] else SAMPLE["Fakes"]


def tally_boxes(filled: bool) -> VGroup:
    g = VGroup()
    for p in SS_PAIRS:
        g.add(Rectangle(width=0.34, height=0.34, stroke_color=col(INK), stroke_width=1.6,
                        fill_color=col(tally_colour(p)), fill_opacity=1.0 if filled else 0.0))
    g.arrange(RIGHT, buff=0.08).move_to(_p3(A["tally_c"]))
    return g


def ff_formula(at=None):
    n_t = sum(bool(p["passes"]) for p in SS_PAIRS)
    n_l = len(SS_PAIRS) - n_t
    return tex_h(rf"f = N_{{T}} / N_{{L}} = {n_t}/{n_l}", 0.22).move_to(_p3(at or A["ff_tex_c"]))


def _pt_bin_labels(edges):
    out = []
    for i in range(len(edges) - 1):
        if edges[i + 1] >= 500:
            out.append(rf">\!{edges[i]:g}")
        else:
            out.append(rf"{edges[i]:g}\!-\!{edges[i + 1]:g}")
    return out


def ff_map() -> VGroup:
    M = FK["ff_map"]
    rows = [DM_TEX[str(d)] for d in M["row_labels_dm"]]
    cols = _pt_bin_labels(M["pt_edges"])
    vg = value_grid(M["values"], rows, cols, fmt="{:.2f}", cell=(0.58, 0.28), label_h=0.125)
    head = tex_h(r"p_T\ [\mathrm{GeV}]", 0.125)
    head.next_to(vg.col_labels, UP, buff=0.06)
    g = VGroup(vg, head)
    g.move_to(_p3(A["ffmap_c"]))
    return g


def ff_dependence():
    return tex_h(r"f = f(\mathrm{period},\ \mathrm{DM},\ N_{\mathrm{jets}},\ p_T)", 0.2).move_to(_p3(A["ffdep_at"]))


def closure_plot(which: str, after: bool) -> VGroup:
    """The same-sign closure in |eta(tau1)| (6 bins) or p_T(tau2) (5 bins, overflow drawn
    80-100): the FF prediction as a pale-slate filled step, the data as slate dots with sqrt(N)
    bars, and a ratio strip underneath (obs / pred). ``after``: the prediction carries the
    correction (pred x before) and the ratio is the frozen closure "after" (= 1)."""
    C = FK["closure"][which]
    e = [float(x) for x in C["edges"]]
    if e[-1] >= 500:
        e[-1] = e[-2] + (e[-2] - e[-3])
        ticks, fmt = e[:-1], "{:.0f}"
    else:
        ticks, fmt = e, "{:.1f}"
    obs = np.asarray(C["obs"], dtype=float)
    pred = np.asarray(C["pred"], dtype=float)
    before = np.asarray(C["before"], dtype=float)
    ratio = np.asarray(C["after"] if after else C["before"], dtype=float)
    pred_drawn = pred * before if after else pred
    top = float(max(obs.max(), pred.max())) * 1.25
    step = 1000.0 if top > 2500 else 500.0
    w, h = A["clos_len"]
    dax = DataAxes([e[0], e[-1], e[1] - e[0]], [0, top, step], w, h, x_ticks=ticks, x_fmt=fmt,
                   y_ticks=list(np.arange(0, top + 1e-9, step)), show_x_labels=False, tick_label_h=0.12)
    dax.move_frame_to(A["clos_main"])
    fill = step_hist(dax, e, pred_drawn, color=SAMPLE["Fakes"], stroke_width=0, fill_opacity=0.95)
    xc = 0.5 * (np.asarray(e[:-1]) + np.asarray(e[1:]))
    dots = VGroup(*[data_dot(dax, xc[i], obs[i], color=SAMPLE["Data"], radius=0.035) for i in range(len(obs))])
    errs = VGroup(*[data_errorbar(dax, xc[i], obs[i], math.sqrt(obs[i]), color=SAMPLE["Data"], stroke_width=1.8)
                    for i in range(len(obs))])
    x_title = r"|\eta(\tau_1)|" if which == "eta" else r"p_T(\tau_2)\ [\mathrm{GeV}]"
    rdax = DataAxes([e[0], e[-1], e[1] - e[0]], [0.8, 1.2, 0.2], w, A["clos_rlen"], x_ticks=ticks, x_fmt=fmt,
                    y_ticks=[0.8, 1.0, 1.2], y_fmt="{:.1f}", tick_label_h=0.12, title_h=0.16, title_buff=0.12,
                    x_title=x_title, y_title=r"\mathrm{obs/pred.}")
    rdax.move_frame_to(A["clos_ratio"])
    rdax.shift(RIGHT * (dax.c2p(e[0], 0)[0] - rdax.c2p(e[0], 0.8)[0]))
    rref = DashedVMobject(Line(rdax.c2p(e[0], 1.0), rdax.c2p(e[-1], 1.0), stroke_color=col(GREY), stroke_width=2.0),
                          num_dashes=30)
    rdots = VGroup(*[data_dot(rdax, xc[i], float(ratio[i]), color=SAMPLE["Data"], radius=0.035)
                     for i in range(len(ratio))])
    sym = r"f(|\eta_{\tau_1}|)" if which == "eta" else r"g(p_T^{\tau_2})"
    lab = tex_h(sym, 0.2, color=darken(SAMPLE["Fakes"], 0.45)).next_to(rdax.frame, RIGHT, buff=0.15)
    g = VGroup(dax, fill, errs, dots, rdax, rref, rdots)
    g.dax, g.fill, g.errs, g.dots, g.rdax, g.rref, g.rdots, g.lab = dax, fill, errs, dots, rdax, rref, rdots, lab
    return g


def osss_parts() -> VGroup:
    lo, hi = FK["osss"]["range_per_category"]
    l1 = tex_h(r"C_{\mathrm{OS/SS}} = C(\mathrm{period},\ N_{\mathrm{jets}},\ D_{\mathrm{BDT}})", 0.2)
    l2 = tex_h(rf"C_{{\mathrm{{OS/SS}}}} \in [{lo:.2f};\ {hi:.3f}]", 0.2)
    l1.move_to(_p3(A["osss_at"][0]))
    l2.move_to(_p3(A["osss_at"][1]))
    return VGroup(l1, l2)


# ---------------------------------------------------------------------------
# corrections, topologies, score plot, category panels, fit
# ---------------------------------------------------------------------------

def corr_rows() -> VGroup:
    tp = CO["taupog"]
    sf = ",\ ".join(f"{tp['id_sf_tight_per_dm'][d]['value']:.3f}" for d in DMS)
    tes = ",\ ".join(f"{tp['tes_per_dm'][d]['value']:.3f}" for d in DMS)
    cst = CO["constants"]
    exprs = [
        (",\ ".join(DM_TEX[d] for d in DMS), GREY),
        (rf"\mathrm{{SF}}_{{\mathrm{{ID}}}} = {sf}", INK),
        (rf"\mathrm{{TES}} = {tes}", INK),
        (rf"\langle \mathrm{{SF}}_{{\mathrm{{ID}}}} \rangle = {MEANS['sf_id_both_legs']:.3f}", INK),
        (rf"\langle \mathrm{{SF}}_{{\mathrm{{trig}}}} \rangle = {MEANS['sf_trig_both_legs']:.3f}", INK),
        (rf"\langle w_{{\mathrm{{PU}}}} \rangle = {MEANS['w_pu']:.3f}", INK),
        (rf"\langle w_{{\mathrm{{L1}}}} \rangle = {MEANS['w_l1']:.3f}", INK),
        (rf"f_{{\tau}}^{{\mathrm{{AR}}}} = {100 * cst['mc_subtraction_ar_fraction']['value']:.1f}\%", INK),
        (rf"f_{{\mathrm{{non\text{{-}}fid}}}} = {100 * cst['nonfiducial_fraction_of_selected_dy']['value']:.0f}\%", INK),
    ]
    g = VGroup()
    for i, (e, cc) in enumerate(exprs):
        lab = tex_h(e, 0.125 if i == 0 else A["corr_h"], color=cc)
        lab.shift(np.array([A["left_col_x"], A["corr_y0"] - A["corr_dy"] * i, 0.0]) - lab.get_left())
        g.add(lab)
    return g


def inputs_list(parked: bool = False) -> VGroup:
    """The 16 BDT input symbols by importance: two columns of eight, or (parked) one narrow
    column at the left edge. Same strings in both, so the two lists Transform 1:1."""
    g = VGroup()
    for i, f in enumerate(BI["features"]):
        if parked:
            P = A["inputs_park"]
            t = tex_h(f["symbol"], P["h"])
            t.shift(np.array([P["x"], P["y0"] - P["dy"] * i, 0.0]) - t.get_left())
        else:
            t = tex_h(f["symbol"], A["inputs_h"])
            t.shift(np.array([A["inputs_cols"][i // 8], A["inputs_y0"] - A["inputs_dy"] * (i % 8), 0.0]) - t.get_left())
        g.add(t)
    return g


def feature_panels() -> dict:
    """2 x 4 small panels: unit-normalised fakes (pale-slate fill) vs signal (red line) for the
    eight leading inputs, plus a two-entry key."""
    FP = A["fpanel"]
    w, h = FP["len"]
    out = {}
    for k, name in enumerate(SHAPE_FEATURES):
        f = FEATURES[name]
        e = np.asarray(f["edges"], dtype=float)
        fk = np.clip(np.asarray(f["fakes"], dtype=float), 0, None)
        sg = np.clip(np.asarray(f["signal_unit"], dtype=float), 0, None)
        top = float(max(fk.max(), sg.max())) * 1.15
        ticks = [float(e[0]), 0.5 * float(e[0] + e[-1]), float(e[-1])]
        dax = DataAxes([e[0], e[-1], e[-1] - e[0]], [0, top, top], w, h, x_ticks=ticks, x_fmt="{:.3g}", y_ticks=[],
                       show_y_labels=False, tick_label_h=0.13, title_h=0.17, title_buff=0.1, x_title=f["symbol"])
        dax.move_frame_to((FP["xs"][k % 4], FP["ys"][k // 4]))
        fill = step_hist(dax, e, fk, color=SAMPLE["Fakes"], stroke_width=0, fill_opacity=0.95)
        line = step_hist(dax, e, sg, color=CHANNEL_LINE["tautau"], stroke_width=2.5)
        out[f"f{k}"] = VGroup(dax, fill, line)
    key = colour_key([(r"\mathrm{Fakes}", SAMPLE["Fakes"]), (r"Z\to\tau\tau", CHANNEL_LINE["tautau"], "line")],
                     label_h=0.18)
    key.rows.arrange(RIGHT, buff=0.4)
    key.move_to(_p3(FP["key_at"]))
    out["fkey"] = key
    return out


def bdt_sketch() -> VGroup:
    """Schematic gradient-boosted classifier: the input bracket, three depth-2 trees (nodes and
    branches in DETECTOR_ACCENT, one highlighted path each in the channel colour), plus signs,
    an arrow into D_BDT."""
    T = A["tree"]
    acc, red = col(DETECTOR_ACCENT), col(CHANNEL_LINE["tautau"])

    def node(p, hl):
        return Circle(radius=T["r"], arc_center=_p3(p), fill_color=lighten(DETECTOR_ACCENT, 0.8), fill_opacity=1.0,
                      stroke_color=red if hl else acc, stroke_width=2.6 if hl else 2.0)

    def edge(p, q, hl):
        return Line(_p3(p), _p3(q), stroke_color=red if hl else acc, stroke_width=3.0 if hl else 2.0)

    trees = VGroup()
    for k, x in enumerate(T["xs"]):
        y_root = T["y"] + T["dy"]
        hl_kid, hl_leaf = (0, 1, 0)[k], (1, 3, 0)[k]
        root = (x, y_root)
        kids = [(x - T["dx"][0], y_root - T["dy"]), (x + T["dx"][0], y_root - T["dy"])]
        leaves = [(kx + s * T["dx"][1], y_root - 2 * T["dy"]) for kx, _ in kids for s in (-1, 1)]
        edges, nodes = VGroup(), VGroup()
        for j, kid in enumerate(kids):
            edges.add(edge(root, kid, j == hl_kid))
        for j, lf in enumerate(leaves):
            edges.add(edge(kids[j // 2], lf, j == hl_leaf))
        nodes.add(node(root, True))
        for j, kid in enumerate(kids):
            nodes.add(node(kid, j == hl_kid))
        for j, lf in enumerate(leaves):
            nodes.add(node(lf, j == hl_leaf))
        trees.add(VGroup(edges, nodes))
    plus = VGroup(*[tex_h("+", 0.3).move_to(np.array([T["xs"][i] + T["plus_dx"], T["y"], 0.0])) for i in range(2)])
    bx, y = T["bracket_x"], T["y"]
    bracket = VGroup(Line([bx, y - 1.5, 0], [bx, y + 1.5, 0], stroke_color=col(INK), stroke_width=2.5),
                     Line([bx, y + 1.5, 0], [bx + 0.22, y + 1.5, 0], stroke_color=col(INK), stroke_width=2.5),
                     Line([bx, y - 1.5, 0], [bx + 0.22, y - 1.5, 0], stroke_color=col(INK), stroke_width=2.5))
    ln_in = Line([bx + 0.35, y, 0], [T["xs"][0] - 1.05, y, 0], stroke_color=col(INK), stroke_width=2.5)
    arrow_in = VGroup(ln_in, arrow_tip_on(ln_in, color=INK, at=1.0, tip_length=0.2))
    ln_out = Line([T["xs"][2] + 1.05, y, 0], [T["out_at"][0] - 0.7, y, 0], stroke_color=col(INK), stroke_width=2.5)
    arrow_out = VGroup(ln_out, arrow_tip_on(ln_out, color=INK, at=1.0, tip_length=0.2))
    out = tex_h(r"D_{\mathrm{BDT}}", 0.28).move_to(_p3(T["out_at"]))
    g = VGroup(bracket, arrow_in, trees, plus, arrow_out, out)
    g.bracket, g.arrow_in, g.trees, g.plus, g.arrow_out, g.out = bracket, arrow_in, trees, plus, arrow_out, out
    return g


def score_parts() -> dict:
    S = BD["sr_score"]
    e = np.asarray(S["edges"], dtype=float)
    y0, y1 = A["score_y"]
    dax = DataAxes([0, 1, 0.25], [y0, y1, 1], 6.4, 2.9, x_fmt="{:.2f}", y_ticks=[k for k in range(1, 5) if y0 <= k <= y1],
                   y_log=True, show_x_labels=False, title_h=0.2, title_buff=0.16, y_title=r"\mathrm{events}")
    dax.move_frame_to(A["plot_main_c"])
    G = S["groups"]
    layers = [("Fakes", np.clip(np.asarray(S["fakes"], dtype=float), 0, None), LAYER_COL["Fakes"]),
              ("rest", np.clip(np.asarray(G["rest"], dtype=float), 0, None), LAYER_COL["rest"]),
              ("DYtautau_nonfid", np.clip(np.asarray(G["DYtautau_nonfid"], dtype=float), 0, None),
               LAYER_COL["DYtautau_nonfid"]),
              ("DYtautau", np.clip(np.asarray(G["DYtautau"], dtype=float), 0, None), LAYER_COL["DYtautau"])]
    stack = stack_hist(dax, e, layers)
    xc = 0.5 * (e[:-1] + e[1:])
    dots = VGroup(*[data_dot(dax, xc[i], float(S["data"][i]), color=SAMPLE["Data"], radius=DOT_R)
                    for i in range(len(xc))])
    cuts = VGroup(*[vref_line(dax, float(c), color=INK, stroke_width=2.0, num_dashes=24)
                    for c in BD["category_edges"][1:-1]])
    auc = tex_h(rf"\mathrm{{AUC}} = {BD['auc_test_mean']:.3f}", 0.2).move_to(dax.c2p(0.5, 10 ** (y1 - 0.3)))
    ratio = ratio_panel(dax, e, np.asarray(S["data"], dtype=float), np.asarray(S["pred"], dtype=float), A["ratio_c"],
                        y_range=(0.5, 1.5, 0.5), y_length=0.9, y_ticks=[0.5, 1.0, 1.5], dot_radius=RDOT_R,
                        x_title=r"D_{\mathrm{BDT}}", y_title=r"\mathrm{data/pred.}", title_h=0.2, title_buff=0.16)
    return {"s_dax": dax, "s_stack": stack, "s_data": dots, "s_ratio": ratio, "s_cuts": cuts, "s_auc": auc}


SCORE_KEYS = ("s_dax", "s_stack", "s_data", "s_ratio", "s_cuts", "s_auc")


def region_layers(fit: str, r: str):
    S = FIT["regions"][fit][r]["samples"]
    rest = np.clip(np.sum([np.asarray(S[n], dtype=float) for n in REST_NAMES], axis=0), 0, None)
    return [("Fakes", np.asarray(S["Fakes"], dtype=float), LAYER_COL["Fakes"]), ("rest", rest, LAYER_COL["rest"]),
            ("DYtautau_nonfid", np.asarray(S["DYtautau_nonfid"], dtype=float), LAYER_COL["DYtautau_nonfid"]),
            ("DYtautau", np.asarray(S["DYtautau"], dtype=float), LAYER_COL["DYtautau"])]


def panel_parts(fit: str = "prefit", placement=None) -> dict:
    """Three narrow category panels at the main plot position: dax, stack, data, (SR0) the
    greyed dropped bins, the category range and the fiducial / fake yields above."""
    placement = placement or A["plot_main"]
    P = {}
    ranges = (r"D_{\mathrm{BDT}} < 0.55", r"0.55 < D_{\mathrm{BDT}} < 0.90", r"D_{\mathrm{BDT}} > 0.90")
    for k, r in enumerate(REGIONS):
        ytop, yt = P_YTOP[r]
        dax = DataAxes([0, 350, 100], [0, ytop, yt[1] - yt[0]], A["panel_w"], A["panel_h"], y_ticks=yt,
                       tick_label_h=0.14, title_h=0.18, title_buff=0.12,
                       x_title=r"m_{\tau\tau}\ [\mathrm{GeV}]" if k == 1 else None)
        dax.move_frame_to((A["panel_x"][k], A["panel_cy"]))
        P[f"p{k}_dax"] = dax
        P[f"p{k}_stack"] = stack_hist(dax, EDGES, region_layers(fit, r))
        D = FIT["regions"][fit][r]["data"]
        P[f"p{k}_data"] = VGroup(*[data_dot(dax, XC[i], float(D[i]), color=SAMPLE["Data"], radius=0.03)
                                   for i in range(NB)])
        if k == 0:
            m_max = float(FIT["dropped_bins"][r]["m_tt_max"])
            lo, hi = dax.c2p(0.0, 0.0), dax.c2p(m_max, ytop)
            P["p0_drop"] = Rectangle(width=hi[0] - lo[0], height=hi[1] - lo[1], fill_color=col(LIGHT_GREY),
                                     fill_opacity=0.6, stroke_width=0).move_to(0.5 * (lo + hi))
        c = BD["categories"][r]
        rng_lab = tex_h(ranges[k], 0.15).move_to(np.array([A["panel_x"][k], A["panel_lab_y"][0], 0.0]))
        yl = VGroup(tex_h(f"{c['DYtautau']:.0f}", 0.16, color=CHANNEL_LINE["tautau"]), tex_h("/", 0.16),
                    tex_h(f"{c['fakes']:.0f}", 0.16, color=darken(SAMPLE["Fakes"], 0.45))).arrange(RIGHT, buff=0.07)
        yl.move_to(np.array([A["panel_x"][k], A["panel_lab_y"][1], 0.0]))
        P[f"p{k}_lab"] = VGroup(rng_lab, yl)
    for m in P.values():
        place(m, placement, PLOT_C)
    return P


PANEL_KEYS = ("p0_dax", "p0_stack", "p0_data", "p0_drop", "p0_lab", "p1_dax", "p1_stack", "p1_data", "p1_lab",
              "p2_dax", "p2_stack", "p2_data", "p2_lab")

NP = {n["name"]: n for n in FIT["nps"]}


def pulls(post: bool) -> VGroup:
    names = list(FIT["nps_shown"])
    p = [NP[n]["pull"] for n in names] if post else [0.0] * len(names)
    c = [NP[n]["constraint"] for n in names] if post else [1.0] * len(names)
    pp = pull_plot([NP_TEX[n] for n in names], p, c, x_length=1.6, row_h=0.34, label_h=0.12, tex=True)
    pp.move_to(_p3(A["pulls_at"]))
    pp.shift(RIGHT * (A["pulls_at"][0] - pp.baseline.get_center()[0]))   # centre the axis, not the names
    return pp


def mu_strings():
    poi = FIT["poi"]
    return f"{poi['value']:.3f}", rf"^{{+{poi['err_up']:.3f}}}_{{-{poi['err_down']:.3f}}}"


def mu_slider() -> VGroup:
    sl = slider(0.8, 1.3, 1.0, err=None, ref=1.0, length=2.6, ticks=[0.8, 1.0, 1.2], fmt="{:.1f}", tick_label_h=0.16)
    sl.move_to(_p3(A["slider_at"]))
    return sl


def asym_marker(sl) -> VGroup:
    poi = FIT["poi"]
    y = sl.axis.get_center()[1]
    cc = col(sl.color)
    bar = Line([sl.x_of(poi["value"] - poi["err_down"]), y, 0], [sl.x_of(poi["value"] + poi["err_up"]), y, 0],
               stroke_color=cc, stroke_width=sl.bar_sw)
    dot = Dot([sl.x_of(poi["value"]), y, 0], radius=sl.marker_r, color=cc)
    return VGroup(bar, dot)


def mu_line(h=0.34, at=None):
    v, e = mu_strings()
    return tex_h(rf"\mu_{{Z}} = {v} {e}", h).move_to(_p3(at or A["mu_line_at"]))


def impact_bars() -> VGroup:
    rows = [(n, float(FIT["grouped_impact"][n]), SLATE) for n in FIT["impact_order"]]
    rows.append(("stat", float(FIT["stat_impact"]), CHANNEL["tautau"]))
    g = VGroup()
    g.bars, g.nums, g.names = VGroup(), VGroup(), VGroup()
    for i, (name, v, cc) in enumerate(rows):
        y = A["imp_y0"] - A["imp_dy"] * i
        lab = text_h(IMPACT_NAME.get(name, name), 0.13)
        lab.shift(np.array([A["imp_name_r"], y, 0.0]) - lab.get_right())
        w = 100.0 * v * A["imp_per_pct"]
        bar = Rectangle(width=w, height=0.18, fill_color=col(cc), fill_opacity=0.85, stroke_width=0)
        bar.shift(np.array([A["imp_bar_x0"], y, 0.0]) - bar.get_left())
        num = text_h(f"{100 * v:.1f}", 0.13).next_to(bar, RIGHT, buff=0.08)
        g.names.add(lab); g.bars.add(bar); g.nums.add(num)
    head = tex_h(r"\Delta\mu_Z / \mu_Z\ [\%]", 0.16).move_to(_p3(A["imp_head"]))
    g.head = head
    g.add(head, g.names, g.bars, g.nums)
    return g


def sigma_fid_line():
    s = FIT["sigma_fid"]
    expr = (rf"\sigma_{{\mathrm{{fid}}}} = {s['value']:.2f} \pm {s['stat']:.2f}_{{\mathrm{{stat}}}}"
            rf" \pm {s['syst']:.2f}_{{\mathrm{{syst}}}}\ \mathrm{{pb}}")
    return tex_h(expr, 0.27).move_to(_p3(A["sig_fid_at"]))


def sigma_panel() -> VGroup:
    s = FIT["sigma_60_120"]
    S = A["sig_ax"]
    lo, hi, x0, L, y0 = S["lo"], S["hi"], S["x0"], S["length"], S["y"]
    ink = col(INK)

    def xs(v):
        return x0 + (float(v) - lo) / (hi - lo) * L

    axis = Line([xs(lo), y0, 0], [xs(hi), y0, 0], stroke_color=ink, stroke_width=2.5)
    ticks, tick_labels = VGroup(), VGroup()
    for v in np.arange(lo, hi + 1, 100.0):
        p = np.array([xs(v), y0, 0.0])
        ticks.add(Line(p, p + DOWN * 0.08, stroke_color=ink, stroke_width=2.5))
        tick_labels.add(text_h(f"{v:.0f}", 0.15).next_to(p + DOWN * 0.08, DOWN, buff=0.08))
    unit_ = tex_h(r"\mathrm{pb}", 0.17).next_to(tick_labels[-1], RIGHT, buff=0.18)
    th = RF["theory"]
    pred = float(th["value"])
    H = S["height"]
    band = None
    if th.get("unc_up") is not None and th.get("unc_down") is not None:      # only a documented one
        lo_b, hi_b = xs(pred - float(th["unc_down"])), xs(pred + float(th["unc_up"]))
        band = Rectangle(width=hi_b - lo_b, height=H, fill_color=col(THEORY), fill_opacity=0.15,
                         stroke_width=0).move_to([0.5 * (lo_b + hi_b), y0 + H / 2, 0.0])
    theory = DashedLine([xs(pred), y0, 0], [xs(pred), y0 + H, 0], dash_length=0.1,
                        stroke_color=col(THEORY), stroke_width=3.0)
    if band is not None:                                                     # the value with its uncertainty
        th_lab = tex_h(rf"{pred:.1f}^{{+{float(th['unc_up']):.0f}}}_{{-{float(th['unc_down']):.0f}}}", 0.17,
                       color=THEORY)
    else:
        th_lab = tex_h(f"{pred:.1f}", 0.17, color=THEORY)
    th_lab.next_to(theory.get_end(), LEFT, buff=0.14)                        # the result label sits above
    cc = col(CHANNEL_LINE["tautau"])
    ypt = y0 + S["rows"][0]
    bar = Line([xs(s["value"] - s["err_down"]), ypt, 0], [xs(s["value"] + s["err_up"]), ypt, 0],
               stroke_color=cc, stroke_width=5.0)
    dot = Dot([xs(s["value"]), ypt, 0], radius=0.1, color=cc)
    refs = VGroup()
    for key, dy in (("cms", S["rows"][1]), ("atlas", S["rows"][2])):
        r = RF[key]
        yr = y0 + dy
        sc = col(SAMPLE["Data"])
        rbar = Line([xs(r["value"] - r["total"]), yr, 0], [xs(r["value"] + r["total"]), yr, 0],
                    stroke_color=sc, stroke_width=4.0)
        rdot = Dot([xs(r["value"]), yr, 0], radius=0.08, color=sc)
        rlab = tex_h(rf"\mathrm{{{r['label']}}}", 0.16, color=SAMPLE["Data"]).next_to(rbar, RIGHT, buff=0.15)
        refs.add(VGroup(rbar, rdot, rlab))
    label = tex_h(rf"\sigma_{{60\text{{--}}120}} = {s['value']:.0f}^{{+{s['err_up']:.0f}}}_{{-{s['err_down']:.0f}}}"
                  rf"\ \mathrm{{pb}}", 0.27)
    label.move_to(_p3(A["sig_label_at"]))
    members = [axis, ticks, tick_labels, unit_] + ([band] if band is not None else []) + \
        [theory, th_lab, bar, dot, refs, label]
    g = VGroup(*members)
    g.axis, g.theory, g.th_lab, g.bar, g.dot, g.refs, g.label, g.band = axis, theory, th_lab, bar, dot, refs, label, band
    g.frame = VGroup(axis, ticks, tick_labels, unit_)
    return g


# ---------------------------------------------------------------------------
# end states (pure) and their z-orders
# ---------------------------------------------------------------------------

def detector_end_state() -> dict:
    """Exact replica of 5-03's last frame (s5_ztautau.ZtautauDetector): the slice, a 1-prong
    tau_h at 40 deg (charge -1, pi0, seed 2) and a 3-prong at 218 deg (charge +1, seed 4)
    with tracker hits and deposits, the two nu lines at +-6 deg, the p_T^miss arrow on their
    bisector; the outer labels were faded by zoom_in. Camera: scale 0.5 on det.c."""
    det = CMSSlice()
    a1, a2 = A["phi_503"]
    d_nu = A["d_nu_503"]
    st = {"det": det}
    for sfx, phi, q, prongs, pi0, seed, dn in (("1", a1, -1, 1, True, 2, +d_nu), ("2", a2, +1, 3, False, 4, -d_nu)):
        sig = signature(det, "tau_h", phi, charge=q, prongs=prongs, pi0=pi0, seed=seed)
        trk, rest = sig[0], VGroup(*sig[1:])
        hits = tracker_hits(det, trk[len(trk) // 2], color=CHANNEL["tautau"])
        st.update({f"trk_{sfx}": trk, f"hits_{sfx}": hits, f"dep_{sfx}": rest,
                   f"nu_{sfx}": signature(det, "nu", phi + dn)})
    phi_met = 0.5 * ((a1 + d_nu) + (a2 - d_nu))
    st["met"] = signature(det, "met", phi_met, r_end=det.radii["hcal"][0] + 0.1)
    return st


ORDER_503 = ("det", "trk_1", "hits_1", "dep_1", "nu_1", "trk_2", "hits_2", "dep_2", "nu_2", "met")


def state_a() -> dict:
    t1, t2 = event_taus()
    labs = pt_labels()
    return {"det": CMSSlice(), "t1": t1, "t2": t2, "met": event_met(), "pt_met": labs["pt_met"],
            "stamp": run_stamp(), "mvis": mass_label("vis"), "mscale": mass_scale(float(EVENT["m_vis"]))}


ORDER_A = ("det", "t1", "t2", "met", "pt_met", "stamp", "mvis", "mscale")


def state_b() -> dict:
    t1, t2 = event_taus()
    st = {"det": CMSSlice(), "t1": t1, "t2": t2, "met": event_met(), "nus": event_nus(), "ell": met_ellipse(),
          "stamp": run_stamp(), "mtt": mass_label("tt"), "grid": heat_grid(parked=True)}
    st.update(shapes_parts())
    return st


ORDER_B = ("det", "t1", "t2", "met", "stamp", "nus", "ell", "grid", "mtt", *SHAPE_KEYS)


def state_c() -> dict:
    P = plot_parts(mc=False, data="bars", ratio=False, key_rows=0, x_labels=True)
    return {"det": parked_slice(), **P}


ORDER_C = ("det", "dax", "stack", "data", "counter")


def state_d() -> dict:
    return {"det": parked_slice(), **plot_parts(fakes=False, key_rows=4)}


ORDER_D = ("det", *PLOT_KEYS)


def state_e() -> dict:
    st = {"det": parked_slice(A["det_cr"]), **plot_parts(A["plot_ghost"], fakes=False, key_rows=4)}
    st.update({"box": cr_box(r"\tau_h^{+}\tau_h^{-}"), "ff_tex": ff_formula(A["ff_tex_end"]), "osss": osss_parts()})
    return st


ORDER_E = ("det", *PLOT_KEYS, "box", "ff_tex", "osss")


def state_f() -> dict:
    P = plot_parts(fakes=True, key_rows=5)
    return {**P, "n_fakes": n_fakes_big(), "pct": pct_big()}


ORDER_F = ("dax", "stack", "data", "counter", "key", "n_fakes", "ratio", "rlabel", "pct")


def state_g() -> dict:
    P = plot_parts(fakes=True, key_rows=5)
    return {**P, "n_fakes": n_fakes_small(), "corr": corr_rows()}


ORDER_G = ("dax", "stack", "data", "counter", "key", "n_fakes", "ratio", "rlabel", "corr")


def state_h() -> dict:
    return {**panel_parts("prefit"), "key": make_key(5)}


ORDER_H = (*PANEL_KEYS, "key")


def state_i() -> dict:
    st = {**panel_parts("postfit", A["plot_left"]), "key": key_i(), "mu_line": mu_line(),
          "sig_fid": sigma_fid_line(), "sig_tot": sigma_panel()}
    return st


ORDER_I = (*PANEL_KEYS, "key", "mu_line", "sig_fid", "sig_tot")


# ---------------------------------------------------------------------------
# the nine scenes
# ---------------------------------------------------------------------------

class TautauEvent(MovingCameraScene):
    """tautau_a1_event / tautau_a2_visible_mass: the camera zooms out of 5-03, the schematic
    tracks fade, the real SR event is drawn, its p_T values converge to m_vis."""

    def construct(self):
        white_background(self)
        prev = detector_end_state()
        add_state(self, prev, ORDER_503)
        det = prev["det"]
        self.camera.frame.scale(A["zoom_503"][0]).move_to(_p3(A["zoom_503"][1]))
        clip_open(self, "tautau_a1_event")
        end = state_a()

        self.play(self.camera.frame.animate.scale(1.0 / A["zoom_503"][0]).move_to(ORIGIN), run_time=1.4,
                  rate_func=EASE)
        if self.camera.frame in self.mobjects:
            self.remove(self.camera.frame)
        self.play(*[FadeOut(prev[k]) for k in ORDER_503[1:]], run_time=0.5)
        labs = pt_labels()
        for key, tau in (("t1", end["t1"]), ("t2", end["t2"])):
            self.play(Create(tau.trk, lag_ratio=0.0), run_time=0.9, rate_func=EASE)
            self.play(FadeIn(tau.hits), run_time=0.3)
            self.play(FadeIn(tau.deposits), run_time=0.45)
            adopt(self, tau, key)
            self.play(FadeIn(labs["pt1" if key == "t1" else "pt2"], shift=RIGHT * 0.25), run_time=0.4)
        met = end["met"]
        self.play(Create(met[0]), run_time=0.7, rate_func=EASE)
        self.play(FadeIn(met[1]), run_time=0.2)
        adopt(self, met, "met")
        self.play(FadeIn(end["pt_met"], shift=RIGHT * 0.25), run_time=0.4)
        self.play(FadeIn(end["stamp"]), run_time=0.4)
        clip_cut(self, "tautau_a2_visible_mass")

        pts = VGroup(labs["pt1"], labs["pt2"])
        self.remove(labs["pt1"], labs["pt2"])
        self.add(pts)
        self.play(ReplacementTransform(pts, end["mvis"]), run_time=0.9, rate_func=EASE)
        ms = end["mscale"]
        self.play(FadeIn(ms.frame), run_time=0.5)
        self.play(GrowFromCenter(ms.dot), run_time=0.4)
        adopt(self, ms, "mscale")
        settle(self, det, end["det"], "det")
        check_order(self, end, ORDER_A)
        self.wait(0.1)


class TautauMass(Scene):
    """tautau_b1_neutrinos / b2_likelihood / b3_shapes: the collinear neutrinos and the
    p_T^miss, the (x1, x2) posterior grid and the formula, the simulated shapes."""

    def construct(self):
        white_background(self)
        prev = state_a()
        add_state(self, prev, ORDER_A)
        clip_open(self, "tautau_b1_neutrinos")
        end = state_b()
        det = prev["det"]

        nus, ell = end["nus"], end["ell"]
        self.play(Create(nus[0][0]), run_time=0.7, rate_func=EASE)
        self.play(FadeIn(nus[0][1]), run_time=0.15)
        self.play(Create(nus[1][0]), run_time=0.5, rate_func=EASE)
        self.play(FadeIn(nus[1][1]), run_time=0.15)
        adopt(self, nus, "nus")
        self.play(FadeIn(ell), run_time=0.5)
        tip = det.c + det.radii["hcal"][0] * unit(float(EVENT["met"]["phi"]))
        self.play(Flash(tip, color=col(SLATE), line_length=0.15, flash_radius=0.3, run_time=0.6))
        clip_cut(self, "tautau_b2_likelihood")

        grid = heat_grid()
        glab = grid_labels(grid)
        self.play(FadeOut(prev["pt_met"]), run_time=0.3)
        self.play(LaggedStart(*[FadeIn(c) for c in grid], lag_ratio=0.012, group=grid), run_time=2.6)
        self.play(FadeIn(glab), run_time=0.3)
        formula = mass_formula()
        self.play(FadeIn(formula, shift=UP * 0.15), run_time=0.5)
        cell = grid[median_cell_index()]
        self.play(Flash(cell.get_center(), color=col(SLATE), line_length=0.15, flash_radius=0.3, run_time=0.6))
        ms = prev["mscale"]
        ms_new = mass_scale(float(EVENT["m_tt"]))
        self.play(FadeOut(prev["mvis"]), FadeIn(end["mtt"]), Transform(ms.dot, ms_new.dot), run_time=0.8,
                  rate_func=EASE)
        settle(self, ms, ms_new, "mscale")
        clip_cut(self, "tautau_b3_shapes")

        self.play(FadeOut(glab), FadeOut(formula), FadeOut(ms_new),
                  grid.animate.scale(A["grid_park"][0]).move_to(_p3(A["grid_park"][1])), run_time=0.9, rate_func=EASE)
        settle(self, grid, end["grid"], "grid")
        self.play(FadeIn(end["shp_dax"]), run_time=0.5)
        self.play(Create(end["shp_vis"]), run_time=0.9, rate_func=rate_functions.linear)
        self.play(Create(end["shp_tt"]), run_time=0.9, rate_func=rate_functions.linear)
        self.play(Create(end["shp_mz"]), FadeIn(end["shp_med"]), FadeIn(end["shp_lab"]), run_time=0.5)
        self.play(LaggedStart(*[FadeIn(p, shift=RIGHT * 0.2) for p in end["perf"]], lag_ratio=0.4,
                              group=end["perf"]), run_time=0.9)
        settle_same(self, prev, end, ("det", "t1", "t2", "met", "stamp"))
        check_order(self, end, ORDER_B)
        self.wait(0.1)


class TautauRain(Scene):
    """tautau_c1_first_entry / c2_rain: the slice parks, the linear m_tautau axes appear,
    the event's m_tautau is the first entry; then the clock, the rain, N = 21,160."""

    def construct(self):
        white_background(self)
        prev = state_b()
        add_state(self, prev, ORDER_B)
        clip_open(self, "tautau_c1_first_entry")
        end = state_c()
        P = plot_parts(mc=False, data="entries", entries=(EVENT["m_tt"],), counter=False, ratio=False, key_rows=0,
                       x_labels=True)

        self.play(*[FadeOut(prev[k]) for k in (*SHAPE_KEYS, "grid")], run_time=0.5)
        park = {"det": parked_slice(), "t1": None, "t2": None, "met": place(event_met(), A["det_park"], DET_CENTER),
                "nus": place(event_nus(), A["det_park"], DET_CENTER), "ell": place(met_ellipse(), A["det_park"], DET_CENTER)}
        park["t1"], park["t2"] = event_taus()
        place(park["t1"], A["det_park"], DET_CENTER)
        place(park["t2"], A["det_park"], DET_CENTER)
        m_spot = P["dax"].c2p(float(EVENT["m_tt"]), 0.9 * A["y_top"])
        self.play(*[ReplacementTransform(prev[k], park[k]) for k in park],
                  prev["mtt"].animate.move_to(m_spot), FadeIn(P["dax"]), run_time=1.8, rate_func=EASE)
        self.add(P["stack"])
        self.bring_to_front(prev["mtt"])
        self.wait(0.2)
        self.play(ReplacementTransform(prev["mtt"], P["data"]), run_time=0.8, rate_func=EASE)
        clip_cut(self, "tautau_c2_rain")

        det, dax, entries = park["det"], P["dax"], P["data"]
        bars = end["data"].copy()
        bars0 = data_bars(dax, visible=False)
        self.remove(entries)
        self.add(bars0, entries)
        clk = clock(A["clock_at"])
        counter = make_counter(len(entries))
        self.play(FadeIn(clk), FadeIn(counter), FadeOut(prev["stamp"]),
                  *[FadeOut(park[k]) for k in ("t1", "t2", "met", "nus", "ell")], run_time=0.45)
        n = ValueTracker(float(len(entries)))
        counter.num.add_updater(lambda m: (m.set_value(float(round(n.get_value()))),
                                           layout_counter(counter.pre, m)))
        order = np.random.default_rng(20260916).permutation(NB)
        drops = hidden_rain(det, dax, EDGES, DATA, n=70, seed=3, color=SAMPLE["Data"])
        self.play(LaggedStart(*[Transform(bars0[i], bars[i]) for i in order], lag_ratio=0.06, group=bars0),
                  LaggedStart(*drops, lag_ratio=0.05),
                  clk.turns.animate(rate_func=rate_functions.ease_in_quad).set_value(3.0),
                  n.animate(rate_func=rate_functions.ease_in_quad).set_value(float(TOTAL)),
                  FadeOut(entries), run_time=3.8)
        counter.num.clear_updaters()
        set_counter(counter, TOTAL)
        self.remove(n, clk.turns)
        self.play(FadeOut(clk), run_time=0.35)
        settle(self, bars0, end["data"], "bars")
        settle(self, counter, end["counter"], "counter")
        settle(self, det, end["det"], "det")
        settle(self, dax, end["dax"], "dax")
        settle(self, P["stack"], end["stack"], "stack")
        check_order(self, end, ORDER_C)
        self.wait(0.1)


class TautauStack(Scene):
    """tautau_d1_simulation / d2_gap: bars become points, the simulation stack slides in
    and explains a third; the ratio panel opens at data / simulation = 2.97."""

    def construct(self):
        white_background(self)
        prev = state_c()
        add_state(self, prev, ORDER_C)
        clip_open(self, "tautau_d1_simulation")
        end = state_d()

        self.play(ReplacementTransform(prev["data"], end["data"]), run_time=1.0, rate_func=EASE)
        self.wait(0.2)
        stack, dax = prev["stack"], prev["dax"]
        base = sr_layers(False)
        for i, poly in enumerate(stack):
            rows = [(n, c if j < i else np.zeros(NB), cc) for j, (n, c, cc) in enumerate(base)]
            poly.become(stack_hist(dax, EDGES, rows)[i])
        self.play(LaggedStart(*[Transform(stack[i], end["stack"][i]) for i in range(len(stack))],
                              lag_ratio=0.15, group=stack), run_time=2.4, rate_func=EASE)
        settle(self, stack, end["stack"], "stack")
        self.play(FadeIn(end["key"]), run_time=0.6)
        clip_cut(self, "tautau_d2_gap")
        ratio = end["ratio"]
        self.play(FadeIn(ratio.dax), FadeIn(ratio.ref), dax.x_labels.animate.set_opacity(0.0),
                  dax.x_title.animate.set_opacity(0.0), run_time=0.6)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in ratio.dots], lag_ratio=0.05, group=ratio.dots),
                  run_time=1.2)
        adopt(self, ratio, "ratio")
        self.wait(0.2)
        self.play(FadeIn(end["rlabel"], shift=LEFT * 0.2), run_time=0.5)
        self.bring_to_front(end["key"])
        settle_same(self, prev, end, ("det", "dax", "counter"))
        check_order(self, end, ORDER_D)
        self.wait(0.1)


class TautauFakes(Scene):
    """tautau_e1_same_sign ... e5_osss: the same-sign control region, a real SS pair, the
    ten-pair tally, the fake-factor map, the closure corrections, C_OS/SS."""

    def construct(self):
        white_background(self)
        prev = state_d()
        add_state(self, prev, ORDER_D)
        clip_open(self, "tautau_e1_same_sign")
        end = state_e()
        box0 = cr_box()

        self.play(*[ReplacementTransform(prev[k], end[k]) for k in ("det", *PLOT_KEYS)], run_time=1.6, rate_func=EASE)
        self.wait(0.2)
        self.play(Create(box0.rect), FadeIn(box0.sym), run_time=0.8)
        adopt(self, box0, "box")
        pair = ss_pair(0)
        c1, c2 = pair.c1, pair.c2
        self.play(Create(c2.trk, lag_ratio=0.0), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(c2.hits), run_time=0.3)
        self.play(FadeIn(c2.deposits), run_time=0.4)
        adopt(self, c2, "tau2")
        self.wait(0.15)
        self.play(Create(c1.jet[0], lag_ratio=0.0), run_time=0.7, rate_func=EASE)
        self.play(FadeIn(VGroup(c1.jet[1], c1.jet[2])), run_time=0.35)
        self.play(Create(c1.trk, lag_ratio=0.0), run_time=0.8, rate_func=EASE)
        self.play(FadeIn(c1.dep), run_time=0.35)
        adopt(self, c1, "tau1")
        labs = pair_labels(pair)
        self.play(FadeIn(labs), run_time=0.4)
        clip_cut(self, "tautau_e2_tally")

        tally = tally_boxes(False)
        self.play(FadeIn(tally), run_time=0.4)
        self.play(tally[0].animate.set_fill(col(tally_colour(SS_PAIRS[0])), opacity=1.0), run_time=0.3)
        self.play(FadeOut(c1), FadeOut(c2), FadeOut(labs), run_time=0.25)
        for i in range(1, len(SS_PAIRS)):
            p = ss_pair(i)
            lines = VGroup(p.c1.lines, p.c2.lines)
            rest = VGroup(p.c1.rest, p.c2.rest)
            self.play(Create(lines, lag_ratio=0.0), FadeIn(rest), run_time=0.45, rate_func=EASE)
            self.play(tally[i].animate.set_fill(col(tally_colour(SS_PAIRS[i])), opacity=1.0), run_time=0.22)
            self.play(FadeOut(lines), FadeOut(rest), run_time=0.22)
        tally_f = tally_boxes(True)
        settle(self, tally, tally_f, "tally")
        tally = tally_f
        f_tex = ff_formula()
        self.play(FadeIn(f_tex, shift=UP * 0.15), run_time=0.8)
        clip_cut(self, "tautau_e3_ff_map")

        fmap, fdep = ff_map(), ff_dependence()
        self.play(FadeIn(fmap), run_time=1.0)
        self.play(FadeIn(fdep, shift=UP * 0.15), run_time=0.6)
        clip_cut(self, "tautau_e4_eta_closure")

        self.play(FadeOut(fmap), FadeOut(fdep), FadeOut(tally), f_tex.animate.move_to(_p3(A["ff_tex_end"])),
                  run_time=0.5, rate_func=EASE)

        def parts(c):
            return [c.dax, c.fill, c.errs, c.dots, c.rdax, c.rref, c.rdots]

        last = None
        for which in ("eta", "pt2"):
            c0, c1 = closure_plot(which, False), closure_plot(which, True)
            self.play(FadeIn(c0.dax), FadeIn(c0.fill), run_time=0.6)
            self.play(FadeIn(c0.errs), FadeIn(c0.dots), run_time=0.4)
            self.play(FadeIn(c0.rdax), FadeIn(c0.rref), run_time=0.4)
            self.play(LaggedStart(*[GrowFromCenter(d) for d in c0.rdots], lag_ratio=0.1, group=c0.rdots), run_time=0.6)
            self.wait(0.5)
            self.play(FadeIn(c0.lab, shift=LEFT * 0.15), run_time=0.4)
            self.play(Transform(c0.rdots, c1.rdots), Transform(c0.fill, c1.fill), run_time=1.0, rate_func=EASE)
            self.wait(0.3)
            if which == "eta":
                self.play(*[FadeOut(m) for m in parts(c0) + [c0.lab]], run_time=0.4)
            last = c0
        clip_cut(self, "tautau_e5_osss")

        self.play(*[FadeOut(m) for m in parts(last) + [last.lab]], run_time=0.4)
        self.play(FadeIn(end["osss"][0], shift=UP * 0.15), run_time=0.7)
        self.play(FadeIn(end["osss"][1], shift=UP * 0.15), run_time=0.6)
        adopt(self, end["osss"], "osss")
        self.play(ReplacementTransform(box0, end["box"]), run_time=0.7, rate_func=EASE)
        settle(self, f_tex, end["ff_tex"], "ff_tex")
        check_order(self, end, ORDER_E)
        self.wait(0.1)


class TautauTransfer(Scene):
    """tautau_f1_apply / f2_template: N_fake = C_OS/SS f x N_AR; the plot returns; the fake
    template enters the bottom of the stack, the ratio drops to 1.02; N_fake = 13592 (64 %)."""

    def construct(self):
        white_background(self)
        prev = state_e()
        add_state(self, prev, ORDER_E)
        clip_open(self, "tautau_f1_apply")
        end = state_f()
        mid = plot_parts(fakes=False, key_rows=4)

        f_tex = prev["ff_tex"]
        self.play(FadeOut(prev["det"]), FadeOut(prev["box"]), FadeOut(prev["osss"]),
                  f_tex.animate.move_to(_p3(A["formula_f"])), run_time=0.8, rate_func=EASE)
        self.play(*[ReplacementTransform(prev[k], mid[k]) for k in PLOT_KEYS], run_time=1.5, rate_func=EASE)
        formula = tex_h(r"N_{\mathrm{fake}} = C_{OS/SS}\, f \times N_{\mathrm{AR}}", 0.24).move_to(_p3(A["formula_f"]))
        self.play(FadeOut(f_tex), FadeIn(formula), run_time=0.6)
        clip_cut(self, "tautau_f2_template")

        tmpl = step_hist(mid["dax"], EDGES, FAKES, color=darken(SAMPLE["Fakes"], 0.55), stroke_width=4.0)
        self.play(Create(tmpl), run_time=1.1, rate_func=rate_functions.linear)
        self.wait(0.2)
        wide = plot_parts(fakes=True, wide=True, key_rows=5)
        self.play(Transform(mid["stack"], wide["stack"]), Transform(mid["ratio"].dots, wide["ratio"].dots),
                  Transform(mid["rlabel"], wide["rlabel"]), ReplacementTransform(mid["key"], end["key"]),
                  FadeOut(tmpl), run_time=1.5, rate_func=EASE)
        self.play(FadeOut(mid["ratio"]), FadeOut(mid["rlabel"]), FadeIn(end["ratio"]), FadeIn(end["rlabel"]),
                  run_time=0.7)
        self.play(ReplacementTransform(formula, end["n_fakes"]), FadeIn(end["pct"], shift=UP * 0.15),
                  run_time=0.9, rate_func=EASE)
        settle_same(self, mid, end, ("dax", "stack", "data", "counter"))
        check_order(self, end, ORDER_F)
        self.wait(0.1)


class TautauCorrections(Scene):
    """tautau_g1_tau_sf / g2_event_weights: the fiducial layer inflates to its raw yield and the
    tau_h scale factors, then the event weights, bring it to the corrected height."""

    def construct(self):
        white_background(self)
        prev = state_f()
        add_state(self, prev, ORDER_F)
        clip_open(self, "tautau_g1_tau_sf")
        end = state_g()
        corr = end["corr"]
        stack, ratio, rlabel = prev["stack"], prev["ratio"], prev["rlabel"]
        self.play(ReplacementTransform(prev["n_fakes"], end["n_fakes"]), FadeOut(prev["pct"]), run_time=0.9,
                  rate_func=EASE)

        def to_stage(k, rt=1.4):
            S = plot_parts(fakes=True, dy_scale=k, key_rows=0, counter=False)
            self.play(Transform(stack, S["stack"]), Transform(ratio.dots, S["ratio"].dots),
                      Transform(rlabel, S["rlabel"]), run_time=rt, rate_func=EASE)

        to_stage(K_RAW, 1.0)
        self.wait(0.2)
        for j in (0, 1, 2, 3, 4):
            self.play(FadeIn(corr[j], shift=RIGHT * 0.25), run_time=0.45, rate_func=EASE)
        to_stage(K_SF)
        clip_cut(self, "tautau_g2_event_weights")
        for j in (5, 6):
            self.play(FadeIn(corr[j], shift=RIGHT * 0.25), run_time=0.45, rate_func=EASE)
        to_stage(1.0)
        for j in (7, 8):
            self.play(FadeIn(corr[j], shift=RIGHT * 0.25), run_time=0.45, rate_func=EASE)
        adopt(self, corr, "corr")
        settle_same(self, prev, end, PLOT_KEYS)
        check_order(self, end, ORDER_G)
        self.wait(0.1)


class TautauBDT(Scene):
    """tautau_h1_inputs ... h5_categories: from a blank frame, the classifier's inputs, the
    fakes-vs-signal shapes per input, a schematic BDT, the D_BDT output with its ratio panel,
    the three categories."""

    def construct(self):
        white_background(self)
        prev = state_g()
        add_state(self, prev, ORDER_G)
        clip_open(self, "tautau_h1_inputs")
        end = state_h()

        self.play(*[FadeOut(prev[k]) for k in ORDER_G], run_time=0.6)          # blank frame
        self.wait(0.2)
        lst = inputs_list()
        self.play(LaggedStart(*[FadeIn(t, shift=RIGHT * 0.2) for t in lst], lag_ratio=0.08, group=lst), run_time=2.4)
        clip_cut(self, "tautau_h2_shapes")

        self.play(Transform(lst, inputs_list(parked=True)), run_time=0.9, rate_func=EASE)
        FP = feature_panels()
        pan = [FP[f"f{k}"] for k in range(len(SHAPE_FEATURES))]
        self.play(LaggedStart(*[FadeIn(p) for p in pan], lag_ratio=0.25), run_time=3.2)
        self.play(FadeIn(FP["fkey"]), run_time=0.4)
        clip_cut(self, "tautau_h3_sketch")

        SK = bdt_sketch()
        self.play(*[FadeOut(p) for p in pan], FadeOut(FP["fkey"]), ReplacementTransform(lst, SK.bracket),
                  run_time=0.8, rate_func=EASE)
        self.play(FadeIn(SK.arrow_in), run_time=0.3)
        self.play(LaggedStart(*[Create(t) for t in SK.trees], lag_ratio=0.3), run_time=1.8)
        self.play(FadeIn(SK.plus), run_time=0.3)
        self.play(FadeIn(SK.arrow_out), FadeIn(SK.out, shift=RIGHT * 0.2), run_time=0.5)
        clip_cut(self, "tautau_h4_score")

        S = score_parts()
        self.play(FadeOut(SK.bracket), FadeOut(SK.arrow_in), *[FadeOut(t) for t in SK.trees], FadeOut(SK.plus),
                  FadeOut(SK.arrow_out), FadeOut(SK.out), run_time=0.5)
        self.play(FadeIn(S["s_dax"]), run_time=0.5)
        self.play(FadeIn(S["s_stack"]), run_time=0.8)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in S["s_data"]], lag_ratio=0.04, group=S["s_data"]),
                  run_time=1.0)
        rp = S["s_ratio"]
        self.play(FadeIn(rp.dax), FadeIn(rp.ref), run_time=0.5)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in rp.dots], lag_ratio=0.04, group=rp.dots), run_time=0.9)
        adopt(self, rp, "s_ratio")
        self.play(Create(S["s_cuts"]), FadeIn(S["s_auc"]), run_time=0.6)
        self.play(FadeIn(end["key"]), run_time=0.4)
        clip_cut(self, "tautau_h5_categories")

        self.play(*[FadeOut(S[k]) for k in SCORE_KEYS], run_time=0.5)
        for k in range(3):
            keys = [f"p{k}_dax", f"p{k}_stack", f"p{k}_data"] + (["p0_drop"] if k == 0 else []) + [f"p{k}_lab"]
            self.play(*[FadeIn(end[key]) for key in keys], run_time=0.6)
        self.bring_to_front(end["key"])
        check_order(self, end, ORDER_H)
        self.wait(0.1)


class TautauFit(Scene):
    """tautau_i1_fit ... i4_sigma_total: the profile-likelihood fit (pulls, mu_Z slider,
    post-fit panels), the grouped impacts, sigma_fid, sigma(60-120) beside the prediction."""

    def construct(self):
        white_background(self)
        prev = state_h()
        add_state(self, prev, ORDER_H)
        clip_open(self, "tautau_i1_fit")
        end = state_i()

        pp0, pp1 = pulls(False), pulls(True)
        sl = mu_slider()
        self.play(FadeIn(pp0), FadeIn(sl), run_time=0.7)
        post = panel_parts("postfit")
        self.play(*[Transform(prev[f"p{k}_stack"], post[f"p{k}_stack"]) for k in range(3)],
                  Transform(pp0.rows, pp1.rows), Transform(sl.marker, asym_marker(sl)), run_time=2.4, rate_func=EASE)
        for k in range(3):
            settle(self, prev[f"p{k}_stack"], post[f"p{k}_stack"], f"p{k}_stack")
        v, e = mu_strings()
        seed = tex_h(rf"\mu_{{Z}} = {v} {e}", 0.18).next_to(sl.marker, UP, buff=A["mu_seed_dy"])
        self.play(FadeIn(seed, shift=UP * 0.15), run_time=0.5)
        clip_cut(self, "tautau_i2_impacts")

        imp = impact_bars()
        self.play(FadeOut(pp0), FadeOut(sl), seed.animate.move_to(_p3(A["mu_seed_i2"])), run_time=0.6)
        self.play(FadeIn(imp.head), FadeIn(imp.names), run_time=0.5)
        self.play(LaggedStart(*[GrowFromEdge(b, LEFT) for b in imp.bars], lag_ratio=0.12, group=imp.bars), run_time=2.0)
        self.play(FadeIn(imp.nums), run_time=0.4)
        adopt(self, imp, "impacts")
        clip_cut(self, "tautau_i3_sigma_fid")

        sig = end["sig_tot"]
        moves = {f"p{k}_{p}": None for k in range(3) for p in ("dax", "data", "lab")}
        anims = [FadeOut(imp)]
        for key in PANEL_KEYS:
            src = post[key] if key.endswith("_stack") else prev[key]
            anims.append(ReplacementTransform(src, end[key]))
        anims += [ReplacementTransform(prev["key"], end["key"]), ReplacementTransform(seed, end["mu_line"]),
                  ReplacementTransform(seed.copy(), end["sig_fid"]), ReplacementTransform(seed.copy(), sig.label)]
        self.play(*anims, run_time=2.0, rate_func=EASE)
        clip_cut(self, "tautau_i4_sigma_total")
        self.play(FadeIn(sig.frame), run_time=0.5)
        self.play(Create(sig.theory), FadeIn(sig.th_lab), *([FadeIn(sig.band)] if sig.band is not None else []),
                  run_time=0.6)
        self.play(GrowFromCenter(sig.bar), GrowFromCenter(sig.dot), run_time=0.6)
        self.play(LaggedStart(*[FadeIn(r, shift=UP * 0.1) for r in sig.refs], lag_ratio=0.4, group=sig.refs),
                  run_time=1.0)
        self.remove(sig.label)
        self.add(sig.label)
        adopt(self, sig, "sig_tot")
        check_order(self, end, ORDER_I)
        self.wait(0.1)
