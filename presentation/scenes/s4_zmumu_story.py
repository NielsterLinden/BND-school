"""Section 4, real data: the Z -> mumu story, chained on 4-02 (``mumu_detector``).

Every number printed here is read from ``presentation/data/zmumu_*.json`` and
formatted in code (module-level anchor asserts stop the render on schema
drift). Eight plain ``Scene``s, each opening on the previous clip's final
frame: pure builders ``state_a() ... state_h()`` + ``ORDER_*`` tuples,
``check_order`` at the end of every clip, seams checked with
``tools/framediff.py``.

    MumuEvent        mumu_a_event        the schematic tracks of 4-02 fade; real SR event
                                         (run 278969) drawn from its muons' phi, charge, p_T;
                                         the p_T values converge to m_mumu; the slice parks
                                         left, log axes appear, first entry
    MumuRain         mumu_b_rain         three more real events, then the clock, the rain and
                                         the 60 data bins growing to their real counts
    MumuStack        mumu_c_stack        the *uncorrected* simulation stack slides in under the
                                         data; the ratio panel opens at data/pred = 0.944
    MumuControl      mumu_d_control      the plot parks; the slice returns inside the dashed
                                         control-region box with a real same-sign pair, the
                                         probe inside a jet
    MumuTenPairs     mumu_e_tenpairs     ten real same-sign pairs: nine anti-isolated probes,
                                         one isolated -> f = N_tight/N_anti = 1/9; the real
                                         6 x 4 fake-factor map
    MumuTransfer     mumu_f_transfer     the fake template enters the bottom of the stack as a
                                         thin wedge; the ratio does not move; 3870 +- 80
    MumuCorrections  mumu_g_corrections  pileup, L1 prefiring, ID / iso scale factors, trigger,
                                         kappa: the prediction walks down, the data never moves
    MumuFit          mumu_h_fit          rebin 60 -> 12, pulls, post-fit; mu_Z, sigma_fid,
                                         sigma(60-120) beside the prediction

Delivered clips (manim sections via clip_open / clip_cut, tools/deliver_chain.py):
    mumu_a1_event  mumu_a2_mass  mumu_a3_first_entry | mumu_b1_more_events  mumu_b2_rain |
    mumu_c1_simulation  mumu_c2_ratio | mumu_d1_control_region  mumu_d2_same_sign_pair |
    mumu_e1_ten_pairs  mumu_e2_fake_factor  mumu_e3_ff_map | mumu_f1_apply  mumu_f2_fake_template |
    mumu_g1_pileup  mumu_g2_prefiring  mumu_g3_lepton_sf |
    mumu_h1_rebin  mumu_h2_fit  mumu_h3_sigma_fid  mumu_h4_sigma_total

Physics honesty: the data points never move, only the prediction and the
ratio do. The log axis starts at 10^1.5 (31.6 events / GeV) so the ~65 / GeV
fake template is a thin wedge at the bottom of the stack, not a band that
would dominate the frame. Track curvature follows kappa_from_pt (exaggerated
scale, honest relative curvature); eta is not drawn (r-phi view). The 12-bin
counts of the fit are per 5 GeV (they sit log10(5) higher; the y title says so).
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
    DOWN, LEFT, RIGHT, UP, Create, DashedLine, DashedVMobject, DecimalNumber, Dot, FadeIn,
    FadeOut, GrowFromCenter, LaggedStart, Line, Rectangle, ReplacementTransform, RoundedRectangle,
    Scene, Succession, SurroundingRectangle, Transform, ValueTracker, VGroup, VMobject, rate_functions,
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

assert SR["data"]["total"] == 10378567
assert len(SR["data"]["counts"]) == 60
assert sum(SR["data"]["counts"]) == SR["data"]["total"]
assert abs(FK["yields"]["sr_fakes"] - 3869.955) < 0.01
assert abs(FIT["sigma_fid"]["value"] - 790.078) < 0.01
assert tuple(SR["stack_order"]) == ("Fakes", "WW", "WZ", "ZZ", "SingleTop", "TTbar", "DYtautau", "DYmumu")
assert abs(SR["fakes"]["total"] - FK["yields"]["sr_fakes"]) < 1e-6
assert len(EV["sr"]["chosen"]) == 4 and len(EV["ss_tp"]["chosen"]) == 10
assert sum(bool(p["passes"]) for p in EV["ss_tp"]["chosen"]) == 1
assert len(FIT["sr_12bin"]["data"]) == 12 and len(FIT["sr_12bin"]["edges"]) == 13

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
SR_EVENTS = EV["sr"]["chosen"]                     # textbook, FSR, forward, boosted
SS_PAIRS = EV["ss_tp"]["chosen"]                   # #0 = the mumu_d display pair, #9 passes

# ---------------------------------------------------------------------------
# one coordinate dictionary
# ---------------------------------------------------------------------------

EASE = rate_functions.ease_in_out_sine
PHI_M, PHI_P, KAPPA_402 = math.radians(35), math.radians(212), 0.28   # 4-02 (s4_zmumu.py)
R_LABEL = R_DET + 0.42                             # channel_common.show_signature
Y_SAFE = 2.42                                      # a label centre above this would enter the title band
DET_PARK = (0.42, (-4.6, -0.7))                    # slice parked left of the plot
DET_CR = (0.84, (-1.2, -0.45))                     # slice inside the control-region box (d, e)
CR_BOX = dict(center=(-1.2, -0.45), width=6.2, height=6.2)
PLOT_MAIN_C = (1.85, 0.85)                         # main plot area centre (natural)
RATIO_C = (1.85, -1.55)
PLOT_C = np.array([1.85, -0.25, 0.0])              # natural centre every plot placement scales about
PLOT_MAIN = (1.0, PLOT_C)
PLOT_GHOST = (0.44, (4.75, 1.24))
PLOT_LEFT = (0.58, (-4.35, -0.35))
Y_EXP = (1.5, 7)                                   # log axis 10^1.5 .. 10^7
RATIO_RANGE = (0.86, 1.06, 0.1)
KEY_LEFT = (-1.25, -2.92)                          # left end of the one-row colour key
COUNTER_AT = (-0.2, 2.12)                          # bottom-left of "N" (natural)
CLOCK_AT = (-0.72, 2.24)
STAMP_AT = (-4.6, -2.3)
LEFT_COL_X = -6.75                                 # left edge of the correction labels (g)
TALLY_C = (4.7, -0.62)
FF_TEX_C = (4.7, -1.28)
FFMAP_C = (4.85, -2.85)
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


def charge_tex(q: int) -> str:
    return r"\mu^{+}" if q > 0 else r"\mu^{-}"


def _p3(p):
    p = np.asarray(p, dtype=float)
    return p if p.shape == (3,) else np.array([p[0], p[1], 0.0])


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


# -- the verified swap and the adoption of individually revealed children ----------

def _look(m):
    """What a mobject draws: per family member with points, (style, bbox)."""
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
    """Swap a live mobject for the builder's (identity for check_order / the next
    clip) after proving they draw the same thing (same styles, same member boxes
    within 6e-3 units; coincident duplicates left by a 5 -> 1 Transform count once)."""
    ok, why = looks_same(live, target)
    if not ok:
        raise AssertionError(f"settle({what}): live mobject does not match the builder: {why}")
    scene.replace(live, target)


def adopt(scene, parent, what=""):
    """Children of ``parent`` revealed one by one sit top-level in ``scene.mobjects``;
    replace that contiguous run by ``parent`` (same draw order)."""
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
    """``rain()`` with each dot invisible until its own flight starts: a short FadeIn
    is prepended to every Succession (LaggedStart begins all children up front, so
    the FadeIn's start state hides the dot; without it all dots sit in the tracker
    as one dark blob from the first frame)."""
    out = []
    for a in rain(*args, **kw):
        move, fade = a.animations
        out.append(Succession(FadeIn(move.mobject, run_time=0.06), move, fade))
    return out


# ---------------------------------------------------------------------------
# plot builders (always at the natural position, then place())
# ---------------------------------------------------------------------------

def sr_layers(stage, fakes):
    rows = []
    for n in LAYERS:
        if n == "Fakes":
            c = FAKES if fakes else np.zeros(60)
        else:
            c = np.asarray(SR["mc"][n][stage], dtype=float) if stage else np.zeros(60)
        rows.append((n, c, SAMPLE[n]))
    return rows


def fit_layers(kind):
    return [(n, np.asarray(B12[kind][n], dtype=float), SAMPLE[n]) for n in LAYERS]


def main_axes(per5=False, x_labels=False) -> DataAxes:
    """The log m_mumu axes. ``x_labels``: tick labels + title on the main plot itself
    (clips a, b, before the ratio panel takes them over). ``per5``: counts per 5 GeV,
    the axis range shifted by log10(5) so the rebinned stack keeps its height on screen
    and the relabelled ticks (10^3 ... 10^7) carry the factor 5."""
    off = math.log10(5.0) if per5 else 0.0
    e0, e1 = Y_EXP[0] + off, Y_EXP[1] + off
    dax = DataAxes([60, 120, 10], [e0, e1, 1], 6.4, 2.9,
                   y_ticks=[k for k in range(1, 9) if e0 < k <= e1 + 1e-9],
                   y_log=True, show_x_labels=x_labels, title_h=0.2, title_buff=0.16,
                   x_title=r"m_{\mu\mu}\ [\mathrm{GeV}]" if x_labels else None,
                   y_title=r"\mathrm{events}\,/\,5\,\mathrm{GeV}" if per5 else r"\mathrm{events}\,/\,\mathrm{GeV}")
    dax.move_frame_to(PLOT_MAIN_C)
    return dax


def entry_dot(dax, m) -> Dot:
    return Dot(dax.c2p(float(m), 1.0), radius=0.05, color=col(SAMPLE["Data"]))


def data_bars(dax, visible=True) -> VGroup:
    g = VGroup(*[data_bar(dax, XC[i], DATA[i], 0.5, color=SAMPLE["Data"], fill_opacity=0.55,
                          stroke_width=1.0) for i in range(60)])
    if not visible:        # the same bars at the floor, fully transparent (Transform source)
        g = VGroup(*[data_bar(dax, XC[i], 0.0, 0.5, color=SAMPLE["Data"], fill_opacity=0.55,
                              stroke_width=1.0).set_opacity(0.0) for i in range(60)])
    return g


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


def make_key(rows: int) -> VGroup:
    entries = [(r"\mathrm{Data}", SAMPLE["Data"], "dot"), (r"Z\to\mu\mu", SAMPLE["DYmumu"]),
               (r"\tau\tau,\ t\bar{t},\ tW,\ VV", SAMPLE["TTbar"])]
    if rows == 4:
        entries.append((r"\mathrm{Fakes}", SAMPLE["Fakes"]))
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


def ratio_value(stage, fakes) -> float:
    t = SR["totals"][stage]
    if fakes:
        return float(t["data_over_pred"])
    return TOTAL / float(t["mc"])


def plot_parts(placement=PLOT_MAIN, *, stage=None, fakes=False, data="dots", bins=60, fit="prefit",
               counter=True, ratio=True, rlabel=True, key_rows=3, entries=(), x_labels=False) -> dict:
    """The plot family of the chain, built at the natural position and then placed
    about PLOT_C: dax, stack (8 layers bottom-up), data (bars / dots / entry dots),
    counter, ratio (panel), rlabel (total data/pred), key."""
    P = {}
    dax = main_axes(per5=bins == 12, x_labels=x_labels)
    P["dax"] = dax
    if bins == 60:
        P["stack"] = stack_hist(dax, EDGES, sr_layers(stage, fakes))
    else:
        P["stack"] = stack_hist(dax, EDGES12, fit_layers(fit))
    if data == "entries":
        P["data"] = VGroup(*[entry_dot(dax, m) for m in entries])
    elif data == "bars":
        P["data"] = data_bars(dax)
    elif data == "bars0":
        P["data"] = data_bars(dax, visible=False)
    else:
        P["data"] = data_dots(dax, bins)
    if counter:
        P["counter"] = make_counter(TOTAL)
    if ratio:
        kw = dict(y_range=RATIO_RANGE, y_length=0.9, y_ticks=[0.9, 1.0], dot_radius=RDOT_R,
                  x_title=r"m_{\mu\mu}\ [\mathrm{GeV}]", y_title=r"\mathrm{data/pred.}",
                  title_h=0.2, title_buff=0.16)
        if bins == 60:
            P["ratio"] = ratio_panel(dax, EDGES, DATA, P["stack"].total, RATIO_C, **kw)
        else:
            P["ratio"] = ratio_panel(dax, EDGES12, np.asarray(B12[f"ratio_{fit}"]), np.ones(12), RATIO_C, **kw)
        if rlabel and bins == 60:
            v = ratio_value(stage, fakes)
            lab = tex_h(f"{v:.3f}", 0.2)
            lab.next_to(P["ratio"].dax.c2p(120, v), RIGHT, buff=0.22)
            P["rlabel"] = lab
    if key_rows:
        P["key"] = make_key(key_rows)
    for m in P.values():
        place(m, placement, PLOT_C)
    return P


PLOT_KEYS = ("dax", "stack", "data", "counter", "ratio", "rlabel", "key")


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


def cr_box() -> VGroup:
    rect = RoundedRectangle(width=CR_BOX["width"], height=CR_BOX["height"], corner_radius=0.25,
                            stroke_color=col(DETECTOR_ACCENT), stroke_width=3.0, fill_opacity=0.0)
    rect.move_to(_p3(CR_BOX["center"]))
    dashed_rect = DashedVMobject(rect, num_dashes=70)
    sym = tex_h(r"\mu^{\pm}\mu^{\pm}", 0.24, color=DETECTOR_ACCENT)
    sym.move_to(rect.get_corner(UP + LEFT) + np.array([0.62, -0.36, 0.0]))
    g = VGroup(dashed_rect, sym)
    g.rect, g.sym = dashed_rect, sym
    return g


def ss_pair(i: int) -> VGroup:
    """Real same-sign tag + probe pair i on a natural slice, placed into the box."""
    pair = pair_from_json(CMSSlice(), SS_PAIRS[i], seed=i)
    tag_phi, probe_trk = phi_end(pair.tag.trk), (pair.probe.mu.trk if pair.in_jet else pair.probe.trk)
    place(pair, DET_CR, DET_CENTER)
    pair.tag_phi, pair.probe_phi = tag_phi, phi_end(probe_trk)
    return pair


def pair_labels(pair) -> VGroup:
    rec = pair.record
    labs = VGroup()
    for who, phi in (("tag", pair.tag_phi), ("probe", pair.probe_phi)):
        lab = mathtex(charge_tex(_q(rec[who])), color=CHANNEL_LINE["mumu"])
        labs.add(lab.move_to(label_spot(phi, DET_CR)))
    return labs


def tally_boxes(filled: bool) -> VGroup:
    g = VGroup()
    for p in SS_PAIRS:
        box = Rectangle(width=0.34, height=0.34, stroke_color=col(INK), stroke_width=1.6,
                        fill_color=col(tally_colour(p)), fill_opacity=1.0 if filled else 0.0)
        g.add(box)
    g.arrange(RIGHT, buff=0.08).move_to(_p3(TALLY_C))
    return g


def tally_colour(p):
    return CHANNEL["mumu"] if p["passes"] else SAMPLE["Fakes"]


def ff_formula():
    n_tight = sum(bool(p["passes"]) for p in SS_PAIRS)
    n_anti = len(SS_PAIRS) - n_tight
    return tex_h(rf"f = N_{{\mathrm{{tight}}}} / N_{{\mathrm{{anti}}}} = {n_tight}/{n_anti}", 0.22).move_to(_p3(FF_TEX_C))


def ff_map() -> VGroup:
    pt, eta = FK["pt_edges"], FK["eta_edges"]
    rows = [rf"{pt[i]:g}\!-\!{pt[i + 1]:g}" for i in range(len(pt) - 1)]
    cols = [rf"{eta[j]:g}\!-\!{eta[j + 1]:g}" for j in range(len(eta) - 1)]
    vg = value_grid(FK["ff"]["nominal"], rows, cols, fmt="{:.2f}", cell=(0.80, 0.28), label_h=0.13)
    head = tex_h(r"p_T \,\backslash\, |\eta|", 0.13)
    head.next_to(vg.row_labels, UP, buff=0.10).align_to(vg.col_labels, DOWN)
    n_col = len(cols)
    barrel = VGroup(*[vg.cells[i * n_col] for i in range(4)])      # |eta| < 0.9, pT < 50: the pairs' cells
    hl = SurroundingRectangle(barrel, color=col(DETECTOR_ACCENT), buff=0.0, stroke_width=3.0)
    g = VGroup(vg, head, hl)
    g.move_to(_p3(FFMAP_C))
    return g


def corr_labels() -> VGroup:
    mf = SR["mean_factors"]["all_mc"]
    sf = CO["scale_factors"]["map_mean"]
    tr = CO["trigger"]
    kap = CO["momentum"]["kappa"]
    exprs = [
        rf"\langle w_{{\mathrm{{PU}}}} \rangle = {mf['pileup']:.3f}",
        rf"\langle w_{{\mathrm{{L1}}}} \rangle = {mf['prefiring']:.3f}",
        rf"\mathrm{{SF}}_{{\mathrm{{ID}}}} = {sf['id']:.3f}",
        rf"\mathrm{{SF}}_{{\mathrm{{iso}}}} = {sf['iso']:.3f}",
        rf"\varepsilon^{{\mathrm{{data}}}}_{{\mathrm{{trig}}}} / \varepsilon^{{\mathrm{{MC}}}}_{{\mathrm{{trig}}}}"
        rf" = {tr['plateau_data']:.3f} / {tr['plateau_mc']:.3f}",
        rf"\kappa_{{|\eta|}} = {min(kap):.3f} \ldots {max(kap):.3f}",
    ]
    g = VGroup()
    for i, e in enumerate(exprs):
        lab = tex_h(e, 0.22)
        lab.shift(np.array([LEFT_COL_X, 2.05 - 0.82 * i, 0.0]) - lab.get_left())
        g.add(lab)
    return g


# -- the fit numbers --------------------------------------------------------------

NP_SHOWN = ("SigModel", "MuonScale", "MuonRes", "PDF", "QCDScale", "Lumi")
NP = {n["name"]: n for n in FIT["nps"]}
PULLS_AT = (-4.35, 0.95)
SLIDER_AT = (-4.35, -1.75)


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


def mu_line(h=0.34, at=(2.9, 2.0)):
    v, e = mu_strings()
    return tex_h(rf"\mu_{{Z}} = {v} {e}", h).move_to(_p3(at))


def sigma_fid_line():
    s = FIT["sigma_fid"]
    expr = (rf"\sigma_{{\mathrm{{fid}}}} = {s['value']:.1f} \pm {s['stat']:.1f}_{{\mathrm{{stat}}}}"
            rf" \pm {s['syst']:.1f}_{{\mathrm{{syst}}}} \pm {s['lumi']:.1f}_{{\mathrm{{lumi}}}}\ \mathrm{{pb}}")
    return tex_h(expr, 0.27).move_to(np.array([2.95, 0.95, 0.0]))


SIG_AX = dict(lo=1850.0, hi=2050.0, x0=0.25, length=5.4, y=-1.95)


def sigma_panel() -> VGroup:
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
    unit = tex_h(r"\mathrm{pb}", 0.17).next_to(tick_labels[-1], RIGHT, buff=0.18)
    pred = float(s["pred"])
    theory = DashedLine([xs(pred), y0, 0], [xs(pred), y0 + 0.95, 0], dash_length=0.1,
                        stroke_color=col(THEORY), stroke_width=3.0)
    th_lab = tex_h(f"{pred:.1f}", 0.17, color=THEORY).next_to(theory.get_end(), UP, buff=0.08)
    ypt = y0 + 0.45
    cc = col(CHANNEL_LINE["mumu"])
    bar = Line([xs(s["value"] - s["total"]), ypt, 0], [xs(s["value"] + s["total"]), ypt, 0],
               stroke_color=cc, stroke_width=5.0)
    dot = Dot([xs(s["value"]), ypt, 0], radius=0.1, color=cc)
    label = tex_h(rf"\sigma_{{60\text{{--}}120}} = {s['value']:.0f} \pm {s['total']:.0f}\ \mathrm{{pb}}", 0.27)
    label.move_to(np.array([2.95, -0.45, 0.0]))
    g = VGroup(axis, ticks, tick_labels, unit, theory, th_lab, bar, dot, label)
    g.axis, g.theory, g.th_lab, g.bar, g.dot, g.label = axis, theory, th_lab, bar, dot, label
    g.frame = VGroup(axis, ticks, tick_labels, unit)
    return g


# ---------------------------------------------------------------------------
# end states (pure) and their z-orders
# ---------------------------------------------------------------------------

def detector_end_state() -> dict:
    """Exact replica of 4-02's last frame (s4_zmumu.ZmumuDetector: CMSSlice +
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
    rec = SR_EVENTS[0]
    P = plot_parts(stage=None, data="entries", entries=(rec["mass"],), counter=False, ratio=False, key_rows=0,
                   x_labels=True)
    return {"det": parked_slice(), "ev": place(event_group(rec), DET_PARK, DET_CENTER),
            "stamp": run_stamp(rec), "dax": P["dax"], "stack": P["stack"], "data": P["data"]}


ORDER_A = ("det", "ev", "stamp", "dax", "stack", "data")


def state_b() -> dict:
    P = plot_parts(stage=None, data="bars", ratio=False, key_rows=0, x_labels=True)
    return {"det": parked_slice(), **P}


ORDER_B = ("det", "dax", "stack", "data", "counter")


def state_c() -> dict:
    return {"det": parked_slice(), **plot_parts(stage="raw", fakes=False, key_rows=3)}


ORDER_C = ("det", *PLOT_KEYS)


def state_d() -> dict:
    st = {"det": parked_slice(DET_CR), **plot_parts(PLOT_GHOST, stage="raw", fakes=False, key_rows=3)}
    pair = ss_pair(0)
    st.update({"box": cr_box(), "tag": pair.tag, "probe": pair.probe, "labs": pair_labels(pair)})
    return st


ORDER_D = ("det", *PLOT_KEYS, "box", "tag", "probe", "labs")


def state_e() -> dict:
    st = {"det": parked_slice(DET_CR), **plot_parts(PLOT_GHOST, stage="raw", fakes=False, key_rows=3)}
    st.update({"box": cr_box(), "tally": tally_boxes(True), "ff_tex": ff_formula(), "ffmap": ff_map()})
    return st


ORDER_E = ("det", *PLOT_KEYS, "box", "tally", "ff_tex", "ffmap")


N_FAKE_AT = (-4.75, 0.35)


def n_fakes_big():
    y = FK["yields"]
    return tex_h(rf"N_{{\mathrm{{fake}}}} = {y['sr_fakes']:.0f} \pm {y['sr_fakes_stat']:.0f}", 0.32).move_to(
        _p3(N_FAKE_AT))


def n_fakes_small() -> VGroup:
    y = FK["yields"]
    lab = tex_h(rf"{y['sr_fakes']:.0f} \pm {y['sr_fakes_stat']:.0f}", 0.2)
    dax = main_axes()            # natural axes: the label sits right of the wedge's end
    lab.next_to(dax.c2p(120, 45.0), RIGHT, buff=0.22)
    return lab


def state_f() -> dict:
    P = plot_parts(stage="raw", fakes=True, key_rows=4)
    return {**P, "n_fakes": n_fakes_big()}


ORDER_F = (*PLOT_KEYS, "n_fakes")


def state_g() -> dict:
    P = plot_parts(stage="nominal", fakes=True, key_rows=4)
    return {**P, "n_fakes": n_fakes_small(), "corr": corr_labels()}


ORDER_G = (*PLOT_KEYS, "n_fakes", "corr")


def state_h() -> dict:
    P = plot_parts(PLOT_LEFT, bins=12, fit="postfit", counter=False, rlabel=False, key_rows=4)
    return {"dax": P["dax"], "stack": P["stack"], "data": P["data"], "ratio": P["ratio"], "key": P["key"],
            "mu_line": mu_line(), "sig_fid": sigma_fid_line(), "sig_tot": sigma_panel()}


ORDER_H = ("dax", "stack", "data", "ratio", "key", "mu_line", "sig_fid", "sig_tot")


def settle_same(scene, prev, end, keys):
    for k in keys:
        settle(scene, prev[k], end[k], k)


# ---------------------------------------------------------------------------
# the eight clips
# ---------------------------------------------------------------------------

class MumuEvent(Scene):
    """mumu_a_event: 4-02's schematic tracks fade, the real event is drawn, its p_T
    values converge to m_mumu, the slice parks, the log axes appear, first entry."""

    def construct(self):
        white_background(self)
        prev = detector_end_state()
        add_state(self, prev, ORDER_402)
        clip_open(self, "mumu_a1_event")
        end = state_a()
        rec = SR_EVENTS[0]

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
        adopt(self, grp, "event 1")
        self.play(*[FadeIn(p, shift=0.25 * (RIGHT if p.get_center()[0] > 0 else LEFT)) for p in pts],
                  run_time=0.5)
        clip_cut(self, "mumu_a2_mass")
        m_tex = mass_tex(rec).move_to(np.array([4.95, 0.2, 0.0]))
        self.remove(*pts)
        self.add(pts)
        self.play(ReplacementTransform(pts, m_tex), run_time=0.9, rate_func=EASE)
        clip_cut(self, "mumu_a3_first_entry")

        m_spot = end["dax"].c2p(float(rec["mass"]), 10 ** 6.7)      # MAIN placement = natural
        self.play(ReplacementTransform(prev["det"], end["det"]), ReplacementTransform(grp, end["ev"]),
                  FadeOut(prev["lab_m"]), FadeOut(prev["lab_p"]), m_tex.animate.move_to(m_spot),
                  FadeIn(end["stamp"]), FadeIn(end["dax"]), run_time=1.8, rate_func=EASE)
        self.add(end["stack"])                      # 8 degenerate layers at the floor (invisible)
        self.bring_to_front(m_tex)
        self.wait(0.2)
        self.play(ReplacementTransform(m_tex, end["data"]), run_time=0.8, rate_func=EASE)
        check_order(self, end, ORDER_A)
        self.wait(0.2)


class MumuRain(Scene):
    """mumu_b_rain: three more real events, each dropping an entry; then the clock,
    the rain and the 60 data bins growing to the real counts (N = 10,378,567)."""

    def construct(self):
        white_background(self)
        prev = state_a()
        add_state(self, prev, ORDER_A)
        clip_open(self, "mumu_b1_more_events")
        end = state_b()
        det, dax, stamp, entries = prev["det"], prev["dax"], prev["stamp"], prev["data"]

        self.play(FadeOut(prev["ev"]), run_time=0.35)
        for rec in SR_EVENTS[1:4]:
            g = place(event_group(rec), DET_PARK, DET_CENTER)
            lines = VGroup(*[mu.trk for mu in g.ev], *[ph[0] for ph in g.photons])
            rest = VGroup(*[VGroup(mu.hits, mu.deposits) for mu in g.ev], *[ph[1] for ph in g.photons])
            self.play(Create(lines, lag_ratio=0.0), Transform(stamp, run_stamp(rec)), run_time=0.7,
                      rate_func=EASE)
            self.play(FadeIn(rest), run_time=0.3)
            dot = Dot(det.get_center(), radius=0.05, color=col(SAMPLE["Data"]))
            self.play(dot.animate.move_to(dax.c2p(float(rec["mass"]), 1.0)), FadeOut(lines), FadeOut(rest),
                      run_time=0.6, rate_func=rate_functions.linear)
            self.remove(dot)
            entries.add(dot)
            if rec is not SR_EVENTS[3]:
                self.wait(0.1)
        clip_cut(self, "mumu_b2_rain")

        bars = end["data"].copy()
        bars0 = data_bars(dax, visible=False)
        self.remove(entries)
        self.add(bars0, entries)
        clk = clock(CLOCK_AT)
        counter = make_counter(len(entries))
        self.play(FadeIn(clk), FadeIn(counter), FadeOut(stamp), run_time=0.45)
        n = ValueTracker(float(len(entries)))
        counter.num.add_updater(lambda m: (m.set_value(float(round(n.get_value()))),
                                           layout_counter(counter.pre, m)))
        order = np.random.default_rng(20260915).permutation(60)
        drops = hidden_rain(det, dax, EDGES, DATA, n=70, seed=3, color=SAMPLE["Data"])
        self.play(LaggedStart(*[Transform(bars0[i], bars[i]) for i in order], lag_ratio=0.02, group=bars0),
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
        settle_same(self, prev, end, ("det", "dax", "stack"))
        check_order(self, end, ORDER_B)
        self.wait(0.2)


class MumuStack(Scene):
    """mumu_c_stack: bars become points; the uncorrected simulation slides in under
    them bottom-up; the ratio panel opens at data/pred = 0.944."""

    def construct(self):
        white_background(self)
        prev = state_b()
        add_state(self, prev, ORDER_B)
        clip_open(self, "mumu_c1_simulation")
        end = state_c()

        self.play(ReplacementTransform(prev["data"], end["data"]), run_time=1.0, rate_func=EASE)
        self.wait(0.2)
        stack = prev["stack"]
        dax = prev["dax"]
        base = sr_layers("raw", False)
        for i, poly in enumerate(stack):            # each (invisible) layer waits at its own base
            rows = [(n, c if j < i else np.zeros(60), cc) for j, (n, c, cc) in enumerate(base)]
            poly.become(stack_hist(dax, EDGES, rows)[i])
        self.play(LaggedStart(*[Transform(stack[i], end["stack"][i]) for i in range(len(stack))],
                              lag_ratio=0.15, group=stack), run_time=2.6, rate_func=EASE)
        settle(self, stack, end["stack"], "stack")
        self.play(FadeIn(end["key"]), run_time=0.6)
        clip_cut(self, "mumu_c2_ratio")
        ratio = end["ratio"]
        self.play(FadeIn(ratio.dax), FadeIn(ratio.ref), dax.x_labels.animate.set_opacity(0.0),
                  dax.x_title.animate.set_opacity(0.0), run_time=0.6)
        self.play(LaggedStart(*[GrowFromCenter(d) for d in ratio.dots], lag_ratio=0.03, group=ratio.dots),
                  run_time=1.3)
        adopt(self, ratio, "ratio")
        self.wait(0.3)
        self.play(FadeIn(end["rlabel"], shift=LEFT * 0.2), run_time=0.5)
        self.bring_to_front(end["key"])
        settle_same(self, prev, end, ("det", "dax", "counter"))
        check_order(self, end, ORDER_C)
        self.wait(0.2)


class MumuControl(Scene):
    """mumu_d_control: the plot parks top right; the slice returns into the dashed
    control-region box with a real same-sign pair, the probe inside a jet."""

    def construct(self):
        white_background(self)
        prev = state_c()
        add_state(self, prev, ORDER_C)
        clip_open(self, "mumu_d1_control_region")
        end = state_d()

        self.play(*[ReplacementTransform(prev[k], end[k]) for k in ("det", *PLOT_KEYS)],
                  run_time=1.6, rate_func=EASE)
        self.wait(0.2)
        box = end["box"]
        self.play(Create(box.rect), FadeIn(box.sym), run_time=0.8)
        adopt(self, box, "box")
        clip_cut(self, "mumu_d2_same_sign_pair")
        tag, probe = end["tag"], end["probe"]
        self.play(Create(tag.trk, lag_ratio=0.0), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(tag.hits), run_time=0.3)
        self.play(FadeIn(tag.deposits), run_time=0.45)
        adopt(self, tag, "tag")
        self.wait(0.2)
        self.play(Create(probe.jet[0], lag_ratio=0.0), run_time=0.7, rate_func=EASE)
        self.play(FadeIn(VGroup(*probe.jet[1:])), run_time=0.4)
        self.play(Create(probe.mu.trk, lag_ratio=0.0), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(VGroup(probe.mu.hits, probe.mu.deposits)), run_time=0.45)
        adopt(self, probe, "probe")
        self.wait(0.3)
        self.play(FadeIn(end["labs"]), run_time=0.4)
        check_order(self, end, ORDER_D)
        self.wait(0.2)


class MumuTenPairs(Scene):
    """mumu_e_tenpairs: ten real same-sign pairs, a ten-box tally (anti-isolated
    grey, isolated gold) -> f = N_tight/N_anti = 1/9; the real fake-factor map."""

    def construct(self):
        white_background(self)
        prev = state_d()
        add_state(self, prev, ORDER_D)
        clip_open(self, "mumu_e1_ten_pairs")
        end = state_e()

        tally = tally_boxes(False)
        self.play(FadeIn(tally), run_time=0.4)
        self.play(tally[0].animate.set_fill(col(tally_colour(SS_PAIRS[0])), opacity=1.0), run_time=0.3)
        self.play(FadeOut(prev["tag"]), FadeOut(prev["probe"]), FadeOut(prev["labs"]), run_time=0.25)
        for i in range(1, len(SS_PAIRS)):
            pair = ss_pair(i)
            tag = pair.tag
            probe_mu = pair.probe.mu if pair.in_jet else pair.probe
            jet = pair.probe.jet if pair.in_jet else VGroup()
            lines = VGroup(*jet[:1], tag.trk, probe_mu.trk)          # jet under the muons (muon_in_jet order)
            rest = VGroup(*jet[1:], tag.hits, tag.deposits, probe_mu.hits, probe_mu.deposits)
            self.play(Create(lines, lag_ratio=0.0), FadeIn(rest), run_time=0.45, rate_func=EASE)
            self.play(tally[i].animate.set_fill(col(tally_colour(SS_PAIRS[i])), opacity=1.0), run_time=0.22)
            self.play(FadeOut(lines), FadeOut(rest), run_time=0.22)
        settle(self, tally, end["tally"], "tally")
        clip_cut(self, "mumu_e2_fake_factor")
        self.play(FadeIn(end["ff_tex"], shift=UP * 0.15), run_time=0.9)
        clip_cut(self, "mumu_e3_ff_map")
        self.play(FadeIn(end["ffmap"]), run_time=1.1)
        settle_same(self, prev, end, ("det", *PLOT_KEYS, "box"))
        check_order(self, end, ORDER_E)
        self.wait(0.2)


class MumuTransfer(Scene):
    """mumu_f_transfer: the plot returns; f applied to the anti-isolated events gives
    the fake template, which enters the bottom of the stack as a thin wedge; the
    ratio does not move; N_fake = 3870 +- 80."""

    def construct(self):
        white_background(self)
        prev = state_e()
        add_state(self, prev, ORDER_E)
        clip_open(self, "mumu_f1_apply")
        end = state_f()
        mid = plot_parts(stage="raw", fakes=False, key_rows=3)

        f_tex = prev["ff_tex"]
        self.play(FadeOut(prev["det"]), FadeOut(prev["box"]), FadeOut(prev["tally"]), FadeOut(prev["ffmap"]),
                  f_tex.animate.move_to(np.array([-4.75, 1.2, 0.0])), run_time=0.8, rate_func=EASE)
        self.play(*[ReplacementTransform(prev[k], mid[k]) for k in PLOT_KEYS], run_time=1.5, rate_func=EASE)
        formula = tex_h(r"N_{\mathrm{fake}} = f \times N_{\mathrm{anti}}", 0.26).move_to(np.array([-4.75, 1.2, 0.0]))
        self.play(ReplacementTransform(f_tex, formula), run_time=0.7, rate_func=EASE)
        clip_cut(self, "mumu_f2_fake_template")
        tmpl = step_hist(mid["dax"], EDGES, FAKES, color=darken(SAMPLE["Fakes"], 0.55), stroke_width=4.0)
        self.play(Create(tmpl), run_time=1.1, rate_func=rate_functions.linear)
        self.wait(0.2)
        self.play(Transform(mid["stack"], end["stack"]), Transform(mid["ratio"].dots, end["ratio"].dots),
                  Transform(mid["rlabel"], end["rlabel"]), ReplacementTransform(mid["key"], end["key"]),
                  FadeOut(tmpl), run_time=1.5, rate_func=EASE)
        self.play(ReplacementTransform(formula, end["n_fakes"]), run_time=0.9, rate_func=EASE)
        settle_same(self, mid, end, ("dax", "stack", "data", "counter", "ratio", "rlabel"))
        check_order(self, end, ORDER_F)
        self.wait(0.2)


class MumuCorrections(Scene):
    """mumu_g_corrections: each correction appears with its real value and the
    prediction steps through the frozen stages raw -> pileup -> prefiring ->
    nominal; the data never move, the ratio climbs from 0.944 to 0.994."""

    def construct(self):
        white_background(self)
        prev = state_f()
        add_state(self, prev, ORDER_F)
        clip_open(self, "mumu_g1_pileup")
        end = state_g()
        corr = end["corr"]
        stack, ratio, rlabel = prev["stack"], prev["ratio"], prev["rlabel"]
        self.play(ReplacementTransform(prev["n_fakes"], end["n_fakes"]), run_time=0.9, rate_func=EASE)
        self.wait(0.2)

        beats = (((0,), "pileup", None), ((1,), "prefiring", "mumu_g2_prefiring"),
                 ((2, 3, 4), "nominal", "mumu_g3_lepton_sf"), ((5,), None, None))
        for b, (idxs, stage, clip) in enumerate(beats):
            if clip:
                clip_cut(self, clip)
            elif b:
                self.wait(0.3)
            for j in idxs:
                self.play(FadeIn(corr[j], shift=RIGHT * 0.25), run_time=0.5, rate_func=EASE)
            if stage is not None:
                S = plot_parts(stage=stage, fakes=True, key_rows=0, counter=False)
                self.play(Transform(stack, S["stack"]), Transform(ratio.dots, S["ratio"].dots),
                          Transform(rlabel, S["rlabel"]), run_time=1.6, rate_func=EASE)
        adopt(self, corr, "corr")
        settle_same(self, prev, end, PLOT_KEYS)
        check_order(self, end, ORDER_G)
        self.wait(0.2)


class MumuFit(Scene):
    """mumu_h_fit: rebin 60 -> 12 (5 GeV); the pulls and mu_Z move as the fit goes
    post-fit and the ratio flattens; mu_Z, sigma_fid and sigma(60-120) end up owning
    the frame beside the prediction."""

    def construct(self):
        white_background(self)
        prev = state_g()
        add_state(self, prev, ORDER_G)
        clip_open(self, "mumu_h1_rebin")
        end = state_h()

        self.play(FadeOut(prev["corr"]), FadeOut(prev["counter"]), FadeOut(prev["rlabel"]),
                  FadeOut(prev["n_fakes"]), run_time=0.5)
        pre12 = plot_parts(bins=12, fit="prefit", counter=False, rlabel=False, key_rows=4)
        anims = [Transform(prev["dax"], pre12["dax"]), Transform(prev["stack"], pre12["stack"])]
        for j in range(12):
            for k in range(5):
                anims.append(Transform(prev["data"][5 * j + k], pre12["data"][j].copy()))
                anims.append(Transform(prev["ratio"].dots[5 * j + k], pre12["ratio"].dots[j].copy()))
        self.play(*anims, run_time=1.7, rate_func=EASE)
        for k in ("dax", "stack", "data", "ratio", "key"):
            settle(self, prev[k], pre12[k], k)
        clip_cut(self, "mumu_h2_fit")

        pp0, pp1 = pulls(False), pulls(True)
        sl = mu_slider(False)
        self.play(FadeIn(pp0), FadeIn(sl), run_time=0.7)
        post12 = plot_parts(bins=12, fit="postfit", counter=False, rlabel=False, key_rows=4)
        poi = FIT["poi"]
        self.play(Transform(pre12["stack"], post12["stack"]), Transform(pre12["ratio"].dots, post12["ratio"].dots),
                  Transform(pp0.rows, pp1.rows),
                  Transform(sl.marker, sl.marker_at(poi["value"], max(poi["err_up"], poi["err_down"]))),
                  run_time=2.6, rate_func=EASE)
        for k in ("dax", "stack", "data", "ratio", "key"):
            settle(self, pre12[k], post12[k], k)
        v, e = mu_strings()
        seed = tex_h(rf"\mu_{{Z}} = {v} {e}", 0.18).next_to(sl.marker, UP, buff=0.62)
        self.play(FadeIn(seed, shift=UP * 0.15), run_time=0.5)
        clip_cut(self, "mumu_h3_sigma_fid")

        sig = end["sig_tot"]
        self.play(*[ReplacementTransform(post12[k], end[k]) for k in ("dax", "stack", "data", "ratio", "key")],
                  FadeOut(pp0), FadeOut(sl),
                  ReplacementTransform(seed, end["mu_line"]),
                  ReplacementTransform(seed.copy(), end["sig_fid"]),
                  ReplacementTransform(seed.copy(), sig.label),
                  run_time=2.0, rate_func=EASE)
        clip_cut(self, "mumu_h4_sigma_total")
        self.play(FadeIn(sig.frame), run_time=0.5)
        self.play(Create(sig.theory), FadeIn(sig.th_lab), run_time=0.6)
        self.play(GrowFromCenter(sig.bar), GrowFromCenter(sig.dot), run_time=0.6)
        self.remove(sig.label)
        self.add(sig.label)
        adopt(self, sig, "sig_tot")
        check_order(self, end, ORDER_H)
        self.wait(0.2)
