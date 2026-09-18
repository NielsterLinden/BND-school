"""Section 3, real data: the Z -> ee chapter, told on the section-2 map (deck owner, 17 Sep 2026).

The chapter opens on the three maps that section 2 ends on (2-21 map_three), grows the ee map and zooms into two
of its nodes: the selection (what we keep, and why) and the fit (a template fit with nuisance parameters).
Result version: the channel's own TRExFitter fit, z-ee delivery of 16 Sep 14:44 (docs/FREEZE.md "Z -> ee"):
mu_Z = 0.942 +- 0.015, sigma(60-120) = 1841 +- 30 pb; the result frame is the one of the tau tau chapter (5-33).

    ZeeSelection (opens on 2-21's last frame)
        ee_map_selection    the ee map grows out of the three maps, the camera flies into the selection node
        ee_sel_event        a real selected event (run 279024): the trigger electron, e+ e- with opposite charge
        ee_sel_pt           electron p_T in real data: the Z electrons peak near m_Z / 2; p_T > 20 GeV
        ee_sel_eta          side view: the tracker ends at |eta| = 2.5; beyond it an electron has no track
        ee_sel_id           electron vs jet: shower width, H/E, isolation (two real candidates)
        ee_sel_wp           the ID levels on real m_ee spectra: the same-sign pairs drop faster than the peak;
                            purity vs efficiency, medium
        ee_sel_mass         the cut list folds into the funnel; the selected pairs fill m_ee: N = 6,320,097
        ee_sel_mc           the simulation through the same selection: the stack, data / pred. = 0.970
        ee_sel_counting     N_data = sigma L eps A + N_bkg -> 1896 pb; the histogram flattens: the shape is lost
    ZeeFit (opens on ZeeSelection's last frame, pure builder ``counting_state``)
        ee_map_fit          zoom out to the ee map, zoom into the fit node: the shape is back
        ee_fit_mu           nu_i = mu s_i + b_i: mu_Z scales Z -> ee only; prediction and ratio move
        ee_fit_nuisance     one nuisance parameter at a time on its real +-1 sigma template (luminosity, electron ID,
                            QCD scale, FSR)
        ee_fit_result       the likelihood; every parameter goes to its post-fit value, the ratio flattens, the pulls;
                            mu_Z = 0.942 +- 0.015
        ee_sigma            sigma = mu_Z sigma_theory = 1841 +- 30 pb beside the prediction (CMS and ATLAS
                            are left for the end of the talk)

Numbers: presentation/data/zee_fit.json (extract_zee_fit.py, the z-ee tarball) and zee_selection.json
(extract_zee_selection.py, the SingleElectron NanoAOD the z-ee team processed, a 10-file sample for shapes and
fractions); asserted below. Physics honesty: the data never move; the prediction moves by the real templates
(linear between nominal and +-1 sigma, display only) and ends on TRExFitter's own post-fit yields. The side view,
the calorimeter panels and the track curvature are schematic; eta, p_T, phi and the ID variables are real.
Reserved areas (06 A1): nothing above y = 2.7 and nothing at x < -5.85, y > 0.22 on a held frame.
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
    DOWN, LEFT, RIGHT, UP, Circle, Create, DashedLine, DashedVMobject, DecimalNumber, Dot, FadeIn, FadeOut, Flash,
    GrowFromCenter, GrowFromEdge, LaggedStart, Line, Polygon, Rectangle, ReplacementTransform, Scene, Succession,
    Transform, ValueTracker, VGroup, always_redraw, rate_functions,
)

from style.bnd_style import *  # noqa: E402,F401,F403
from style.bnd_style import _cap_scale  # noqa: E402
import s2_map as M  # noqa: E402
from s2_pipeline import camera, icon_funnel, masked_play  # noqa: E402
from s2_trigger import cross, tick  # noqa: E402

EASE = rate_functions.ease_in_out_sine
ACC = DETECTOR_ACCENT
GREEN_LINE = CHANNEL_LINE["ee"]

# ---------------------------------------------------------------------------
# frozen data + anchors (a mismatch stops the render; never adjusted here)
# ---------------------------------------------------------------------------
FIT = load_data("zee_fit")
SEL = load_data("zee_selection")

EDGES = np.asarray(FIT["edges"], dtype=float)
XC = 0.5 * (EDGES[:-1] + EDGES[1:])
DATA = np.asarray(FIT["data"], dtype=float)
ORDER_S = list(FIT["stack_order"])                       # bottom-up
NOM = {s: np.asarray(FIT["prefit"]["samples"][s], dtype=float) for s in ORDER_S}
POST = {s: np.asarray(FIT["postfit"]["samples"][s], dtype=float) for s in ORDER_S}
TPL = FIT["templates"]
NPD = {n["name"]: n for n in FIT["nps"]}
GAM = np.asarray([g["value"] for g in FIT["gammas"]], dtype=float)
MU, SIG, CNT, TH = FIT["mu"], FIT["sigma"], FIT["counting"], FIT["theory"]
EV = SEL["event"]
FAKE = SEL["fake_candidate"]
WPS = list(SEL["wp_names"])
M_EDGES = np.asarray(SEL["m_edges"], dtype=float)
PT_EDGES = np.asarray(SEL["pt_edges"], dtype=float)
PRE_RATIO = float(DATA.sum() / sum(NOM[s].sum() for s in ORDER_S))

assert int(DATA.sum()) == 6320097 and len(DATA) == 30 and EDGES[0] == 60 and EDGES[-1] == 120
assert ORDER_S == ["WJets", "TTbar", "Diboson", "DYtautau", "DYee"]
assert round(MU["value"], 3) == 0.942 and round(MU["err_up"], 3) == 0.015 and round(MU["err_down"], 3) == 0.015
assert round(SIG["value"]) == 1841 and round(SIG["err"]) == 30 and round(SIG["ref"], 1) == 1954.1
assert round(CNT["sigma"]) == 1896 and round(CNT["n_bkg"]) == 44357 and round(CNT["eff_acc"], 3) == 0.202
assert round(PRE_RATIO, 3) == 0.970
assert (round(TH["unc_up"]), round(TH["unc_down"])) == (15, 21)
assert EV["run"] == 279024 and EV["event"] == 104459341 and round(EV["m_ee"], 1) == 91.4
assert [round(e["pt"], 1) for e in EV["electrons"]] == [45.5, 41.4] and [e["charge"] for e in EV["electrons"]] == [-1, 1]
assert round(FAKE["hoe"], 2) == 0.58 and round(FAKE["relIso"], 2) == 2.24
assert WPS == ["none", "veto", "loose", "medium", "tight"]
assert [round(SEL["wp"][w]["eff_rel_none"], 2) for w in WPS] == [1.0, 0.83, 0.77, 0.64, 0.5]
assert SEL["totals"]["window"] == 6320097 and int(np.argmax(SEL["pt"]["in_window"])) == 8
assert round(NPD["ElectronID"]["pull"], 2) == 0.59 and round(NPD["L1Prefiring"]["pull"], 2) == 2.02

# ---------------------------------------------------------------------------
# layout (scene units)
# ---------------------------------------------------------------------------
L = dict(
    det=(1.0, (-1.9, -0.25)),                       # the event view: the inner detector of a 2.4x slice
    det_park=(0.28, (-5.0, -3.0)),
    lad_x=3.45, lad_tick_x=3.05, lad_y=(2.15, 1.5, 0.85, 0.2, -0.45, -1.1), lad_h=0.24,
    funnel=(5.0, 0.45), funnel_w=2.0,
    pt_c=(-1.35, 0.55), pt_wh=(5.6, 2.7),
    sv_c=(-2.1, 0.3), sv_u=0.8,                     # side view centre, scene units per metre
    id_c=((-3.35, 0.75), (0.75, 0.75)), id_cell=0.24,
    wp_c=(-2.85, 0.8), wp_wh=(4.1, 2.5), to_c=(1.35, 0.8), to_wh=(2.1, 2.5),
    sel_main=(-1.3, 0.95), sel_wh=(6.0, 2.4), sel_ratio=(-1.3, -1.25),
    formula_x=4.7,
    fit_main=(1.85, 1.0), fit_wh=(6.4, 2.35), fit_ratio=(1.85, -1.2),
    key_y=-2.72, key_left_sel=-4.0, key_left_fit=-0.85,
    mu_sl=(-3.85, 1.2), model_at=(-3.85, 2.2), pulls_x=-3.75, pulls_top=0.0,
)
Y_EXP = (2.5, 6.5)
RATIO_RANGE, RATIO_TICKS = (0.84, 1.12, 0.1), (0.9, 1.0, 1.1)
DOT_R, RDOT_R = 0.04, 0.034
COLOURS = {"DYee": SAMPLE["DYee"], "DYtautau": SAMPLE["DYtautau"], "TTbar": SAMPLE["TTbar"],
           "Diboson": SAMPLE["WW"], "WJets": SAMPLE["Fakes"]}
NP_ANIM = ("Lumi", "ElectronID", "QCDScale", "PS_FSR")
NP_MORE = ("Pileup", "L1Prefiring", "ElectronRECO", "XS_DYtautau")
NP_ROWS = NP_ANIM + NP_MORE
CUT_ROWS = (r"\mathrm{HLT}:\ p_T^{\,e} > 27\ \mathrm{GeV}", r"e^{+}e^{-}", r"p_T > 20\ \mathrm{GeV}", r"|\eta| < 2.5",
            r"\mathrm{ID}:\ \mathrm{medium}", r"60 < m_{ee} < 120\ \mathrm{GeV}")


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def P3(p):
    p = np.asarray(p, dtype=float)
    return p if p.shape == (3,) else np.array([p[0], p[1], 0.0])


def tex_h(expr: str, h: float, color=INK):
    return mathtex(expr, color=color).scale(_cap_scale("tex", h))


def text_h(s: str, h: float, color=INK):
    return text(s, color=color).scale(_cap_scale("text", h))


def tex_int(n) -> str:
    return f"{int(round(n)):,}".replace(",", "{,}")


def clear_updaters(mobs):
    for m in mobs:
        for f in m.get_family():
            f.clear_updaters()


def live_mobjects(scene, keep=()):
    ids = {id(k) for k in keep}
    return [m for m in scene.mobjects if id(m) not in ids and not isinstance(m, ValueTracker)
            and m.family_members_with_points()]


def reset_to(scene, state: dict, order) -> None:
    scene.remove(*scene.mobjects)
    add_state(scene, state, order)


# ---------------------------------------------------------------------------
# the map (s2_map geometry, ee colours)
# ---------------------------------------------------------------------------

def ee_map() -> dict:
    return M.map_state(M.ALL, "ee")


def zoom_into(scene, st: dict, about, detail, run_time=1.5):
    """The camera flies into a node of the map: the map leaves through the zoom, ``detail`` grows out of the node."""
    anims = [Transform(st[k], camera(st[k].copy(), M.ZOOM, about, M.DETAIL).set_opacity(0.0), remover=True)
             for k in M.ORDER_MAP]
    anims.append(ReplacementTransform(camera(detail.copy(), 1.0 / M.ZOOM, M.DETAIL, about).set_opacity(0.0), detail))
    masked_play(scene, *anims, run_time=run_time, rate_func=EASE)


def zoom_out_to_map(scene, about, run_time=1.5) -> dict:
    """Everything on screen shrinks into the node ``about``; the ee map comes back around it."""
    content = live_mobjects(scene)
    clear_updaters(content)
    end = ee_map()
    anims = [ReplacementTransform(camera(end[k].copy(), M.ZOOM, about, M.DETAIL).set_opacity(0.0), end[k])
             for k in M.ORDER_MAP]
    g = VGroup(*content)
    anims.append(Transform(g, camera(g.copy(), 1.0 / M.ZOOM, M.DETAIL, about).set_opacity(0.0), remover=True))
    masked_play(scene, *anims, run_time=run_time, rate_func=EASE)
    reset_to(scene, end, M.ORDER_MAP)
    return end


# ---------------------------------------------------------------------------
# the cut list
# ---------------------------------------------------------------------------

def cut_row(i: int) -> VGroup:
    lab = tex_h(CUT_ROWS[i], L["lad_h"])
    lab.shift(np.array([L["lad_x"], L["lad_y"][i], 0.0]) - lab.get_left())
    tk = tick((L["lad_tick_x"], L["lad_y"][i]), s=0.13)
    g = VGroup(lab, tk)
    g.lab, g.tick = lab, tk
    return g


def add_cut(scene, i: int, run_time=0.7) -> VGroup:
    r = cut_row(i)
    scene.play(FadeIn(r.lab, shift=LEFT * 0.2), run_time=0.45 * run_time / 0.7)
    scene.play(Create(r.tick), run_time=0.35)
    return r


def big_funnel(center=None, width=None) -> VGroup:
    f = icon_funnel(GREEN_LINE, CHANNEL["ee"])
    f.scale((width or L["funnel_w"]) / f.width).move_to(P3(center or L["funnel"]))
    return f


def funnel_spout(f) -> np.ndarray:
    return f.get_bottom() + UP * 0.05


# ---------------------------------------------------------------------------
# 3-04: the real event
# ---------------------------------------------------------------------------

R_INNER = 2.4 * R_DET          # the event view: a slice this large, drawn only up to the HCAL (outer radius 2.7)
INNER = ("beampipe", "tracker", "tob", "ecal", "band", "hcal")


def event_display() -> dict:
    """The real event on the inner detector (tracker, ECAL, HCAL) of a 2.4x slice at L['det'].
    Parts: det (the drawn layers), tracks, hits, deps, labs, stamp, ring."""
    geo = CMSSlice(center=DET_CENTER, radius=R_INNER)
    det = VGroup(*[g for name, g in geo.parts if name in INNER])
    tracks, hits, deps, labs = VGroup(), VGroup(), VGroup(), VGroup()
    for e in EV["electrons"]:
        sig = signature(geo, "e", float(e["phi"]), charge=int(e["charge"]), kappa=kappa_from_pt(e["pt"]), sw=5.0)
        trk = sig[0]
        tracks.add(trk)
        hits.add(tracker_hits(geo, trk, color=CHANNEL["ee"], radius=0.035))
        deps.add(VGroup(*sig[1:]))
    for m in (det, tracks, hits, deps):
        place(m, L["det"], DET_CENTER)
    spots = {-1: (-5.0, -2.55), 1: (1.4, 2.1)}
    for e in EV["electrons"]:
        q = int(e["charge"])
        sym = tex_h(r"e^{-}" if q < 0 else r"e^{+}", 0.3, color=GREEN_LINE)
        pt = tex_h(rf"p_T = {e['pt']:.1f}\ \mathrm{{GeV}}", 0.2)
        grp = VGroup(sym, pt).arrange(DOWN, buff=0.12, aligned_edge=LEFT).move_to(P3(spots[q]))
        grp.sym, grp.pt = sym, pt
        labs.add(grp)
    stamp = text_h(f"Run {EV['run']}   Event {EV['event']}", 0.13, color=GREY).move_to(P3((-1.9, -3.2)))
    lead_dep = deps[0]
    ring = Circle(radius=0.55, stroke_color=col(ACC), stroke_width=4.0).move_to(lead_dep.get_center())
    return {"det": det, "tracks": tracks, "hits": hits, "deps": deps, "labs": labs, "stamp": stamp, "ring": ring}


# ---------------------------------------------------------------------------
# 3-05: electron p_T
# ---------------------------------------------------------------------------

def pt_axes() -> DataAxes:
    top = 1.1 * max(np.asarray(SEL["pt"]["in_window"]) + np.asarray(SEL["pt"]["outside"]))
    dax = DataAxes([0, 100, 20], [0, top, top / 4], *L["pt_wh"], x_ticks=[0, 20, 40, 60, 80, 100], y_ticks=[],
                   show_y_labels=False, tick_label_h=0.18, x_title=r"p_T\ [\mathrm{GeV}]", title_h=0.22, title_buff=0.14)
    dax.move_frame_to(P3(L["pt_c"]))
    dax.add(tex_h(r"\mathrm{electrons}", 0.2).rotate(np.pi / 2).next_to(dax.frame, LEFT, buff=0.2))
    return dax


def pt_stack(dax, frac=1.0) -> VGroup:
    out = frac * np.asarray(SEL["pt"]["outside"], dtype=float)
    inn = frac * np.asarray(SEL["pt"]["in_window"], dtype=float)
    return stack_hist(dax, PT_EDGES, [("out", out, GREY), ("in", inn, CHANNEL["ee"])])


# ---------------------------------------------------------------------------
# 3-06: side view (r-z), CMS proportions in metres, schematic
# ---------------------------------------------------------------------------
TRK_R, TRK_Z = 1.1, 2.8
EB_R = (1.29, 1.75)
EE_Z, EE_R = (3.15, 3.9), (0.31, 1.62)
ETA_TRK = 2.5


def sv(z, r):
    u = L["sv_u"]
    c = P3(L["sv_c"])
    return c + np.array([z * u, r * u, 0.0])


def side_view() -> dict:
    ink = col(INK)
    th25 = eta_theta(ETA_TRK)
    r_end = TRK_Z * math.tan(th25)
    beam = Line(sv(-4.25, 0), sv(4.25, 0), stroke_color=col(GREY), stroke_width=2.0)
    trk = VGroup()
    for up in (1, -1):
        pts = [sv(0, 0), sv(TRK_Z, up * r_end), sv(TRK_Z, up * TRK_R), sv(-TRK_Z, up * TRK_R), sv(-TRK_Z, up * r_end)]
        trk.add(Polygon(*pts, fill_color=col(DETECTOR["tracker"]["fill"]), fill_opacity=1.0,
                        stroke_color=col(DETECTOR["tracker"]["stroke"]), stroke_width=1.5))
    ecal = VGroup()
    style = dict(fill_color=col(DETECTOR["ecal"]["fill"]), fill_opacity=1.0, stroke_color=col(DETECTOR["ecal"]["stroke"]),
                 stroke_width=1.5)
    zb = EB_R[0] / math.tan(eta_theta(1.479))
    for up in (1, -1):
        ecal.add(Polygon(sv(-zb, up * EB_R[0]), sv(zb, up * EB_R[0]), sv(zb, up * EB_R[1]), sv(-zb, up * EB_R[1]), **style))
        for side in (1, -1):
            ecal.add(Polygon(sv(side * EE_Z[0], up * EE_R[0]), sv(side * EE_Z[1], up * EE_R[0]),
                             sv(side * EE_Z[1], up * EE_R[1]), sv(side * EE_Z[0], up * EE_R[1]), **style))
    lines = VGroup()
    z_to = 4.15
    for side in (1, -1):
        for up in (1, -1):
            lines.add(DashedLine(sv(0, 0), sv(side * z_to, up * z_to * math.tan(th25)), dash_length=0.1,
                                 stroke_color=col(ACC), stroke_width=3.0))
    lab = tex_h(r"|\eta| = 2.5", 0.22, color=ACC).next_to(sv(z_to, z_to * math.tan(th25)), RIGHT, buff=0.08)
    ip = Dot(sv(0, 0), radius=0.05, color=ink)
    return {"beam": beam, "trk": trk, "ecal": ecal, "eta_lines": lines, "eta_lab": lab, "ip": ip}


def sv_particle(eta, up, color, dashed=False, with_hits=True) -> VGroup:
    """A straight particle in the side view from the vertex to the ECAL (barrel or endcap), its deposit and hits."""
    th = eta_theta(abs(eta))
    side = 1 if eta >= 0 else -1
    tan = math.tan(th)
    if abs(eta) < 1.479:
        r0, r1 = EB_R
        p0, p1 = (side * r0 / tan, up * r0), (side * r1 / tan, up * r1)
    else:
        z0, z1 = EE_Z
        p0, p1 = (side * z0, up * z0 * tan), (side * z1, up * z1 * tan)
    cc = col(color)
    if dashed:
        line = DashedLine(sv(0, 0), sv(*p0), dash_length=0.08, stroke_color=cc, stroke_width=3.5)
    else:
        line = Line(sv(0, 0), sv(*p0), stroke_color=cc, stroke_width=4.0)
    d = np.array([math.cos(th) * side, math.sin(th) * up])
    n = np.array([-d[1], d[0]]) * 0.045
    a, b = np.array(p0), np.array(p1)
    dep = Polygon(sv(*(a + n)), sv(*(b + n)), sv(*(b - n)), sv(*(a - n)), fill_color=cc, fill_opacity=0.95,
                  stroke_color=darken(cc, 0.35), stroke_width=1.2)
    hits = VGroup()
    if with_hits:
        for r in (0.25, 0.45, 0.7, 0.95):
            z = side * r / tan
            if abs(z) < TRK_Z:
                hits.add(Dot(sv(z, up * r), radius=0.035, color=col(DETECTOR["tracker"]["stroke"])))
    g = VGroup(line, hits, dep)
    g.line, g.hits, g.dep = line, hits, dep
    return g


# ---------------------------------------------------------------------------
# 3-07: electron vs jet (eta-phi calorimeter picture, schematic geometry, real variables)
# ---------------------------------------------------------------------------

def id_panel(kind: str, center) -> dict:
    c = P3(center)
    n, w = 9, L["id_cell"]
    rng = np.random.default_rng(7 if kind == "e" else 11)
    E = np.zeros((n, n))
    if kind == "e":
        E[4, 4] = 1.0
        for (i, j), v in {(3, 4): 0.35, (5, 4): 0.3, (4, 3): 0.22, (4, 5): 0.25, (3, 3): 0.08, (5, 5): 0.1}.items():
            E[i, j] = v
    else:
        for i in range(n):
            for j in range(n):
                r2 = (i - 4.2) ** 2 + (j - 3.8) ** 2
                E[i, j] = 0.75 * math.exp(-r2 / 5.0) * rng.uniform(0.5, 1.0)
        E[E < 0.12] = 0.0
        for i, j, v in ((1, 6, 0.3), (7, 7, 0.25), (6, 1, 0.2)):
            E[i, j] = v
    colour = CHANNEL["ee"] if kind == "e" else SLATE
    cells = VGroup()
    for i in range(n):
        for j in range(n):
            p = c + np.array([(j - 4) * w, (4 - i) * w, 0.0])
            a = float(E[i, j])
            fill = mix(WHITE, col(colour).to_hex(), 0.06 + 0.85 * a) if a > 0 else WHITE
            cells.add(Rectangle(width=w, height=w, fill_color=col(fill), fill_opacity=1.0,
                                stroke_color=col(LIGHT_GREY), stroke_width=1.0).move_to(p))
    cone = DashedVMobject(Circle(radius=1.0, arc_center=c, stroke_color=col(GREY), stroke_width=2.5), num_dashes=30)
    trk_pts = [(0.0, 0.0)] if kind == "e" else [(-0.1, 0.12), (0.45, -0.35), (-0.55, 0.5), (0.62, 0.52), (0.1, -0.7)]
    trks = VGroup(*[Dot(c + np.array([x, y, 0.0]), radius=0.06, color=col(INK)) for x, y in trk_pts])
    # H / E: the energy behind the cluster
    rec = EV["electrons"][0] if kind == "e" else FAKE
    base_y = c[1] - n * w / 2
    bx = c[0] + n * w / 2 + 0.38
    hE = 1.3
    eb = Rectangle(width=0.22, height=hE, fill_color=col(colour), fill_opacity=0.9, stroke_width=0).move_to(
        np.array([bx, base_y + hE / 2, 0.0]))
    hh = max(0.02, hE * float(rec["hoe"]))
    hb = Rectangle(width=0.22, height=hh, fill_color=col(DEPOSIT["hcal"]), fill_opacity=0.95, stroke_width=0).move_to(
        np.array([bx + 0.3, base_y + hh / 2, 0.0]))
    base = Line(np.array([bx - 0.2, base_y, 0]), np.array([bx + 0.5, base_y, 0]), stroke_color=col(INK), stroke_width=2.0)
    blab = VGroup(tex_h("E", 0.15).next_to(eb, DOWN, buff=0.1), tex_h("H", 0.15).next_to(hb, DOWN, buff=0.1))
    blab[1].set_y(blab[0].get_y())
    bars = VGroup(base, eb, hb, blab)
    sym = tex_h(r"e" if kind == "e" else r"\mathrm{jet}", 0.3, color=GREEN_LINE if kind == "e" else GREY)
    sym.move_to(c + np.array([0.0, n * w / 2 + 0.3, 0.0]))
    vals = VGroup(
        tex_h(rf"\sigma_{{i\eta i\eta}} = {rec['sieie']:.4f}", 0.19),
        tex_h(rf"H/E = {rec['hoe']:.2f}", 0.19),
        tex_h(rf"I_{{\mathrm{{rel}}}} = {rec['relIso']:.2f}", 0.19),
    ).arrange(DOWN, buff=0.16, aligned_edge=LEFT)
    vals.next_to(np.array([c[0] + 0.2, base_y - 0.25, 0.0]), DOWN, buff=0.0)
    return {"cells": cells, "cone": cone, "trks": trks, "bars": bars, "sym": sym, "vals": vals}


# ---------------------------------------------------------------------------
# 3-08: working points
# ---------------------------------------------------------------------------

def wp_axes() -> DataAxes:
    dax = DataAxes([60, 120, 10], [0.5, 5.5, 1], *L["wp_wh"], x_ticks=[70, 80, 90, 100, 110], y_ticks=[1, 2, 3, 4, 5],
                   y_log=True, tick_label_h=0.16, x_title=r"m_{ee}\ [\mathrm{GeV}]", title_h=0.2, title_buff=0.12)
    return dax.move_frame_to(P3(L["wp_c"]))


def wp_hists(dax, w) -> VGroup:
    os_ = np.asarray(SEL["wp"][w]["os"], dtype=float)
    ss = np.asarray(SEL["wp"][w]["ss"], dtype=float)
    g = VGroup(step_hist(dax, M_EDGES, os_, color=CHANNEL["ee"], stroke_width=2.5, fill_opacity=0.3),
               step_hist(dax, M_EDGES, ss, color=GREY, stroke_width=2.5, fill_opacity=0.45))
    return g


def wp_key(dax) -> VGroup:
    rows = VGroup(
        VGroup(Rectangle(width=0.3, height=0.18, fill_color=col(CHANNEL["ee"]), fill_opacity=0.5,
                         stroke_color=col(CHANNEL["ee"]), stroke_width=2), tex_h(r"e^{+}e^{-}", 0.17)).arrange(RIGHT, buff=0.12),
        VGroup(Rectangle(width=0.3, height=0.18, fill_color=col(GREY), fill_opacity=0.6, stroke_color=col(GREY),
                         stroke_width=2), tex_h(r"e^{\pm}e^{\pm}", 0.17)).arrange(RIGHT, buff=0.12),
    ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
    rows.next_to(dax.c2p(62.0, 10 ** 5.4), DOWN + RIGHT, buff=0.0)
    return rows


def wp_name(w) -> VGroup:
    s = {"none": "no ID", "veto": "veto", "loose": "loose", "medium": "medium", "tight": "tight"}[w]
    return text_h(s, 0.2, color=GREEN_LINE if w == "medium" else INK).move_to(P3((L["wp_c"][0] + 1.3, 2.4)))


def purity(w) -> float:
    return 1.0 - float(SEL["wp"][w]["ss_over_os"])


def to_axes() -> DataAxes:
    dax = DataAxes([0.4, 1.0, 0.2], [0.88, 1.0, 0.04], *L["to_wh"], x_ticks=[0.4, 0.6, 0.8, 1.0], y_ticks=[0.9, 0.94, 0.98],
                   x_fmt="{:.1f}", y_fmt="{:.2f}", tick_label_h=0.15)
    dax.move_frame_to(P3(L["to_c"]))
    dax.add(text_h("efficiency", 0.15).next_to(dax.x_labels, DOWN, buff=0.14))
    dax.add(text_h("purity", 0.15).rotate(np.pi / 2).next_to(dax.y_labels, LEFT, buff=0.14))
    return dax


def to_point(dax, w, lit=False) -> VGroup:
    x, y = float(SEL["wp"][w]["eff_rel_none"]), purity(w)
    d = Dot(dax.c2p(x, y), radius=0.07 if lit else 0.055, color=col(GREEN_LINE if lit else INK))
    return d


# ---------------------------------------------------------------------------
# the m_ee plots (selection layout and fit layout)
# ---------------------------------------------------------------------------

def main_axes(kind: str, x_labels: bool) -> DataAxes:
    c, wh = (L["sel_main"], L["sel_wh"]) if kind == "sel" else (L["fit_main"], L["fit_wh"])
    dax = DataAxes([60, 120, 10], [*Y_EXP, 1], *wh, x_ticks=[70, 80, 90, 100, 110], y_ticks=[3, 4, 5, 6], y_log=True,
                   show_x_labels=x_labels, tick_label_h=0.17, x_title=r"m_{ee}\ [\mathrm{GeV}]" if x_labels else None,
                   y_title=r"\mathrm{events}\,/\,2\,\mathrm{GeV}", title_h=0.19, title_buff=0.14)
    return dax.move_frame_to(P3(c))


def ratio_axes(dax, kind: str) -> DataAxes:
    c = L["sel_ratio"] if kind == "sel" else L["fit_ratio"]
    rax = DataAxes(dax.x_range, list(RATIO_RANGE), dax.x_length, 0.85, x_ticks=dax.x_ticks, y_ticks=list(RATIO_TICKS),
                   y_fmt="{:.1f}", tick_label_h=0.15, x_title=r"m_{ee}\ [\mathrm{GeV}]", y_title=r"\mathrm{data/pred.}",
                   title_h=0.17, title_buff=0.12)
    rax.move_frame_to(P3(c))
    rax.shift(RIGHT * (dax.c2p(60, dax.y_base)[0] - rax.c2p(60, rax.y_base)[0]))
    ref = DashedVMobject(Line(rax.c2p(60, 1.0), rax.c2p(120, 1.0), stroke_color=col(GREY), stroke_width=2.0), num_dashes=40)
    g = VGroup(rax, ref)
    g.rax, g.ref = rax, ref
    return g


def ratio_dots(rax: DataAxes, total) -> VGroup:
    r = DATA / np.asarray(total, dtype=float)
    e = np.sqrt(DATA) / np.asarray(total, dtype=float)
    return VGroup(*[VGroup(data_errorbar(rax, XC[i], r[i], e[i], color=SAMPLE["Data"], stroke_width=1.8),
                           data_dot(rax, XC[i], r[i], color=SAMPLE["Data"], radius=RDOT_R)) for i in range(30)])


def layers_of(counts: dict):
    return [(s, np.asarray(counts[s], dtype=float), COLOURS[s]) for s in ORDER_S]


def model(mu=1.0, theta=None, gam=None) -> dict:
    """The prediction per sample at (mu, theta, gamma): nominal + |t| (template - nominal), mu on Z -> ee only.
    Display interpolation between the fit's own templates; the end state is TRExFitter's post-fit yield."""
    theta = theta or {}
    out = {}
    for s in ORDER_S:
        v = NOM[s].copy()
        for n, t in theta.items():
            if abs(t) < 1e-12 or n not in TPL[s]:
                continue
            tpl = np.asarray(TPL[s][n]["up" if t > 0 else "down"], dtype=float)
            v = v + abs(t) * (tpl - NOM[s])
        if s == "DYee":
            v = v * mu
        if gam is not None:
            v = v * gam
        out[s] = v
    return out


def total_of(counts: dict):
    return np.sum([counts[s] for s in ORDER_S], axis=0)


def data_dots(dax) -> VGroup:
    return VGroup(*[data_dot(dax, XC[i], DATA[i], color=SAMPLE["Data"], radius=DOT_R) for i in range(30)])


def data_bars(dax, frac=1.0) -> VGroup:
    return VGroup(*[data_bar(dax, XC[i], max(DATA[i] * frac, 1.0), 1.0, color=SAMPLE["Data"], fill_opacity=0.55,
                             stroke_width=1.0) for i in range(30)])


def ee_key(left_x) -> VGroup:
    entries = [(r"\mathrm{Data}", SAMPLE["Data"], "dot"), (r"Z\to ee", COLOURS["DYee"]),
               (r"\tau\tau,\ t\bar{t},\ VV,\ W\!+\!\mathrm{jets}", COLOURS["TTbar"])]
    key = colour_key(entries, label_h=0.18)
    row = key.rows[2]
    sw = row[1]
    names = ("DYtautau", "TTbar", "Diboson", "WJets")
    w, h = sw.width, sw.height
    stripes = VGroup(*[Rectangle(width=w / len(names), height=h, stroke_width=0, fill_color=col(COLOURS[n]),
                                 fill_opacity=0.95) for n in names]).arrange(RIGHT, buff=0).move_to(sw)
    frame = Rectangle(width=w, height=h, fill_opacity=0.0, stroke_color=darken(COLOURS["TTbar"], 0.25),
                      stroke_width=1.0).move_to(sw)
    row.submobjects[1] = VGroup(stripes, frame)
    key.rows.arrange(RIGHT, buff=0.4)
    key.rows.shift(np.array([left_x, L["key_y"], 0.0]) - key.rows.get_left())
    return key


def plot_state(kind: str, counts: dict | None, flat=False) -> dict:
    """dax, stack, data (dots), ratio (axes + ref), rdots, rlabel, key: the m_ee plot of a layout."""
    dax = main_axes(kind, x_labels=False)
    if counts is None:
        counts = NOM
    tot = total_of(counts)
    if flat:
        counts = {s: np.full(30, counts[s].mean()) for s in ORDER_S}
    st = {"dax": dax, "stack": stack_hist(dax, EDGES, layers_of(counts))}
    if flat:
        st["data"] = VGroup(*[data_dot(dax, XC[i], DATA.mean(), color=SAMPLE["Data"], radius=DOT_R) for i in range(30)])
    else:
        st["data"] = data_dots(dax)
    ra = ratio_axes(dax, kind)
    st["ratio"] = ra
    if flat:
        r = DATA.sum() / tot.sum()
        st["rdots"] = VGroup(*[VGroup(data_errorbar(ra.rax, XC[i], r, 0.0, color=SAMPLE["Data"], stroke_width=1.8),
                                      data_dot(ra.rax, XC[i], r, color=SAMPLE["Data"], radius=RDOT_R)) for i in range(30)])
    else:
        st["rdots"] = ratio_dots(ra.rax, tot)
    st["key"] = ee_key(L["key_left_sel"] if kind == "sel" else L["key_left_fit"])
    return st


PLOT_ORDER = ("dax", "stack", "data", "ratio", "rdots", "key")


# ---------------------------------------------------------------------------
# 3-11 end state (the seam between the two scenes)
# ---------------------------------------------------------------------------

def formula_lines() -> dict:
    x = L["formula_x"]
    l1 = tex_h(r"N_{\mathrm{data}} = \sigma \cdot L \cdot \varepsilon A + N_{\mathrm{bkg}}", 0.27).move_to(P3((x, 2.05)))
    l2 = tex_h(r"\sigma = \frac{N_{\mathrm{data}} - N_{\mathrm{bkg}}}{\varepsilon A \cdot L}", 0.27).move_to(P3((x, 1.0)))
    l3 = tex_h(rf"= \frac{{{tex_int(CNT['n_data'])} - {tex_int(CNT['n_bkg'])}}}"
               rf"{{{CNT['eff_acc']:.3f} \times {CNT['lumi_pb'] / 1000:.2f}\ \mathrm{{fb}}^{{-1}}}}", 0.24).move_to(P3((x, -0.1)))
    l4 = tex_h(rf"= {CNT['sigma']:.0f}\ \mathrm{{pb}}", 0.3).move_to(P3((x, -1.0)))
    l4.shift(RIGHT * (l3.get_left()[0] - l4.get_left()[0]))
    return {"f1": l1, "f2": l2, "f3": l3, "f4": l4}


def counting_state() -> dict:
    st = plot_state("sel", NOM, flat=True)
    lab = tex_h(f"{PRE_RATIO:.3f}", 0.19).next_to(st["ratio"].rax.c2p(120, PRE_RATIO), RIGHT, buff=0.18)
    st["rlabel"] = lab
    st.update(formula_lines())
    return st


ORDER_COUNT = PLOT_ORDER + ("rlabel", "f1", "f2", "f3", "f4")


# ---------------------------------------------------------------------------
# the fit panels
# ---------------------------------------------------------------------------

def mu_slider() -> Slider:
    sl = slider(0.8, 1.2, 1.0, ref=1.0, length=2.6, ticks=[0.8, 0.9, 1.0, 1.1, 1.2], fmt="{:.1f}", tick_label_h=0.15,
                marker_r=0.08)
    sl.shift(P3(L["mu_sl"]) - sl.axis.get_center())
    lab = tex_h(r"\mu_{Z}", 0.26).next_to(sl.axis, LEFT, buff=0.25)
    g = VGroup(sl, lab)
    g.sl, g.lab = sl, lab
    return g


def pulls_frame(n_rows: int, post=False):
    names = [NPD[n]["label"] for n in NP_ROWS[:n_rows]]
    p = [NPD[n]["pull"] if post else 0.0 for n in NP_ROWS[:n_rows]]
    c = [NPD[n]["constraint"] if post else None for n in NP_ROWS[:n_rows]]
    pp = pull_plot(names, p, c, x_range=(-3.0, 3.0), x_length=2.8, row_h=0.36, label_h=0.16, label_buff=0.18)
    pp.shift(np.array([L["pulls_x"], L["pulls_top"], 0.0]) - np.array([pp.baseline.get_center()[0], pp.band.get_top()[1], 0.0]))
    return pp


def pull_dot(pp, row: int, value: float, err=None) -> VGroup:
    y = pp.rows[row][2].get_center()[1]
    x = pp.x_of(value)
    dot = Dot(np.array([x, y, 0.0]), radius=0.055, color=col(INK))
    if err:
        bar = Line(np.array([pp.x_of(value - err), y, 0]), np.array([pp.x_of(value + err), y, 0]), stroke_color=col(INK),
                   stroke_width=2.5)
    else:
        bar = Line(np.array([x, y, 0]), np.array([x, y, 0]), stroke_color=col(INK), stroke_width=0)
    return VGroup(bar, dot)


def model_tex(with_theta: bool):
    s = r"\nu_i = \mu\, s_i(\vec\theta) + b_i(\vec\theta)" if with_theta else r"\nu_i = \mu\, s_i + b_i"
    return tex_h(s, 0.26).move_to(P3(L["model_at"]))


def likelihood_tex():
    return tex_h(r"L(\mu, \vec\theta) = \prod_{i} \mathrm{Pois}\big(n_i \,\big|\, \nu_i(\mu, \vec\theta)\big)\,"
                 r"\prod_{j} G(\theta_j)", 0.24).move_to(P3((L["fit_main"][0], -3.62)))


# ---------------------------------------------------------------------------
# 3-16: the result frame (the tau tau chapter's, 5-33)
# ---------------------------------------------------------------------------
SIGF = dict(lo=1750.0, hi=2250.0, x0=-2.6, length=8.0, y=-2.9, height=4.9, rows=dict(ee=4.25))


def sigma_panel() -> VGroup:
    S = SIGF
    lo, hi, x0, Lx, y0, H = S["lo"], S["hi"], S["x0"], S["length"], S["y"], S["height"]
    ink = col(INK)

    def xs(v):
        return x0 + (float(v) - lo) / (hi - lo) * Lx

    axis = Line([xs(lo), y0, 0], [xs(hi), y0, 0], stroke_color=ink, stroke_width=2.5)
    ticks, tick_labels = VGroup(), VGroup()
    for v in np.arange(lo, hi + 1, 50.0):
        p = np.array([xs(v), y0, 0.0])
        big = abs(v % 100) < 1e-9
        ticks.add(Line(p, p + DOWN * (0.1 if big else 0.06), stroke_color=ink, stroke_width=2.5))
        if big:
            tick_labels.add(text_h(f"{v:.0f}", 0.15).next_to(p + DOWN * 0.1, DOWN, buff=0.08))
    unit_ = tex_h(r"\mathrm{pb}", 0.17).next_to(tick_labels[-1], RIGHT, buff=0.3)
    pred = float(TH["value"])
    lo_b, hi_b = xs(pred - float(TH["unc_down"])), xs(pred + float(TH["unc_up"]))
    band = Rectangle(width=hi_b - lo_b, height=H, fill_color=col(THEORY), fill_opacity=0.15, stroke_width=0).move_to(
        [0.5 * (lo_b + hi_b), y0 + H / 2, 0.0])
    theory = DashedLine([xs(pred), y0, 0], [xs(pred), y0 + H, 0], dash_length=0.1, stroke_color=col(THEORY), stroke_width=3.0)
    th_lab = tex_h(rf"{pred:.1f}^{{+{float(TH['unc_up']):.0f}}}_{{-{float(TH['unc_down']):.0f}}}", 0.17, color=THEORY)
    th_lab.next_to(theory.get_end(), LEFT, buff=0.14)
    g_ee = col(GREEN_LINE)
    ye = y0 + S["rows"]["ee"]
    point = VGroup(Line([xs(SIG["value"] - SIG["err"]), ye, 0], [xs(SIG["value"] + SIG["err"]), ye, 0], stroke_color=g_ee,
                        stroke_width=5.0),
                   Dot([xs(SIG["value"]), ye, 0], radius=0.1, color=g_ee))
    label = tex_h(rf"\sigma_{{60\text{{--}}120}} = {SIG['value']:.0f} \pm {SIG['err']:.0f}\ \mathrm{{pb}}", 0.3).move_to(
        P3((2.9, 2.3)))
    g = VGroup(axis, ticks, tick_labels, unit_, band, theory, th_lab, point, label)
    g.frame, g.band, g.theory, g.th_lab, g.point, g.label = \
        VGroup(axis, ticks, tick_labels, unit_), band, theory, th_lab, point, label
    return g


def mu_line(at=(-3.7, 2.3), h=0.28):
    return tex_h(rf"\mu_Z = {MU['value']:.3f} \pm {MU['err_up']:.3f}", h).move_to(P3(at))


# ===========================================================================
# the scenes
# ===========================================================================

class ZeeSelection(Scene):
    """ee_map_selection ... ee_sel_counting."""

    def construct(self):
        white_background(self)
        three = M.state_three()
        add_state(self, three, M.ORDER_THREE)
        clip_open(self, "ee_map_selection")

        # -- the ee map grows out of the three --------------------------------------------------------
        st = ee_map()
        anims = []
        for k in M.ORDER_MAP:
            anims.append(Transform(three[f"ee:{k}"], st[k]))
            for ch in ("mumu", "tautau"):
                anims.append(FadeOut(three[f"{ch}:{k}"]))
        anims += [FadeOut(three[f"lab_{ch}"]) for ch in M.CHANNELS]
        masked_play(self, *anims, run_time=1.6, rate_func=EASE)
        reset_to(self, st, M.ORDER_MAP)
        self.wait(0.3)
        ev = event_display()
        zoom_into(self, st, M.M_SEL, ev["det"], run_time=1.6)
        clip_cut(self, "ee_sel_event")

        # -- 3-04: a real selected event -----------------------------------------------------------------
        self.play(Create(ev["tracks"], lag_ratio=0.0), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(ev["hits"]), FadeIn(ev["deps"]), run_time=0.45)
        self.play(FadeIn(ev["labs"][0].sym), FadeIn(ev["labs"][1].sym), FadeIn(ev["stamp"]), run_time=0.4)
        self.play(FadeIn(ev["labs"][0].pt, shift=UP * 0.1), FadeIn(ev["labs"][1].pt, shift=UP * 0.1), run_time=0.5)
        # the trigger: the leading electron fired it
        self.play(GrowFromCenter(ev["ring"]), run_time=0.5)
        self.play(Flash(ev["ring"].get_center(), color=col(ACC), line_length=0.16, flash_radius=0.45, num_lines=10),
                  run_time=0.5)
        rows = [add_cut(self, 0)]
        # two electrons, opposite charge: the tracks bend opposite ways
        glow = VGroup(*[t.copy().set_stroke(col(ACC), width=9, opacity=0.6) for t in ev["tracks"]])
        self.play(FadeIn(glow), ev["labs"][0].sym.animate.scale(1.25), ev["labs"][1].sym.animate.scale(1.25), run_time=0.45)
        self.play(FadeOut(glow), ev["labs"][0].sym.animate.scale(1 / 1.25), ev["labs"][1].sym.animate.scale(1 / 1.25),
                  run_time=0.45)
        rows.append(add_cut(self, 1))
        clip_cut(self, "ee_sel_pt")

        # -- 3-05: p_T --------------------------------------------------------------------------------------
        display = VGroup(ev["det"], ev["tracks"], ev["hits"], ev["deps"], ev["ring"])
        s0, c0 = L["det"]
        s1, c1 = L["det_park"]
        tgt = display.copy()
        tgt.scale(s1 / s0, about_point=P3(c0)).shift(P3(c1) - P3(c0))
        dax = pt_axes()
        pt_marks = VGroup()
        for e in EV["electrons"]:
            p = dax.c2p(float(e["pt"]), 0.0)
            pt_marks.add(Line(p + DOWN * 0.02, p + UP * 0.35, stroke_color=col(GREEN_LINE), stroke_width=5.0))
        self.play(Transform(display, tgt), FadeOut(ev["stamp"]), FadeIn(dax),
                  *[ReplacementTransform(ev["labs"][i].pt, pt_marks[i]) for i in range(2)],
                  FadeOut(ev["labs"][0].sym), FadeOut(ev["labs"][1].sym), run_time=1.3, rate_func=EASE)
        stack0, stack1 = pt_stack(dax, 0.0), pt_stack(dax, 1.0)
        self.add(stack0)
        self.bring_to_front(pt_marks)
        self.play(Transform(stack0, stack1), FadeOut(pt_marks), run_time=1.3, rate_func=EASE)
        half = vref_line(dax, 91.19 / 2, color=GREY)
        half_lab = tex_h(r"m_Z/2", 0.2, color=GREY).next_to(dax.c2p(91.19 / 2, dax.y_top), UP, buff=0.08)
        self.play(Create(half), FadeIn(half_lab), run_time=0.6)
        veil = Rectangle(width=abs(dax.c2p(20, 0)[0] - dax.c2p(0, 0)[0]), height=dax.frame.height, stroke_width=0,
                         fill_color=col(WHITE), fill_opacity=0.5)
        veil.move_to(dax.frame.get_center()).align_to(dax.frame, LEFT)
        cut = DashedLine(dax.c2p(20, 0), dax.c2p(20, dax.y_top), dash_length=0.1, stroke_color=col(ACC), stroke_width=3.5)
        self.play(FadeIn(veil), Create(cut), run_time=0.7)
        rows.append(add_cut(self, 2))
        clip_cut(self, "ee_sel_eta")

        # -- 3-06: |eta| < 2.5 -----------------------------------------------------------------------------
        self.play(FadeOut(VGroup(dax, stack0, half, half_lab, veil, cut)), run_time=0.6)
        svw = side_view()
        self.play(FadeIn(svw["beam"]), FadeIn(svw["trk"]), FadeIn(svw["ecal"]), FadeIn(svw["ip"]), run_time=0.8)
        parts = []
        for e in EV["electrons"]:
            up = 1 if math.sin(float(e["phi"])) >= 0 else -1
            parts.append(sv_particle(float(e["eta"]), up, CHANNEL["ee"]))
        self.play(*[Create(p.line) for p in parts], run_time=0.7, rate_func=EASE)
        self.play(*[FadeIn(p.hits) for p in parts], *[FadeIn(p.dep) for p in parts], run_time=0.5)
        self.play(*[Create(l) for l in svw["eta_lines"]], FadeIn(svw["eta_lab"]), run_time=0.8)
        ghost = sv_particle(2.8, -1, GREY, dashed=True, with_hits=False)
        self.play(Create(ghost.line), run_time=0.6)
        self.play(FadeIn(ghost.dep), run_time=0.35)
        q = tex_h(r"e\,/\,\gamma\,?", 0.24, color=GREY).next_to(sv(EE_Z[1], -1.0), RIGHT, buff=0.15)
        self.play(FadeIn(q, shift=DOWN * 0.1), run_time=0.45)
        rows.append(add_cut(self, 3))
        clip_cut(self, "ee_sel_id")

        # -- 3-07: electron vs jet -------------------------------------------------------------------------
        side = VGroup(*svw.values(), *parts, ghost, q)
        self.play(FadeOut(side), run_time=0.6)
        pe, pj = id_panel("e", L["id_c"][0]), id_panel("jet", L["id_c"][1])
        for pn in (pe, pj):
            self.play(FadeIn(pn["cells"]), FadeIn(pn["sym"]), run_time=0.55)
            self.play(FadeIn(pn["trks"]), Create(pn["cone"]), run_time=0.5)
            self.play(FadeIn(pn["bars"]), run_time=0.4)
            self.play(LaggedStart(*[FadeIn(v, shift=UP * 0.08) for v in pn["vals"]], lag_ratio=0.35), run_time=0.9)
        ok = tick(pe["sym"].get_right() + RIGHT * 0.35, s=0.16, color=GREEN_LINE)
        no = cross(pj["sym"].get_right() + RIGHT * 0.35, s=0.12, color=GREY)
        self.play(Create(ok), Create(no), run_time=0.5)
        clip_cut(self, "ee_sel_wp")

        # -- 3-08: working points --------------------------------------------------------------------------
        panels = VGroup(*pe.values(), *pj.values(), ok, no)
        self.play(FadeOut(panels), run_time=0.6)
        wax = wp_axes()
        hist = wp_hists(wax, "none")
        key = wp_key(wax)
        name = wp_name("none")
        tax = to_axes()
        self.play(FadeIn(wax), FadeIn(key), FadeIn(tax), run_time=0.6)
        self.play(FadeIn(hist), FadeIn(name), run_time=0.6)
        pts = VGroup(to_point(tax, "none"))
        self.play(GrowFromCenter(pts[0]), run_time=0.3)
        curve = VGroup()
        for a, b in zip(WPS, WPS[1:]):
            seg = Line(tax.c2p(SEL["wp"][a]["eff_rel_none"], purity(a)), tax.c2p(SEL["wp"][b]["eff_rel_none"], purity(b)),
                       stroke_color=col(GREY), stroke_width=2.0)
            p = to_point(tax, b)
            self.play(Transform(hist, wp_hists(wax, b)), Transform(name, wp_name(b)), Create(seg), GrowFromCenter(p),
                      run_time=0.9, rate_func=EASE)
            curve.add(seg)
            pts.add(p)
            self.bring_to_front(pts)
        self.play(Transform(hist, wp_hists(wax, "medium")), Transform(name, wp_name("medium")),
                  Transform(pts[3], to_point(tax, "medium", lit=True)), run_time=0.8, rate_func=EASE)
        ring = Circle(radius=0.15, stroke_color=col(GREEN_LINE), stroke_width=3.0).move_to(pts[3])
        self.play(Create(ring), run_time=0.4)
        rows.append(add_cut(self, 4))
        clip_cut(self, "ee_sel_mass")

        # -- 3-09: the funnel fills the m_ee histogram -------------------------------------------------------
        self.play(FadeOut(VGroup(wax, hist, key, name, tax, pts, curve, ring)), FadeOut(display), run_time=0.6)
        rows.append(add_cut(self, 5))
        fun = big_funnel()
        self.play(*[r.animate.scale(0.15).move_to(fun.get_top() + DOWN * 0.25).set_opacity(0.0) for r in rows],
                  FadeIn(fun, scale=0.7), run_time=1.0, rate_func=EASE)
        self.remove(*rows)
        dax = main_axes("sel", x_labels=True)
        bars0 = data_bars(dax, 0.0)
        bars0.set_opacity(0.0)
        self.play(FadeIn(dax), run_time=0.5)
        self.add(bars0)
        rng = np.random.default_rng(5)
        spout = funnel_spout(fun)
        p = np.log10(np.maximum(DATA, 10.0)) - Y_EXP[0]
        p = p / p.sum()
        drops = []
        for _ in range(70):
            j = int(rng.choice(30, p=p))
            x = rng.uniform(EDGES[j], EDGES[j + 1])
            d = Dot(spout, radius=0.045, color=col(SAMPLE["Data"]))
            drops.append(Succession(d.animate(run_time=0.55, rate_func=rate_functions.ease_in_quad).move_to(dax.c2p(x, DATA[j])),
                                    FadeOut(d, run_time=0.15)))
        counter = DecimalNumber(0, num_decimal_places=0, group_with_commas=True, color=col(INK)).scale(_cap_scale("tex", 0.24))
        pre = tex_h(r"N =", 0.24)
        pre.move_to(P3((-3.6, 2.45)))
        counter.next_to(pre, RIGHT, buff=0.14)
        tracker = ValueTracker(0.0)
        counter.add_updater(lambda m: m.set_value(tracker.get_value()).next_to(pre, RIGHT, buff=0.14))
        self.add(pre, counter)
        self.play(LaggedStart(*drops, lag_ratio=0.03), Transform(bars0, data_bars(dax, 1.0)),
                  tracker.animate.set_value(float(DATA.sum())), run_time=2.6, rate_func=rate_functions.linear)
        counter.clear_updaters()
        self.remove(tracker)
        clip_cut(self, "ee_sel_mc")

        # -- 3-10: the simulation through the same selection ---------------------------------------------------
        sim = M.sim_file(fun.get_center() + RIGHT * 1.35 + UP * 1.0, 0.5, 0.66)
        self.play(FadeIn(sim, shift=LEFT * 0.3), run_time=0.5)
        self.play(sim.animate.scale(0.3).move_to(fun.get_top() + DOWN * 0.2).set_opacity(0.0), run_time=0.7, rate_func=EASE)
        self.remove(sim)
        S = plot_state("sel", NOM)
        grow = [GrowFromEdge(S["stack"][i], DOWN) for i in range(len(ORDER_S))]
        self.play(LaggedStart(*grow, lag_ratio=0.12), Transform(bars0, S["data"]), run_time=1.3)
        self.remove(bars0)
        self.add(S["data"])
        self.bring_to_front(S["data"])
        self.play(FadeOut(dax.x_labels), FadeOut(dax.x_title), FadeIn(S["ratio"]), FadeIn(S["key"]), run_time=0.7)
        self.remove(dax)
        self.add(S["dax"])
        self.bring_to_back(S["dax"])
        self.play(LaggedStart(*[FadeIn(d) for d in S["rdots"]], lag_ratio=0.03), run_time=0.9)
        C = counting_state()
        self.play(FadeIn(C["rlabel"], shift=LEFT * 0.1), run_time=0.4)
        clip_cut(self, "ee_sel_counting")

        # -- 3-11: counting, and what it throws away -----------------------------------------------------------
        self.play(FadeOut(fun), FadeOut(pre), FadeOut(counter), run_time=0.5)
        self.play(FadeIn(C["f1"], shift=UP * 0.1), run_time=0.7)
        self.wait(0.3)
        self.play(FadeIn(C["f2"], shift=UP * 0.1), run_time=0.7)
        self.play(FadeIn(C["f3"], shift=UP * 0.1), run_time=0.8)
        self.play(FadeIn(C["f4"], shift=UP * 0.1), run_time=0.6)
        self.wait(0.4)
        self.play(Transform(S["stack"], C["stack"]), Transform(S["data"], C["data"]), Transform(S["rdots"], C["rdots"]),
                  run_time=1.6, rate_func=EASE)
        end = counting_state()
        reset_to(self, end, ORDER_COUNT)
        self.wait(0.1)


class ZeeFit(Scene):
    """ee_map_fit ... ee_sigma."""

    def construct(self):
        white_background(self)
        prev = counting_state()
        add_state(self, prev, ORDER_COUNT)
        clip_open(self, "ee_map_fit")

        # -- 3-12: back to the map, into the fit --------------------------------------------------------------
        st = zoom_out_to_map(self, M.M_SEL, run_time=1.5)
        self.wait(0.35)
        F = plot_state("fit", NOM)
        detail = VGroup(*[F[k] for k in PLOT_ORDER])
        zoom_into(self, st, M.M_FIT, detail, run_time=1.6)
        reset_to(self, F, PLOT_ORDER)
        clip_cut(self, "ee_fit_mu")

        # -- 3-13: mu scales Z -> ee only ------------------------------------------------------------------------
        mu = ValueTracker(1.0)
        th = {n: ValueTracker(0.0) for n in NP_ROWS}
        gt = ValueTracker(0.0)                               # 0 = gammas at 1, 1 = gammas at their post-fit value

        def current():
            g = 1.0 + gt.get_value() * (GAM - 1.0)
            return model(mu.get_value(), {n: th[n].get_value() for n in NP_ROWS}, g)

        dax, rax = F["dax"], F["ratio"].rax
        self.remove(F["stack"], F["rdots"])
        stack = always_redraw(lambda: stack_hist(dax, EDGES, layers_of(current())))
        rdots = always_redraw(lambda: ratio_dots(rax, total_of(current())))
        self.add(stack, rdots)
        self.bring_to_front(F["data"])
        mtex = model_tex(False)
        ms = mu_slider()
        self.play(FadeIn(mtex, shift=DOWN * 0.1), FadeIn(ms), run_time=0.7)
        ms.sl.marker.add_updater(lambda m: m.become(ms.sl.marker_at(mu.get_value())))
        for v, t in ((1.15, 1.2), (0.85, 1.6), (1.0, 0.9)):
            self.play(mu.animate.set_value(v), run_time=t, rate_func=EASE)
        clip_cut(self, "ee_fit_nuisance")

        # -- 3-14: one nuisance parameter at a time ----------------------------------------------------------------
        mtex2 = model_tex(True)
        pp = pulls_frame(len(NP_ROWS))
        shown = VGroup(pp.band, pp.zero, pp.baseline, pp.ticks, pp.tick_labels)
        self.play(Transform(mtex, mtex2), FadeIn(shown), run_time=0.7)
        dots = {}
        for i, n in enumerate(NP_ROWS):
            d = pull_dot(pp, i, 0.0)
            d.add_updater(lambda m, i=i, n=n: m.become(pull_dot(pp, i, th[n].get_value())))
            dots[n] = d
        for i, n in enumerate(NP_ANIM):
            lab = pp.rows[i][0]
            self.play(FadeIn(lab, shift=RIGHT * 0.1), FadeIn(dots[n]), run_time=0.45)
            glow = lab.copy().set_color(col(ACC))
            self.add(glow)
            self.play(th[n].animate.set_value(1.0), run_time=0.8, rate_func=EASE)
            self.play(th[n].animate.set_value(-1.0), run_time=1.1, rate_func=EASE)
            self.play(th[n].animate.set_value(0.0), FadeOut(glow), run_time=0.7, rate_func=EASE)
        clip_cut(self, "ee_fit_result")

        # -- 3-15: the fit ------------------------------------------------------------------------------------------
        lik = likelihood_tex()
        self.play(FadeIn(lik, shift=UP * 0.1), run_time=0.8)
        self.play(*[FadeIn(pp.rows[i][0], shift=RIGHT * 0.1) for i in range(len(NP_ANIM), len(NP_ROWS))],
                  *[FadeIn(dots[n]) for n in NP_MORE], run_time=0.6)
        self.wait(0.3)
        anims = [mu.animate.set_value(MU["value"]), gt.animate.set_value(1.0)]
        anims += [th[n].animate.set_value(NPD[n]["pull"]) for n in NP_ROWS]
        self.play(*anims, run_time=2.8, rate_func=EASE)
        for m in (stack, rdots, ms.sl.marker, *dots.values()):
            m.clear_updaters()
        post = plot_state("fit", POST)
        self.play(Transform(stack, post["stack"]), Transform(rdots, post["rdots"]),
                  Transform(ms.sl.marker, ms.sl.marker_at(MU["value"], MU["err_up"])),
                  *[Transform(dots[n], pull_dot(pp, i, NPD[n]["pull"], NPD[n]["constraint"])) for i, n in enumerate(NP_ROWS)],
                  run_time=0.9, rate_func=EASE)
        self.remove(mu, gt, *th.values())
        mline = tex_h(rf"\mu_Z = {MU['value']:.3f} \pm {MU['err_up']:.3f}", 0.24).next_to(ms.sl.marker, UP, buff=0.35)
        self.play(FadeIn(mline, shift=UP * 0.1), FadeOut(mtex), run_time=0.6)
        clip_cut(self, "ee_sigma")

        # -- 3-16: sigma ------------------------------------------------------------------------------------------
        others = [m for m in live_mobjects(self) if m is not mline]
        clear_updaters(others)
        top = mu_line()
        sig = sigma_panel()
        self.play(*[FadeOut(m) for m in others], ReplacementTransform(mline, top), run_time=0.9, rate_func=EASE)
        self.play(FadeIn(sig.frame), run_time=0.5)
        self.play(FadeIn(sig.band), Create(sig.theory), FadeIn(sig.th_lab), run_time=0.6)
        self.play(ReplacementTransform(top.copy(), sig.label), GrowFromCenter(sig.point), run_time=0.9)
        reset_to(self, {"mu": top, "sig": sig}, ("mu", "sig"))
        self.wait(0.1)
