"""Section 2: from detector readout to NanoAOD (online reconstruction, trigger, skim). Schematic.

One scene, ``TriggerChain``, cut into one-idea clips (clip_open / clip_cut). A spine along the
bottom stays on screen the whole time; each stage opens a detail panel above it:

    trig_a_collisions   the slice, bunch crossings flash, events stream into L1 at 40 MHz
    trig_b_l1_primitives  L1 step 1: coarse calorimeter towers and muon segments light up
    trig_c_l1_objects   L1 step 2: the lit blocks become candidate objects (mu, e/gamma, jet)
    trig_d_l1_decision  L1 step 3: global trigger, a menu of ~440 seeds, some fire -> accept; O(us)
    trig_e_l1_rate      panel closes; 100 kHz leave L1 for the HLT
    trig_f_hlt_paths    HLT: paths of filters run on the event, each gives one trigger bit
    trig_g_hlt_rate     O(100 ms) per event, 1 kHz leave the HLT
    trig_h_raw          the RAW event: detector data + L1 result + HLT bits + some objects
    trig_i_nanoaod      the skim: detector data dropped, the rest squeezed into NanoAOD

Numbers on screen are the ones in the user's outline (17 Sep 2026): 40 MHz, 100 kHz, 1 kHz,
~440 L1 seeds (Run 2), O(us) / O(100 ms). Everything else is schematic; words are limited to
the names of the stages and data formats (L1, HLT, RAW, NanoAOD).
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    DOWN, LEFT, RIGHT, UP, AnimationGroup, Create, FadeIn, FadeOut, Flash, LaggedStart, Line,
    Polygon, Rectangle, ReplacementTransform, RoundedRectangle, Scene, Square, Succession,
    TAU, Transform, VGroup, Dot, rate_functions,
)

from style.bnd_style import *  # noqa: E402,F401,F403

EASE = rate_functions.ease_in_out_sine
ACC = DETECTOR_ACCENT

# ---------------------------------------------------------------------------
# anchors (title band y > 2.7 and top-left block x < -5.85, y > 0.22 stay empty)
# ---------------------------------------------------------------------------
SPINE_Y = -1.55
A = {
    "slice": (-5.1, SPINE_Y), "slice_r": 1.0,
    "l1": (-2.3, SPINE_Y), "hlt": (0.7, SPINE_Y), "raw": (3.5, SPINE_Y), "nano": (5.9, SPINE_Y),
    "box_wh": (1.2, 1.2),
    "rate_y": SPINE_Y + 0.36, "clock_y": SPINE_Y - 1.05,
    "panel": (1.6, 1.15), "panel_wh": (10.2, 2.7),
}


def P(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    return p if p.shape == (3,) else np.array([p[0], p[1], 0.0])


_TEX_N: list = []


def tex(expr: str, h: float = 0.30, color=INK, **kw):
    """MathTex scaled so that a capital N is ``h`` high (uniform symbol size)."""
    if not _TEX_N:
        _TEX_N.append(mathtex("N").height)
    return mathtex(expr, color=color, **kw).scale(h / _TEX_N[0])


def arrow(a, b, color=INK, sw=3.0) -> VGroup:
    ln = Line(P(a), P(b), stroke_color=col(color), stroke_width=sw)
    return VGroup(ln, arrow_tip_on(ln, color=color, at=1.0, tip_length=0.18))


def stage_box(center, label) -> VGroup:
    w, h = A["box_wh"]
    b = RoundedRectangle(width=w, height=h, corner_radius=0.12, stroke_color=col(ACC), stroke_width=2.5,
                         fill_color=col(WHITE), fill_opacity=1.0).move_to(P(center))
    g = VGroup(b, tex(label, h=0.34).move_to(P(center)))
    g.box = b
    return g


def file_glyph(w, h, fill=LIGHT_GREY, stroke=INK, sw=2.0) -> VGroup:
    e = 0.22 * min(w, h)
    body = Polygon(P((-w / 2, -h / 2)), P((w / 2, -h / 2)), P((w / 2, h / 2 - e)), P((w / 2 - e, h / 2)),
                   P((-w / 2, h / 2)), fill_color=col(fill), fill_opacity=1.0, stroke_color=col(stroke),
                   stroke_width=sw)
    fold = Polygon(P((w / 2 - e, h / 2)), P((w / 2 - e, h / 2 - e)), P((w / 2, h / 2 - e)),
                   fill_color=lighten(stroke, 0.45), fill_opacity=1.0, stroke_color=col(stroke), stroke_width=sw)
    g = VGroup(body, fold)
    g.body = body
    return g


def node_file(center, w, h) -> VGroup:
    f = file_glyph(w, h)
    lines = VGroup(*[Line(P((-0.3 * w, y * h)), P((0.25 * w, y * h)), stroke_color=col(INK), stroke_width=1.4)
                     for y in (0.15, -0.05, -0.25)])
    return VGroup(f, lines).move_to(P(center))


def panel() -> RoundedRectangle:
    w, h = A["panel_wh"]
    return RoundedRectangle(width=w, height=h, corner_radius=0.18, stroke_color=col(ACC), stroke_width=2.5,
                            fill_color=lighten(ACC, 0.92), fill_opacity=1.0).move_to(P(A["panel"]))


def callout(node_center, half=0.6) -> VGroup:
    """Two thin lines from the node's top corners to the panel's bottom edge (a zoom callout)."""
    px, py = A["panel"]
    pw, ph = A["panel_wh"]
    x, y = node_center
    top = y + A["box_wh"][1] / 2
    xl, xr = max(px - pw / 2 + 0.2, x - 1.6), min(px + pw / 2 - 0.2, x + 1.6)
    return VGroup(*[Line(P((x + s * half, top)), P((xe, py - ph / 2)), stroke_color=col(ACC), stroke_width=1.8)
                    for s, xe in ((-1, xl), (1, xr))])


def stream(src, dst, n, run_time, seed=0, color=INK, spread=0.18) -> AnimationGroup:
    """``n`` event dots flying from ``src`` to ``dst`` (then fading): the event rate."""
    rng = np.random.default_rng(seed)
    anims = []
    for i in range(n):
        dy = rng.uniform(-spread, spread)
        d = Dot(P(src) + P((0, dy)), radius=0.045, color=col(color))
        dur = 0.55
        anims.append(Succession(FadeIn(d, run_time=0.05),
                                d.animate(run_time=dur, rate_func=rate_functions.linear).move_to(P(dst) + P((0, dy))),
                                FadeOut(d, run_time=0.08)))
    return LaggedStart(*anims, lag_ratio=1.0 / max(1, n) * 3.0, run_time=run_time)


def latency(center, expr) -> VGroup:
    c = clock(P(center) + LEFT * 0.55, radius=0.24)
    lab = tex(expr, h=0.24).next_to(c, RIGHT, buff=0.14)
    g = VGroup(c, lab)
    g.clock, g.lab = c, lab
    return g


# ---------------------------------------------------------------------------
# panel contents
# ---------------------------------------------------------------------------
# L1: primitives (coarse towers + muon segments), objects, global decision
L1_X = (-2.2, 1.3, 4.6)
GRID_COLS, GRID_ROWS, CELL = 7, 3, 0.30
LIT_CALO = ((1, 2), (2, 2), (2, 1), (5, 0))        # (col, row) towers above threshold
LIT_MU = (3, 4)                                    # muon-segment columns


def l1_grid():
    cx, cy = L1_X[0], A["panel"][1] + 0.25
    towers, segs = VGroup(), VGroup()
    x0 = cx - (GRID_COLS - 1) * CELL / 2
    for r in range(GRID_ROWS):
        for c in range(GRID_COLS):
            towers.add(Square(side_length=CELL * 0.9, stroke_color=col(DETECTOR["hcal"]["stroke"]),
                              stroke_width=1.2, fill_color=col(CMS["hcal"]), fill_opacity=1.0)
                       .move_to(P((x0 + c * CELL, cy + (1 - r) * CELL))))
    for c in range(GRID_COLS):
        segs.add(Rectangle(width=CELL * 0.9, height=0.12, stroke_color=col(DETECTOR["muon"]["stroke"]),
                           stroke_width=1.2, fill_color=col(CMS["muon"]), fill_opacity=1.0)
                 .move_to(P((x0 + c * CELL, cy - 1.5 * CELL - 0.18))))
    g = VGroup(towers, segs)
    g.towers, g.segs = towers, segs
    return g


def l1_objects():
    x = L1_X[1]
    y0 = A["panel"][1]
    out = VGroup()
    for expr, dy in ((r"e/\gamma", 0.75), (r"\mathrm{jet}", 0.0), (r"\mu", -0.75)):
        ring = RoundedRectangle(width=1.05, height=0.52, corner_radius=0.2, stroke_color=col(ACC), stroke_width=2.2,
                                fill_color=col(WHITE), fill_opacity=1.0).move_to(P((x, y0 + dy)))
        out.add(VGroup(ring, tex(expr, h=0.26).move_to(ring)))
    return out


MENU_COLS, MENU_ROWS = 20, 22                       # 440 seeds
MENU_CELL = 0.085


def l1_menu():
    cx, cy = L1_X[2] - 0.3, A["panel"][1]
    rng = np.random.default_rng(7)
    g = VGroup()
    x0 = cx - (MENU_COLS - 1) * MENU_CELL / 2
    y0 = cy + (MENU_ROWS - 1) * MENU_CELL / 2
    for r in range(MENU_ROWS):
        for c in range(MENU_COLS):
            # most seeds are physics; a band of calibration / monitoring seeds at the bottom, drawn in a grey tint
            fill = LIGHT_GREY if r < MENU_ROWS - 4 else lighten(GREY, 0.3)
            g.add(Square(side_length=MENU_CELL * 0.78, stroke_width=0, fill_color=col(fill), fill_opacity=1.0)
                  .move_to(P((x0 + c * MENU_CELL, y0 - r * MENU_CELL))))
    fired = sorted(rng.choice(MENU_ROWS * MENU_COLS - 80, 9, replace=False).tolist())
    g.fired = fired
    return g


# HLT: paths of filters
HLT_ROWS = (2.05, 1.5, 0.95, 0.4)
HLT_FX = (-1.6, -0.4, 0.8, 2.0, 3.2)
HLT_FAIL_AT = (2, None, 0, None)                    # index of the failing filter per path (None = accepted)
HLT_BIT_X = 4.7


def hlt_paths():
    rows = VGroup()
    for y in HLT_ROWS:
        start = Dot(P((-2.6, y)), radius=0.06, color=col(INK))
        link = Line(P((-2.6, y)), P((HLT_FX[-1], y)), stroke_color=col(GREY), stroke_width=2.0)
        filt = VGroup(*[Square(side_length=0.34, stroke_color=col(ACC), stroke_width=2.2, fill_color=col(WHITE),
                               fill_opacity=1.0).move_to(P((x, y))) for x in HLT_FX])
        tail = arrow((HLT_FX[-1] + 0.2, y), (HLT_BIT_X - 0.3, y), color=GREY, sw=2.0)
        r = VGroup(link, start, filt, tail)
        r.filt = filt
        rows.add(r)
    return rows


def hlt_bits():
    return VGroup(*[tex("0" if f is not None else "1", h=0.30, color=INK if f is None else GREY)
                    .move_to(P((HLT_BIT_X, y))) for f, y in zip(HLT_FAIL_AT, HLT_ROWS)])


def cross(center, s=0.12, color=GREY) -> VGroup:
    c = P(center)
    return VGroup(Line(c + P((-s, -s)), c + P((s, s)), stroke_color=col(color), stroke_width=3.0),
                  Line(c + P((-s, s)), c + P((s, -s)), stroke_color=col(color), stroke_width=3.0))


def tick(center, s=0.22, color=ACC) -> VGroup:
    c = P(center)
    return VGroup(Line(c + P((-s, 0.0)), c + P((-0.3 * s, -0.7 * s)), stroke_color=col(color), stroke_width=6),
                  Line(c + P((-0.3 * s, -0.7 * s)), c + P((s, 0.8 * s)), stroke_color=col(color), stroke_width=6))


# RAW event: a file with four rows
RAW_C, RAW_WH = (0.9, 1.15), (5.4, 2.45)
RAW_ROW_Y = (1.9, 1.38, 0.86, 0.34)
RAW_LAB_X = -1.1


def raw_rows():
    x0 = RAW_LAB_X + 0.85
    # 1 detector data: a mini slice and its readout words
    det = mini_slice(0.12, P((RAW_LAB_X, RAW_ROW_Y[0])))
    rng = np.random.default_rng(3)
    words = VGroup(*[Rectangle(width=0.07, height=float(rng.uniform(0.08, 0.36)), stroke_width=0,
                               fill_color=col(DEPOSIT["hcal"] if i % 3 else DEPOSIT["muon"]), fill_opacity=1.0)
                     .move_to(P((x0 + 0.11 * i, RAW_ROW_Y[0]))) for i in range(30)])
    r1 = VGroup(det, words)
    # 2 L1 result
    l1 = VGroup(tex(r"\mathrm{L1}", h=0.24).move_to(P((RAW_LAB_X, RAW_ROW_Y[1]))),
                *[Square(side_length=0.2, stroke_color=col(ACC), stroke_width=1.5,
                         fill_color=col(ACC if i in (2, 5, 9) else WHITE), fill_opacity=1.0)
                  .move_to(P((x0 + 0.1 + 0.28 * i, RAW_ROW_Y[1]))) for i in range(12)])
    # 3 HLT bits (the four of the HLT panel land here, more follow)
    hl = tex(r"\mathrm{HLT}", h=0.24).move_to(P((RAW_LAB_X, RAW_ROW_Y[2])))
    bits = VGroup(*[tex(b, h=0.26, color=INK if b == "1" else GREY).move_to(P((x0 + 0.1 + 0.28 * i, RAW_ROW_Y[2])))
                    for i, b in enumerate("010100100010")])
    r3 = VGroup(hl, bits)
    # 4 some reconstructed objects
    objs = VGroup(*[tex(e, h=0.26).move_to(P((x0 + 0.25 + 0.75 * i, RAW_ROW_Y[3])))
                    for i, e in enumerate((r"\mu", r"e", r"\gamma", r"\mathrm{jet}", r"p_T^{\mathrm{miss}}"))])
    r4 = VGroup(Dot(P((RAW_LAB_X, RAW_ROW_Y[3])), radius=0.0, fill_opacity=0.0), objs)
    rows = VGroup(r1, l1, r3, r4)
    rows.bits = bits
    return rows


def trigger_spine() -> dict:
    """The spine along the bottom (pure builder): slice -> L1 -> HLT -> RAW -> NanoAOD.
    Keys in build order = ``ORDER_TRIG``. Section 2's map (scenes/s2_map.py) draws this
    same spine at 0.6x as its experiment row and zooms onto it."""
    det = mini_slice(A["slice_r"] / R_DET, P(A["slice"]))
    l1 = stage_box(A["l1"], r"\mathrm{L1}")
    hlt = stage_box(A["hlt"], r"\mathrm{HLT}")
    w = A["box_wh"][0]
    raw_icon = node_file(A["raw"], 0.78, 1.0)
    nano_icon = node_file(A["nano"], 0.42, 0.56)
    raw_lab = tex(r"\mathrm{RAW}", h=0.24).next_to(raw_icon, DOWN, buff=0.18)
    nano_lab = tex(r"\mathrm{NanoAOD}", h=0.24).next_to(nano_icon, DOWN, buff=0.18)
    sx, lx, hx, rx, nx = A["slice"][0], A["l1"][0], A["hlt"][0], A["raw"][0], A["nano"][0]
    y = SPINE_Y
    arr = [arrow((sx + A["slice_r"] + 0.08, y), (lx - w / 2 - 0.08, y)),
           arrow((lx + w / 2 + 0.08, y), (hx - w / 2 - 0.08, y)),
           arrow((hx + w / 2 + 0.08, y), (rx - 0.5, y)),
           arrow((rx + 0.5, y), (nx - 0.32, y))]
    return {"det": det, "arr0": arr[0], "l1": l1, "arr1": arr[1], "hlt": hlt, "arr2": arr[2],
            "raw_icon": raw_icon, "raw_lab": raw_lab, "arr3": arr[3], "nano_icon": nano_icon, "nano_lab": nano_lab}


ORDER_TRIG = ("det", "arr0", "l1", "arr1", "hlt", "arr2", "raw_icon", "raw_lab", "arr3", "nano_icon", "nano_lab")


class TriggerChain(Scene):
    """Standalone: opens on a white frame and builds the spine (``enter``), ends on the NanoAOD
    skim (``leave``). ``s2_map.TriggerInMap`` overrides both hooks to zoom in from and out to the map."""

    def enter(self, st: dict) -> None:
        white_background(self)
        clip_open(self, "trig_a_collisions", hold=0.1)
        arr = [st[f"arr{i}"] for i in range(4)]
        self.play(FadeIn(st["det"]), run_time=0.6)
        self.play(LaggedStart(Create(arr[0]), FadeIn(st["l1"]), Create(arr[1]), FadeIn(st["hlt"]), Create(arr[2]),
                              FadeIn(st["raw_icon"]), FadeIn(st["raw_lab"]), Create(arr[3]), FadeIn(st["nano_icon"]),
                              FadeIn(st["nano_lab"]), lag_ratio=0.25), run_time=1.6)

    def leave(self) -> None:
        self.wait(0.1)

    def construct(self):
        # -- spine ------------------------------------------------------------
        st = trigger_spine()
        self.enter(st)
        det, l1, hlt, raw_icon, nano_icon = st["det"], st["l1"], st["hlt"], st["raw_icon"], st["nano_icon"]
        arr = [st[f"arr{i}"] for i in range(4)]
        mids = [0.5 * (a[0].get_start() + a[0].get_end()) for a in arr]
        rates = [tex(e, h=0.24).move_to(m + UP * 0.38)
                 for e, m in zip((r"40\,\mathrm{MHz}", r"100\,\mathrm{kHz}", r"1\,\mathrm{kHz}"), mids)]

        # bunch crossings: flashes with a few tracks, faster and faster
        c0 = det.c
        rng = np.random.default_rng(11)
        r_trk = det.radii["tob"][1]
        for j, dur in enumerate((0.5, 0.35, 0.25, 0.18, 0.14, 0.12)):
            phis = rng.uniform(0, TAU, 3)
            trks = VGroup(*[track(det, float(p), float(rng.choice([-1, 1]) * rng.uniform(1.2, 3.0)), r_end=r_trk,
                                  color=INK, sw=2.5) for p in phis])
            self.play(Flash(c0, color=col(INK), line_length=0.08, flash_radius=0.14, num_lines=8, run_time=dur),
                      Create(trks, lag_ratio=0.0, run_time=dur), rate_func=rate_functions.linear)
            self.remove(trks)
        a0s, a0e = arr[0][0].get_start(), arr[0][0].get_end()
        self.play(stream(a0s, a0e, 24, 1.8, seed=1), FadeIn(rates[0], shift=UP * 0.1, run_time=0.5))

        # -- L1 -----------------------------------------------------------------
        clip_cut(self, "trig_b_l1_primitives")
        pan = panel()
        call = callout(A["l1"])
        grid = l1_grid()
        self.play(FadeIn(pan, scale=0.9), Create(call), l1.box.animate.set_fill(lighten(ACC, 0.8)), run_time=0.7)
        self.play(FadeIn(grid, lag_ratio=0.02), run_time=0.7)
        lit = VGroup(*[grid.towers[r * GRID_COLS + c] for c, r in LIT_CALO])
        lit_mu = VGroup(*[grid.segs[c] for c in LIT_MU])
        self.play(LaggedStart(*[t.animate.set_fill(col(DEPOSIT["hcal"])) for t in lit],
                              *[s.animate.set_fill(col(DEPOSIT["muon"])) for s in lit_mu], lag_ratio=0.2),
                  run_time=1.0)

        clip_cut(self, "trig_c_l1_objects")
        objs = l1_objects()
        a1 = arrow((L1_X[0] + 1.35, A["panel"][1]), (L1_X[1] - 0.75, A["panel"][1]), color=ACC)
        self.play(Create(a1), run_time=0.4)
        src = (VGroup(lit[0], lit[1]), VGroup(lit[2], lit[3]), lit_mu)
        self.play(LaggedStart(*[ReplacementTransform(s.copy(), o) for s, o in zip(src, objs)], lag_ratio=0.3),
                  run_time=1.3, rate_func=EASE)

        clip_cut(self, "trig_d_l1_decision")
        menu = l1_menu()
        a2 = arrow((L1_X[1] + 0.65, A["panel"][1]), (menu.get_left()[0] - 0.15, A["panel"][1]), color=ACC)
        n440 = tex("440", h=0.26).next_to(menu, RIGHT, buff=0.25).shift(UP * 0.6)
        self.play(Create(a2), FadeIn(menu, lag_ratio=0.002), run_time=1.0)
        self.play(FadeIn(n440, shift=LEFT * 0.1), run_time=0.4)
        lat1 = latency((A["l1"][0], A["clock_y"]), r"\mathcal{O}(\mu\mathrm{s})")
        self.play(FadeIn(lat1), run_time=0.3)
        ok = tick(n440.get_center() + DOWN * 1.2)
        self.play(LaggedStart(*[menu[i].animate.set_fill(col(ACC)).scale(2.2) for i in menu.fired], lag_ratio=0.1),
                  lat1.clock.turns.animate.set_value(1.0), run_time=1.2)
        self.remove(lat1.clock.turns)
        self.play(Create(ok), run_time=0.4)

        clip_cut(self, "trig_e_l1_rate")
        l1_detail = VGroup(grid, a1, objs, a2, menu, n440, ok)
        self.play(FadeOut(l1_detail), FadeOut(call), FadeOut(pan), l1.box.animate.set_fill(col(WHITE)), run_time=0.7)
        a1s, a1e = arr[1][0].get_start(), arr[1][0].get_end()
        self.play(stream(a0s, a0e, 18, 1.6, seed=2), stream(a1s, a1e, 5, 1.6, seed=3),
                  FadeIn(rates[1], shift=UP * 0.1, run_time=0.5))

        # -- HLT ----------------------------------------------------------------
        clip_cut(self, "trig_f_hlt_paths")
        pan = panel()
        call = callout(A["hlt"])
        paths = hlt_paths()
        self.play(FadeIn(pan, scale=0.9), Create(call), hlt.box.animate.set_fill(lighten(ACC, 0.8)), run_time=0.7)
        self.play(LaggedStart(*[FadeIn(r, shift=RIGHT * 0.15) for r in paths], lag_ratio=0.15), run_time=0.9)
        tokens = [Dot(P((-2.6, y)), radius=0.09, color=col(INK)) for y in HLT_ROWS]
        self.add(*tokens)
        alive = [True] * len(HLT_ROWS)
        marks = VGroup()
        for k, fx in enumerate(HLT_FX):
            anims = []
            for i, (y, fail) in enumerate(zip(HLT_ROWS, HLT_FAIL_AT)):
                if not alive[i]:
                    continue
                sq = paths[i].filt[k]
                anims.append(tokens[i].animate.move_to(P((fx, y))))
                if fail == k:
                    m = cross((fx, y))
                    marks.add(m)
                    anims += [sq.animate.set_fill(col(LIGHT_GREY)).set_stroke(col(GREY)), FadeIn(m)]
                    alive[i] = False
                else:
                    anims.append(sq.animate.set_fill(lighten(ACC, 0.55)))
            self.play(*anims, run_time=0.45 if k == 0 else 0.32, rate_func=EASE)
            for i, (fail) in enumerate(HLT_FAIL_AT):
                if fail == k:
                    self.play(FadeOut(tokens[i]), run_time=0.12)
        bits = hlt_bits()
        live = [tokens[i] for i, f in enumerate(HLT_FAIL_AT) if f is None]
        self.play(*[t.animate.move_to(P((HLT_BIT_X - 0.35, t.get_y()))) for t in live], run_time=0.35)
        self.play(*[FadeOut(t) for t in live], LaggedStart(*[FadeIn(b, scale=1.3) for b in bits], lag_ratio=0.15),
                  run_time=0.6)

        clip_cut(self, "trig_g_hlt_rate")
        lat2 = latency((A["hlt"][0], A["clock_y"]), r"\mathcal{O}(100\,\mathrm{ms})")
        self.play(FadeIn(lat2), run_time=0.3)
        a2s, a2e = arr[2][0].get_start(), arr[2][0].get_end()
        self.play(lat2.clock.turns.animate.set_value(2.0), stream(a0s, a0e, 18, 2.0, seed=4),
                  stream(a1s, a1e, 6, 2.0, seed=5), stream(a2s, a2e, 2, 2.0, seed=6),
                  FadeIn(rates[2], shift=UP * 0.1, run_time=0.5))
        self.remove(lat2.clock.turns)

        # -- RAW ----------------------------------------------------------------
        clip_cut(self, "trig_h_raw")
        call2 = callout(A["raw"], half=0.45)
        raw = file_glyph(*RAW_WH).move_to(P(RAW_C))
        rows = raw_rows()
        self.play(FadeOut(paths), FadeOut(marks), ReplacementTransform(call, call2),
                  hlt.box.animate.set_fill(col(WHITE)), raw_icon[0].body.animate.set_fill(lighten(ACC, 0.8)),
                  run_time=0.7)
        self.play(FadeIn(raw, scale=0.95), run_time=0.5)
        self.play(FadeIn(rows[0]), run_time=0.6)
        self.play(FadeIn(rows[1], shift=RIGHT * 0.1), run_time=0.5)
        self.play(FadeIn(rows[2][0]), ReplacementTransform(bits, VGroup(*rows.bits[:4])), run_time=0.8, rate_func=EASE)
        self.play(FadeIn(VGroup(*rows.bits[4:]), lag_ratio=0.1), run_time=0.5)
        self.play(FadeIn(rows[3], shift=RIGHT * 0.1, lag_ratio=0.1), run_time=0.7)

        clip_cut(self, "trig_i_nanoaod")
        # the skim: the detector readout is dropped, the rest squeezed into a small file
        call3 = callout(A["nano"], half=0.25)
        self.play(ReplacementTransform(call2, call3), raw_icon[0].body.animate.set_fill(col(LIGHT_GREY)),
                  nano_icon[0].body.animate.set_fill(lighten(ACC, 0.8)), run_time=0.6)
        raw_all = VGroup(raw, rows)
        self.play(raw_all.animate.scale(0.62).move_to(P((-1.3, A["panel"][1]))), run_time=0.8, rate_func=EASE)
        self.play(rows[0].animate.set_opacity(0.2), run_time=0.5)
        skim = arrow((0.55, A["panel"][1]), (2.35, A["panel"][1]), color=ACC)
        rest = VGroup(rows[1], rows[2], rows[3])
        body = rest.copy().set_opacity(1.0).scale(0.85)
        nano = file_glyph(body.width + 0.25, body.height + 0.3).move_to(P((3.9, A["panel"][1])))
        body.move_to(nano.get_center() + LEFT * 0.04)
        self.play(Create(skim), FadeIn(nano, scale=0.9), run_time=0.5)
        self.play(ReplacementTransform(rest.copy(), body), run_time=1.1, rate_func=EASE)
        a3s, a3e = arr[3][0].get_start(), arr[3][0].get_end()
        self.play(stream(a0s, a0e, 12, 1.6, seed=7), stream(a1s, a1e, 4, 1.6, seed=8),
                  stream(a2s, a2e, 2, 1.6, seed=9), stream(a3s, a3e, 2, 1.6, seed=10))
        self.leave()
