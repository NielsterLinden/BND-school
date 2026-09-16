"""Manim helpers for the BND-school Z cross-section talk (ManimCE 0.20.1, cairo).

Every scene imports colours, font and the drawing helpers from here and never
hardcodes a hex or a font. The palette itself lives in the manim-free
``style/palette.py`` so matplotlib scripts share the same values.

    from style.bnd_style import *   # in scenes (see the bootstrap in CLAUDE.md)

Colour code: six chapter colours (``CHAPTER``) for everything that is ours,
the official CMS logo colours (``CMS``) for the detector drawing. The detector
geometry *is* the logo's: ``CMSLogo`` (the square, a quarter of the layers)
and ``CMSSlice`` (the full transverse view) share ``LOGO_R``, so the section-2
build animation can turn the one into the other.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from manim import (
    DOWN, LEFT, RIGHT, UP, TAU, ORIGIN, BOLD,
    AnnularSector, Annulus, Axes, Circle, DashedLine, DashedVMobject, Dot, Ellipse, FadeOut,
    Line, LogBase, ManimColor, MathTex, ParametricFunction, Polygon, Rectangle,
    RoundedRectangle, Sector, Succession, Text, ValueTracker, VGroup, VMobject,
    always_redraw, rate_functions,
)

from style.palette import (  # noqa: F401  (re-exported)
    CHAPTER, CHAPTER_NAME, PURPLE, CYAN, GREEN, GOLD, RED, SLATE, BLACK, WHITE,
    mix, tint, shade,
    BG, INK, GREY, LIGHT_GREY, CHANNEL, CHANNEL_LINE, TRACK_OUTLINE, COMBINED, THEORY,
    DETECTOR_ACCENT, HIGHLIGHT, SAMPLE, PARTICLE, CMS, DETECTOR, DEPOSIT, FLAG, FONT,
)

# ---------------------------------------------------------------------------
# colours / text
# ---------------------------------------------------------------------------

def col(x) -> ManimColor:
    """Anything -> ManimColor (hex string, ManimColor, or a palette key)."""
    if isinstance(x, ManimColor):
        return x
    for table in (CHANNEL_LINE, SAMPLE, PARTICLE, DEPOSIT, CMS):
        if isinstance(x, str) and x in table:
            return ManimColor(table[x])
    return ManimColor(x)


def lighten(color, amount: float) -> ManimColor:
    """Move ``color`` toward white by ``amount`` in [0, 1] (same mix as
    ``palette.tint`` so Manim and matplotlib tints agree)."""
    return ManimColor(tint(col(color).to_hex(), amount))


def darken(color, amount: float) -> ManimColor:
    """Move ``color`` toward black by ``amount`` in [0, 1] (``palette.shade``)."""
    return ManimColor(shade(col(color).to_hex(), amount))


def text(s: str, color=INK, **kw) -> Text:
    """Plain text in the deck font. Reserve for tick labels / values / symbols;
    narrative text belongs in PowerPoint (see CLAUDE.md)."""
    return Text(s, font=FONT, color=col(color), **kw)


def mathtex(expr: str, color=INK, **kw) -> MathTex:
    return MathTex(expr, color=col(color), **kw)


def white_background(scene) -> None:
    scene.camera.background_color = ManimColor(BG)


# ---------------------------------------------------------------------------
# frozen data
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def load_data(name: str):
    """Load presentation/data/<name>.json (dict/list) or .csv (list[dict]).
    Numbers are frozen there by the channel code; scenes never recompute."""
    pj, pc = DATA_DIR / f"{name}.json", DATA_DIR / f"{name}.csv"
    if pj.exists():
        return json.loads(pj.read_text())
    if pc.exists():
        with pc.open() as fh:
            return list(csv.DictReader(fh))
    raise FileNotFoundError(f"no {name}.json/.csv in {DATA_DIR}")


# ---------------------------------------------------------------------------
# Feynman primitives (plain geometry, no package)
# ---------------------------------------------------------------------------

def _p3(p):
    p = np.asarray(p, dtype=float)
    return p if p.shape == (3,) else np.array([p[0], p[1], 0.0])


def _e(phi):
    return np.array([math.cos(phi), math.sin(phi), 0.0])


unit = _e          # public name: unit vector at angle ``phi`` (radians)


def fline(p0, p1, color=INK, sw: float = 4.0) -> Line:
    """Fermion line without an arrow."""
    return Line(_p3(p0), _p3(p1), stroke_color=col(color), stroke_width=sw)


def fermion(p_from, p_to, color=INK, sw: float = 4.0, tip_length: float = 0.22) -> VGroup:
    """Fermion line with a mid-point arrow encoding fermion flow. Two segments,
    so ``Create`` pops the arrowed half first; use ``fline`` for a smooth draw
    and add the arrow afterwards if that matters."""
    a, b = _p3(p_from), _p3(p_to)
    mid = (a + b) / 2
    seg_a = Line(a, mid, stroke_color=col(color), stroke_width=sw)
    seg_a.add_tip(tip_length=tip_length, tip_width=tip_length * 0.85)
    seg_b = Line(mid, b, stroke_color=col(color), stroke_width=sw)
    return VGroup(seg_a, seg_b)


def arrow_tip_on(line: Line, color=INK, tip_length: float = 0.22, at: float = 0.5) -> VGroup:
    """A standalone arrowhead placed at fraction ``at`` along ``line`` (so a
    line can be Created smoothly and the tip faded in afterwards)."""
    a, b = line.get_start(), line.get_end()
    p = a + (b - a) * at
    stub = Line(p - (b - a) * 1e-3, p, stroke_width=0)
    stub.add_tip(tip_length=tip_length, tip_width=tip_length * 0.85)
    tip = stub.tip
    tip.set_color(col(color))
    return VGroup(tip)


def dashed(p0, p1, color=INK, sw: float = 3.5, dash_length: float = 0.14) -> DashedLine:
    """Dashed propagator / undetected particle (neutrino, scalar)."""
    return DashedLine(_p3(p0), _p3(p1), dash_length=dash_length,
                      stroke_color=col(color), stroke_width=sw)


def wavy(p0, p1, color=INK, amplitude: float = 0.10, wavelength: float = 0.55,
         sw: float = 3.5, n: int | None = None) -> VMobject:
    """Photon / Z / W: sine offset perpendicular to the line. ``wavelength`` is
    held fixed so long and short bosons read as the same particle."""
    a, b = _p3(p0), _p3(p1)
    d = b - a
    L = float(np.linalg.norm(d))
    u = d / L
    nperp = np.array([-u[1], u[0], 0.0])
    waves = max(1, round(L / wavelength))
    n = n or max(120, 40 * waves)
    pts = [a + u * (t * L) + nperp * amplitude * math.sin(TAU * waves * t)
           for t in np.linspace(0.0, 1.0, n)]
    vm = VMobject(stroke_color=col(color), stroke_width=sw)
    vm.set_points_smoothly(pts)
    return vm


def gluon(p0, p1, color=INK, n_loops: int | None = None, radius: float = 0.15,
          sw: float = 3.5) -> ParametricFunction:
    """Gluon as a coiled spring (prolate cycloid)."""
    a, b = _p3(p0), _p3(p1)
    L = float(np.linalg.norm(b - a))
    n_loops = n_loops or max(3, round(L / 0.45))
    pitch = L / (2 * np.pi * n_loops)

    def f(t):
        th = 2 * np.pi * n_loops * t
        return np.array([pitch * th - radius * np.sin(th), radius * np.cos(th), 0.0])

    curve = ParametricFunction(f, t_range=[0, 1, 0.004],
                               color=col(color), stroke_width=sw)
    curve.put_start_and_end_on(a, b)
    return curve


def vertex_dot(p, color=INK, radius: float = 0.06) -> Dot:
    return Dot(_p3(p), radius=radius, color=col(color))


def shower_tree(root, direction: float, depth: int = 3, length: float = 1.2,
                spread: float = 0.30, shrink: float = 0.72, color=None,
                photon_color=None, sw: float = 3.0, seed: int = 3) -> VGroup:
    """Electromagnetic cascade in Feynman style: an electron radiates a photon
    (e -> e gamma), the photon converts (gamma -> e+ e-), and so on for
    ``depth`` generations, every generation shorter by ``shrink``. Returns a
    VGroup whose submobjects are the generations (each a VGroup of lines) so
    ``LaggedStart(*[Create(g) for g in tree])`` grows it outward. Straight
    electron lines in the ee colour, wavy photons in the photon colour."""
    rng = np.random.default_rng(seed)
    cc = col(PARTICLE["e"] if color is None else color)
    pc = col(PARTICLE["photon"] if photon_color is None else photon_color)
    gens = [VGroup() for _ in range(depth)]
    # (point, angle, kind, generation)
    frontier = [(_p3(root), direction, "e", 0)]
    while frontier:
        p, ang, kind, g = frontier.pop(0)
        if g >= depth:
            continue
        L = length * shrink ** g
        s = sw * shrink ** (0.5 * g)
        jit = rng.uniform(-0.25, 0.25) * spread
        if kind == "e":
            a_e, a_g = ang + spread * 0.55 + jit, ang - spread * 0.9 + jit
            q_e, q_g = p + L * _e(a_e), p + 0.8 * L * _e(a_g)
            gens[g].add(fline(p, q_e, color=cc, sw=s),
                        wavy(p, q_g, color=pc, sw=s * 0.85, amplitude=0.05, wavelength=0.28))
            frontier += [(q_e, a_e, "e", g + 1), (q_g, a_g, "gamma", g + 1)]
        else:
            a1, a2 = ang + spread * 0.7 + jit, ang - spread * 0.7 + jit
            q1, q2 = p + L * _e(a1), p + L * _e(a2)
            gens[g].add(fline(p, q1, color=cc, sw=s), fline(p, q2, color=cc, sw=s))
            frontier += [(q1, a1, "e", g + 1), (q2, a2, "e", g + 1)]
    return VGroup(*gens)


def hadron_blob(p, width: float = 0.5, height: float = 0.32, color=None, angle: float = 0.0) -> Ellipse:
    """Hadronisation: a soft grey blob where quarks turn into hadrons."""
    cc = col(PARTICLE["hadron"] if color is None else color)
    return Ellipse(width=width, height=height, fill_color=cc, fill_opacity=0.35,
                   stroke_color=cc, stroke_width=1.5).rotate(angle).move_to(_p3(p))


def pion_lines(p, angles, length: float = 1.3, color=None, sw: float = 3.5,
               tip: bool = True) -> VGroup:
    """Straight hadron lines fanning out of ``p`` at the given ``angles``
    (radians), with an arrowhead at 60 % (fermion-flow style, no arrow for
    neutral pions when ``tip=False``)."""
    cc = col(PARTICLE["pion"] if color is None else color)
    g = VGroup()
    for a in angles:
        ln = fline(p, _p3(p) + length * _e(a), color=cc, sw=sw)
        g.add(VGroup(ln, arrow_tip_on(ln, color=cc, at=0.6, tip_length=0.18)) if tip else VGroup(ln))
    return g


# ---------------------------------------------------------------------------
# CMS: the logo and the transverse slice share one geometry (logo units,
# outer muon radius = 0.996 of the logo square side). Radii from Izaak
# Neutelings' CMS_logo.tex; the muon system is subdivided into four stations
# and three yoke rings only in the "detail" layer of the slice.
# ---------------------------------------------------------------------------

LOGO_R = {"tib": 0.059, "tob": 0.1145, "ecal": 0.220, "band": 0.231,
          "hcal": 0.385, "magnet": 0.578, "muon": 0.996}
# the seven flat rings of the logo, inside out: (name, r_in, r_out, DETECTOR key)
LOGO_RINGS = [
    ("tracker",  0.0,             LOGO_R["tib"],    "tracker"),
    ("tob",      LOGO_R["tib"],   LOGO_R["tob"],    "tob"),
    ("ecal",     LOGO_R["tob"],   LOGO_R["ecal"],   "ecal"),
    ("band",     LOGO_R["ecal"],  LOGO_R["band"],   "band"),
    ("hcal",     LOGO_R["band"],  LOGO_R["hcal"],   "hcal"),
    ("solenoid", LOGO_R["hcal"],  LOGO_R["magnet"], "solenoid"),
    ("muon",     LOGO_R["magnet"], LOGO_R["muon"],  "muon"),
]
STATION_FRAC = [(0.578, 0.675), (0.740, 0.815), (0.865, 0.930), (0.960, 0.996)]   # MB1..MB4
YOKE_FRAC = [(0.675, 0.740), (0.815, 0.865), (0.930, 0.960)]
PIXEL_FRAC = (0.018, 0.028, 0.040)
STRIP_FRAC = (0.055, 0.072, 0.090, 0.108)
BEAMPIPE_FRAC = 0.010
N_ECAL, N_HCAL, N_SECTORS = 60, 36, 12
R_DET = 2.95                          # delivered outer radius (scene units)
DET_CENTER = np.array([0.0, -0.25, 0.0])   # top edge at y = 2.70, under the title band
_C15 = math.cos(math.radians(15.0))   # apothem / circumradius of a dodecagon

# the four muons of the logo: TikZ ``to[out,in,looseness]`` segments, endpoints
# in polar (deg, r) logo units, sampled in ``logo_muon_paths``
_LOGO_MUONS = [
    ((0.003, 0.0), [(55.8, 167.5, 0.95, (14.94, 0.45)), (-12.5, 170, 1.04, (11.65, 0.51)),
                    (-10, -134, 0.85, (10.4, 0.766)), (46, -132, 1.2, (13.16, 0.815))]),
    ((0.003, 0.0), [(56.8, 178, 0.933, (28.79, 0.45)), (-2, -179, 0.953, (24.98, 0.51)),
                    (1, -116, 0.845, (25.33, 0.954))]),
    ((0.002, 0.0), [(61.5, -158.5, 0.875, (41.45, 0.47)), (21.5, -167, 0.9, (37.0, 0.57)),
                    (13.0, -131, 0.77, (32.85, 1.075))]),
    ((0.006, 0.0), [(73.0, -147, 0.85, (49.33, 0.50)), (33.0, -149, 1.0, (46.74, 0.58)),
                    (31.0, -131, 0.80, (42.849, 1.2114))]),
]


def _bezier(p0, p1, p2, p3, n=30):
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3


def logo_muon_paths() -> list[np.ndarray]:
    """The four logo muon curves as point arrays in logo units (origin = the
    logo's bottom-left corner, side = 1)."""
    paths = []
    for start, segs in _LOGO_MUONS:
        p = np.array([start[0], start[1], 0.0])
        pts = [p]
        for out_deg, in_deg, loose, (deg, r) in segs:
            q = r * _e(math.radians(deg))
            d = float(np.linalg.norm(q - p)) * 0.3915 * loose
            c1 = p + d * _e(math.radians(out_deg))
            c2 = q + d * _e(math.radians(in_deg))
            pts += list(_bezier(p, c1, c2, q)[1:])
            p = q
        paths.append(np.array(pts))
    return paths


def _ring_style(key, sw=0.6, fo=1.0):
    return dict(fill_color=col(DETECTOR[key]["fill"]), fill_opacity=fo,
                stroke_color=col(DETECTOR[key]["stroke"]), stroke_width=sw)


def _ring(center, r_in, r_out, angle=TAU, start_angle=0.0, **style):
    if r_in <= 0:
        return Sector(radius=r_out, angle=angle, start_angle=start_angle, arc_center=center, **style) \
            if angle < TAU - 1e-6 else Circle(radius=r_out, arc_center=center, **style)
    if angle >= TAU - 1e-6:
        return Annulus(inner_radius=r_in, outer_radius=r_out, arc_center=center, **style)
    return AnnularSector(inner_radius=r_in, outer_radius=r_out, angle=angle,
                         start_angle=start_angle, arc_center=center, **style)


def logo_rings(center, side: float, angle: float = TAU / 4, start_angle: float = 0.0,
               sw: float = 0.6) -> VGroup:
    """The seven flat logo layers (inside out) as sectors of ``angle`` about
    ``center``, radii = ``LOGO_R`` x ``side``. ``angle=TAU`` gives the full
    rings that are the base layer of ``CMSSlice`` (pixel-identical), so the
    logo -> slice animation can hand over to a real slice at the end."""
    c = _p3(center)
    g = VGroup()
    for name, r0, r1, key in LOGO_RINGS:
        g.add(_ring(c, r0 * side, r1 * side, angle, start_angle, **_ring_style(key, sw=sw)))
    return g


def _dodeca_sector(center, r_in, r_out, k, gap=0.0, **style) -> Polygon:
    """Trapezoid k (of 12) of a dodecagonal annulus with apothems r_in, r_out;
    ``gap`` insets the box along both flat sides (muon chambers)."""
    a0, a1 = math.radians(15 + 30 * k), math.radians(45 + 30 * k)
    A, B = center + r_in / _C15 * _e(a0), center + r_in / _C15 * _e(a1)
    C, D = center + r_out / _C15 * _e(a1), center + r_out / _C15 * _e(a0)
    if gap > 0:
        u = (B - A) / np.linalg.norm(B - A)
        A, B, C, D = A + gap * u, B - gap * u, C - gap * u, D + gap * u
    return Polygon(A, B, C, D, **style)


class CMSLogo(VGroup):
    """The CMS logo: an off-white rounded square whose bottom-left corner is
    the beam spot, a quarter of the detector layers, four red muons bending in
    the field, and the wordmark. ``corner`` = bottom-left corner, ``side`` =
    square side (the outer muon radius is 0.996 side).

    Attributes: .bg, .frame, .rings (VGroup of 7 quarter layers), .muons,
    .wordmark, .corner, .side."""

    def __init__(self, corner=ORIGIN, side: float = 2.0, wordmark: bool = True, **kw):
        super().__init__(**kw)
        self.corner, self.side = _p3(corner), side
        c = self.corner
        centre = c + side / 2 * (RIGHT + UP)
        self.bg = RoundedRectangle(width=side, height=side, corner_radius=0.02 * side,
                                   fill_color=col(CMS["bg"]), fill_opacity=1.0,
                                   stroke_width=0).move_to(centre)
        self.rings = logo_rings(c, side, TAU / 4, 0.0)
        # clip to the square: the muon layer (r = 0.996 side) pokes out of the
        # square only near the axes, so the sectors already fit; the corner
        # region beyond r = side shows the off-white background as in the logo
        self.muons = logo_muons(c, side)
        self.frame = RoundedRectangle(width=side, height=side, corner_radius=0.02 * side,
                                      fill_opacity=0, stroke_color=col(INK),
                                      stroke_width=2.0).move_to(centre)
        self.add(self.bg, self.rings, self.muons)
        if wordmark:
            self.wordmark = text("CMS", color=CMS["font"], weight=BOLD)
            self.wordmark.scale_to_fit_height(0.135 * side)
            self.wordmark.move_to(c + np.array([0.058 * side, 0.888 * side, 0.0]), aligned_edge=UP + LEFT)
            self.add(self.wordmark)
        else:
            self.wordmark = None
        self.add(self.frame)


def logo_muons(corner, side: float, color=None, core_sw: float = 5.0, rim_sw: float = 8.0) -> VGroup:
    """The four red logo muons with their white rim (TikZ ``double``), scaled
    to ``side`` from ``corner``. Submobjects: (rim, core) per muon."""
    cc = col(CMS["muon_track"] if color is None else color)
    c = _p3(corner)
    g = VGroup()
    for pts in logo_muon_paths():
        P = [c + side * p for p in pts]
        rim = VMobject(stroke_color=col(WHITE), stroke_width=rim_sw).set_points_smoothly(P)
        core = VMobject(stroke_color=cc, stroke_width=core_sw).set_points_smoothly(P)
        g.add(VGroup(rim, core))
    return g


class CMSSlice(VGroup):
    """Schematic CMS barrel, transverse view, official logo colours and logo
    radii. Two layers of drawing: ``.base`` = the seven flat rings of the
    logo (what the logo animation produces), ``.detail`` = beam pipe, tracker
    rings, ECAL cells, HCAL towers, muon stations and yoke (what makes it a
    detector).

    Attributes (all in scene coordinates):
      .parts     ordered list of (name, VGroup) inside out, for build animations:
                 beampipe, tracker, tob, ecal, band, hcal, solenoid, muon,
                 mb1, yoke1, mb2, yoke2, mb3, yoke3, mb4
      .layer     name -> VGroup (also 'stations', 'yoke' = all of them)
      .radii     name -> (r_in, r_out); .stations, .yokes lists of (r_in, r_out)
      .ecal_cells / .hcal_towers  lists of AnnularSector, index = phi bin
      .point_at(r, phi), .mid_radius(name), .cell_index(layer, phi), .cells(layer)
    No subsystem names on screen (animation-only rule)."""

    def __init__(self, center=DET_CENTER, radius: float = R_DET, **kw):
        super().__init__(**kw)
        self.c = _p3(center)
        self.s = s = radius / LOGO_R["muon"]      # logo side that gives this outer radius
        self.radius = radius
        self.radii: dict[str, tuple[float, float]] = {}
        self.layer: dict[str, VGroup] = {}
        self.parts: list[tuple[str, VGroup]] = []
        c = self.c

        # -- base: the seven logo rings, full circle ---------------------------
        self.base = logo_rings(c, s, TAU, 0.0)
        base_by_name = {name: ring for (name, *_), ring in zip(LOGO_RINGS, self.base)}

        # -- detail ---------------------------------------------------------
        bp = VGroup(Circle(radius=BEAMPIPE_FRAC * s, arc_center=c, **_ring_style("beampipe", sw=1.5)))
        self._add_part("beampipe", bp, (0.0, BEAMPIPE_FRAC * s))
        tr_rings = VGroup(*[Circle(radius=r * s, arc_center=c, fill_opacity=0,
                                   stroke_color=col(DETECTOR["tracker"]["stroke"]),
                                   stroke_width=1.6 if r in PIXEL_FRAC else 1.1)
                            for r in PIXEL_FRAC + STRIP_FRAC
                            if r < LOGO_R["tib"]])   # rings under the tob annulus would be hidden in the composite
        self._add_part("tracker", VGroup(base_by_name["tracker"], tr_rings), (0.0, LOGO_R["tib"] * s))
        self._add_part("tob", VGroup(base_by_name["tob"]), (LOGO_R["tib"] * s, LOGO_R["tob"] * s))
        self.ecal_cells = self._segmented("ecal", (LOGO_R["tob"] * s, LOGO_R["ecal"] * s), N_ECAL,
                                          _ring_style("ecal", sw=0.7), base_by_name["ecal"])
        self._add_part("band", VGroup(base_by_name["band"]), (LOGO_R["ecal"] * s, LOGO_R["band"] * s))
        self.hcal_towers = self._segmented("hcal", (LOGO_R["band"] * s, LOGO_R["hcal"] * s), N_HCAL,
                                           _ring_style("hcal", sw=0.9), base_by_name["hcal"])
        self._add_part("solenoid", VGroup(base_by_name["solenoid"]),
                       (LOGO_R["hcal"] * s, LOGO_R["magnet"] * s))
        self._add_part("muon", VGroup(base_by_name["muon"]), (LOGO_R["magnet"] * s, LOGO_R["muon"] * s))
        # stations and yokes: dodecagonal, corners kept inside the blue circle
        self.stations = [(a * s, b * s) for a, b in STATION_FRAC]
        self.yokes = [(a * s, b * s) for a, b in YOKE_FRAC]
        st_style = dict(fill_color=col(shade(CMS["muon"], 0.10)), fill_opacity=1.0,
                        stroke_color=col(DETECTOR["muon"]["stroke"]), stroke_width=1.0)
        yk_style = _ring_style("yoke", sw=0.8)
        stations_all, yokes_all = VGroup(), VGroup()
        self.station_groups: list[VGroup] = []
        for i, (a, b) in enumerate(self.stations):
            g = VGroup(*[_dodeca_sector(c, a * _C15, b * _C15, k, gap=0.02 * s, **st_style)
                         for k in range(N_SECTORS)])
            self.station_groups.append(g); stations_all.add(g)
            self._add_part(f"mb{i + 1}", g, (a, b))
            if i < len(self.yokes):
                ya, yb = self.yokes[i]
                y = VGroup(*[_dodeca_sector(c, ya * _C15, yb * _C15, k, **yk_style)
                             for k in range(N_SECTORS)])
                yokes_all.add(y)
                self._add_part(f"yoke{i + 1}", y, (ya, yb))
        self.layer["stations"], self.layer["yoke"] = stations_all, yokes_all
        self.outer_radius = radius
        # everything that is not a base ring
        self.detail = VGroup(*[m for _, grp in self.parts for m in grp if m not in self.base])

    def _add_part(self, name, group, radii):
        self.parts.append((name, group))
        self.layer[name] = group
        self.radii[name] = radii
        self.add(group)

    def _segmented(self, name, radii, n, style, base_ring):
        r0, r1 = radii
        cells = [AnnularSector(inner_radius=r0, outer_radius=r1, angle=TAU / n,
                               start_angle=TAU * i / n, arc_center=self.c, **style)
                 for i in range(n)]
        self._add_part(name, VGroup(base_ring, *cells), (r0, r1))
        return cells

    # -- geometry helpers ---------------------------------------------------
    def point_at(self, radius: float, phi: float) -> np.ndarray:
        return self.c + radius * _e(phi)

    def mid_radius(self, name: str) -> float:
        r0, r1 = self.radii[name]
        return 0.5 * (r0 + r1)

    def cell_index(self, layer: str, phi: float) -> int:
        n = N_ECAL if layer == "ecal" else N_HCAL
        return int(math.floor((phi % TAU) / TAU * n)) % n

    def cells(self, layer: str):
        return self.ecal_cells if layer == "ecal" else self.hcal_towers


def field_kappa(det: CMSSlice, kappa: float, yoke_ratio: float = -0.6):
    """Curvature as a function of radius: ``kappa`` inside the solenoid, zero
    in the coil, reversed (return field) in the muon system."""
    r_sol_in, r_sol_out = det.radii["solenoid"]

    def k(r):
        if r < r_sol_in:
            return kappa
        if r < r_sol_out:
            return 0.0
        return yoke_ratio * kappa

    return k


def curved_track(center, r_end: float, angle0: float, kappa, color=INK,
                 sw: float = 3.5, n: int = 200, r_start: float = 0.0) -> VMobject:
    """A charged track in a solenoid field: circular arc of curvature ``kappa``
    (1/scene-unit, sign = charge; or a callable kappa(r), see ``field_kappa``)
    leaving ``center`` at ``angle0`` and ending at radius ``r_end``.
    Neutral particles: ``kappa=0``. The sampled points are kept in ``.pts``."""
    center = _p3(center)
    kf = kappa if callable(kappa) else (lambda r: kappa)
    pts, s, ds = [], 0.0, 0.01
    p, th = np.array([0.0, 0.0, 0.0]), angle0
    while np.linalg.norm(p) < r_end and s < 50:
        r = float(np.linalg.norm(p))
        if r >= r_start:
            pts.append(center + p)
        p = p + ds * _e(th)
        th += kf(r) * ds
        s += ds
    vm = VMobject(stroke_color=col(color), stroke_width=sw)
    vm.set_points_as_corners(pts[:: max(1, len(pts) // n)] + [pts[-1]])
    vm.pts = pts
    return vm


def track(det: CMSSlice, phi: float, kappa: float, r_end: float | None = None,
          color=INK, sw: float = 3.5, in_field: bool = True, r_start: float = 0.0,
          outline=None, outline_sw: float = 2.5):
    """Charged-particle track from the beam spot. ``r_end`` defaults to the
    outer radius; ``in_field`` applies the solenoid/yoke field map.
    ``outline`` (a colour, e.g. ``TRACK_OUTLINE["mumu"]``) draws a thin dark
    rim under the track, the logo's double-line style, so pale colours read
    on white; the result is then a VGroup(rim, core) that still carries ``.pts``."""
    r_end = det.outer_radius if r_end is None else r_end
    k = field_kappa(det, kappa) if in_field else kappa
    core = curved_track(det.c, r_end, phi, k, color=color, sw=sw, r_start=r_start)
    if outline is None:
        return core
    rim = curved_track(det.c, r_end, phi, k, color=outline, sw=sw + outline_sw, r_start=r_start)
    g = VGroup(rim, core)
    g.pts = core.pts
    return g


def neutral_track(det: CMSSlice, phi: float, r_end: float, color=INK, sw: float = 3.0,
                  dash_length: float = 0.12, r_start: float = 0.0) -> DashedLine:
    """Photon / neutrino: straight dashed line (no tracker hits)."""
    return DashedLine(det.point_at(r_start, phi), det.point_at(r_end, phi),
                      dash_length=dash_length, stroke_color=col(color), stroke_width=sw)


def calo_hit(det: CMSSlice, layer: str, phi: float, width: int = 1, color=None,
             opacity: float = 0.95, extent: float = 1.0, side_opacity: float = 0.45,
             stroke_width: float = 0.8, stroke_color=None) -> VGroup:
    """Energy deposit: the cell at ``phi`` (plus ``width-1`` neighbours, half
    on each side, drawn fainter) filled in the subsystem's DEPOSIT colour.
    ``extent`` < 1 fills only the inner fraction of the cell (a smaller deposit)."""
    cc = col(DEPOSIT[layer] if color is None else color)
    sc = cc if stroke_color is None else col(stroke_color)
    n = N_ECAL if layer == "ecal" else N_HCAL
    r0, r1 = det.radii[layer]
    i0 = det.cell_index(layer, phi)
    g = VGroup()
    half = (width - 1) // 2
    for d in range(-half, width - half):
        i = (i0 + d) % n
        op = opacity if d == 0 else side_opacity
        g.add(AnnularSector(inner_radius=r0, outer_radius=r0 + extent * (r1 - r0),
                            angle=TAU / n, start_angle=TAU * i / n, arc_center=det.c,
                            fill_color=cc, fill_opacity=op, stroke_color=sc,
                            stroke_width=stroke_width))
    return g


def muon_hits(det: CMSSlice, trk, color=None, width: float = 0.075) -> VGroup:
    """Hit boxes in every muon station the track crosses (DT-style segments):
    a filled rectangle of ``width`` along the local track direction, in the
    track colour with a darker outline so it stands out on the track itself."""
    cc = col(DEPOSIT["muon"] if color is None else color)
    g = VGroup()
    for a, b in det.stations:
        inside = [p for p in trk.pts if a + 0.015 <= np.linalg.norm(p - det.c) <= b - 0.015]
        if len(inside) >= 2:
            p0, p1 = inside[0], inside[-1]
            u = (p1 - p0) / np.linalg.norm(p1 - p0)
            nrm = np.array([-u[1], u[0], 0.0]) * width / 2
            g.add(Polygon(p0 + nrm, p1 + nrm, p1 - nrm, p0 - nrm, fill_color=cc,
                          fill_opacity=0.9, stroke_color=darken(cc, 0.35), stroke_width=1.2))
    return g


def tracker_hits(det: CMSSlice, trk, color=None, radius: float = 0.022) -> VGroup:
    """Dots where the track crosses the pixel and strip layers."""
    cc = col(DETECTOR["tracker"]["stroke"] if color is None else color)
    g = VGroup()
    rings = [r * det.s for r in PIXEL_FRAC + STRIP_FRAC]
    rs = np.array([np.linalg.norm(p - det.c) for p in trk.pts])
    for r in rings:
        idx = np.where((rs[:-1] < r) & (rs[1:] >= r))[0]
        if len(idx):
            g.add(Dot(trk.pts[int(idx[0])], radius=radius, color=cc))
    return g


def _phi_of(det, p):
    d = p - det.c
    return math.atan2(d[1], d[0])


def signature(det: CMSSlice, kind: str, phi: float, charge: int = -1, kappa: float = 0.35,
              color=None, sw: float = 3.5, seed: int = 0, prongs: int = 1, pi0: bool = True,
              brem: bool = False, r_end: float | None = None) -> VGroup:
    """Standard detector signatures on a CMSSlice, as one VGroup whose
    submobjects appear in physical order (track, then deposits/hits):
      'e'      bent track to the ECAL, EM cluster in the ECAL (ee green);
               ``brem=True`` adds a soft photon cell next to it
      'mu'     outlined track through everything, MIP dots in both calorimeters,
               stubs in the four muon stations (mumu gold)
      'gamma'  dashed line, ECAL cluster, no track
      'tau_h'  ``prongs`` (1 or 3) narrow tracks into the HCAL, ``pi0`` adds the
               photon pair in the ECAL (tautau red): a narrow jet
      'jet'    a spray of soft tracks, wide ECAL + HCAL deposits
      'nu'     dashed grey line leaving the detector
      'met'    missing transverse momentum: a slate dashed arrow from the centre
    ``charge`` sets the bending direction; ``color`` overrides the default."""
    rng = np.random.default_rng(seed)
    g = VGroup()
    if kind == "e":
        cc = col(CHANNEL["ee"] if color is None else color)
        t = track(det, phi, charge * kappa, r_end=det.radii["ecal"][0], color=cc, sw=sw,
                  outline=TRACK_OUTLINE["ee"])
        phi_end = _phi_of(det, t.pts[-1])
        g.add(t, calo_hit(det, "ecal", phi_end, width=3, color=cc, side_opacity=0.6,
                          stroke_color=TRACK_OUTLINE["ee"], stroke_width=1.4))
        if brem:   # bremsstrahlung photon: straight, on the outside of the bend
            phi_g = phi - charge * 0.045
            g.add(neutral_track(det, phi_g, det.radii["ecal"][0],
                                color=PARTICLE["photon"], sw=sw * 0.5, r_start=det.radii["tob"][0] * 0.6),
                  calo_hit(det, "ecal", phi_g, width=1,
                           color=PARTICLE["photon"], opacity=0.6, extent=0.7))
    elif kind == "mu":
        cc = col(CHANNEL["mumu"] if color is None else color)
        t = track(det, phi, charge * kappa, color=cc, sw=sw, outline=TRACK_OUTLINE["mumu"])
        mips = VGroup()
        for layer in ("ecal", "hcal"):
            r0 = det.radii[layer][0]
            p_in = next(p for p in t.pts if np.linalg.norm(p - det.c) >= r0)
            mips.add(calo_hit(det, layer, _phi_of(det, p_in), width=1, color=cc,
                              opacity=0.9, extent=0.35, stroke_color=TRACK_OUTLINE["mumu"],
                              stroke_width=1.2))
        g.add(t, mips, muon_hits(det, t, color=cc))
    elif kind == "gamma":
        cc = col(PARTICLE["photon"] if color is None else color)
        g.add(neutral_track(det, phi, det.radii["ecal"][0], color=cc, sw=sw * 0.8),
              calo_hit(det, "ecal", phi, width=3, color=cc))
    elif kind == "tau_h":
        cc = col(CHANNEL["tautau"] if color is None else color)
        tracks = VGroup()
        dphis = (0.0,) if prongs == 1 else (-0.07, 0.0, 0.07)
        for i, dphi in enumerate(dphis):
            q = charge if prongs == 1 else (charge, -charge, charge)[i]
            tracks.add(track(det, phi + dphi, q * kappa * (1.0 + rng.uniform(-0.3, 0.6)),
                             r_end=det.radii["hcal"][0], color=cc, sw=sw * 0.8,
                             outline=TRACK_OUTLINE["tautau"], outline_sw=1.8))
        phi_end = _phi_of(det, tracks[len(dphis) // 2].pts[-1])
        deposits = VGroup(calo_hit(det, "hcal", phi_end, width=1 if prongs == 1 else 3,
                                   color=cc, extent=0.85, side_opacity=0.5,
                                   stroke_color=TRACK_OUTLINE["tautau"], stroke_width=1.2))
        if pi0:    # pi0 -> gamma gamma: neutral, so straight along the tau axis
            deposits.add(calo_hit(det, "ecal", phi, width=2, color=cc,
                                  opacity=0.85, side_opacity=0.85, extent=0.9,
                                  stroke_color=TRACK_OUTLINE["tautau"], stroke_width=1.2))
        g.add(tracks, deposits)
    elif kind == "jet":
        cc = col(PARTICLE["hadron"] if color is None else color)
        sprays = VGroup()
        for _ in range(6):
            dphi = rng.uniform(-0.16, 0.16)
            r_end = det.radii["ecal"][0] if rng.uniform() < 0.4 else det.radii["hcal"][0]
            sprays.add(track(det, phi + dphi, rng.choice([-1, 1]) * rng.uniform(0.2, 0.7),
                             r_end=r_end, color=cc, sw=sw * 0.6))
        g.add(sprays, calo_hit(det, "ecal", phi, width=5, color=cc, extent=0.7, side_opacity=0.5),
              calo_hit(det, "hcal", phi, width=3, color=cc, extent=0.9, side_opacity=0.5))
    elif kind == "nu":
        cc = col(PARTICLE["nu"] if color is None else color)
        line = neutral_track(det, phi, det.outer_radius + 0.25, color=cc, sw=sw * 0.8)
        g.add(line, arrow_tip_on(Line(det.c, det.point_at(det.outer_radius + 0.25, phi)),
                                 color=cc, at=1.0, tip_length=0.18))
    elif kind == "met":
        cc = col(PARTICLE["met"] if color is None else color)
        r_end = det.radii["hcal"][1] if r_end is None else r_end
        line = neutral_track(det, phi, r_end, color=cc, sw=sw * 1.3, dash_length=0.18)
        g.add(line, arrow_tip_on(Line(det.c, det.point_at(r_end, phi)), color=cc, at=1.0,
                                 tip_length=0.28))
    else:
        raise ValueError(f"unknown signature kind {kind!r}")
    return g


# ---------------------------------------------------------------------------
# real events on the slice (section 4): the display pieces of the 4-02 reveal
# ---------------------------------------------------------------------------

def mini_slice(scale: float, center) -> CMSSlice:
    """A CMSSlice drawn directly at ``scale`` x the delivered size, centred at
    ``center``: the same points (to 1e-16) as ``place(CMSSlice(), (scale,
    center), DET_CENTER)``, but with ``.c``, ``.radii``, ``.s`` and
    ``.point_at`` correct — those are plain attributes that ``scale()`` /
    ``shift()`` never update. Use it when signatures are built *on* the small
    slice (section-2 icons, the transient pairs of mumu_e); chained frames
    keep building on a natural-scale slice and ``place()`` the result."""
    return CMSSlice(center=center, radius=R_DET * scale)


def kappa_from_pt(pt: float, k45: float = 0.28, lo: float = 0.10, hi: float = 0.55) -> float:
    """Track curvature (1/scene-unit) of a muon with transverse momentum
    ``pt`` [GeV]: the 1/p_T law of a solenoid (kappa = 0.3 B / p_T),
    calibrated so that 45 GeV — a typical Z muon — gives ``k45 = 0.28``, the
    curvature of the two schematic muons of clip 4-02 (``s4_zmumu.py``), and
    clipped to [lo, hi] so a soft muon never curls up and a hard one still
    visibly bends. The scale is exaggerated on purpose (a 45 GeV muon in
    3.8 T has a 40 m radius); only the *relative* curvature of two tracks is
    honest. Sign = charge * kappa, exactly as ``signature()`` applies it
    (``track(det, phi, charge * kappa)``)."""
    return float(min(hi, max(lo, k45 * (45.0 / float(pt)))))


def muon_pieces(det: CMSSlice, phi: float, charge: int, kappa: float, sw: float = 3.5) -> VGroup:
    """One reconstructed muon exactly as the 4-02 reveal draws it
    (``channel_common.show_signature``): VGroup(track, tracker hits,
    deposits) with deposits = the MIP cells in both calorimeters + the four
    station stubs, all from ``signature(det, "mu", ...)``. eta is not drawn
    (r-phi view). Attributes ``.trk``, ``.hits``, ``.deposits`` for a staged
    reveal (``Create(trk)``, ``FadeIn(hits)``, ``FadeIn(deposits)``)."""
    sig = signature(det, "mu", phi, charge=charge, kappa=kappa, sw=sw)
    trk, deposits = sig[0], VGroup(*sig[1:])
    hits = tracker_hits(det, trk, color=CHANNEL["mumu"])
    g = VGroup(trk, hits, deposits)
    g.trk, g.hits, g.deposits = trk, hits, deposits
    return g


def _muon_records(ev) -> list:
    if isinstance(ev, (list, tuple)):
        return list(ev)
    for key in ("mu", "muons"):
        if key in ev:
            return list(ev[key])
    raise KeyError("event record has neither 'mu' nor 'muons'")


def _charge(m) -> int:
    return int(np.sign(float(m["charge"])))


def event_from_json(det: CMSSlice, ev) -> VGroup:
    """A real event from ``zmumu_events.json`` drawn on ``det``: one
    ``muon_pieces`` per muon record (``ev["mu"]`` or ``ev["muons"]``, or a
    plain list of records; each needs ``phi`` [rad], ``charge`` [+-1], ``pt``
    [GeV]) at the real phi, bending by ``charge * kappa_from_pt(pt)``; every
    other field (eta, iso, idx, ...) is carried, not drawn. Build it on the
    slice at natural scale and ``place()`` it together with the detector.
    Returns VGroup(muon_0, muon_1, ...) in record order (p_T-ordered in the
    JSON) with ``.muons`` = the records and ``.record`` = ``ev``."""
    muons = _muon_records(ev)
    g = VGroup(*[muon_pieces(det, float(m["phi"]), _charge(m), kappa_from_pt(m["pt"]))
                 for m in muons])
    g.muons, g.record = muons, ev
    return g


def muon_in_jet(det: CMSSlice, phi: float, charge: int, kappa: float, seed: int = 0) -> VGroup:
    """A non-prompt muon: ``signature("jet")`` centred where the muon track
    enters the ECAL (so the jet's calorimeter cluster sits on the muon), then
    the muon itself (``muon_pieces``) on top, so its MIP cells and station
    stubs paint over the jet's deposits. Returns VGroup(jet, mu) with
    ``.jet``, ``.mu``."""
    mu = muon_pieces(det, phi, charge, kappa)
    r0 = det.radii["ecal"][0]
    p_in = next(p for p in mu.trk.pts if np.linalg.norm(p - det.c) >= r0)
    jet = signature(det, "jet", _phi_of(det, p_in), seed=seed)
    g = VGroup(jet, mu)
    g.jet, g.mu = jet, mu
    return g


def pair_from_json(det: CMSSlice, pair, seed: int = 0) -> VGroup:
    """A real same-sign tag + probe pair from ``zmumu_events.json``
    (``pair["tag"]`` / ``pair["probe"]`` records with ``phi``, ``charge``,
    ``pt``; ``pair["passes"]``; ``pair["jets"]`` with ``contains_probe``).
    The tag is a plain ``muon_pieces``; the probe is a ``muon_in_jet`` when
    ``pair["passes"]`` is False or a jet with ``contains_probe`` exists
    (anti-isolated: the muon sits in a jet), else a plain ``muon_pieces``.
    Both tracks bend the same way (same sign). Returns VGroup(tag, probe)
    with ``.tag``, ``.probe``, ``.in_jet``, ``.passes``, ``.record``."""
    tag, probe = pair["tag"], pair["probe"]
    passes = pair.get("passes")
    jets = pair.get("jets") or []
    in_jet = (passes is False) or any(bool(j.get("contains_probe")) for j in jets)
    t = muon_pieces(det, float(tag["phi"]), _charge(tag), kappa_from_pt(tag["pt"]))
    if in_jet:
        p = muon_in_jet(det, float(probe["phi"]), _charge(probe), kappa_from_pt(probe["pt"]), seed=seed)
    else:
        p = muon_pieces(det, float(probe["phi"]), _charge(probe), kappa_from_pt(probe["pt"]))
    g = VGroup(t, p)
    g.tag, g.probe, g.in_jet, g.passes, g.record = t, p, in_jet, passes, pair
    return g


# ---------------------------------------------------------------------------
# data plots: axes with our own ticks, and the primitives on top of them
# ---------------------------------------------------------------------------

_CAP_H: dict = {}     # cache: height of a capital letter / a MathTex digit at scale 1


def _cap_scale(kind: str, h: float) -> float:
    """Scale factor that makes a capital letter (``kind='text'``, deck font) or
    a MathTex digit (``kind='tex'``) ``h`` high, so labels with and without
    descenders / superscripts come out at one uniform size."""
    if kind not in _CAP_H:
        _CAP_H[kind] = (text("H") if kind == "text" else mathtex("1")).height
    return h / _CAP_H[kind]


def _fmt_exp(k) -> str:
    return "%d" % round(k) if abs(float(k) - round(k)) < 1e-9 else "%g" % k


class DataAxes(VGroup):
    """Plot frame with manual ticks and numeric labels in the deck font.
    Position it FIRST (``move_to``), then build content through ``.c2p``.

    Linear y (default): ``y_range = [y0, y1, step]`` in data units, ticks every
    ``step`` unless ``y_ticks`` is given; ``.y_base = y0``, ``.y_top = y1``;
    ``c2p`` is the plain Axes mapping (no clamping unless ``y_floor`` is set).

    Logarithmic y (``y_log=True``): ``y_range = [e0, e1, 1]`` in **exponents**
    (the axis spans 10**e0 .. 10**e1) and the Axes get
    ``y_axis_config={"scaling": LogBase()}``. ``c2p(x, y)`` still takes real
    values (counts), never exponents, and clamps ``y`` to ``.y_floor``
    (default ``.y_base = 10**e0``), so zero / empty bins sit on the axis
    instead of raising ``log(0)``. Ticks default to every decade; ``y_ticks``
    is then a list of exponents; labels are ``10^{k}`` in MathTex sized so the
    digits match the x labels (``y_fmt`` is ignored). ``.y_top = 10**e1``.

    ``y_floor`` (both kinds): the value every lower ``y`` is drawn at; ``None``
    = no clamping on a linear axis, ``y_base`` on a log axis.

    Attributes: ``.ax`` (manim Axes), ``.x_range/.y_range`` (as given),
    ``.x_ticks/.y_ticks`` (tick values; exponents on log), ``.x_ticks_v``,
    ``.y_ticks_v``, ``.x_labels``, ``.y_labels``, ``.frame`` (invisible
    Rectangle over the plot area), ``.x_title/.y_title`` (or None),
    ``.y_log``, ``.y_base``, ``.y_top``, ``.y_floor``, ``.x_length``,
    ``.y_length``, ``.x_fmt/.y_fmt``, ``.tick_label_h``, ``.axis_color``,
    ``.stroke_width``, ``.tick_len`` (so a companion panel can copy the style).
    ``move_frame_to(center)`` positions the plot area itself (``move_to``
    positions the bounding box incl. labels and titles)."""

    def __init__(self, x_range, y_range, x_length, y_length, x_ticks=None, y_ticks=None,
                 x_fmt="{:.0f}", y_fmt="{:.0f}", tick_label_h=0.22, axis_color=INK,
                 stroke_width=2.5, tick_len=0.09, show_x_labels=True, show_y_labels=True,
                 x_title: str | None = None, y_title: str | None = None, title_h=0.26,
                 title_buff=0.22, y_log: bool = False, y_floor: float | None = None, **kw):
        super().__init__(**kw)
        ink = col(axis_color)
        self.x_range, self.y_range = list(x_range), list(y_range)
        self.x_length, self.y_length = x_length, y_length
        self.x_fmt, self.y_fmt, self.tick_label_h = x_fmt, y_fmt, tick_label_h
        self.axis_color, self.stroke_width, self.tick_len = axis_color, stroke_width, tick_len
        self.y_log = bool(y_log)
        if self.y_log:
            self.y_base, self.y_top = 10.0 ** y_range[0], 10.0 ** y_range[1]
            self.y_floor = self.y_base if y_floor is None else float(y_floor)
            axes_kw = {"y_axis_config": {"scaling": LogBase()}}
        else:
            self.y_base, self.y_top = y_range[0], y_range[1]
            self.y_floor = None if y_floor is None else float(y_floor)
            axes_kw = {}
        self.ax = Axes(x_range=list(x_range), y_range=list(y_range),
                       x_length=x_length, y_length=y_length, tips=False,
                       axis_config={"stroke_color": ink, "stroke_width": stroke_width,
                                    "include_ticks": False}, **axes_kw)
        self.add(self.ax)
        if x_ticks is None:
            x_ticks = list(np.arange(x_range[0], x_range[1] + 1e-9, x_range[2]))
        if y_ticks is None:
            y_ticks = list(np.arange(y_range[0], y_range[1] + 1e-9, y_range[2]))
        self.x_ticks, self.y_ticks = list(x_ticks), list(y_ticks)
        x0 = x_range[0]
        y_axis_v = self.y_base                               # data value of the x axis line
        y_mid_v = self._yval((y_range[0] + y_range[1]) / 2)  # data value of the frame centre
        self.x_ticks_v, self.x_labels = VGroup(), VGroup()
        for xv in x_ticks:
            base = self.ax.c2p(xv, y_axis_v)
            self.x_ticks_v.add(Line(base, base + DOWN * tick_len, stroke_color=ink,
                                    stroke_width=stroke_width))
            if show_x_labels:
                lab = text(x_fmt.format(xv), color=ink).scale_to_fit_height(tick_label_h)
                lab.next_to(base, DOWN, buff=0.16)
                self.x_labels.add(lab)
        self.y_ticks_v, self.y_labels = VGroup(), VGroup()
        for yv in y_ticks:
            base = self.ax.c2p(x0, self._yval(yv))
            self.y_ticks_v.add(Line(base, base + LEFT * tick_len, stroke_color=ink,
                                    stroke_width=stroke_width))
            if show_y_labels:
                if self.y_log:      # decade label 10^k, digits as high as the x labels
                    lab = mathtex(r"10^{%s}" % _fmt_exp(yv), color=ink)
                    lab.scale(tick_label_h / lab[0][0].height)
                    lab.next_to(base, LEFT, buff=0.16)
                    lab.shift(UP * (base[1] - lab[0][0].get_center()[1]))   # centre the "10" on the tick
                else:
                    lab = text(y_fmt.format(yv), color=ink).scale_to_fit_height(tick_label_h)
                    lab.next_to(base, LEFT, buff=0.16)
                self.y_labels.add(lab)
        self.frame = Rectangle(width=x_length, height=y_length, stroke_width=0,
                               fill_opacity=0.0).move_to(
            self.ax.c2p((x_range[0] + x_range[1]) / 2, y_mid_v))
        self.add(self.x_ticks_v, self.y_ticks_v, self.x_labels, self.y_labels, self.frame)
        self.x_title = self.y_title = None
        if x_title is not None:      # axis titles are data labels (allowed)
            xt = mathtex(x_title, color=ink).scale_to_fit_height(title_h)
            xt.next_to(self.x_labels if len(self.x_labels) else self.x_ticks_v, DOWN, buff=title_buff)
            xt.set_x(self.ax.c2p((x_range[0] + x_range[1]) / 2, y_axis_v)[0])
            self.add(xt); self.x_title = xt
        if y_title is not None:
            yt = mathtex(y_title, color=ink).scale_to_fit_height(title_h).rotate(np.pi / 2)
            yt.next_to(self.y_labels if len(self.y_labels) else self.y_ticks_v, LEFT, buff=title_buff)
            yt.set_y(self.ax.c2p(x0, y_mid_v)[1])
            self.add(yt); self.y_title = yt

    def _yval(self, e):
        """Axis coordinate (exponent on log) -> data value."""
        return 10.0 ** e if self.y_log else e

    def move_frame_to(self, center):
        """Shift so the *plot area* (``.frame``) is centred at ``center``.
        ``move_to`` centres the bounding box, which includes labels and
        titles, so the axes themselves land off-centre; use this when a
        constant names the frame centre. Returns self."""
        self.shift(_p3(center) - self.frame.get_center())
        return self

    def c2p(self, x, y):
        """Data (x, y) -> scene point. ``y`` is a real value on both kinds of
        axis; it is clamped to ``.y_floor`` when that is set (always on log)."""
        if self.y_floor is not None:
            y = np.maximum(y, self.y_floor) if np.ndim(y) else max(float(y), self.y_floor)
        return self.ax.c2p(x, y)


def data_bar(dax: DataAxes, x_center, value, half_width, color=INK, base=None,
             fill_opacity=0.85, stroke_width=1.5) -> Rectangle:
    cc = col(color)
    base = dax.y_base if base is None else base
    p_base, p_top = dax.c2p(x_center, base), dax.c2p(x_center, value)
    p_l, p_r = dax.c2p(x_center - half_width, base), dax.c2p(x_center + half_width, base)
    bar = Rectangle(width=abs(p_r[0] - p_l[0]), height=max(abs(p_top[1] - p_base[1]), 1e-3),
                    stroke_color=cc, stroke_width=stroke_width, fill_color=cc,
                    fill_opacity=fill_opacity)
    bar.move_to([p_base[0], (p_base[1] + p_top[1]) / 2, 0.0])
    return bar


def data_trace(dax: DataAxes, xs, ys, color=INK, stroke_width=3.0, dashed_=False) -> VMobject:
    pts = [dax.c2p(float(x), float(y)) for x, y in zip(xs, ys)]
    line = VMobject(stroke_color=col(color), stroke_width=stroke_width)
    line.set_points_as_corners(pts)
    return DashedVMobject(line, num_dashes=max(12, len(pts) // 3)) if dashed_ else line


def step_hist(dax: DataAxes, edges, counts, color=INK, stroke_width=3.0, fill_opacity=0.0):
    """A histogram as a step outline (every edge twice: the polyline IS the binning).
    ``fill_opacity > 0`` returns a closed polygon down to the axis instead."""
    xs, ys = [], []
    for i, n in enumerate(counts):
        xs += [edges[i], edges[i + 1]]
        ys += [n, n]
    if fill_opacity > 0:
        y0 = dax.y_base
        pts = [dax.c2p(xs[0], y0)] + [dax.c2p(x, y) for x, y in zip(xs, ys)] + [dax.c2p(xs[-1], y0)]
        return Polygon(*pts, stroke_color=col(color), stroke_width=stroke_width,
                       fill_color=col(color), fill_opacity=fill_opacity)
    return data_trace(dax, xs, ys, color=color, stroke_width=stroke_width)


def data_band(dax: DataAxes, xs, y_lo, y_hi, color=INK, opacity=0.15) -> Polygon:
    top = [dax.c2p(float(x), float(y)) for x, y in zip(xs, y_hi)]
    bot = [dax.c2p(float(x), float(y)) for x, y in zip(reversed(list(xs)), reversed(list(y_lo)))]
    return Polygon(*top, *bot, stroke_width=0, fill_color=col(color), fill_opacity=opacity)


def data_dot(dax: DataAxes, x, y, color=INK, radius=0.05, **kw) -> Dot:
    return Dot(dax.c2p(float(x), float(y)), radius=radius, color=col(color), **kw)


def data_errorbar(dax: DataAxes, x, mean, err, cap=0.0, color=INK, stroke_width=2.5) -> VGroup:
    """Vertical mean +/- err bar. ``cap`` in data-x units (0 = no caps, HEP style)."""
    cc = col(color)
    lo, hi = dax.c2p(x, mean - err), dax.c2p(x, mean + err)
    g = VGroup(Line(lo, hi, stroke_color=cc, stroke_width=stroke_width))
    if cap > 0:
        d = abs(dax.c2p(x + cap, mean)[0] - dax.c2p(x - cap, mean)[0]) / 2
        g.add(Line(hi + LEFT * d, hi + RIGHT * d, stroke_color=cc, stroke_width=stroke_width),
              Line(lo + LEFT * d, lo + RIGHT * d, stroke_color=cc, stroke_width=stroke_width))
    return g


def vref_line(dax: DataAxes, x, color=GREY, stroke_width=2.0, num_dashes=28) -> DashedVMobject:
    """Dashed vertical reference line (e.g. m_Z) spanning the axis height."""
    line = Line(dax.c2p(x, dax.y_base), dax.c2p(x, dax.y_top),
                stroke_color=col(color), stroke_width=stroke_width)
    return DashedVMobject(line, num_dashes=num_dashes)


def stack_hist(dax: DataAxes, edges, layers, stroke_width: float = 0, fill_opacity: float = 0.95) -> VGroup:
    """Stacked histogram. ``layers = [(name, counts, colour), ...]`` **bottom-up**
    (the first tuple is the bottom of the stack; ``counts`` has len(edges)-1
    entries). Layer i is ONE ``Polygon`` between the cumulative step of the
    layers below it and the cumulative step including it, every edge twice
    (the outline IS the binning), built through ``dax.c2p`` so it is log-safe:
    the bottom of layer 0 and every empty bin are clamped to ``dax.y_floor``.
    An all-zero layer is a degenerate, invisible polygon at the floor with the
    same number of points, so two stacks made from the same layer list are
    always ``Transform``-able 1:1 (8 layers -> 8 layers) whatever the numbers.
    ``stroke_width`` (default 0: no hairlines) outlines each polygon in its
    own colour; degenerate layers always get stroke 0.
    Attributes: ``.names`` (list), ``.cum`` (array (n_layers+1, n_bins): row 0
    zeros, row i = sum of layers < i, last row = the total), ``.counts``
    (n_layers, n_bins), ``.edges``, ``.layer`` (name -> polygon), ``.total``
    (= ``.cum[-1]``). They describe the construction values (a Transform does
    not update them)."""
    edges = np.asarray(edges, dtype=float)
    nb = len(edges) - 1
    names = [lay[0] for lay in layers]
    counts = (np.array([np.asarray(lay[1], dtype=float) for lay in layers])
              if len(layers) else np.zeros((0, nb)))
    if len(layers) and counts.shape[1] != nb:
        raise ValueError(f"stack_hist: {counts.shape[1]} counts per layer for {nb} bins")
    cum = np.vstack([np.zeros(nb), np.cumsum(counts, axis=0)]) if len(layers) else np.zeros((1, nb))

    def step_pts(row):
        pts = []
        for i, v in enumerate(row):
            pts += [dax.c2p(edges[i], v), dax.c2p(edges[i + 1], v)]
        return pts

    g = VGroup()
    g.layer = {}
    for i, (name, _, colour) in enumerate(layers):
        cc = col(colour)
        degenerate = not bool(np.any(counts[i] > 0))
        top, bot = step_pts(cum[i + 1]), step_pts(cum[i])
        poly = Polygon(*top, *bot[::-1], stroke_color=cc,
                       stroke_width=0 if degenerate else stroke_width,
                       fill_color=cc, fill_opacity=fill_opacity)
        g.add(poly)
        g.layer[name] = poly
    g.names, g.cum, g.counts, g.edges = names, cum, counts, edges
    g.total = cum[-1]
    return g


def ratio_panel(dax_main: DataAxes, edges, num, den, center, y_range=(0.9, 1.1, 0.1),
                y_length: float = 0.9, num_err=None, dot_radius: float = 0.04,
                dot_color=SAMPLE["Data"], **axes_kw) -> VGroup:
    """Data / prediction panel under a main plot. A *linear* ``DataAxes`` with
    the main plot's ``x_range``, ``x_length``, ``x_ticks`` and tick style
    (``y_fmt="{:.1f}"``), ``move_to(center)`` FIRST and then shifted in x so
    its x axis lies exactly under the main plot's (bins line up; the y of
    ``center`` is kept, its x only approximately). Then a dashed GREY line at
    y = 1 and one ``data_dot`` per bin at (bin centre, num/den); bins with
    den == 0 are skipped (``.bins`` lists the drawn ones), so the dot count
    is the number of populated bins. ``num_err`` (per-bin absolute error on
    ``num``) adds a ``data_errorbar`` (err/den) under each dot. Extra
    ``axes_kw`` go to ``DataAxes`` (``x_title``, ``show_x_labels``, ...).
    Attributes: ``.dax``, ``.ref`` (the y=1 line), ``.dots`` (VGroup: one Dot
    per drawn bin, so ``Transform(rp.dots, rp2.dots)`` is 60 -> 60), ``.errs``
    (VGroup, empty without num_err), ``.ratio`` (array, nan where den == 0),
    ``.bins``."""
    num, den = np.asarray(num, dtype=float), np.asarray(den, dtype=float)
    edges = np.asarray(edges, dtype=float)
    kw = dict(x_ticks=dax_main.x_ticks, x_fmt=dax_main.x_fmt, y_fmt="{:.1f}",
              tick_label_h=dax_main.tick_label_h, axis_color=dax_main.axis_color,
              stroke_width=dax_main.stroke_width, tick_len=dax_main.tick_len)
    kw.update(axes_kw)
    dax = DataAxes(dax_main.x_range, list(y_range), dax_main.x_length, y_length, **kw)
    dax.move_to(_p3(center))
    x0, x1 = dax_main.x_range[0], dax_main.x_range[1]
    dax.shift(RIGHT * (dax_main.c2p(x0, dax_main.y_base)[0] - dax.c2p(x0, dax.y_base)[0]))
    ref = DashedVMobject(Line(dax.c2p(x0, 1.0), dax.c2p(x1, 1.0),
                              stroke_color=col(GREY), stroke_width=2.0), num_dashes=40)
    ok = den > 0
    ratio = np.full(len(den), np.nan)
    ratio[ok] = num[ok] / den[ok]
    xc = 0.5 * (edges[:-1] + edges[1:])
    dots, errs, bins = VGroup(), VGroup(), []
    for i in np.where(ok)[0]:
        dots.add(data_dot(dax, xc[i], ratio[i], color=dot_color, radius=dot_radius))
        if num_err is not None:
            errs.add(data_errorbar(dax, xc[i], ratio[i], float(num_err[i]) / den[i],
                                   color=dot_color, stroke_width=1.8))
        bins.append(int(i))
    g = VGroup(dax, ref, errs, dots)
    g.dax, g.ref, g.dots, g.errs, g.ratio, g.bins = dax, ref, dots, errs, ratio, bins
    return g


def colour_key(entries, swatch=(0.34, 0.22), label_h: float = 0.24, buff: float = 0.14,
               row_buff: float = 0.12, tex: bool = True) -> VGroup:
    """The one colour key a chained story may carry (at most 4 rows).
    ``entries = [(label, colour[, kind]), ...]`` top to bottom; ``kind`` is
    ``"fill"`` (square swatch: a stacked sample, the default), ``"dot"`` (data
    marker: dot with a short vertical bar), ``"line"`` or ``"dashed"`` (a short
    line: a fit / theory line). Labels are sample or process *symbols*
    (MathTex when ``tex``, else plain deck-font Text at one uniform size),
    never narrative. Rows are left-aligned around ORIGIN; move it afterwards.
    Attributes: ``.rows`` (VGroup of rows, each VGroup(anchor, swatch, label);
    use ``group=key.rows`` in a LaggedStart over them), ``.swatches``,
    ``.labels`` (lists)."""
    if len(entries) > 4:
        raise ValueError("colour_key: at most 4 rows (CLAUDE.md: one small key)")
    w, h = swatch
    rows, swatches, labels = VGroup(), [], []
    for e in entries:
        label, colour = e[0], e[1]
        kind = e[2] if len(e) > 2 else "fill"
        cc = col(colour)
        anchor = Rectangle(width=w, height=h, stroke_width=0, fill_opacity=0.0)   # alignment column
        if kind == "fill":
            sw_ = Rectangle(width=w, height=h, fill_color=cc, fill_opacity=0.95,
                            stroke_color=darken(cc, 0.25), stroke_width=1.0)
        elif kind == "dot":
            sw_ = VGroup(Line(DOWN * h / 2, UP * h / 2, stroke_color=cc, stroke_width=2.0),
                         Dot(ORIGIN, radius=0.05, color=cc))
        elif kind == "line":
            sw_ = Line(LEFT * w / 2, RIGHT * w / 2, stroke_color=cc, stroke_width=3.0)
        elif kind == "dashed":
            sw_ = DashedLine(LEFT * w / 2, RIGHT * w / 2, dash_length=0.07,
                             stroke_color=cc, stroke_width=3.0)
        else:
            raise ValueError(f"colour_key: unknown kind {kind!r}")
        sw_.move_to(anchor)
        lab = mathtex(label, color=INK) if tex else text(label, color=INK)
        lab.scale(_cap_scale("tex" if tex else "text", label_h))
        lab.next_to(anchor, RIGHT, buff=buff)
        rows.add(VGroup(anchor, sw_, lab))
        swatches.append(sw_); labels.append(lab)
    rows.arrange(DOWN, aligned_edge=LEFT, buff=row_buff).move_to(ORIGIN)
    g = VGroup(rows)
    g.rows, g.swatches, g.labels = rows, swatches, labels
    return g


def pull_plot(names, pulls, constraints, x_range=(-2.0, 2.0), x_length: float = 3.2,
              row_h: float = 0.42, label_h: float = 0.2, label_buff: float = 0.2,
              tex: bool = False, color=INK, band_color=LIGHT_GREY) -> VGroup:
    """Nuisance-parameter pulls (theta_hat - theta_0)/Delta theta with their
    post-fit constraints, one row per NP top to bottom in the order given:
    the NP name as the fit spells it (plain Text; MathTex when ``tex``, e.g.
    schematic ``\\theta_1``), a horizontal bar of half-length ``constraint``
    and a dot at ``pull``. ``constraints`` entries may be None (no bar).
    Behind the rows: the +-1 band (``band_color``), a dashed line at 0, a
    baseline with integer ticks and tick labels. Built around ORIGIN.
    Attributes: ``.rows`` (VGroup of VGroup(label, bar, dot); use
    ``group=pp.rows`` for a LaggedStart), ``.band``, ``.zero``, ``.baseline``,
    ``.ticks``, ``.tick_labels``, ``.labels`` (list), ``.x_of(v)`` -> scene x
    of pull value v (reads the baseline, valid after move / shift / scale)."""
    x0, x1 = float(x_range[0]), float(x_range[1])
    n = len(names)
    H = n * row_h
    ink = col(color)
    left, right = -x_length / 2, x_length / 2

    def xs(v):
        return left + (float(v) - x0) / (x1 - x0) * x_length

    band = Rectangle(width=xs(1.0) - xs(-1.0), height=H, fill_color=col(band_color),
                     fill_opacity=0.7, stroke_width=0).move_to([xs(0.0), 0.0, 0.0])
    zero = DashedLine([xs(0.0), -H / 2, 0], [xs(0.0), H / 2, 0], dash_length=0.1,
                      stroke_color=ink, stroke_width=1.5)
    baseline = Line([left, -H / 2, 0], [right, -H / 2, 0], stroke_color=ink, stroke_width=2.0)
    ticks, tick_labels = VGroup(), VGroup()
    for k in range(math.ceil(x0), math.floor(x1) + 1):
        p = np.array([xs(k), -H / 2, 0.0])
        ticks.add(Line(p, p + DOWN * 0.07, stroke_color=ink, stroke_width=2.0))
        tick_labels.add(text(f"{k:g}", color=ink).scale(_cap_scale("text", 0.16))
                        .next_to(p + DOWN * 0.07, DOWN, buff=0.08))
    rows, labels = VGroup(), []
    for i, (name, p, c) in enumerate(zip(names, pulls, constraints)):
        y = H / 2 - row_h * (i + 0.5)
        lab = (mathtex(name, color=ink) if tex else text(name, color=ink))
        lab.scale(_cap_scale("tex" if tex else "text", label_h))
        lab.next_to([left, y, 0], LEFT, buff=label_buff)
        bar = Line([xs(p - c), y, 0], [xs(p + c), y, 0], stroke_color=ink,
                   stroke_width=2.5 if c else 0) if c is not None else \
            Line([xs(p), y, 0], [xs(p), y, 0], stroke_color=ink, stroke_width=0)
        dot = Dot([xs(p), y, 0], radius=0.055, color=ink)
        rows.add(VGroup(lab, bar, dot))
        labels.append(lab)
    g = VGroup(band, zero, baseline, ticks, tick_labels, rows)
    g.band, g.zero, g.baseline, g.ticks, g.tick_labels, g.rows, g.labels = \
        band, zero, baseline, ticks, tick_labels, rows, labels
    g.x_range = (x0, x1)

    def x_of(v):
        a, b = g.baseline.get_start()[0], g.baseline.get_end()[0]
        return a + (float(v) - x0) / (x1 - x0) * (b - a)

    g.x_of = x_of
    return g


def value_grid(values, row_labels, col_labels, vmax=None, fmt: str = "{:.2f}",
               cell=(0.70, 0.36), label_h: float = 0.16, color=SLATE, vmin: float = 0.0,
               tex: bool = True) -> VGroup:
    """A table of numbers as a heat map (e.g. the fake-factor map, rows = p_T
    bins, columns = |eta| bins). Cell (i, j) is filled with a tint of ``color``
    proportional to (v - vmin)/(vmax - vmin) — ``vmax`` is taken from the data
    when None, never assumed — and the value is printed with ``fmt`` (Text:
    ink on light cells, white on dark). Row labels sit left of the first
    column, column labels above the first row (MathTex when ``tex``: bin
    ranges / symbols). NaN cells are left blank. Built around ORIGIN.
    Attributes: ``.cells`` (VGroup, row-major Rectangles), ``.texts`` (VGroup,
    same order, blank cells skipped), ``.row_labels``, ``.col_labels``
    (VGroups), ``.values`` (array), ``.vmin``, ``.vmax``."""
    V = np.asarray(values, dtype=float)
    nr, nc = V.shape
    vmax = float(np.nanmax(V)) if vmax is None else float(vmax)
    span = vmax - vmin if vmax > vmin else 1.0
    w, h = cell
    cells, texts = VGroup(), VGroup()
    for i in range(nr):
        for j in range(nc):
            c = np.array([(j - (nc - 1) / 2) * w, ((nr - 1) / 2 - i) * h, 0.0])
            v = V[i, j]
            a = 0.0 if np.isnan(v) else float(np.clip((v - vmin) / span, 0.0, 1.0))
            fill = WHITE if np.isnan(v) else mix(WHITE, col(color).to_hex(), 0.08 + 0.80 * a)
            cells.add(Rectangle(width=w, height=h, fill_color=col(fill), fill_opacity=1.0,
                                stroke_color=col(WHITE), stroke_width=1.5).move_to(c))
            if not np.isnan(v):
                texts.add(text(fmt.format(v), color=WHITE if a > 0.55 else INK)
                          .scale(_cap_scale("text", label_h)).move_to(c))
    rl, cl = VGroup(), VGroup()
    for i, s in enumerate(row_labels):
        lab = mathtex(s, color=INK) if tex else text(s, color=INK)
        lab.scale(_cap_scale("tex" if tex else "text", label_h))
        lab.next_to(cells[i * nc], LEFT, buff=0.14)
        rl.add(lab)
    for j, s in enumerate(col_labels):
        lab = mathtex(s, color=INK) if tex else text(s, color=INK)
        lab.scale(_cap_scale("tex" if tex else "text", label_h))
        lab.next_to(cells[j], UP, buff=0.10)
        cl.add(lab)
    g = VGroup(cells, texts, rl, cl)
    g.cells, g.texts, g.row_labels, g.col_labels = cells, texts, rl, cl
    g.values, g.vmin, g.vmax = V, vmin, vmax
    return g


class Slider(VGroup):
    """See ``slider()``."""

    def x_of(self, v) -> float:
        """Scene x of value ``v`` on the axis (reads the axis line, so it stays
        valid after move / shift / scale)."""
        a, b = self.axis.get_start()[0], self.axis.get_end()[0]
        return a + (float(v) - self.lo) / (self.hi - self.lo) * (b - a)

    def marker_at(self, value, err=None) -> VGroup:
        """A fresh marker (VGroup(bar, dot)) at ``value`` +- ``err`` on the
        axis as it is now, for ``Transform(sl.marker, sl.marker_at(v, e))``."""
        y = self.axis.get_center()[1]
        cc = col(self.color)
        e = 0.0 if not err else float(err)
        bar = Line([self.x_of(value - e), y, 0], [self.x_of(value + e), y, 0],
                   stroke_color=cc, stroke_width=self.bar_sw if e else 0)
        dot = Dot([self.x_of(value), y, 0], radius=self.marker_r, color=cc)
        g = VGroup(bar, dot)
        g.bar, g.dot = bar, dot
        return g


def slider(lo, hi, value, err=None, ref=None, length: float = 3.2, color=INK, ref_color=THEORY,
           ticks=None, fmt: str = "{:.2f}", tick_label_h: float = 0.2, marker_r: float = 0.08,
           bar_sw: float = 4.0, ref_h: float = 0.5) -> Slider:
    """A one-parameter meter (mu_Z, epsilon_data / epsilon_MC, ...): a
    horizontal axis from ``lo`` to ``hi`` with ticks at ``ticks`` (default:
    lo, ref, hi) labelled with the numbers (``fmt``), a dashed vertical line
    rising ``ref_h`` above the axis at ``ref`` in ``ref_color`` (purple = the
    expectation / theory; it stays clear of the tick labels below), and ``.marker`` = VGroup(bar, dot) at ``value`` +- ``err`` (bar
    of stroke 0 when err is None / 0, so markers always Transform 1:1).
    Built around ORIGIN; move / scale afterwards. ``x_of(v)`` gives the scene
    x of a value, ``marker_at(v, e)`` a fresh marker for
    ``Transform(sl.marker, sl.marker_at(...))`` (the marker slides, the bar
    shrinks). Attributes: ``.axis`` (Line), ``.ticks``, ``.labels``, ``.ref``
    (or None), ``.marker``, ``.bar``, ``.dot``, ``.lo``, ``.hi``, ``.value``,
    ``.err``, ``.ref_value``."""
    sl = Slider()
    sl.lo, sl.hi, sl.value, sl.err, sl.ref_value = float(lo), float(hi), float(value), err, ref
    sl.color, sl.marker_r, sl.bar_sw = color, marker_r, bar_sw
    ink = col(color)
    sl.axis = Line(LEFT * length / 2, RIGHT * length / 2, stroke_color=ink, stroke_width=2.5)
    sl.add(sl.axis)
    if ticks is None:
        ticks = [lo, hi] if ref is None else sorted({float(lo), float(ref), float(hi)})
    sl.ticks, sl.labels = VGroup(), VGroup()
    for tv in ticks:
        p = np.array([sl.x_of(tv), 0.0, 0.0])
        sl.ticks.add(Line(p, p + DOWN * 0.09, stroke_color=ink, stroke_width=2.5))
        sl.labels.add(text(fmt.format(tv), color=ink).scale_to_fit_height(tick_label_h)
                      .next_to(p + DOWN * 0.09, DOWN, buff=0.1))
    sl.add(sl.ticks, sl.labels)
    sl.ref = None
    if ref is not None:
        sl.ref = DashedLine([sl.x_of(ref), -0.06, 0], [sl.x_of(ref), ref_h, 0],
                            dash_length=0.1, stroke_color=col(ref_color), stroke_width=2.5)
        sl.add(sl.ref)
    sl.marker = sl.marker_at(value, err)
    sl.bar, sl.dot = sl.marker.bar, sl.marker.dot
    sl.add(sl.marker)
    return sl


# ---------------------------------------------------------------------------
# line shapes for explainer histograms (schematic, seeded)
# ---------------------------------------------------------------------------

def breit_wigner(m, m0=91.1876, gamma=2.4952):
    return 1.0 / ((m * m - m0 * m0) ** 2 + m0 * m0 * gamma * gamma)


def dscb(t, a_l, n_l, a_r, n_r):
    """Double-sided Crystal Ball, unit height at t=0."""
    if t < -a_l:
        return (n_l / a_l) ** n_l * math.exp(-0.5 * a_l * a_l) * (n_l / a_l - a_l - t) ** (-n_l)
    if t > a_r:
        return (n_r / a_r) ** n_r * math.exp(-0.5 * a_r * a_r) * (n_r / a_r - a_r + t) ** (-n_r)
    return math.exp(-0.5 * t * t)


def schematic_zpeak(edges, peak=91.19, sigma=2.4, tail=(0.9, 3.0), height=1.0,
                    seed=11, jitter=0.05, floor=0.01):
    """Seeded bin heights of a Z peak with a radiative left tail. Same seed ->
    same jitter, so a 'before/after' pair differs only where the physics does."""
    rng = np.random.default_rng(seed)
    centres = 0.5 * (np.asarray(edges[:-1]) + np.asarray(edges[1:]))
    h = np.array([height * dscb((x - peak) / sigma, *tail, 1.6, 5.0) for x in centres])
    h *= 1.0 + jitter * rng.standard_normal(len(h))
    h[h < floor * height] = 0.0
    return h


# ---------------------------------------------------------------------------
# story primitives: clock, event rain, the pipeline spine
# ---------------------------------------------------------------------------

def clock(center, radius: float = 0.32, color=INK, sw: float = 2.5) -> VGroup:
    """A small clock face whose hand follows ``.turns`` (a ValueTracker in
    full turns, clockwise from 12): the hand is an ``always_redraw``, so
    ``clk.turns.animate.set_value(3)`` in a play spins it. That puts the
    tracker into ``scene.mobjects``: ``scene.remove(clk.turns)`` after the
    play (recipes trap 7). Attributes: ``.face``, ``.ticks``, ``.hand``,
    ``.pivot``, ``.turns``."""
    c = _p3(center)
    ink = col(color)
    g = VGroup()
    g.face = Circle(radius=radius, arc_center=c, stroke_color=ink, stroke_width=sw,
                    fill_color=col(WHITE), fill_opacity=1.0)
    g.ticks = VGroup(*[Line(c + 0.82 * radius * _e(a), c + radius * _e(a), stroke_color=ink,
                            stroke_width=sw * 0.8) for a in (0, TAU / 4, TAU / 2, 3 * TAU / 4)])
    g.turns = ValueTracker(0.0)
    g.hand = always_redraw(lambda: Line(
        c, c + 0.72 * radius * _e(TAU / 4 - TAU * g.turns.get_value()),
        stroke_color=ink, stroke_width=sw * 1.2))
    g.pivot = Dot(c, radius=0.035, color=ink)
    g.add(g.face, g.ticks, g.hand, g.pivot)
    return g


def rain(det: CMSSlice, dax: DataAxes, edges, counts, n: int = 60, seed: int = 3, color=INK,
         radius: float = 0.045, flight: float = 0.55, fade: float = 0.2,
         r_max: float | None = None) -> list:
    """Event dots streaming from the detector into the histogram (recipes
    section 4): ``n`` dots start at seeded random points inside the tracker of
    ``det`` *as drawn now* (centre and radius read from the mobject, so a
    ``place()``d slice works; ``r_max`` defaults to the tracker's outer
    radius) and fly (``flight`` s, ease-in) to a bin sampled with
    ``p = counts / counts.sum()``, landing on that bin's final height
    (``dax.c2p`` clamps empties on log), then fade (``fade`` s). Returns one
    ``Succession(move, FadeOut)`` per dot for
    ``LaggedStart(*rain(...), lag_ratio=0.06)``; the dot of animation ``a`` is
    ``a.animations[0].mobject``. Each Succession leaves an empty Group
    placeholder in the scene afterwards, which ``check_order`` skips."""
    rng = np.random.default_rng(seed)
    counts = np.asarray(counts, dtype=float)
    edges = np.asarray(edges, dtype=float)
    p = counts / counts.sum() if counts.sum() > 0 else np.full(len(counts), 1.0 / len(counts))
    c = det.get_center()
    R = det.width / 2
    r_max = R * LOGO_R["tob"] / LOGO_R["muon"] if r_max is None else r_max
    anims = []
    for _ in range(n):
        r, ph = r_max * math.sqrt(rng.uniform()), rng.uniform(0, TAU)
        j = int(rng.choice(len(counts), p=p))
        x = rng.uniform(edges[j], edges[j + 1])
        dot = Dot(c + r * _e(ph), radius=radius, color=col(color))
        dst = dax.c2p(x, counts[j])
        anims.append(Succession(
            dot.animate(run_time=flight, rate_func=rate_functions.ease_in_quad).move_to(dst),
            FadeOut(dot, run_time=fade)))
    return anims


def spine(node_icons, xs, y: float = 0.0, box=(1.24, 1.24), done=(), accent=DETECTOR_ACCENT,
          arrow_sw: float = 3.0, box_sw: float = 2.5, fit: float = 0.80) -> VGroup:
    """The analysis-pipeline strip of section 2: one rounded box per node at
    (xs[i], y), ``node_icons[i]`` shrunk (never enlarged) to fit ``fit`` x the
    box and centred in it, and an arrow between consecutive boxes. Nodes whose
    index is in ``done`` are tinted with ``accent`` (the stages already told);
    the others are white. Build order = z-order: boxes, arrows, icons.
    Attributes: ``.boxes`` (VGroup), ``.arrows`` (VGroup of VGroup(line, tip)),
    ``.icons`` (VGroup: the objects passed in), ``.nodes`` (list of
    VGroup(box, icon) *views* — not added, for ``about_point`` / ``get_center``),
    ``.centers`` (array (n, 3)), ``.done`` (tuple)."""
    w, h = box
    n = len(node_icons)
    if len(xs) != n:
        raise ValueError("spine: one x per icon")
    acc = col(accent)
    boxes, icons, arrows, nodes = VGroup(), VGroup(), VGroup(), []
    centers = np.array([[float(x), float(y), 0.0] for x in xs])
    for i, (icon, c) in enumerate(zip(node_icons, centers)):
        fill = lighten(acc, 0.80) if i in done else col(WHITE)
        b = RoundedRectangle(width=w, height=h, corner_radius=0.12, stroke_color=acc,
                             stroke_width=box_sw, fill_color=fill, fill_opacity=1.0).move_to(c)
        f = min(1.0, fit * w / max(icon.width, 1e-6), fit * h / max(icon.height, 1e-6))
        if f < 1.0:
            icon.scale(f)
        icon.move_to(c)
        boxes.add(b); icons.add(icon); nodes.append(VGroup(b, icon))
    for i in range(n - 1):
        a = centers[i] + RIGHT * (w / 2 + 0.08)
        b = centers[i + 1] + LEFT * (w / 2 + 0.08)
        ln = Line(a, b, stroke_color=col(INK), stroke_width=arrow_sw)
        arrows.add(VGroup(ln, arrow_tip_on(ln, color=INK, at=1.0, tip_length=0.2)))
    g = VGroup(boxes, arrows, icons)
    g.boxes, g.arrows, g.icons, g.nodes, g.centers, g.done = boxes, arrows, icons, nodes, centers, tuple(done)
    return g


# ---------------------------------------------------------------------------
# chaining helpers
# ---------------------------------------------------------------------------

def place(mobj, placement, natural_center=ORIGIN):
    """Affine placement shared by chained clips: ``placement = (scale, center)``
    (scale about ``natural_center``, then move that point to ``center``; both
    centres may be 2- or 3-vectors). Returns ``mobj``."""
    s, center = placement
    natural_center = _p3(natural_center)
    mobj.scale(s, about_point=natural_center)
    mobj.shift(_p3(center) - natural_center)
    return mobj


def add_state(scene, state: dict, order) -> None:
    """Open a chained clip on the previous clip's final frame: add the
    builder's mobjects ``state[k]`` for ``k in order`` (build order =
    z-order)."""
    for k in order:
        scene.add(state[k])


def _is_placeholder(m) -> bool:
    """Empty Mobject left by ``wait()`` or by a finished Succession/FadeOut."""
    return len(m.get_family()) == 1 and not m.has_points()


def check_order(scene, state: dict, order) -> None:
    """Z-order guard at the end of a chained clip (recipes section 5, rule 1):
    the live ``scene.mobjects`` — skipping wait placeholders, empty groups and
    ValueTrackers — must be exactly ``[state[k] for k in order]``, by identity.
    Raises AssertionError with the names it found otherwise, so manim exits
    before combining the partial files."""
    live = [m for m in scene.mobjects
            if not _is_placeholder(m) and not isinstance(m, ValueTracker)]
    want = [state[k] for k in order]
    if len(live) != len(want) or any(a is not b for a, b in zip(live, want)):
        by_id = {id(v): k for k, v in state.items()}
        got = [by_id.get(id(m), f"<{type(m).__name__}>") for m in live]
        raise AssertionError(f"scene z-order {got} != ORDER {list(order)}")


def clip_open(scene, name: str, hold: float = 0.3) -> None:
    """Open delivered clip ``name`` (a manim section; tools/deliver_chain.py writes
    every section as its own MP4) on a short hold of the current frame."""
    scene.next_section(name)
    scene.wait(hold)


def clip_cut(scene, name: str, tail: float = 0.1, hold: float = 0.3) -> None:
    """End the running clip on its last change (``tail``) and open clip ``name``
    on that same frame (``hold``): the two clips share the frame exactly. A name
    repeated at the end of one scene and the start of the next is one clip
    (joined by tools/deliver_chain.py)."""
    scene.wait(tail)
    clip_open(scene, name, hold)


__all__ = [
    "CHAPTER", "CHAPTER_NAME", "PURPLE", "CYAN", "GREEN", "GOLD", "RED", "SLATE", "BLACK", "WHITE",
    "mix", "tint", "shade",
    "BG", "INK", "GREY", "LIGHT_GREY", "CHANNEL", "CHANNEL_LINE", "TRACK_OUTLINE", "COMBINED",
    "THEORY", "DETECTOR_ACCENT", "HIGHLIGHT", "SAMPLE", "PARTICLE", "CMS", "DETECTOR", "DEPOSIT", "FLAG",
    "FONT",
    "col", "lighten", "darken", "text", "mathtex", "white_background", "load_data",
    "unit", "fline", "fermion", "arrow_tip_on", "dashed", "wavy", "gluon", "vertex_dot",
    "shower_tree", "hadron_blob", "pion_lines",
    "LOGO_R", "LOGO_RINGS", "STATION_FRAC", "YOKE_FRAC", "N_ECAL", "N_HCAL", "N_SECTORS",
    "R_DET", "DET_CENTER", "logo_muon_paths", "logo_rings", "logo_muons", "CMSLogo", "CMSSlice",
    "field_kappa", "curved_track", "track", "neutral_track", "calo_hit", "muon_hits",
    "tracker_hits", "signature",
    "DataAxes", "data_bar", "data_trace", "step_hist", "data_band", "data_dot",
    "data_errorbar", "vref_line", "breit_wigner", "dscb", "schematic_zpeak", "place",
    "mini_slice", "kappa_from_pt", "muon_pieces", "event_from_json", "muon_in_jet",
    "pair_from_json",
    "stack_hist", "ratio_panel", "colour_key", "pull_plot", "value_grid", "Slider", "slider",
    "clock", "rain", "spine", "add_state", "check_order", "clip_open", "clip_cut",
]
