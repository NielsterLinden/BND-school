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
    AnnularSector, Annulus, Axes, Circle, DashedLine, DashedVMobject, Dot, Ellipse, Line,
    ManimColor, MathTex, ParametricFunction, Polygon, Rectangle, RoundedRectangle, Sector,
    Text, VGroup, VMobject,
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
# data plots: axes with our own ticks, and the primitives on top of them
# ---------------------------------------------------------------------------

class DataAxes(VGroup):
    """Plot frame with manual ticks and numeric labels in the deck font.
    Position it FIRST (``move_to``), then build content through ``.c2p``."""

    def __init__(self, x_range, y_range, x_length, y_length, x_ticks=None, y_ticks=None,
                 x_fmt="{:.0f}", y_fmt="{:.0f}", tick_label_h=0.22, axis_color=INK,
                 stroke_width=2.5, tick_len=0.09, show_x_labels=True, show_y_labels=True,
                 x_title: str | None = None, y_title: str | None = None, title_h=0.26,
                 title_buff=0.22, **kw):
        super().__init__(**kw)
        ink = col(axis_color)
        self.x_range, self.y_range = list(x_range), list(y_range)
        self.ax = Axes(x_range=list(x_range), y_range=list(y_range),
                       x_length=x_length, y_length=y_length, tips=False,
                       axis_config={"stroke_color": ink, "stroke_width": stroke_width,
                                    "include_ticks": False})
        self.add(self.ax)
        if x_ticks is None:
            x_ticks = list(np.arange(x_range[0], x_range[1] + 1e-9, x_range[2]))
        if y_ticks is None:
            y_ticks = list(np.arange(y_range[0], y_range[1] + 1e-9, y_range[2]))
        x0, y0 = x_range[0], y_range[0]
        self.x_ticks_v, self.x_labels = VGroup(), VGroup()
        for xv in x_ticks:
            base = self.ax.c2p(xv, y0)
            self.x_ticks_v.add(Line(base, base + DOWN * tick_len, stroke_color=ink,
                                    stroke_width=stroke_width))
            if show_x_labels:
                lab = text(x_fmt.format(xv), color=ink).scale_to_fit_height(tick_label_h)
                lab.next_to(base, DOWN, buff=0.16)
                self.x_labels.add(lab)
        self.y_ticks_v, self.y_labels = VGroup(), VGroup()
        for yv in y_ticks:
            base = self.ax.c2p(x0, yv)
            self.y_ticks_v.add(Line(base, base + LEFT * tick_len, stroke_color=ink,
                                    stroke_width=stroke_width))
            if show_y_labels:
                lab = text(y_fmt.format(yv), color=ink).scale_to_fit_height(tick_label_h)
                lab.next_to(base, LEFT, buff=0.16)
                self.y_labels.add(lab)
        self.frame = Rectangle(width=x_length, height=y_length, stroke_width=0,
                               fill_opacity=0.0).move_to(
            self.ax.c2p((x_range[0] + x_range[1]) / 2, (y_range[0] + y_range[1]) / 2))
        self.add(self.x_ticks_v, self.y_ticks_v, self.x_labels, self.y_labels, self.frame)
        self.x_title = self.y_title = None
        if x_title is not None:      # axis titles are data labels (allowed)
            xt = mathtex(x_title, color=ink).scale_to_fit_height(title_h)
            xt.next_to(self.x_labels if len(self.x_labels) else self.x_ticks_v, DOWN, buff=title_buff)
            xt.set_x(self.ax.c2p((x_range[0] + x_range[1]) / 2, y0)[0])
            self.add(xt); self.x_title = xt
        if y_title is not None:
            yt = mathtex(y_title, color=ink).scale_to_fit_height(title_h).rotate(np.pi / 2)
            yt.next_to(self.y_labels if len(self.y_labels) else self.y_ticks_v, LEFT, buff=title_buff)
            yt.set_y(self.ax.c2p(x0, (y_range[0] + y_range[1]) / 2)[1])
            self.add(yt); self.y_title = yt

    def c2p(self, x, y):
        return self.ax.c2p(x, y)


def data_bar(dax: DataAxes, x_center, value, half_width, color=INK, base=None,
             fill_opacity=0.85, stroke_width=1.5) -> Rectangle:
    cc = col(color)
    base = dax.y_range[0] if base is None else base
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
        y0 = dax.y_range[0]
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
    line = Line(dax.c2p(x, dax.y_range[0]), dax.c2p(x, dax.y_range[1]),
                stroke_color=col(color), stroke_width=stroke_width)
    return DashedVMobject(line, num_dashes=num_dashes)


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
# chaining helper
# ---------------------------------------------------------------------------

def place(mobj, placement, natural_center=ORIGIN):
    """Affine placement shared by chained clips: ``placement = (scale, center)``."""
    s, center = placement
    mobj.scale(s, about_point=natural_center)
    mobj.shift(np.asarray(center, dtype=float) - np.asarray(natural_center, dtype=float))
    return mobj


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
]
