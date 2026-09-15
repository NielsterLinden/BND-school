"""Section 7: the closing slide.

    MicDrop   the CMS slice is already there. A hand holding a microphone comes
              in from the top; the fingers open, the microphone falls into the
              centre of the detector and bursts into rays in the colours of the
              three BND-school flags (NL, BE, DE). The rays fly outward off the
              frame and the detector breaks into wedges that fly off with them;
              on the clean slide "Thank you!" / "Any questions?" appears, centred
              (the one clip where words are allowed: the closing slide).

Colours: the hand is line art (off-white fill, ink stroke), the microphone is
slate; the rays use ``FLAG`` from the palette (outro only).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    BOLD, DEGREES, DOWN, ORIGIN, TAU, UP, AnimationGroup, Annulus, AnnularSector, ApplyFunction, Arc,
    Circle, Dot, FadeIn, FadeOut, Line, Polygon, RoundedRectangle, Scene, Sector, Succession, Transform,
    VGroup, ValueTracker, Wait, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403

EASE = rate_functions.ease_in_out_sine
HAND_FILL = CMS["bg"]                      # off-white line-art hand
HAND_STROKE = INK

# microphone (head down, handle up), placed by the centre of its head
MIC_HEAD_R = 0.24
MIC_HANDLE_W, MIC_HANDLE_H = 0.20, 0.72
HEAD0 = np.array([0.0, 2.88, 0.0])         # head touching the detector rim (top at 2.70); the hand stays inside the frame
HEAD1 = DET_CENTER.copy()                  # where it lands

# rays: NL / BE / DE in turn, in several bursts at random angles, each ray grows
# out of the centre and flies off; bursts overlap, so it stays chaotic
RAY_ORDER = ["NL", "BE", "DE"]
BURSTS = [(0.05, 11), (0.40, 9), (0.75, 10), (1.10, 8), (1.40, 6)]   # (start time after impact, rays)
FLY = 11.0                                 # how far rays and detector wedges travel (off frame)
N_WEDGES = 12                              # the detector breaks into this many wedges
N_SPARKS = 60
SEED = 11


def microphone(head_center) -> VGroup:
    c = np.asarray(head_center, dtype=float)
    head = Circle(radius=MIC_HEAD_R, arc_center=c, fill_color=col(tint(SLATE, 0.30)),
                  fill_opacity=1.0, stroke_color=col(INK), stroke_width=2.0)
    mesh = VGroup()
    for k in (-0.5, 0.0, 0.5):
        d = k * MIC_HEAD_R
        half = np.sqrt(MIC_HEAD_R ** 2 - d ** 2) * 0.96
        mesh.add(Line(c + [-half, d, 0], c + [half, d, 0]), Line(c + [d, -half, 0], c + [d, half, 0]))
    mesh.set_stroke(col(tint(SLATE, 0.70)), width=1.2)
    collar = RoundedRectangle(width=MIC_HANDLE_W * 1.35, height=0.10, corner_radius=0.03,
                              fill_color=col(GREY), fill_opacity=1.0,
                              stroke_color=col(INK), stroke_width=1.5)
    collar.move_to(c + [0.0, MIC_HEAD_R + 0.03, 0.0])
    handle = RoundedRectangle(width=MIC_HANDLE_W, height=MIC_HANDLE_H, corner_radius=0.06,
                              fill_color=col(SLATE), fill_opacity=1.0,
                              stroke_color=col(INK), stroke_width=1.5)
    handle.move_to(c + [0.0, MIC_HEAD_R + 0.08 + MIC_HANDLE_H / 2, 0.0])
    return VGroup(handle, collar, head, mesh)


def _pad(w, h, center, angle=0.0, r=0.06) -> RoundedRectangle:
    p = RoundedRectangle(width=w, height=h, corner_radius=r, fill_color=col(HAND_FILL),
                         fill_opacity=1.0, stroke_color=col(HAND_STROKE), stroke_width=2.0)
    return p.rotate(angle).move_to(np.asarray(center, dtype=float))


def hand(handle_center):
    """A right hand gripping the handle from the right, forearm leaving the frame
    at the top right. Returns (behind, front, fingers, thumb, knuckles): ``behind``
    is drawn under the microphone (forearm, palm, thumb), ``front`` over it."""
    hc = np.asarray(handle_center, dtype=float)
    palm_c = hc + [0.46, 0.02, 0.0]
    forearm = _pad(0.62, 1.6, palm_c + [0.55, 0.95, 0.0], angle=-32 * DEGREES, r=0.12)
    palm = _pad(0.58, 0.72, palm_c, r=0.14)
    knuckle_x = palm_c[0] - 0.22
    fingers = VGroup()
    for i, dy in enumerate((0.27, 0.11, -0.05, -0.21)):
        f = _pad(0.62, 0.135, [knuckle_x - 0.31, palm_c[1] + dy, 0.0], r=0.06)
        fingers.add(f)
    thumb = _pad(0.50, 0.14, [knuckle_x - 0.15, palm_c[1] - 0.40, 0.0], angle=-18 * DEGREES, r=0.06)
    knuckles = [np.array([knuckle_x, palm_c[1] + dy, 0.0]) for dy in (0.27, 0.11, -0.05, -0.21)]
    thumb_root = np.array([knuckle_x + 0.08, palm_c[1] - 0.33, 0.0])
    behind = VGroup(forearm, palm, thumb)
    return behind, fingers, thumb, knuckles, thumb_root


def flag_ray(country: str, theta: float, center, length: float, half_w: float) -> VGroup:
    """A tapered ray from ``center`` in direction ``theta`` with three stripes
    (the flag's colours, side by side) and a thin ink outline so white reads."""
    c = np.asarray(center, dtype=float)
    ys = np.linspace(-half_w, half_w, 4)
    g = VGroup()
    for k, colour in enumerate(FLAG[country]):
        g.add(Polygon([0, 0, 0], [length, ys[k], 0], [length, ys[k + 1], 0],
                      fill_color=col(colour), fill_opacity=1.0, stroke_width=0))
    g.add(Polygon([0, 0, 0], [length, -half_w, 0], [length, half_w, 0],
                  fill_opacity=0, stroke_color=col(INK), stroke_width=1.4))
    return g.rotate(theta, about_point=[0, 0, 0]).shift(c)


def _dir(a: float) -> np.ndarray:
    return np.array([np.cos(a), np.sin(a), 0.0])


def star_burst(center, r_out: float, r_in: float, n: int, rng, **style) -> Polygon:
    """A jagged explosion star (alternating outer / inner radii, jittered)."""
    c = np.asarray(center, dtype=float)
    pts = []
    for k in range(2 * n):
        r = (r_out if k % 2 == 0 else r_in) * rng.uniform(0.75, 1.25)
        pts.append(c + r * _dir(k * TAU / (2 * n) + rng.uniform(-0.08, 0.08)))
    return Polygon(*pts, **style)


def delayed(delay: float, anim):
    return Succession(Wait(delay), anim) if delay > 0 else anim


def shattered(det: CMSSlice) -> list[VGroup]:
    """The same drawing as ``det``, cut into ``N_WEDGES`` angular wedges: every
    full ring becomes sectors / arcs, every segmented cell goes to the wedge
    its centre lies in. Returns the wedges (each drawn inside out)."""
    c = det.c
    dphi = TAU / N_WEDGES
    phi0 = TAU / 24                      # wedge edges on the dodecagon's edges (15 + 30 k degrees)
    # the ECAL and HCAL base rings are fully covered by their cells / towers: leave
    # them out so a cell straddling a wedge edge is never half-hidden by a base sector
    covered = {det.base[i] for i, (name, *_) in enumerate(LOGO_RINGS) if name in ("ecal", "hcal")}
    wedges = [VGroup() for _ in range(N_WEDGES)]
    for leaf in det.family_members_with_points():
        if leaf in covered:
            continue
        if isinstance(leaf, Annulus):
            for k in range(N_WEDGES):
                wedges[k].add(AnnularSector(inner_radius=leaf.inner_radius, outer_radius=leaf.outer_radius,
                                            angle=dphi, start_angle=phi0 + k * dphi, arc_center=c).match_style(leaf))
        elif isinstance(leaf, Circle) and np.linalg.norm(leaf.get_center() - c) < 1e-3:
            for k in range(N_WEDGES):
                piece = (Sector(radius=leaf.radius, angle=dphi, start_angle=phi0 + k * dphi, arc_center=c)
                         if leaf.get_fill_opacity() > 0 else
                         Arc(radius=leaf.radius, start_angle=phi0 + k * dphi, angle=dphi, arc_center=c))
                wedges[k].add(piece.match_style(leaf))
        else:
            d = leaf.get_center() - c
            k = int(np.floor(((np.arctan2(d[1], d[0]) - phi0) % TAU) / dphi)) % N_WEDGES
            wedges[k].add(leaf.copy())
    return wedges


class MicDrop(Scene):

    def construct(self):
        white_background(self)
        det = CMSSlice()
        self.add(det)
        self.wait(0.3)

        # the hand with the microphone, from the top
        mic = microphone(HEAD0)
        handle_c = HEAD0 + [0.0, MIC_HEAD_R + 0.08 + MIC_HANDLE_H * 0.40, 0.0]
        behind, fingers, thumb, knuckles, thumb_root = hand(handle_c)
        hand_all = VGroup(behind, fingers)
        self.play(FadeIn(hand_all, shift=0.6 * UP), FadeIn(mic, shift=0.6 * UP), run_time=0.9, rate_func=EASE)
        self.add(behind, mic, fingers)      # thumb behind the handle, fingers in front
        self.wait(0.5)

        # the hand opens
        self.play(*[f.animate.rotate(70 * DEGREES, about_point=k) for f, k in zip(fingers, knuckles)],
                  thumb.animate.rotate(-45 * DEGREES, about_point=thumb_root),
                  run_time=0.55, rate_func=EASE)

        # the drop
        self.play(mic.animate.shift(HEAD1 - HEAD0), run_time=0.85, rate_func=rate_functions.ease_in_quad)

        # -- the impact: the detector cracks into wedges, flash, shock waves, sparks, camera shake
        rng = np.random.default_rng(SEED)
        wedges = shattered(det)
        self.remove(det); self.add(*wedges)
        wedge_dirs = [_dir(TAU / 24 + (k + 0.5) * TAU / N_WEDGES) for k in range(N_WEDGES)]

        star2 = star_burst(HEAD1, 1.05, 0.45, 14, rng, fill_color=col(INK), fill_opacity=1.0, stroke_width=0)
        star1 = star_burst(HEAD1, 1.45, 0.55, 12, rng, fill_color=col(FLAG["BE"][1]), fill_opacity=1.0,
                           stroke_color=col(HIGHLIGHT), stroke_width=3.0)
        stars = VGroup(star2, star1).scale(1 / 15, about_point=HEAD1)

        def shock_ring(delay, r_max, sw, t):
            ring = Circle(radius=0.06, arc_center=HEAD1, fill_opacity=0, stroke_color=col(INK), stroke_width=sw)
            return ring, delayed(delay, ApplyFunction(
                lambda m, f=r_max / 0.06: m.scale(f, about_point=HEAD1).set_stroke(width=0.5, opacity=0.0),
                ring, run_time=t, rate_func=rate_functions.ease_out_quad))

        spark_cols = [INK, HIGHLIGHT] + FLAG["NL"][::2] + FLAG["BE"][1:] + FLAG["DE"][1:]
        sparks, spark_anims = VGroup(), []
        for _ in range(N_SPARKS):
            a = rng.uniform(0, TAU)
            ln = rng.uniform(0.12, 0.35)
            sp = Line(HEAD1 + 0.15 * _dir(a), HEAD1 + (0.15 + ln) * _dir(a),
                      stroke_color=col(spark_cols[rng.integers(len(spark_cols))]), stroke_width=rng.uniform(2, 5))
            sparks.add(sp)
            v = rng.uniform(2.5, 6.5) * _dir(a)
            spark_anims.append(delayed(rng.uniform(0, 0.25), ApplyFunction(
                lambda m, v=v: m.shift(v).set_stroke(opacity=0.0), sp,
                run_time=rng.uniform(0.45, 0.9), rate_func=rate_functions.ease_out_quad)))

        shake = ValueTracker(0.0)
        shaker = Dot(radius=0).set_opacity(0)

        def _shake(m, dt):
            amp = shake.get_value()
            self.camera.frame_center = np.array([rng.normal(0, amp), rng.normal(0, amp), 0.0])

        shaker.add_updater(_shake)

        rings = [shock_ring(0.0, 4.0, 10.0, 0.7), shock_ring(0.10, 3.2, 6.0, 0.6),
                 shock_ring(0.40, 4.5, 7.0, 0.7), shock_ring(0.80, 4.5, 5.0, 0.7)]
        self.add(*[r for r, _ in rings], sparks, stars, shaker)
        self.remove(mic)
        shake.set_value(0.16)
        crack = [ApplyFunction(lambda m, v=0.14 * d: m.shift(v), w, run_time=0.2,
                               rate_func=rate_functions.ease_out_quad) for w, d in zip(wedges, wedge_dirs)]
        impact = AnimationGroup(
            Succession(ApplyFunction(lambda m: m.scale(15, about_point=HEAD1), stars, run_time=0.16,
                                     rate_func=rate_functions.ease_out_cubic),
                       FadeOut(stars, run_time=0.22)),
            *[a for _, a in rings], *spark_anims, *crack,
            FadeOut(hand_all, shift=1.5 * UP, run_time=0.6),
            shake.animate(run_time=0.7, rate_func=rate_functions.ease_out_quad).set_value(0.0),
        )

        # -- the rays: bursts at random angles, each grows then flies off; the wedges go with them
        ray_anims = []
        i = 0
        for t0, n in BURSTS:
            for _ in range(n):
                a = rng.uniform(0, TAU)
                r = flag_ray(RAY_ORDER[i % 3], a, HEAD1, rng.uniform(3.4, 5.8), rng.uniform(0.14, 0.40))
                i += 1
                full = r.copy()
                r.scale(1e-3, about_point=HEAD1)          # invisible until its burst
                self.add(r)
                ray_anims.append(Succession(
                    Wait(t0 + rng.uniform(0, 0.18)),
                    Transform(r, full, run_time=rng.uniform(0.22, 0.4), rate_func=rate_functions.ease_out_cubic),
                    Wait(rng.uniform(0, 0.15)),
                    ApplyFunction(lambda m, v=FLY * _dir(a): m.shift(v), r,
                                  run_time=rng.uniform(0.45, 0.8), rate_func=rate_functions.ease_in_quad)))
        wedge_anims = [delayed(rng.uniform(0.35, 1.05), ApplyFunction(
            lambda m, v=FLY * d, ang=rng.uniform(-150, 150) * DEGREES: m.shift(v).rotate(ang), w,
            run_time=rng.uniform(0.55, 0.9), rate_func=rate_functions.ease_in_quad))
            for w, d in zip(wedges, wedge_dirs)]
        self.play(AnimationGroup(impact, *ray_anims, *wedge_anims))
        shaker.remove_updater(_shake)
        self.camera.frame_center = ORIGIN.copy()
        self.remove(shaker, sparks, stars, *[r for r, _ in rings], *wedges,
                    *[m for a in ray_anims for m in a.mobject])

        # the clean slide
        thanks = text("Thank you!", weight=BOLD).scale_to_fit_height(0.95)
        quest = text("Any questions?").scale_to_fit_height(0.48)
        words = VGroup(thanks, quest).arrange(DOWN, buff=0.55).move_to(ORIGIN)
        self.play(FadeIn(thanks, shift=0.3 * UP), run_time=0.7, rate_func=EASE)
        self.play(FadeIn(quest, shift=0.3 * UP), run_time=0.6, rate_func=EASE)
        self.wait(0.1)
