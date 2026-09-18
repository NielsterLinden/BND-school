"""Section 7: the closing clip.

    MicDropResults   opens on comb_published's last frame (the six-row result plot). The CMS
                     point of the published row grows into the CMS detector, our own combined
                     point hands over the microphone (a hand holding it comes in, the fingers
                     open), and the drop takes the whole slide apart: rays in the colours of the
                     three BND-school flags (NL, BE, DE), the detector breaks into wedges.

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
    Circle, Dot, FadeIn, FadeOut, Line, Polygon, Rectangle, RoundedRectangle, Scene, Sector, Succession,
    Transform, VGroup, VMobject, ValueTracker, Wait, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s6_combination_story as s6  # noqa: E402  (MicDropResults opens on the last frame of comb_published)

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
FLY = 11.0                                 # how far rays and detector wedges travel (off frame)
N_WEDGES = 12                              # the detector breaks into this many wedges
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


def shattered(det: CMSSlice, n_wedges: int = N_WEDGES) -> list[VGroup]:
    """The same drawing as ``det``, cut into ``n_wedges`` angular wedges: every
    full ring becomes sectors / arcs, every segmented cell goes to the wedge
    its centre lies in. Returns the wedges (each drawn inside out)."""
    c = det.c
    N_WEDGES = n_wedges
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


# ---------------------------------------------------------------------------
# the closing clip of the talk: the results plot itself turns into the mic drop
# ---------------------------------------------------------------------------

N_WEDGES_X = 18        # many pieces, sparks and rays: "violent, chaotic, extreme"
N_SPARKS_X = 170
N_DEBRIS_X = 48
N_CRACKS_X = 16
BURSTS_X = [(0.00, 14), (0.20, 12), (0.42, 14), (0.65, 12), (0.88, 12), (1.12, 10)]


def _crack(center, rng) -> VMobject:
    """A jagged crack running outward from ``center``."""
    a = rng.uniform(0, TAU)
    pts, r, r_max = [np.asarray(center, dtype=float)], 0.2, rng.uniform(3.2, 7.0)
    while r < r_max:
        r += rng.uniform(0.35, 0.85)
        pts.append(center + r * _dir(a + rng.uniform(-0.28, 0.28)))
    m = VMobject(stroke_color=col(INK), stroke_width=rng.uniform(2.5, 6.0))
    m.set_points_as_corners(pts)
    return m


class MicDropResults(Scene):
    """The closing clip: opens on comb_published's last frame (the six-row result plot).

    The CMS point of the published row grows into the CMS detector, our own combined point
    hands over the microphone, and the drop takes the whole slide apart."""

    def construct(self):
        white_background(self)
        st = s6.summary_state(6)
        s6.restack(self, st, s6.ORDER_SUM)
        self.wait(0.45)

        cms_marker, comb_marker = st["row_cms"][4], st["row_comb"][4]
        keep = (cms_marker, comb_marker)
        fade = VGroup(*[m for k in s6.ORDER_SUM if k in st
                        for m in ([st[k]] if k not in ("row_cms", "row_comb") else [x for x in st[k] if x not in keep])])

        # 1. the CMS result opens out into the CMS detector
        det = CMSSlice()
        det0 = mini_slice(0.028, cms_marker.get_center())
        self.add(det0)
        self.bring_to_front(*keep)
        self.play(FadeOut(fade), FadeOut(cms_marker, scale=0.2), Transform(det0, det), run_time=1.5, rate_func=EASE)
        self.remove(det0)
        self.add(det, comb_marker)
        self.wait(0.2)

        # 2. our own combined point hands over the microphone
        mic = microphone(HEAD0)
        handle_c = HEAD0 + [0.0, MIC_HEAD_R + 0.08 + MIC_HANDLE_H * 0.40, 0.0]
        behind, fingers, thumb, knuckles, thumb_root = hand(handle_c)
        hand_all = VGroup(behind, fingers)
        mic0 = mic.copy().scale(0.05).move_to(comb_marker.get_center())
        self.add(mic0)
        self.play(FadeOut(comb_marker, scale=0.3), Transform(mic0, mic), run_time=1.0,
                  rate_func=rate_functions.ease_out_cubic)
        self.remove(mic0)
        self.add(mic)
        self.play(FadeIn(hand_all, shift=0.5 * UP), run_time=0.45, rate_func=EASE)
        self.add(behind, mic, fingers)          # thumb behind the handle, fingers in front
        self.wait(0.3)

        # 3. the hand opens, the microphone falls harder
        self.play(*[f.animate.rotate(78 * DEGREES, about_point=k) for f, k in zip(fingers, knuckles)],
                  thumb.animate.rotate(-52 * DEGREES, about_point=thumb_root), run_time=0.42, rate_func=EASE)
        self.play(mic.animate.shift(HEAD1 - HEAD0), run_time=0.5, rate_func=rate_functions.ease_in_quad)

        # 4. impact: flash, stars, shock rings, cracks, sparks, debris, the detector in pieces
        rng = np.random.default_rng(SEED + 3)
        wedges = shattered(det, N_WEDGES_X)
        self.remove(det)
        self.add(*wedges)
        wedge_dirs = [_dir(TAU / 24 + (k + 0.5) * TAU / N_WEDGES_X) for k in range(N_WEDGES_X)]

        stars = VGroup(star_burst(HEAD1, 1.30, 0.50, 16, rng, fill_color=col(INK), fill_opacity=1.0, stroke_width=0),
                       star_burst(HEAD1, 1.80, 0.62, 13, rng, fill_color=col(FLAG["BE"][1]), fill_opacity=1.0,
                                  stroke_color=col(HIGHLIGHT), stroke_width=4.0),
                       star_burst(HEAD1, 2.30, 0.70, 11, rng, fill_color=col(FLAG["NL"][0]), fill_opacity=0.55,
                                  stroke_width=0)).scale(1 / 18, about_point=HEAD1)
        flash = Rectangle(width=16, height=10, stroke_width=0, fill_color=col(HIGHLIGHT), fill_opacity=0.5).move_to(ORIGIN)

        def shock_ring(delay, r_max, sw, t):
            ring = Circle(radius=0.06, arc_center=HEAD1, fill_opacity=0, stroke_color=col(INK), stroke_width=sw)
            return ring, delayed(delay, ApplyFunction(
                lambda m, f=r_max / 0.06: m.scale(f, about_point=HEAD1).set_stroke(width=0.5, opacity=0.0),
                ring, run_time=t, rate_func=rate_functions.ease_out_quad))

        rings = [shock_ring(0.0, 5.5, 14.0, 0.7), shock_ring(0.08, 4.2, 9.0, 0.6), shock_ring(0.26, 6.0, 8.0, 0.7),
                 shock_ring(0.50, 5.0, 7.0, 0.7), shock_ring(0.80, 6.5, 6.0, 0.8), shock_ring(1.15, 5.5, 4.0, 0.8)]

        cracks, crack_anims = VGroup(), []
        for _ in range(N_CRACKS_X):
            cr = _crack(HEAD1, rng)
            full = cr.copy()
            cr.scale(1e-3, about_point=HEAD1)
            cracks.add(cr)
            crack_anims.append(Succession(
                Transform(cr, full, run_time=0.16, rate_func=rate_functions.ease_out_cubic),
                ApplyFunction(lambda m: m.set_stroke(opacity=0.0), cr, run_time=0.5)))

        spark_cols = [INK, HIGHLIGHT] + FLAG["NL"][::2] + FLAG["BE"][1:] + FLAG["DE"][1:]
        sparks, spark_anims = VGroup(), []
        for _ in range(N_SPARKS_X):
            a = rng.uniform(0, TAU)
            ln = rng.uniform(0.12, 0.45)
            sp = Line(HEAD1 + 0.15 * _dir(a), HEAD1 + (0.15 + ln) * _dir(a),
                      stroke_color=col(spark_cols[rng.integers(len(spark_cols))]), stroke_width=rng.uniform(2, 6))
            sparks.add(sp)
            v = rng.uniform(3.5, 9.5) * _dir(a)
            spark_anims.append(delayed(rng.uniform(0, 0.35), ApplyFunction(
                lambda m, v=v: m.shift(v).set_stroke(opacity=0.0), sp,
                run_time=rng.uniform(0.4, 0.9), rate_func=rate_functions.ease_out_quad)))

        debris, debris_anims = VGroup(), []
        for _ in range(N_DEBRIS_X):
            a, r0 = rng.uniform(0, TAU), rng.uniform(0.2, 1.3)
            c0 = HEAD1 + r0 * _dir(a)
            s = rng.uniform(0.05, 0.19)
            colour = spark_cols[rng.integers(len(spark_cols))]
            tri = Polygon(*[c0 + s * _dir(b + rng.uniform(-0.35, 0.35)) for b in (0.0, TAU / 3, 2 * TAU / 3)],
                          fill_color=col(colour), fill_opacity=1.0, stroke_color=col(INK), stroke_width=1.0)
            debris.add(tri)
            debris_anims.append(delayed(rng.uniform(0, 0.2), ApplyFunction(
                lambda m, v=rng.uniform(4.0, 9.0) * _dir(a), ang=rng.uniform(-8, 8): m.shift(v).rotate(ang),
                tri, run_time=rng.uniform(0.6, 1.0), rate_func=rate_functions.ease_out_quad)))

        # the microphone itself is destroyed: its four pieces fly out spinning
        mic_anims = [delayed(rng.uniform(0, 0.12), ApplyFunction(
            lambda m, v=rng.uniform(3.0, 7.0) * _dir(rng.uniform(0, TAU)), ang=rng.uniform(-10, 10): m.shift(v).rotate(ang),
            piece, run_time=rng.uniform(0.5, 0.9), rate_func=rate_functions.ease_out_quad)) for piece in mic]

        shake = ValueTracker(0.0)
        shaker = Dot(radius=0).set_opacity(0)

        def _shake(m, dt):
            amp = shake.get_value()
            self.camera.frame_center = np.array([rng.normal(0, amp), rng.normal(0, amp), 0.0])

        shaker.add_updater(_shake)

        self.add(*[r for r, _ in rings], cracks, sparks, debris, stars, flash, shaker)
        shake.set_value(0.34)
        crack_out = [ApplyFunction(lambda m, v=0.22 * d: m.shift(v), w, run_time=0.18,
                                   rate_func=rate_functions.ease_out_quad) for w, d in zip(wedges, wedge_dirs)]
        impact = AnimationGroup(
            Succession(ApplyFunction(lambda m: m.scale(18, about_point=HEAD1), stars, run_time=0.14,
                                     rate_func=rate_functions.ease_out_cubic),
                       FadeOut(stars, run_time=0.25)),
            ApplyFunction(lambda m: m.set_fill(opacity=0.0), flash, run_time=0.32,
                          rate_func=rate_functions.ease_out_quad),
            *[a for _, a in rings], *crack_anims, *spark_anims, *debris_anims, *mic_anims, *crack_out,
            FadeOut(hand_all, shift=2.2 * UP, run_time=0.5),
            Succession(Wait(0.18), shake.animate(run_time=1.0, rate_func=rate_functions.ease_out_quad).set_value(0.0)),
        )

        # 5. the rays: more bursts, at random angles, each grows and flies off with the wedges
        ray_anims, i = [], 0
        for t0, n in BURSTS_X:
            for _ in range(n):
                a = rng.uniform(0, TAU)
                r = flag_ray(RAY_ORDER[i % 3], a, HEAD1, rng.uniform(3.8, 6.8), rng.uniform(0.14, 0.48))
                i += 1
                full = r.copy()
                r.scale(1e-3, about_point=HEAD1)
                self.add(r)
                ray_anims.append(Succession(
                    Wait(t0 + rng.uniform(0, 0.16)),
                    Transform(r, full, run_time=rng.uniform(0.18, 0.34), rate_func=rate_functions.ease_out_cubic),
                    Wait(rng.uniform(0, 0.12)),
                    ApplyFunction(lambda m, v=FLY * _dir(a): m.shift(v), r,
                                  run_time=rng.uniform(0.4, 0.75), rate_func=rate_functions.ease_in_quad)))
        wedge_anims = [delayed(rng.uniform(0.25, 0.9), ApplyFunction(
            lambda m, v=FLY * d, ang=rng.uniform(-240, 240) * DEGREES: m.shift(v).rotate(ang), w,
            run_time=rng.uniform(0.5, 0.85), rate_func=rate_functions.ease_in_quad))
            for w, d in zip(wedges, wedge_dirs)]
        self.play(AnimationGroup(impact, *ray_anims, *wedge_anims))
        shaker.remove_updater(_shake)
        self.camera.frame_center = ORIGIN.copy()
        self.clear()

        # the clean slide
        thanks = text("Thank you!", weight=BOLD).scale_to_fit_height(0.95)
        quest = text("Any questions?").scale_to_fit_height(0.48)
        VGroup(thanks, quest).arrange(DOWN, buff=0.55).move_to(ORIGIN)
        self.play(FadeIn(thanks, scale=1.7), run_time=0.55, rate_func=rate_functions.ease_out_back)
        self.play(FadeIn(quest, shift=0.3 * UP), run_time=0.55, rate_func=EASE)
        self.wait(0.1)
