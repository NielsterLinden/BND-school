# Manim recipes for particle-physics talks

A handover note for an agent starting a new Manim presentation project from scratch.
Everything here was learned the hard way on an LHCb momentum-calibration talk
(~30 animation clips, ManimCE 0.20.1, cairo renderer, rendered headless on lxplus).

You do **not** need a Feynman-diagram package, a detector-geometry package, or a 3D
engine. All four things below are built from plain `Line`, `Polygon`, `VMobject` and
`Dot`. That is the main lesson: physics pictures in Manim are geometry plus a couple
of small helper functions, kept in one shared module.

Read the sections you need. Section 0 and section 6 apply to every project.

---

## 0. Setup that actually works

**Renderer.** Use the default **cairo** renderer. Never pass `--renderer=opengl`:
it needs a display and will fail on any batch or remote node.

**Render command.** `manim render --quality l -o <name> scenes/<file>.py <SceneClass>`
for iteration, `--quality h` for delivery. Write `--quality h`, not `-qh`; the short
form is not accepted in 0.20.x.

**TeX.** `MathTex` needs a real TeX installation with `standalone.cls` and `dvisvgm`.
A conda `texlive-core` package is binaries only and will not work. On a CERN machine
put CVMFS TeX Live on PATH ahead of everything:

```bash
export PATH="<env>/bin:/cvmfs/sft.cern.ch/lcg/external/texlive/latest/bin/x86_64-linux:$PATH"
```

Elsewhere, install a full TeX Live. If `latex --version` works and `standalone.cls`
exists, `MathTex` will work.

**Inspect every frame you render.** Never assume a scene looks right.

```bash
ffmpeg -sseof -0.15 -i clip.mp4 -frames:v 1 final.png   # last frame
ffmpeg -ss 3.2   -i clip.mp4 -frames:v 1 t32.png        # a chosen time
```

`-sseof -0.05` decodes nothing at 15 fps (the `--quality l` frame rate). Use `-0.15`
or more.

**Render serially.** Two concurrent `manim render` calls that compile the same
`MathTex` string race on the same cached SVG file, and one of them dies with
`FileNotFoundError`.

**One style module, no exceptions.** Put every colour, font and plot helper in
`style/`. Scenes import from it and never hardcode a hex or a font name. Split the
palette into a **manim-free** `style/palette.py` so matplotlib scripts producing the
static figures can import the identical hex values. Scene files import from
`style/lhcb_style.py`, which re-exports the palette and adds the Manim helpers.

**White background** reads better on a projector than Manim's default dark:

```python
self.camera.background_color = ManimColor("#FFFFFF")
```

**Scene bootstrap** used by every scene file:

```python
from __future__ import annotations
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # project root on path
from manim import ...          # noqa: E402
from style.lhcb_style import ...  # noqa: E402
```

Note the `sys.path.insert` resolves through symlinks to the real tree. If you ever
overlay a scene file on another checkout, copy it, do not symlink it.

---

## 1. Feynman diagrams

There is no Feynman package and you do not want one. A tree-level diagram is four
things: straight fermion lines, a dashed or wavy boson line, labels at the line ends,
and a fixed dictionary of vertex coordinates.

### The coordinate dictionary

Put every vertex and every label anchor in one dict, in scene units. This is the
single thing you tune when the diagram looks wrong, and it makes a diagram
reproducible across chained clips.

```python
_FEY = {
    "b0": (-5.0, 3.35),   "wk": (-1.75, 3.35),   # W emission point, on the b line
    "u0": (-5.0, 1.1),    "u1": (1.5, 1.1),      # spectator quark
    "wv": (-1.0, 2.1),                           # W -> c sbar vertex
    "cbar_end": (1.5, 3.35),
    "c_end":    (1.5, 2.9),
    "sbar_end": (1.5, 1.55),
}
```

### The three line kinds

```python
from manim import Line, DashedLine, VMobject, ManimColor, TAU
import numpy as np, math

def fline(p0, p1, sw=2.6):                       # fermion: plain straight line
    return Line(np.array([*p0, 0.0]), np.array([*p1, 0.0]),
                stroke_color=ManimColor(INK), stroke_width=sw)

def wline(p0, p1):                               # massive boson: dashed
    return DashedLine(np.array([*p0, 0.0]), np.array([*p1, 0.0]),
                      dash_length=0.12, stroke_color=ManimColor(INK), stroke_width=2.4)

def wavy(p0, p1, amplitude=0.09, waves=4, sw=2.4, n=140):
    """Photon / gluon: a sine offset perpendicular to the line direction."""
    p0, p1 = np.array([*p0, 0.0]), np.array([*p1, 0.0])
    d = p1 - p0
    L = float(np.linalg.norm(d)); u = d / L
    nperp = np.array([-u[1], u[0], 0.0])
    pts = [p0 + u * (t * L) + nperp * amplitude * math.sin(TAU * waves * t)
           for t in np.linspace(0.0, 1.0, n)]
    vm = VMobject(stroke_color=ManimColor(INK), stroke_width=sw)
    vm.set_points_smoothly(pts)          # smooth, not set_points_as_corners
    return vm
```

Tuning notes that saved rounds of re-rendering:

- `amplitude` around 0.09 and 4 to 6 waves reads well over a 3 to 5 unit span.
  Scale `waves` with the length so the wavelength stays constant, otherwise a long
  photon looks like a different particle from a short one.
- `n=140` sample points is enough for a short photon. A long one needs more or the
  smoothing rounds off the crests.
- A gluon is the same helper with a larger amplitude and tighter waves, or a coil
  drawn as a parametric helix projected to 2D. A sine is usually enough for a talk.

### Reading order beats decoration

Physics people read a quark-line diagram left to right. Three rules that made ours
legible:

1. **Keep the spectator line perfectly straight and horizontal.** In `B+ -> J/psi K+`
   the `b` line and the `cbar` continuation are drawn as two segments at the *same*
   y, so they read as one line with the `W` emitted at its midpoint.
2. **Never cross lines that do not have to cross.** The `c` from the W vertex is
   placed just below the `cbar` (0.45 units apart), not above it, so no line crossing
   appears. Crossing lines read as a different topology.
3. **Name the bound state by placing its label between the two quark labels**, offset
   to the right. No brace mobject needed. `J/psi` sits between the `cbar` and `c` end
   labels; `K+` sits between `sbar` and `u`. The eye pairs them for free.

### Build order is z-order

Animate the pieces in exactly the order you would `add` them in a static builder.
Manim's z-order follows add order, and if a later clip rebuilds the final frame from
a static builder, any mismatch shows up as a visible seam. Practical sequence that
works for a tree decay:

1. fade in the initial-state labels (`B+`, `bbar`, `u`),
2. `Create` the two straight fermion lines,
3. `Create` the boson line and the continuation line together,
4. fade in the boson label and the first daughter label,
5. `Create` the second-generation lines,
6. fade in the remaining quark labels,
7. fade in the bound-state labels last, with a slight `scale=0.8` so they arrive.

### Annihilation and decay in a later clip

To show `ccbar -> gamma -> mu+ mu-` without redrawing the first diagram: keep every
line and label from the previous clip untouched, fade out only the `J/psi` label, then
`Create` two short new stubs converging into an annihilation vertex, `Flash` the
vertex, `Create` the wavy photon, then the two lepton legs. Do not re-animate what is
already on screen.

`Flash` on a vertex point is fine. Do not `Flash` or pulse text or formulas; it reads
as an error rather than as emphasis.

---

## 2. A 3D detector

The detector is a 2D side-view drawing that gets **extruded**. That single decision
makes everything else possible: the same geometry constants serve the flat schematic,
the 3D solids, and the track model.

### Step 1: the 2D schematic, as one class with named parts

Define every subsystem boundary as a module-level constant in scene units, with a
stated scale (for example `0.64 scene units per metre`, with the vertex detector
deliberately oversized for readability). Then build one `VGroup` subclass whose parts
are attributes:

```python
class DetectorSchematic(VGroup):
    def __init__(self, show_beamline=True, **kw):
        super().__init__(**kw)
        self.beamline    = DashedLine(...)
        self.velo        = VGroup(self.velo_vessel, self.velo_plates)
        self.magnet      = VGroup(self.magnet_wash, self.magnet_coils)
        self.t_stations  = VGroup(...)
        ...
        self.add(self.rich1, self.rich2, self.velo, self.ut, self.magnet,
                 self.t_stations, self.ecal, self.hcal, self.muon)
```

Named parts let a scene build the detector up piece by piece (`FadeIn(det.velo)`,
`DrawBorderThenFill(det.magnet_coils)`) without redrawing anything.

Details that make it look like a real detector schematic:

- Give each subsystem a **very light tint** of its canonical colour (roughly 75 %
  toward white). The detector is a backdrop. Tracks and highlights must lead.
- Draw the RICH-style vessels **first** so everything crossing them paints on top.
- Make the vertex detector a **hollow outline** with thin sensor plates flanking an
  open gap, so the collision point can be drawn inside it.
- Draw the dipole magnet as **two coil polygons whose aperture widens downstream**,
  plus a pale filled polygon between them marking the field region. The wash is what
  makes "the track bends here" obvious without a single word.

### Step 2: extrusion

One function turns any x-y polygon into a solid:

```python
def prism_z(pts, half_d, *, fill, fill_opacity=1.0, stroke_color=INK,
            stroke_width=1.8, face_opacity=None, stroke_opacity=1.0):
    """An x-y polygon extruded symmetrically in z. Faces: back (-z), sides, front (+z)."""
    P = [np.asarray(p, dtype=float) for p in pts]
    off = np.array([0.0, 0.0, float(half_d)])
    fo = ((fill_opacity,) * 3) if face_opacity is None else face_opacity
    kw = dict(stroke_color=ManimColor(stroke_color), stroke_width=stroke_width,
              stroke_opacity=stroke_opacity, fill_color=fill)
    back  = [p - off for p in P]
    front = [p + off for p in P]
    n = len(P)
    faces  = [Polygon(*back, fill_opacity=fo[0], **kw)]
    faces += [Polygon(back[i], back[(i+1) % n], front[(i+1) % n], front[i],
                      fill_opacity=fo[1], **kw) for i in range(n)]
    faces.append(Polygon(*front, fill_opacity=fo[2], **kw))
    return VGroup(*faces)
```

Key properties, all of them used:

- **`half_d -> 0` collapses every face onto the flat polygon.** So a `depth` parameter
  scaling all half-depths gives you `depth=0.001` (visually identical to the 2D
  drawing) and `depth=1.0` (full solid). A `Transform` between the two is a clean
  linear inflate, and it is how a flat schematic becomes 3D on screen.
- **Paint order inside one convex single-colour prism never matters.** Across
  elements it does: order them **far to near for the final view**. If the camera will
  end up looking upstream, build the downstream systems first.
- **Face opacities stack.** Four plate prisms in a row put eight translucent faces
  over the centre of the image; at per-face opacity `o` the centre reads
  `1 - (1 - o)^8`. Per-face 0.12 looked like a solid 0.64 wall. Ghosting values that
  work are around **0.03 to 0.10 per face**, not 0.3.
- Give each named part a key so later scenes can address it:
  `grp.named_parts = dict(zip(names, parts))`.

### Step 3: the turn

Use a `ThreeDScene`, but **rotate the object, not the camera**.

```python
self.play(Rotate(g, angle=PI/2, axis=UP, about_point=ORIGIN),
          run_time=3.2, rate_func=rate_functions.ease_in_out_sine)
```

Why this and not `set_camera_orientation` / `move_camera`:

- The cairo `ThreeDCamera` at its default orientation is the **identity** rotation,
  and its perspective factor `fd / (fd - z)` is exactly 1 at `z = 0`. So a
  `ThreeDScene` whose opening state is flat at `z = 0` renders **pixel-identical** to
  a plain 2D `Scene`. You can chain a 3D clip onto a 2D clip with no visible seam.
- Sweeping `phi`/`theta` away from that start pose hits spherical-camera
  degeneracies and the motion looks wrong. Rotating the object is equivalent to
  walking the camera around, with none of that.

Choose the turn direction by what you want in front. Turning so the downstream end
comes toward the viewer puts the last tracking station face-on. Turning the other way
puts the interaction point in front and you look *down the beam through* the detector.
The second one is far more striking, and it needs the next trick.

### Step 4: see-through

For the "look down the beam" view, make the things you look through translucent as
the depth inflates. Typical end-state per-face opacities: vessels 0.10, tracker plane
0.12, field volume 0.22. One catch: a near-white fill has no colour left to show at
0.22 opacity, so the field volume "reads as nothing". Give the see-through version a
**real colour cast** (for example lighten the brand blue by 0.6 instead of using the
near-white fill). A `Transform` morphs colour and depth together, so the swap is free.

Grade sizes to beat perspective. Stations further from the viewer shrink; scale each
successive station by about 8 % so their rims peek past the nearest one and the view
reads as nested layers rather than one square.

---

## 3. Particles flying through the detector

### The track model

One function, three regimes, slope continuous at both magnet boundaries:

```python
def track_y(x, slope, kappa):
    """Straight -> constant-curvature arc -> straight. kappa > 0 bends up."""
    if x <= MAG_X0:
        return slope * (x - IP_X)
    y_in = slope * (MAG_X0 - IP_X)
    if x <= MAG_X1:
        dx = x - MAG_X0
        return y_in + slope * dx + 0.5 * kappa * dx * dx
    L = MAG_X1 - MAG_X0
    y_out = y_in + slope * L + 0.5 * kappa * L * L
    return y_out + (slope + kappa * L) * (x - MAG_X1)
```

Physics encoded, all of it visible to an audience without narration:

- **sign of `kappa` = sign of the charge.** A positive and a negative track from the
  same vertex splay apart in the magnet. Define one canonical `(slope, kappa)` pair
  for `+` and `-` and reuse it across every clip so the two tracks are the same two
  tracks each time.
- **`|kappa|` is inversely the momentum.** A stiff (high momentum) track barely bends.
  Make the momentum differences visible: a factor of two in `kappa` is obvious, 10 %
  is not.
- The quadratic only applies inside the field region. Outside it the track is
  straight, which is what makes the magnet read as the thing that did it.

Turn it into a mobject as a polyline, not a `ParametricFunction`:

```python
def track_path(slope, kappa, color, stroke_width=3.5, x_start=IP_X, x_end=X_END, n=240):
    pts = [np.array([float(x), track_y(float(x), slope, kappa), 0.0])
           for x in np.linspace(x_start, x_end, n)]
    vm = VMobject(stroke_color=color, stroke_width=stroke_width)
    vm.set_points_as_corners(pts)
    return vm
```

`set_points_as_corners` with `n=240` is smooth on screen and, unlike smoothing,
never overshoots at the magnet boundaries. `Create` on it draws left to right, which
is exactly the flight direction. `run_time` around 2 s per track, with
`rate_func=ease_in_out_sine`, reads as a particle travelling rather than a line
appearing.

### Where the tracks stop matters

Stop each track where its species actually stops. Muons go through everything and end
inside the muon stations. Hadrons stop in the hadron calorimeter. Photons are neutral,
so they do not bend and leave no hits: draw them straight, tangent to the track at the
emission point, running on through the field with no kinks. An audience of physicists
notices the difference immediately, and it does the explaining for you.

Past the last tracking station the particle is no longer measured. Split the track
there: solid coloured line up to the last tracker, then a **light grey dashed
extrapolation**. Split by arc length so the drawing hand-over is seamless:

```python
seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
cum = np.concatenate([[0.0], np.cumsum(seg)])
frac = float(np.interp(split_x, xs, cum) / cum[-1])
draw = Succession(Create(blue, run_time=frac * T,       rate_func=linear),
                  Create(grey, run_time=(1 - frac) * T, rate_func=linear))
```

### 3D tracks

Extend the same model with a third coordinate. The bend stays in `y` (the magnet
bends in one plane), the other transverse coordinate is linear in `x`:

```python
def fly_world(x_nat, y_nat, w_nat):
    """Natural side-view coords (beam x, bend y, horizontal w) -> world after the turn."""
    return np.array([w_nat, y_nat, -x_nat + DOLLY_DZ])
```

Pick the horizontal slopes so the tracks **fan out and never cross on screen**. Two
3D lines that meet in projection become unreadable. This is worth a couple of
iterations: render, look, adjust the slopes.

### Flying in: a dolly, not a camera move

Translate the whole (already rotated) object toward the camera along `+z`. The
perspective factor `fd / (fd - z)` gives true parallax for free.

```python
self.play(*[m.animate.shift(DOLLY_DZ * OUT) for m in g.named_parts.values()],
          run_time=3.4, rate_func=rate_functions.ease_in_out_sine)
```

Three things to get right:

- **`fd = 20` is the focal distance. Keep every point well short of `z = 20`** or the
  projection explodes. A dolly of about 7 with content spanning a few units is safe.
- **Part the obstacles.** The magnet coils would fill the view, so they slide
  vertically out of the aperture during the same play (add an extra `UP` / `DOWN`
  shift on those parts only).
- **Ghost what you fly into.** Vertex detector, vessels and upstream tracker drop to
  per-face opacity 0.03 to 0.04 during the dolly. Remember the stacking rule above.

Keep the per-part deviations in one dict so the animation and any static rebuild of
the final frame read from the same source:

```python
_FLY_SPEC = {                       # name: (extra shift, fill opacity, stroke opacity)
    "coil_up": (COIL_PART * UP, None, None),
    "coil_dn": (COIL_PART * DOWN, None, None),
    "wash":    (None, 0.10, None),
    "velo_plates": (None, 0.03, 0.12),
}
```

### A decay, end to end

```python
# 1. the parent appears at the collision point
self.play(Flash(ip, color=AMBER, line_length=0.22, flash_radius=0.45, run_time=0.6),
          FadeIn(bdot), FadeIn(lab_b), run_time=0.6)
self.wait(0.35)

# 2. the daughters fly, each with a dot riding the drawing tip
flights = []
for trk in trks:
    pd = Dot(point=trk.get_start(), radius=0.055, color=BLUE)
    flights.append(AnimationGroup(Create(trk), MoveAlongPath(pd, trk.copy()),
                                  run_time=1.6, rate_func=rate_functions.linear))
self.play(FadeOut(bdot), FadeOut(lab_b), LaggedStart(*flights, lag_ratio=0.22))

# 3. the hits, as a SEPARATE play
self.play(*[Flash(trk.get_end(), color=AMBER, line_length=0.14, flash_radius=0.24)
            for trk in trks], run_time=0.35)
```

Two traps in that block, both cost a render round:

- **`MoveAlongPath` needs `trk.copy()`**, not `trk` itself, when `trk` is also being
  `Create`d in the same group.
- **Never put `Flash` inside a `Succession`.** Its spoke lines are added to the scene
  when the play starts and hang there fully drawn until their slot comes up. Flashes
  go in their own `play` after the flights.
- `rate_func=linear` on the flights. Easing makes particles look like they accelerate
  and decelerate, which they do not.

---

## 4. Invariant mass plots

### Build the axes yourself

Manim's `Axes` tick labels do not match a thesis look. Wrap `Axes` in a small class
that draws its own tick marks and numeric labels in the deck font, and exposes `.c2p`
for data-to-scene coordinates. Then every plot primitive takes the axes object plus
data coordinates:

```python
dax = DataAxes(x_range, y_range, x_length, y_length, x_ticks=[...], x_fmt="{:.0f}")
dax.move_to(PANEL_CENTER)          # position FIRST
bar = data_bar(dax, x_center, value, half_width, color=...)   # then build content
```

Positioning before building content matters: the helpers convert through `.c2p` at
construction time.

Give the class `data_bar`, `data_trace`, `data_band`, `data_dot`, `data_errorbar`.
That is the whole toolkit a physics talk needs.

### The line shape

For an explainer histogram you want a real mass line shape, not a Gaussian. Use a
**double-sided Crystal Ball** with a heavier left tail. The left tail is the physics:
final-state radiation carries energy away, so events pile up below the true mass.

```python
def dscb(t, a_l, n_l, a_r, n_r):
    if t < -a_l:
        return ((n_l/a_l)**n_l * math.exp(-0.5*a_l*a_l) * (n_l/a_l - a_l - t)**(-n_l))
    if t > a_r:
        return ((n_r/a_r)**n_r * math.exp(-0.5*a_r*a_r) * (n_r/a_r - a_r + t)**(-n_r))
    return math.exp(-0.5 * t * t)

SYM_TAIL = (1.6, 5.0)   # right side, and the radiation-free reference
RAD_TAIL = (0.9, 3.0)   # left side with the radiative tail
```

Generate bin heights deterministically with a **seeded** RNG, one draw per bin in bin
order. Then the radiative and non-radiative versions of the same histogram differ
only by the tail excess, bin by bin, and you can animate one into the other honestly:

```python
rng = np.random.default_rng(11)
h = HIST_H * dscb((x - PEAK_X)/SIG, *left_tail, *SYM_TAIL)
h *= 1.0 + 0.05 * float(rng.standard_normal())
```

Drop bins below a small threshold so the axis does not end in a fringe of hairlines.

### Filling the histogram is the story

The sequence that works, over about 8 seconds:

1. One event flies through the detector, three momenta peel off the tracks as symbols,
   they converge into `m^2 = (p1 + p2 + p3)^2`, and **one bar** grows.
2. Three more events, faster, with visibly different kinematics, each adding a bar.
3. A clock appears and spins up while event dots stream from all over the detector
   into the plot, and all the remaining bars grow in random order.

```python
self.play(LaggedStart(*bar_anims, lag_ratio=0.02),
          LaggedStart(*dot_anims, lag_ratio=0.06),
          turns.animate(rate_func=rate_functions.ease_in_quad).set_value(3.0),
          run_time=2.6)
```

Sample the destination bin of each streaming dot from the histogram itself
(`p=heights/heights.sum()`) so the rain visibly concentrates on the peak. Use
`GrowFromEdge(bar, DOWN)` for new bars and `Transform` for bars that already exist.

### Real data

For results plots, do not recompute anything at render time. Freeze the numbers into
thin artifacts (`data/*.json`, `*.csv`, `*.npz`) produced by the analysis code, and
load them in the scene. Rendering must never import a heavy analysis stack.

Draw a frozen histogram as a **step outline**, every bin edge twice, so the polyline
*is* the binning and nobody can accuse you of smoothing:

```python
xs, ys = [], []
for i, n in enumerate(counts):
    xs += [edges[i], edges[i+1]]
    ys += [n, n]
return data_trace(dax, xs, ys, color=color, stroke_width=3.0)
```

Show a "before" distribution as a filled area under its own step, and each corrected
version as a coloured step on top. Three overlaid steps is the maximum that stays
readable.

### Marking the numbers

- **Reference mass**: a dashed vertical line, plus (if you have a target precision) a
  pale band showing how wide that tolerance is on this axis. A band that turns out to
  be thinner than the dashed line itself is a very effective silent statement.
- **Width**: draw it where it is measured, as a horizontal arrow spanning
  `peak +/- 1.1774 * sigma` (the Gaussian FWHM). If you overlay several, **stack them
  at different heights**, otherwise widths that differ by 3 % sit on top of each other.
- **Peak versus mean.** Decide which one your story is about and never mix them. A
  radiative tail moves the **mean** a lot and the **peak** very little. If the method
  is mean-seeking, mark the mean with a solid marker that slides, and leave the peak
  and the bars alone. Nothing else on screen should move while that marker slides.

### Numbers should end up owning the frame

The pattern that landed best: each number first appears **small, on the thing that
measured it** (the peak value beside the reference line, the width value at the right
cap of its own arrow). Then the whole plot shrinks and parks at the far left while the
numbers fly right and grow to reading size. `Transform(seed_copy, final_number)` does
it in one play, together with `panel.animate.scale(...).move_to(...)`.

### Words

Keep narrative text out of the clips if the deck is assembled in PowerPoint. Let the
speaker's slide carry titles and captions. What the clip may draw itself:

- variables, units and particle symbols that label a depicted object,
- numeric tick labels and value labels, which are the data,
- one small colour key, if a colour code is load-bearing across several chained clips
  and the slide software cannot follow it.

---

## 5. Chaining clips

Long animations are easier to deliver and to re-cut as several short clips, each
opening exactly where the previous one ended. This needs discipline, and it pays off.

**Rule 1: every clip has a pure builder function** that constructs its final frame
with no animation.

```python
def clip_reco_final_state() -> dict:
    ...
    return {"det": det, "trks": trks, "dax": dax, "bars": bars}
```

The next clip starts with `self.add(*[st[k] for k in ORDER])` and animates from there.
The builder and the scene must produce the *same* mobjects in the *same* order. Keep a
module-level `ORDER` tuple and assert against it at the end of the scene:

```python
assert [name_of(m) for m in self.mobjects if m.has_points()] == list(ORDER)
```

A failed assert exits manim before it combines the partial files, which leaves the
previous mp4 in place. Check for `Rendered <Scene>` in the log.

**Rule 2: put any placement in one affine helper**, so both clips place things
identically:

```python
def place(mobj, place):          # place = (scale, center)
    s, center = place
    mobj.scale(s, about_point=NATURAL_CENTER)
    mobj.shift(center - NATURAL_CENTER)
    return mobj
```

**Rule 3: verify the joint with a pixel diff.** A small tool comparing the last frame
of clip N against the first frame of clip N+1 catches everything: a missing mobject, a
0.1 unit shift, a wrong z-order. Compute mean absolute difference, max, the fraction of
pixels over a threshold, and a "blob" count (pixels whose whole 5x5 neighbourhood is
over threshold). H.264 ringing on glyph edges and grid lines erodes to a blob count of
zero, while a genuinely misplaced object scores in the hundreds.

Calibration from real joints at 1080p: `mean|d|` 0.2 to 0.6, `max|d|` up to about 60
on hairlines, over-threshold fraction under 0.1 %, blob 0. Fail the check at
`mean|d| > 1.0`, fraction `> 0.5 %`, or blob `> 50`.

**Rule 4: end a clip on its last change.** No held final frame, no trailing `wait`
beyond about 0.2 s. The slide software rests on the last frame. Deliver a
`<clip>_final.png` alongside each mp4.

---

## 6. Trap list (ManimCE 0.20.x, cairo)

Every one of these cost at least one render round.

1. **The last frame is written just before the final interpolation.** `finish()` calls
   `interpolate(1)` *after* the last rendered frame, so a clip ending on a `play`
   rests at `alpha = (n-1)/n`: a fade-out is still about 7 % visible. End with
   `self.wait(0.1)` to render the true final state.

2. **`Indicate` and `.animate(rate_func=there_and_back)` leave the mobject at peak
   scale.** `finish()` snaps the mobject to its method target regardless of the rate
   function, so every later frame shows it scaled. For an emphasis pulse, animate a
   throwaway `m.copy()` ghost and `self.remove(ghost)` right after the play. The real
   mobject is never animated.

3. **`LaggedStart` / `AnimationGroup` over children of an already-added `VGroup`
   re-stacks the z-order.** With `group=None` it wraps the children in a new `Group`,
   `play` adds that on top, and the parent gets split out of `self.mobjects`. Symptom:
   bars covering a line they were drawn under. Fix: pass `group=parent_vgroup`.

4. **Every `self.wait()` leaves an empty placeholder `Mobject` in `self.mobjects`.**
   A final z-order guard must skip mobjects with
   `len(m.get_family()) == 1 and not m.has_points()`.

5. **The cairo static-frame cache silently drops side-effect changes.**
   `get_moving_mobjects` returns everything from the first animation target or
   updater onward; earlier mobjects are rendered once into a cached static image for
   the whole play. So a mobject repainted as a side effect of *another* mobject's
   updater is never redrawn, and only sometimes, which makes it maddening to debug.
   Put the updater on the mobject it repaints.

6. **`AnimationGroup.begin()` begins every child up front.** A per-child side effect
   (un-hiding something, say) must fire at that child's first lazy `interpolate`, not
   in `begin()`.

7. **`tracker.animate.set_value()` puts the `ValueTracker` into the scene** and it
   carries one point, so a leaf-count guard flags it as a stray. `self.remove(tracker)`
   after the play.

8. **`Scene.remove(group)` does not remove children re-parented by a `LaggedStart`.**
   Remove the children.

9. **Never pass a plain Python list as `FadeIn(..., shift=...)`.** `list * -1` gives
   `[]` and the broadcast crashes. Worse, the crashed render *hangs*: the writer
   thread waits on its queue forever and it looks like a stall. Use `UP * 0.3` or a
   numpy array.

10. **`MathTex.get_parts_by_tex` (plural) does not exist** in some builds;
    `Mobject.__getattr__` intercepts it and throws a confusing `TypeError`. The
    singular `get_part_by_tex` returns only the first match. Also,
    `substrings_to_isolate` does not split `.submobjects`, so never scan those to find
    a part.

11. **`ImageMobject.set_opacity` scales `orig_alpha_pixel_array`.** If you mutate
    `pixel_array` yourself, keep that reference in sync or fades will undo your work.
    Mutating `pixel_array` per frame *is* supported: the camera re-reads it every
    frame with no cache, which is how a heatmap can recolour through dozens of frames.

12. **ManimCE only.** `Create`, not `ShowCreation`. Confirm any snippet you borrow
    from the web is ManimCE and not 3b1b's `manimgl`; the APIs diverged years ago and
    the error messages do not say so.

---

## 7. Checklist for a new clip

1. Decide the single thing the clip says. If two things move, cut one.
2. Write the builder for the final frame first.
3. Render at `--quality l`, grab the final frame, look at it.
4. Check: nothing clipped at 1920x1080, no unintended overlaps, motion reads.
5. Run the frame diff against the previous clip's final frame.
6. Render at `--quality h`, export the final frame PNG, deliver both.

## 9. The logo → slice build and outlined tracks (2026-09-15)

**Logo → slice** (`scenes/s2_cms_logo.py`). The detector geometry is the CMS
logo's (`LOGO_R`), so the logo is literally a quarter of the slice. Four
`ValueTracker`s (`side`, corner `cx`/`cy`, sector `angle`, box opacity) drive
`always_redraw(lambda: logo_rings(corner, side, angle))`: (a) the logo at
side 2.2 centred on the slide; (b) *expand*: side → `R_DET/0.996`, corner →
`DET_CENTER`, the off-white box and the frame fade to nothing — the square
grows toward the top right because its corner is the beam spot; (c) *round
up*: angle π/2 → 2π sweeps every layer counter-clockwise from 3 o'clock; the
four red logo muons (`logo_muons`, redrawn with the same trackers) fade during
the sweep; (d) the dynamic rings are swapped for `CMSSlice().base` (identical
pixels) and `CMSSlice().detail` fades in with a `LaggedStart`. Traps: `Create`
on a (rim, core) pair needs `lag_ratio=0` or the rim draws first; `Create(m)`
adds `m` itself to the scene, so remove the *submobjects* afterwards, not the
parent VGroup; the `Text` wordmark must be redrawn with an opacity tracker if
it is to ride along with the box.

**Outlined tracks.** `track(det, phi, kappa, outline=TRACK_OUTLINE[ch])`
returns VGroup(rim, core) with `.pts` kept on the group; the pale μμ gold
reads on white this way (the logo's `double` muon style). Channel deposits get
the same rim through `calo_hit(..., stroke_color=TRACK_OUTLINE[ch])`.

**Chained channel clips.** `scenes/channel_common.py`: `ChannelDiagram`
(the Drell-Yan diagram shifted by `DIAGRAM_SHIFT`, labels beside the legs),
`crossfade_to_detector`, `show_signature`, `zoom_in` (a `MovingCameraScene`
closes in on the inner detector at the end of the ee and ττ event displays —
the logo tracker is small). Each `*Detector` scene rebuilds the `*Process`
end state through `end_state()` and `tools/framediff.py` checks the seam.
