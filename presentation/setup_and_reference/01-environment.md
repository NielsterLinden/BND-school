# 01 — Environment and render

## The env (already built, shared with Thesis-Code)

Isolated conda-forge env at `/data/atlas/users/nterlind/venvs/presentation`
(python 3.11, ffmpeg, pango, pycairo from conda-forge; `manim==0.20.1` and
`manim-slides==5.6.0` from pip). Rebuild recipe if it ever disappears:
`Thesis-Code/presentation/environment.yml` + `requirements.txt`.

TeX for `MathTex`: the user's TeX Live 2026 (`standalone.cls` and `dvisvgm`
verified present). CVMFS TeX Live is the fallback on a CERN node.

```bash
export PATH="/data/atlas/users/nterlind/venvs/presentation/bin:/data/atlas/users/nterlind/texlive/2026/bin/x86_64-linux:$PATH"
cd /project/atlas/Users/nterlind/BND-school/presentation
manim --version        # Manim Community v0.20.1
```

Verified on `stbc-i2` (28 cores, 250 GB) on 2026-09-14: `manim`, `ffmpeg`,
`latex`, `dvisvgm`, `kpsewhich standalone.cls` all resolve; the Drell-Yan
scene renders at `-q l` and `-q h`.

**Never** mix with `source setup.sh` (LCG ROOT/python 3.13) in the same shell.
Analysis code freezes numbers into `presentation/data/`; the Manim env only reads them.

## Render

```bash
python tools/render.py <section 1-6> <name> scenes/<file>.py <Class> [-q l|h] [--frames 1.5,3,end]
```

Under the hood: `manim render --quality l|h --media_dir work/<name>/media -o <name> scenes/<file>.py <Class>`.
- `--quality h` = 1920×1080 at 60 fps; `--quality l` = 854×480 at 15 fps (iterate).
- Write `--quality h`, never `-qh` (the short form collides with the manim-slides
  `-h` help flag while manim-slides is installed in the env).
- Default cairo renderer; `--renderer=opengl` needs a display and fails here.
- Render **serially**: concurrent renders that compile the same `MathTex`
  string race on the cached SVG and one dies with `FileNotFoundError`.
- A failed scene (exception or assert) leaves the previous MP4 in place; the
  script checks for `Rendered` in the log and stops otherwise.
- Frame grab for inspection: `tools/grab.sh clip.mp4 out.png [t|end]`
  (`-sseof -0.15`; `-0.05` decodes nothing at 15 fps).

## Inspect before you deliver

Look at every frame the render script writes (Read the PNG). Check: nothing
clipped at the 16:9 border, no overlaps, motion reads, the last frame is the
intended held state, the top ~15 % of the frame is free for the user's title and the
top-left block (3 × 9 cm: x < −5.85, y > 0.22) is empty (`tools/keepout.py` prints it).

## Known traps (full list: MANIM_PHYSICS_RECIPES.md §6)

1. The last rendered frame precedes the final `interpolate(1)`: end every scene
   with `self.wait(0.1)`.
2. `Indicate` / `there_and_back` leave the mobject at peak scale; pulse a copy.
3. `LaggedStart` over children of an added `VGroup` re-stacks z-order: pass `group=`.
4. `Flash` never inside a `Succession`; own `play` after the flights.
5. `MoveAlongPath(dot, trk.copy())` when `trk` is being `Create`d in the same group.
6. `FadeIn(..., shift=UP*0.3)` with numpy, never a python list (hangs the writer).
7. ManimCE API only (`Create`, not `ShowCreation`).
