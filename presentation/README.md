# presentation — the Manim clips of the BND-school Z cross-section talk

The talk itself is `2026.09.18 -- Presentation BND -- OpenData.pptx` in the repository root (PowerPoint, the clips
embedded). This folder is how its 99 animated clips were made: `clips/` (the delivered MP4s, `CLIPLIST.tsv`),
`scenes/` (one Manim scene file per chapter), `style/` (palette and every drawing primitive), `data/` (the frozen
numbers the clips show, with the extractors that wrote them), `tools/` (render, deliver, seam and keep-out checks).

How it works, adapted from the Thesis-Code (ATLAS four-top MSc defense) and the LHCb momentum-calibration decks.
Read in this order:

| file | what |
|------|------|
| `docs/01-environment.md` | the conda env, TeX, render command, what breaks |
| `docs/02-deliverables-and-naming.md` | MP4-only delivery, section folders, numbering, versioning |
| `docs/03-sections-storyboard.md` | the seven sections and every delivered clip: what moves, and the source of what it shows |
| `docs/04-agents.md` | the subagents, what each may touch, the per-clip workflow |
| `docs/05-colour-schema.md` | the standard colour schema: six chapter colours + CMS logo colours, semantic roles, detector tints |
| `docs/06-chapter-anchors.md` | what every chapter shares (title band, colours, symbols, number format) and what is free |
| `docs/07-manim-physics-recipes.md` | drawing recipes: Feynman lines, detector, tracks, mass plots, chaining, trap list |
| `docs/briefs/*.md` | the approved brief of each chapter (story, numbers and their version, clip list) |
| `CLAUDE.md` | the short contract every agent reads first; the agent definitions are in `.claude/agents/` |

## The system in one paragraph

Plain ManimCE `Scene`s (v0.20.1, cairo renderer, headless) render continuous
MP4 clips on a white background. Every colour, font and drawing primitive comes
from one style module (`style/bnd_style.py`, palette in the manim-free
`style/palette.py`). Clips contain animation and physics symbols only; the user
writes every title and caption in PowerPoint. Result plots read frozen JSON/CSV
in `data/`, never the analysis stack. `tools/render.py` renders at low quality
into `work/` for inspection (always look at the PNG frames it writes), then at
high quality delivers a numbered, versioned MP4 into `clips/<section>/`. Longer
stories are cut into chained clips that open on the previous clip's final frame,
verified with `tools/framediff.py`.

## What was taken from the two earlier decks

- **Thesis-Code/presentation** (ATLAS, 60+ scenes): the assembly model
  (PowerPoint by the user, MP4 clips by the agents), the animation-only rule,
  the style-module discipline, the `DataAxes` + `data_*` plotting primitives,
  the Feynman `gluon`/`fermion`/dashed helpers, the transverse detector slice
  class, the agent split (render loop / plot recreator / fidelity checker /
  style warden / assembler / explainer), the render command and its `-qh` trap.
- **LHCb deck** (`docs/07-manim-physics-recipes.md`): the wavy-boson helper, the vertex
  coordinate dictionary, build-order = z-order, the track model, the invariant
  mass plot sequence (one event → one bar → the rain), the double-sided Crystal
  Ball with radiative tail, clip chaining with builder functions and pixel
  diffs, ending on the last change, and the 12-item trap list.

## What is different here

- No `.pptx` shape handoff and no `_final.png` in the delivery folder: the user
  asked for MP4s only. Frames still get written to `work/` for our own checks.
- Seven fixed sections with numbered clip folders (see `docs/02-…`).
- Three-channel colour code (ee green, μμ gold, ττ red; the
  whole palette is the six chapter colours, `docs/05-colour-schema.md`) is the one
  load-bearing colour convention; the fit-input sample names from
  `docs/CONVENTIONS.md` each have a colour in the palette.
- The CMS slice replaces the ATLAS slice (five barrel layers, solenoid between
  HCAL and the muon system), and tracks curve as arcs in the solenoid field.
- The Drell-Yan diagram (q q̄ → Z/γ* → ℓ⁺ℓ⁻, `scenes/s1_drell_yan.py`, with `DrellYanEE/MuMu/TauTau` variants) was the
  pipeline proof; it is now the base class of the process clips rather than a clip of its own.
