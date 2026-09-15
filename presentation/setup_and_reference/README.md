# setup_and_reference — how the Manim clips for the BND-school talk are made

Findings from the Thesis-Code (ATLAS four-top MSc defense) and the LHCb
momentum-calibration decks, adapted to this project. Read in this order:

| file | what |
|------|------|
| `01-environment.md` | the conda env, TeX, render command, what breaks |
| `02-deliverables-and-naming.md` | MP4-only delivery, section folders, numbering, versioning |
| `03-sections-storyboard.md` | the six sections, candidate clips per section (draft, to be decided with the user) |
| `04-agents.md` | the subagents, what each may touch, the per-clip workflow |
| `05-colour-schema.md` | the standard colour schema: NL/BE/DE flag colours, semantic roles, detector tints |
| `agents/*.md` | the agent definitions (symlinked into `.claude/agents/`) |
| `../MANIM_PHYSICS_RECIPES.md` | drawing recipes: Feynman lines, detector, tracks, mass plots, chaining, trap list |
| `../CLAUDE.md` | the short contract every agent reads first |

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
- **LHCb deck** (`MANIM_PHYSICS_RECIPES.md`): the wavy-boson helper, the vertex
  coordinate dictionary, build-order = z-order, the track model, the invariant
  mass plot sequence (one event → one bar → the rain), the double-sided Crystal
  Ball with radiative tail, clip chaining with builder functions and pixel
  diffs, ending on the last change, and the 12-item trap list.

## What is different here

- No `.pptx` shape handoff and no `_final.png` in the delivery folder: the user
  asked for MP4s only. Frames still get written to `work/` for our own checks.
- Six fixed sections with numbered clip folders (see `02-…`).
- Three-channel colour code (ee Dutch blue, μμ Dutch red, ττ German gold; the
  whole palette is the NL/BE/DE flags, `05-colour-schema.md`) is the one
  load-bearing colour convention; the fit-input sample names from
  `fitting/CONVENTIONS.md` each have a colour in the palette.
- The CMS slice replaces the ATLAS slice (five barrel layers, solenoid between
  HCAL and the muon system), and tracks curve as arcs in the solenoid field.
- A first clip exists as a pipeline proof: `1_theory/1-01_drell_yan_v1.mp4`
  (q q̄ → Z/γ* → ℓ⁺ℓ⁻), with `DrellYanEE/MuMu/TauTau` variants in the same file.
