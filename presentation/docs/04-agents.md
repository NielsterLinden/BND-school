# 04 — Agents

Definitions live in `.claude/agents/*.md` (repository root), where Claude Code loads them. Invoke
with `@<name>` or by asking for the task. All of them own `presentation/`
only and never edit analysis code, fit inputs, or channel docs.

| agent | job | may edit | Bash |
|-------|-----|----------|------|
| `manim-render-loop` | render at -q l, look at frames, fix layout/timing/overlap, repeat | `scenes/**`, `style/bnd_style.py` | render + ffmpeg |
| `physics-explainer` | build mechanism / method animations (Feynman, detector, tag-and-probe, fit) grounded in the channel docs; schematic values only | `scenes/**`, new primitives in `style/bnd_style.py` | render |
| `plot-recreator` | freeze channel numbers into `data/*.json` (in the LCG env) and draw result plots natively, asserting against handoff numbers | `data/**`, `scenes/**`, plot helpers in `style/` | LCG python for extraction; manim for render |
| `physics-fidelity-checker` | read-only critic: depicted physics, symbols and numbers vs handoffs / docs / CONVENTIONS | nothing | none |
| `style-warden` | no hardcoded hex/font, palette use, channel colour code, no narrative text | `scenes/**` (style fixes only) | none |
| `clip-deliverer` | final -q h render, numbering/versioning via `tools/render.py`, registry, report | `tools/**`, `clips/CLIPLIST.tsv`, `docs/03-…` | render, ffmpeg |
| `channel-animation-director` | **entry point for a channel team**: long interview → chapter brief → orchestrates the agents below and above → drafts to critique → delivery | `docs/briefs/<chapter>.md`, the team's rows of `03-…` | ffmpeg (preview) |
| `channel-story-builder` | implements an approved brief as a chained clip set, renders drafts once at -q l, applies critique | `scenes/s<S>_<chapter>_story.py`, add-only helpers in `style/bnd_style.py` | render |
| `deck-continuity-checker` | read-only: the chapter still fits the talk (fixed anchors of `06-…`, entry/exit frames, symbols, add-only shared files), judged against the brief | nothing | none |

## Channel chapters (Z→ee, Z→ττ, …)

Teams start from `.claude/prompts/presentation_chapter_animations.md` (paste it into
`claude --agent channel-animation-director` or a normal session). The director
interviews the team, writes `briefs/<chapter>.md`, and runs:
`plot-recreator` (numbers) → `channel-story-builder` (drafts) → preview → critique
rounds → `physics-fidelity-checker` + `style-warden` + `deck-continuity-checker` →
`clip-deliverer`. Shared anchors vs free choices: `06-chapter-anchors.md`; worked
example: `briefs/zmumu.md`.

## Per-clip workflow

1. User (or the storyboard) names the clip and the single idea it shows.
2. `physics-explainer` or `plot-recreator` writes the scene (builder first,
   then the animation). For data plots, the numbers are frozen first.
3. `manim-render-loop` iterates at `-q l` until the frames are right.
4. `physics-fidelity-checker` and `style-warden` review (read-only / trivial fixes).
5. `clip-deliverer` renders at `-q h`, delivers `clips/<S>_…/<S>-<NN>_<name>_v<K>.mp4`,
   updates the storyboard table, reports the path.

One agent at a time on the render env (serial renders).
