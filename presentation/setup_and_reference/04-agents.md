# 04 — Agents

Definitions live in `presentation/setup_and_reference/agents/*.md` and are
symlinked into `/project/atlas/Users/nterlind/BND-school/.claude/agents/` so
Claude Code loads them (if a symlink is not picked up, copy the file). Invoke
with `@<name>` or by asking for the task. All of them own `presentation/`
only and never edit analysis code, fit inputs, or channel docs.

| agent | job | may edit | Bash |
|-------|-----|----------|------|
| `manim-render-loop` | render at -q l, look at frames, fix layout/timing/overlap, repeat | `scenes/**`, `style/bnd_style.py` | render + ffmpeg |
| `physics-explainer` | build mechanism / method animations (Feynman, detector, tag-and-probe, fit) grounded in the channel docs; schematic values only | `scenes/**`, new primitives in `style/bnd_style.py` | render |
| `plot-recreator` | freeze channel numbers into `data/*.json` (in the LCG env) and draw result plots natively, asserting against handoff numbers | `data/**`, `scenes/**`, plot helpers in `style/` | LCG python for extraction; manim for render |
| `physics-fidelity-checker` | read-only critic: depicted physics, symbols and numbers vs handoffs / docs / CONVENTIONS | nothing | none |
| `style-warden` | no hardcoded hex/font, palette use, channel colour code, no narrative text | `scenes/**` (style fixes only) | none |
| `clip-deliverer` | final -q h render, numbering/versioning via `tools/render.py`, registry, report | `tools/**`, `clips/CLIPLIST.tsv`, `setup_and_reference/03-…` | render, ffmpeg |

## Per-clip workflow

1. User (or the storyboard) names the clip and the single idea it shows.
2. `physics-explainer` or `plot-recreator` writes the scene (builder first,
   then the animation). For data plots, the numbers are frozen first.
3. `manim-render-loop` iterates at `-q l` until the frames are right.
4. `physics-fidelity-checker` and `style-warden` review (read-only / trivial fixes).
5. `clip-deliverer` renders at `-q h`, delivers `clips/<S>_…/<S>-<NN>_<name>_v<K>.mp4`,
   updates the storyboard table, reports the path.

One agent at a time on the render env (serial renders).
