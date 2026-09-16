**Template prompt: build the animated chapter of a channel (Z→ee, Z→ττ, …) for the BND-school talk**

The agents live in `presentation/setup_and_reference/agents/` and are loaded by Claude Code
from `.claude/agents/` (symlinks). Pull first; start Claude Code **in the repository root**
(a session started before the pull does not see new agents).

Two ways to start, pick one:

1. `claude --agent channel-animation-director`, then paste the block below.
2. A normal `claude` session in the repository root, paste the block below as your first message.

Fill in the `<…>`; leave a line out if you don't know yet (the director will ask).

---

```text
I'm <your name / role> from the Z→<ee | ττ | …> group. I want to make the animated chapter for
our channel in the BND-school Z cross-section talk.

Work as the channel-animation-director: read presentation/setup_and_reference/agents/channel-animation-director.md
completely and follow it (if you are already running as that agent, just follow your instructions).
Before building anything, interview me with many questions (AskUserQuestion, in rounds) until the
chapter brief is approved. Then build it with the project's agents:
  - plot-recreator            freezes our numbers into presentation/data/<chapter>_*.json (LCG env)
  - channel-story-builder     writes scenes/s<S>_<chapter>_story.py and renders the drafts at -q l
  - physics-fidelity-checker, style-warden, deck-continuity-checker   review once, at delivery
  - clip-deliverer            renders -q h into presentation/clips/ when I approve
Keep the shared anchors of presentation/setup_and_reference/06-chapter-anchors.md; everything
else is our choice. Hand me the draft preview quickly so I can critique it; don't spend long
checking things yourself.

Chapter: <ee | tautau | …>        Section: <3_zee | 5_ztautau | …>
Time budget for the animated part: <e.g. 3 minutes>
What we want the audience to remember: <one sentence, or "help me choose">
Things we already know we want: <e.g. the visible-mass story, a real event display, the fake factor>
Things we don't want: <e.g. no tag-and-probe, no numbers before the result>
Result version to show: <e.g. v3 Tight nominal, or "ask me">
Git: <commit locally when approved | leave uncommitted>; push only when I say so.
```

---

What you get back: a brief in `presentation/setup_and_reference/briefs/<chapter>.md`, draft clips in
`presentation/work/<clip>/`, a preview `presentation/work/preview/<chapter>_preview_480p.mp4`, and after
your approval the numbered MP4s in `presentation/clips/<S>_<section>/`.

Reading list if you want to know the rules first: `presentation/CLAUDE.md`,
`presentation/setup_and_reference/06-chapter-anchors.md`, the worked example
`presentation/setup_and_reference/briefs/zmumu.md`.
