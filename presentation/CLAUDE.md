# CLAUDE.md — presentation/ (Manim clips for the BND-school Z cross-section talk)

Read the repository guide `.claude/CLAUDE.md` first, then this. Overview: `presentation/README.md`; conventions:
`presentation/docs/01 … 07`; chapter briefs: `presentation/docs/briefs/`; agent definitions: `.claude/agents/`.

> **Final state (18 Sep 2026).** The talk was given with `2026.09.18 -- Presentation BND -- OpenData.pptx` (repository
> root). `clips/` holds exactly the 99 clips that deck uses, one file per clip (its latest version), listed in
> `clips/CLIPLIST.tsv` and described in `docs/03-sections-storyboard.md`. Earlier cuts, superseded versions and scenes
> that fed no delivered clip were removed; they are in the git history (tag `pre-cleanup-2026-09-18`).
Drawing know-how: `presentation/docs/07-manim-physics-recipes.md`.
A channel team building its own chapter starts from `.claude/prompts/presentation_chapter_animations.md`
(agent `channel-animation-director`); what all chapters share and what is free:
`docs/06-chapter-anchors.md`.

## What this directory produces

**MP4 clips only.** The user builds the 30-minute deck in PowerPoint and inserts
the clips. No manim-slides, no HTML/PDF deck, no speaker notes, no `.pptx` handoff,
no PNGs in the delivery folder. Six sections plus the outro:

| S | folder            | topic                      |
|---|-------------------|----------------------------|
| 1 | `1_theory`        | theory and relevance       |
| 2 | `2_cms_methods`   | CMS and methods            |
| 3 | `3_zee`           | Z → ee                     |
| 4 | `4_zmumu`         | Z → μμ                     |
| 5 | `5_ztautau`       | Z → ττ                     |
| 6 | `6_combination`   | combination and conclusion |
| 7 | `7_outro`         | closing slide (mic drop)   |

Deliverable: `presentation/clips/<S>_<section>/<S>-<NN>_<name>_v<K>.mp4`
(1920×1080, 60 fps, H.264, white background). `NN` = clip number within the
section from `clips/CLIPLIST.tsv` (creation order, never reused), `K` = version.
`tools/render.py` assigns both; never name a delivered file by hand.
A re-render delivers the next version `_v<K+1>` next to the old one; remove the superseded file once the new one
is in the deck (the history keeps it).

## Environment (every Bash session)

```bash
export PATH="/data/atlas/users/nterlind/venvs/presentation/bin:/data/atlas/users/nterlind/texlive/2026/bin/x86_64-linux:$PATH"
cd /project/atlas/Users/nterlind/BND-school/presentation
```

Do **not** `source fitting/setup.sh` (the LCG stack) in the same shell: the analysis and
the Manim env are separate. Cairo renderer only (never `--renderer=opengl`).
Render serially (two renders compiling the same `MathTex` race on the cache).

## Render loop

```bash
python tools/render.py 4 mumu_process scenes/s4_zmumu.py ZmumuProcess            # -q l, into work/
python tools/render.py 4 mumu_process scenes/s4_zmumu.py ZmumuProcess -q h       # deliver to clips/
python tools/framediff.py work/mumu_process/mumu_process_final.png work/mumu_into_detector/mumu_into_detector_t0.png  # chain seam
python tools/keepout.py clips/5_ztautau/*.mp4          # nothing in the title band / top-left block
# a chain cut into one-idea clips (clip_open / clip_cut = manim sections): renders the scenes in
# order, joins equal section names across scenes, numbers + delivers each clip, prints all seams
python tools/deliver_chain.py 4 scenes/s4_zmumu_story.py MumuEvent MumuFill MumuStack MumuTagProbe \
    MumuTagProbeScan MumuCorrections MumuFit -q h \
    --seam-before clips/4_zmumu/4-02_mumu_into_detector_v1.mp4
```
Seams at 480p fail on mean |d| alone (blob 0) from H.264 noise even between two cuts of
the same rendered frame; judge seams at 1080p.

Scenes (clip → file and class: `clips/CLIPLIST.tsv`): `s1_theory_story.py` (section 1), `s2_cms_logo.py` (logo →
slice), `s2_map.py` (the two-row map of section 2), `s3_zee_story.py` (the ee chapter on the map: selection and fit),
`s4_zmumu.py` + `s4_zmumu_story.py`, `s5_ztautau.py` + `s5_ztautau_story.py`, `s6_combination_story.py`,
`s7_outro.py` (the results plot turns into the mic drop; flag-coloured rays: `FLAG` in the palette, outro only).
Libraries without a clip of their own: `s1_drell_yan.py` (the Drell-Yan diagram, base of the process clips),
`channel_common.py` (process clip → detector clip), `s2_pipeline.py` (icons, camera, `masked_play`) and
`s2_trigger.py` (the trigger chain drawn inside the map).

Always look at the frames `tools/render.py` writes to `work/<name>/` (Read the
PNG) before delivering. Iterate at `-q l`; deliver at `-q h`.

## Rules

- **Animation only.** No titles, captions, bullets or narrative text in a clip.
  Exempt: particle/process symbols (`q`, `Z/γ*`, `μ⁺`, `τ_h`, `p_T`, `η`, `m_ℓℓ`),
  axis tick labels and numeric values that *are* the data. When in doubt, leave
  text out.
- **Style module only.** Every colour and font from `style/bnd_style.py`
  (palette in `style/palette.py`); never hardcode a hex. **The colour schema is
  the six chapter colours** (`CHAPTER`: 1 theory purple `#702677`, 2 detector
  cyan `#00B3C3`, 3 ee green `#059F77`, 4 μμ gold `#E3D88B`, 5 ττ red `#E6223C`,
  6 combination slate `#2B2D42`) and their tints/shades (`palette.tint/shade`).
  The channel colours (ee green, μμ gold, ττ red; thin strokes and glyphs via
  `CHANNEL_LINE`) are the colour code of the whole talk. The CMS detector alone
  is drawn in the official CMS logo colours (`CMS`, `DETECTOR`), logo radii,
  centred at `DET_CENTER=(0,-0.25)`, outer radius `R_DET=2.95`.
  See `docs/05-colour-schema.md`.
- **Frozen numbers.** Result plots read `presentation/data/*.json|csv` written
  by the channel code; a scene never imports uproot/awkward/ROOT or recomputes.
  Every printed number traces to a channel `handoff.md` / docs page.
- **End on the last change.** `self.wait(0.1)` at the end, no long held frame;
  PowerPoint rests on the last frame.
- Keep the top ~15 % of the frame (y > 2.7) clear: the user's title band.
- Keep the **top-left block** clear too: 3 cm wide × 9 cm tall on the slide
  (x < −5.85 and y > 0.22; `CORNER_X_MAX`, `CORNER_Y_MIN`, `in_keepout` in
  `style/bnd_style.py`): the deck's chapter identifier. `tools/keepout.py` scans an MP4 for
  both regions; `render.py` / `deliver_chain.py` print its line after every render. Exempt
  (full-frame by design): the zoom ending of `tautau_detector` (and `tautau_a1_event`, which opens on it),
  the outro rays.
- Never commit `work/` or `media/` (see `.gitignore`); commit scenes, style,
  tools, data JSON, docs and the delivered MP4s in `clips/`.
