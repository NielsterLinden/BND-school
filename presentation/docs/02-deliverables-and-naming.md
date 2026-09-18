# 02 — Deliverables and naming

## What the user gets

**Only MP4 files**, grouped by section, numbered, versioned:

```
presentation/clips/
  CLIPLIST.tsv                      registry: number  section  name  scene_file  scene_class
  1_theory/        1-02_sm_to_process_v1.mp4   1-03_<name>_v<K>.mp4 ...
  2_cms_methods/   2-01_...
  3_zee/           3-03_...
  4_zmumu/         4-01_...
  5_ztautau/       5-02_...
  6_combination/   6-02_...
  7_outro/         7-02_...
```

- `<S>` section 1..7. `<NN>` clip number inside the section, assigned by
  `tools/render.py` in creation order and recorded in `CLIPLIST.tsv`; it
  identifies the clip, it is **not** the slide order (the user orders in
  PowerPoint). Numbers are never reused or renumbered.
  *One exception, made by the deck owner on 17 Sep 2026:* section 4 (Z→μμ) was renumbered once for its
  final pass, 4-01…4-21 in play order (`briefs/zmumu.md` §7); sections 2 and 5 were re-cut the same way.
  Since 18 Sep 2026 the folders hold the final clips only, one file per clip: superseded versions, earlier cuts
  and clips the talk does not use were removed (git history, tag `pre-cleanup-2026-09-18`), so gaps in the
  numbering (no 1-01, 1-05, 1-06, 3-01, 3-02, 5-01, 6-01, 7-01) are retired clips.
- `_v<K>` version. A re-render of the same clip gets the next version; old
  versions stay until the user deletes them, so a PowerPoint that links `_v2`
  never breaks.
- `<name>` snake_case, short, the *thing shown* (`drell_yan`, `cms_slice_muon`,
  `mass_peak_fill_ee`), never a slide title.

Format: 1920×1080, 60 fps, H.264 (yuv420p), white background, no audio.
Insert in PowerPoint with Start: Automatically, "Rewind after Playing" off,
so the slide rests on the last frame.

No `_final.png`, no `.pptx` shapes in `clips/` (the user asked for MP4s only).
The inspection frames live in `work/<name>/` and are git-ignored.

## Transfer to the laptop

`clips/` is committed (MP4s included, ~43 MB on 16 Sep 2026): `git pull` on the
laptop, or copy with `rsync`/`scp` from `stbc-i*`.

## Chained clips

A long story is several clips where clip N+1 opens on clip N's final frame
(same builder function, same `ORDER`, checked with
`python tools/framediff.py work/<a>/<a>_final.png work/<b>/<b>_t0.png`).
Name them with a shared prefix and a letter: `mass_peak_a_one_event`,
`mass_peak_b_rain`, `mass_peak_c_numbers`.

## Side-by-side clips

Clips meant to play simultaneously on one slide (e.g. the three channel
Feynman diagrams) share one layout and identical `run_time`s so they stay in
sync when started together. Note the pairing in the scene docstring.
