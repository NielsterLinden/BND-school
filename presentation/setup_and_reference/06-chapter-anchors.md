# 06 — Chapter anchors: what every channel chapter shares, and what is free

Written 2026-09-16 for the teams building their own chapter (Z→ee in `3_zee`,
Z→ττ in `5_ztautau`, or any later chapter). The Z→μμ chapter (`4_zmumu`,
`scenes/s4_zmumu_story.py`) and the method chain (`2_cms_methods`,
`scenes/s2_pipeline.py`) are the worked examples, **not templates to copy**.
A chapter may tell a completely different story. It keeps the anchors below
so the audience always knows where they are.

Three levels:

| level | meaning | who may change it |
|-------|---------|-------------------|
| **Fixed** | the talk breaks visually or physically without it | the deck owner only (ask first) |
| **Default** | what the other chapters do; reuse unless your story needs otherwise | your team, recorded in your brief |
| **Idea** | a menu to pick from | your team |

---

## A. Fixed anchors (the whole talk)

### A1. Frame and format
- 1920×1080, 60 fps, white background (`white_background(self)`), MP4 only, delivered by `tools/render.py`.
- **y > 2.7 stays empty**: that is the user's PowerPoint title band.
- **The top-left block stays empty**: 3 cm wide × 9 cm tall on the slide, x < −5.85 and
  y > 0.22 (the 33.867 × 19.05 cm slide, 1 cm = 0.42 scene units): the deck's chapter identifier
  sits there (deck owner, 16 Sep 2026). Same status as the title band — every clip, every frame,
  including labels, parked plots, lists and keys. Constants `TITLE_BAND_Y`, `CORNER_X_MAX`,
  `CORNER_Y_MIN` and `in_keepout(mobject)` in `style/bnd_style.py`; check the rendered MP4 with
  `python tools/keepout.py <mp4>` (`render.py` / `deliver_chain.py` print it after every render).
  Layout recipe: long left-hand columns (value lists, corrections, pulls, impact names) start
  at x ≥ −5.8 above y = 0.22 or live entirely below it; a camera zoom picks its centre so the
  block stays white. Known exceptions (`keepout.EXEMPT`): a zoomed event display fills the title
  band (3-02 `ee_detector`, 5-03 `tautau_detector` and 5-04 `tautau_a1_event`, which opens on it;
  accepted before this rule) and the outro rays.
- **No narrative text** in a clip. Allowed: physics symbols (`e^\pm`, `\tau_h`, `p_T`,
  `m_{\ell\ell}`, `\mu_Z`), axis ticks, numbers that *are* the data.
- A clip ends on its last change (`self.wait(0.1)`); PowerPoint rests on the last frame.
- Chained clips open on the previous clip's exact final frame (pure builder per end state,
  `ORDER` tuple, `add_state` / `check_order`; seam checked with `tools/framediff.py` at delivery).

### A2. Colour roles (`style/palette.py`, `05-colour-schema.md`)
Every colour comes from `style.bnd_style`; never a hex in a scene.

| role | colour | constant |
|------|--------|----------|
| Z→ee: electrons, signal fill, the chapter's own accents | green `#059F77` | `CHANNEL["ee"]`, `CHANNEL_LINE["ee"]`, `PARTICLE["e"]` |
| Z→μμ | gold `#E3D88B` (lines: shade 0.30) | `CHANNEL["mumu"]`, `CHANNEL_LINE["mumu"]` |
| Z→ττ: τ_h, signal fill | red `#E6223C` | `CHANNEL["tautau"]`, `PARTICLE["tau"]` |
| theory / prediction line or band | purple | `THEORY` |
| method / detector-that-is-not-the-drawing (efficiencies, sliders, "a measurement is happening") | cyan | `DETECTOR_ACCENT` |
| data points, axes, ink | slate | `SAMPLE["Data"]`, `INK` |
| combined result | slate | `COMBINED` |
| backgrounds | tints of slate / cyan / purple | `SAMPLE[...]` (never your channel colour) |
| "look here" flash | red; **in ττ clips slate** (red is the channel) | `HIGHLIGHT` / `SLATE` |
| the CMS detector drawing | official logo colours | `CMSSlice`, `DETECTOR`, `DEPOSIT` |

**Never use another channel's colour for your own content** (no gold in the ee chapter,
no green in the ττ chapter) except where the three channels appear together.
A colour that is not a chapter colour or a `tint`/`shade`/`lighten`/`darken` of one
needs the deck owner's approval and a palette entry (see A5 for known gaps).

### A3. Physics symbols and names (`fitting/CONVENTIONS.md`)
- Samples `DYee, DYmumu, DYtautau, TTbar, SingleTop, WW, WZ, ZZ, Fakes` (your analysis may
  group them differently, e.g. "Diboson"; draw with the convention names and colours, or
  ask for a palette entry).
- POI `\mu_{Z}`; regions `ee_SR, mumu_SR, tautau_SR` (never printed as words in a clip).
- L = 16393.381 pb⁻¹; σ(Z/γ*→ℓℓ, m > 50) = 6077.22/3 pb per flavour; window 60 < m_ℓℓ < 120 GeV.
- Mass symbols: `m_{ee}`, `m_{\mu\mu}`, `m_{\tau\tau}` (say which one: visible, collinear, MMC...).
- Result format as printed by μμ: `\sigma = 1931 \pm 30\ \mathrm{pb}` (frozen result, 17 Sep 2026); split uncertainties with
  `_{\mathrm{stat}}`, `_{\mathrm{syst}}`, `_{\mathrm{lumi}}`; asymmetric ones as `^{+222}_{-194}`.

### A4. Frozen numbers
- A printed number comes from `presentation/data/<chapter>_*.json`, written by an
  `extract_<chapter>_*.py` in the **LCG** env, never computed in a scene.
- Each JSON carries provenance (source file, selection, version/tag, date) and the scene
  asserts its anchor values at import, so a silently changed input stops the render.
- **Numbers move** (z-tautau v2.1 → v3 changed μ_Z from 1.205 to 1.071; the combination
  still quotes v2.1). A chapter shows the version its team names in the brief, and the
  brief records which one.

### A5. Detector
- The slice is `CMSSlice()` at `DET_CENTER = (0, -0.25)`, `R_DET = 2.95`; tracks via
  `track(..., outline=TRACK_OUTLINE[ch])` and `signature(det, kind, phi)`.
- Channel leptons light their deposits in their chapter colour; everything else in `DEPOSIT`.
- Known palette gaps (raise them, don't invent a hex): no `SAMPLE["WJets"]` (ee and ττ use
  W+jets), no colour for a non-fiducial `DYtautau` part, no "fail" / "anti-ID" role.

---

## B. Default anchors (the story skeleton every chapter already has)

### B1. How the talk hands you the audience
```
section 2: pipe_a_map … pipe_h_three   spine of 7 nodes, schematic; ends with three strips
                                       (ee green / μμ gold / ττ red)
section 3–5, per channel:  <ch>_process  →  <ch>_detector          (delivered: 3-01/3-02,
                           Feynman act       event display act       4-01/4-02, 5-01…5-03)
                           → your chain (this is what you build)
section 6: three_to_one                three channel points slide into one combined point
                                       next to the purple theory line
```
- **Entry**, default: chain on your detector clip's final frame (3-02 `ee_detector`,
  5-03 `tautau_detector`), as μμ chains on 4-02. Alternative: open fresh on a white frame.
- **Exit**, default: your last clip ends on **the chapter's result as a point with error
  bars beside the `THEORY` prediction line** (or, if there is no result yet, on the plot that
  the result will come from). That is the object section 6 collects.

### B2. The pipeline spine as a locator
Section 2 teaches seven nodes: **detector → files (data + simulation) → selection →
corrections (tag-and-probe) → backgrounds (control region) → comparison (stack + ratio) →
fit → σ**. The audience knows this map. A chapter:
- picks the 2–4 nodes where it differs from the others and spends its time there;
- may skip nodes entirely (ee has no fakes step today; ττ's story is the mass and the fakes);
- may recall the spine as a small strip with the current node lit in the channel colour
  (`spine(..., done=..., accent=CHANNEL_LINE[ch])`) — optional, μμ does not.

### B3. Layout defaults (from `s4_zmumu_story.py`; reuse so chapters look like siblings)
| thing | default |
|-------|---------|
| detector parked beside a plot | scale 0.42 at (−4.6, −0.7) |
| main plot | centred on (1.85, −0.25); main panel centre (1.85, 0.85), ratio panel (1.85, −1.55) |
| mass axis | 60–120 GeV; log y when the backgrounds must be visible |
| colour key | one row, left end at (−1.25, −2.92) |
| control region | the slice (scale 0.84) inside a dashed box |
| result | σ axis with the prediction line; data point in the channel line colour |

### B4. Recurring motifs (the audience has seen them; reuse = instant recognition)
- **One event → one entry → the rain** (`rain`, `clock`): a real or schematic event in the
  slice becomes one count in the mass plot, then many.
- **Stack + data + ratio** (`stack_hist`, `data_dot`, `ratio_panel`, `colour_key`).
- **The data never moves**; only the prediction and the ratio move (physics honesty).
- **Control region box** → template **slides into the stack**.
- **Slider** with a `THEORY` reference (`slider`) for μ_Z or a scale factor.
- **Pull plot** (`pull_plot`) for nuisance parameters.
- **Value grid** (`value_grid`) for a fake-factor or scale-factor map.

### B5. Rules from the deck owner's critique of the ττ chapter (2026-09-16; apply to every chapter)
- **Shapes that explain a method show one process only.** A "before / after" of a reconstruction (m_vis → m_ττ, a
  calibration, a resolution) is drawn from the signal simulation of one sample with one clean peak; never a stitched
  or mixed template whose sub-samples add shoulders the audience cannot place.
- **Physics objects are named, never coded.** τ_h decay modes are `1-prong`, `1-prong+π⁰`, `3-prong`, `3-prong+π⁰`
  (not DM0/DM1/DM10/DM11); the same for any working point, era or category that has a code in the analysis, and for
  nuisance-parameter labels in a pull plot (`\tau_h\ \mathrm{ID}\ (1\text{-}\mathrm{prong})`, not `TauID_DM0`;
  a category by its cut, `D_{\mathrm{BDT}} > 0.90`, not `c2`).
- **A correction is shown after its disagreement.** Before a scale factor, closure correction or reweighting appears
  (`f(|η|)`, `g(p_T)`, `SF`), show the data-vs-prediction plot it fixes (prediction filled, data points, a ratio strip);
  then the ratio becomes the correction and the bars slide to 1.
- **Multi-dimensional dependences are written, not tabulated.** One representative map (`value_grid`) at most; the full
  dependence as a schematic function, `f = f(\mathrm{period}, \mathrm{DM}, N_{\mathrm{jets}}, p_T)`,
  `C = C(\mathrm{period}, N_{\mathrm{jets}}, D_{\mathrm{BDT}})`. Ranges as `x \in [a;\ b]`, not `a \ldots b`.
- **A classifier gets its own scene, from a blank frame**: the list of inputs (mass variables visibly absent when the
  score must be mass-agnostic) → unit-normalised distributions of the dominant background vs the signal, input by
  input → a schematic of the classifier → its output `D_{\mathrm{BDT}}` with data, a ratio panel and a y floor high
  enough that the shapes fill the panel → the categories it defines.
- **Uncertainty groups are descriptive.** No fit-internal jargon in the impact ranking: "Gammas" is
  "Template stat. (γ)", and any group name a non-expert cannot read is spelled out.
- **The result frame is a comparison.** The σ axis is large, the result label sits above it, and it carries the
  prediction with its uncertainty band and the value printed with that uncertainty, under the
  chapter's red/gold/green point. **No published CMS / ATLAS points in a channel chapter** (deck owner, 17 Sep 2026):
  they are revealed only at the end of the talk (`combination/result.md` has them with references). The repository documents no uncertainty on the aMC@NLO/FEWZ σ(60–120) reference;
  the band is the relative uncertainty of the published NNLO+NNLL NNPDF3.1 prediction of the same quantity,
  1940 +15 −21 pb (CMS-SMP-20-004, arXiv:2408.03744, Table 5), frozen with its source by the extractor
  (`extract_ztautau_reference.py` → `theory.unc_up/unc_down`) — reuse it, never type a band into a scene.

---

## C. Ideas (pick freely, record in the brief)

### C1. Colour ideas inside the palette
- **Chapter ramp**: one hue, several roles: `shade(ch, 0.30)` strokes and glyphs,
  `CHANNEL[ch]` fills, `tint(ch, 0.45)` secondary fills, `tint(ch, 0.80)` boxes and highlight
  areas. ee green is dark enough for text; ττ red too; μμ gold needs the shade.
- **Before / after**: the "before" state as a `GREY` ghost, the "after" in full colour.
- **Pass / fail**: pass in the channel colour, fail as `LIGHT_GREY` fill with an `INK` rim.
- **Measured vs simulated**: data `SLATE`, simulation in `THEORY` tints
  (`lighten(THEORY, 0.78)` files, as in the pipeline chain).
- **Energy in the ECAL (ee)**: a green ramp `tint(GREEN, 0.8 → 0.0)` per cell for the
  deposited energy; brem photons `PARTICLE["photon"]`.
- **Invisible things (ττ)**: neutrinos dashed `PARTICLE["nu"]`, `p_T^miss` a slate dashed
  arrow; "what we reconstruct" solid red vs "what was there" dashed `tint(RED, 0.55)`.
- **Dominant fakes (ττ: most of the SR)**: `SAMPLE["Fakes"]` is a pale slate on purpose,
  so a big fake background reads as "data-driven, not signal", not as a problem.
- **Uncertainty**: a translucent band in the colour of the thing it belongs to
  (`data_band(..., opacity=0.15)`); systematic vs statistical as two opacities.

### C2. Story ideas per channel (not decided; the team chooses)
- **ee**: bremsstrahlung and conversions, the ECAL cluster as the electron, the trigger
  (`HLT_Ele27_WPTight_Gsf`) and its turn-on, electron ID working points, the ττ→ee
  feed-down (LHE truth split) and W+jets, the energy scale; with no result yet, end on
  the post-selection plot or the fit setup.
- **ττ**: τ_h decay modes and DeepTau, why the visible mass sits low (the neutrinos) and
  how the MET-corrected / MMC mass fixes it, the fake factor (jet → τ_h is the dominant
  background), the BDT categories, TauPOG corrections, a larger uncertainty than μμ and
  why.

---

## D. Where chapter files go

| what | path |
|------|------|
| the interview result | `setup_and_reference/briefs/<chapter>.md` (from `briefs/TEMPLATE.md`) |
| scenes | `scenes/s<S>_<chapter>_story.py` (one file per chain; own module-level anchor dict) |
| frozen data | `data/<chapter>_*.json` + `data/extract_<chapter>_*.py` |
| storyboard rows | `03-sections-storyboard.md`, your section only |
| drafts | `work/<clip>/`, previews in `work/preview/<chapter>_preview_480p.mp4` |
| delivered clips | `clips/<S>_<section>/`, only through `tools/render.py -q h` |

Shared files (`style/bnd_style.py`, `style/palette.py`, `tools/`) serve every chapter:
add, don't change. A new primitive gets a new name; an existing signature or default
never changes, because delivered clips depend on it.
