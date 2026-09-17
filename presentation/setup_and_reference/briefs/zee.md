# Chapter brief — Z→ee (3_zee)

Status: first delivery (2026-09-17, built on the deck owner's storyline of the same day; critique pending).
Anchors: `../06-chapter-anchors.md`. Scene: `scenes/s3_zee_story.py`.

## 1. Scope
- Section folder / number: `3_zee` (3); clips 3-03 … 3-16 after the kept openers 3-01 `ee_process`, 3-02 `ee_detector`.
  **Clip number = play order.**
- Time budget: ≈ 5 min in the talk, 14 clips, ≈ 86 s of animation; one idea per clip, the speaker talks over held frames.
- Audience: PhD students, mostly outside CMS.
- **The one sentence**: we keep the events that look like a Z → ee (trigger, two opposite-charge electrons, p_T, |η|,
  electron ID); counting them would give a cross section but throws the shape away, so we fit templates whose
  nuisance parameters move the prediction, and get σ(60–120) = 1841 ± 30 pb.

## 2. Story
- The chapter is told **on the section-2 map**. It opens on the three maps 2-21 ends on; the ee map grows, the camera
  flies into the **selection** node (added to the map on 17 Sep for this chapter), back out, and into the **fit** node.
  The corrections node (pile-up, prefiring, electron reco/ID scale factors) is passed: they are already in the simulation.
- Selection, in the deck owner's order: trigger + two electrons with opposite charge → p_T > 20 GeV → |η| < 2.5 →
  identification and isolation working points (purity vs efficiency, medium) → the selected pairs fill m_ee → the
  simulation predicts signal and backgrounds → counting N_data = σ·L·εA + N_bkg, which throws the shape away.
- Why the cuts, physically: p_T — Z electrons carry about m_Z/2 (the Jacobian peak of the real spectrum), soft
  electrons mostly come from non-Z pairs, the trigger needs 27 GeV; |η| — the tracker ends at |η| = 2.5, beyond it the
  ECAL endcap still sees the shower but there is no track, so an electron cannot be told from a photon.
- Fit: a multi-parameter template fit on real inputs: μ_Z scales Z→ee only; each nuisance parameter moves the prediction
  by its real ±1σ template (luminosity, electron ID, QCD scale, FSR), the ratio panel shows it; the likelihood; all
  parameters go to their post-fit values (TRExFitter's post-fit yields), the pulls, μ_Z; σ = μ_Z σ_theory.
- Left out: the ECAL-gap / electron-ID artefact and the missing trigger scale factor (combLieke/README.md "The ee
  channel") are not shown; the ElectronID pull (+0.59 ± 0.08) is on screen as fitted.
- Tone: real data from the event on; the data never move.

## 3. Numbers on screen
Result version: **z-ee delivery of 16 Sep 2026 14:44** (`z-ee/Zee_fit.tar.gz`, unchanged by the freeze, `FREEZE.md`
"Z → ee"): the channel's own TRExFitter fit (deck owner's choice, 17 Sep: each chapter shows its own fit; the
combination's ee line 2094 ⁺¹²²₋₁₁₄ pb is for section 6).

| symbol | printed | source | frozen in |
|--------|---------|--------|-----------|
| event | Run 279024 Event 104459341; e⁻ p_T = 45.5 GeV (η 0.84), e⁺ 41.4 GeV (η 0.05), m_ee = 91.4 GeV | SingleElectron NanoAOD, z-ee selection | `zee_selection.json` `event` |
| trigger | HLT: p_T^e > 27 GeV | z-ee/handoff.md (`HLT_Ele27_WPTight_Gsf`) | `zee_selection.json` `cuts` |
| p_T spectrum | shape only, peak 40–45 GeV, cut at 20 | 10-file data sample | `zee_selection.json` `pt` |
| |η| | 2.5 (tracker); ECAL barrel 1.479, endcap 3.0 (geometry only) | CMS detector | `zee_selection.json` `cuts` |
| ID candidates | e: σ_iηiη 0.0087, H/E 0.00, I_rel 0.00; jet: 0.0152, 0.58, 2.24 | data sample | `zee_selection.json` `event`, `fake_candidate` |
| working points | m_ee of e⁺e⁻ / e^±e^± pairs, none/veto/loose/medium/tight; efficiency 1, 0.83, 0.77, 0.64, 0.50; purity = 1 − N_SS/N_OS 0.90 … 0.99 | data sample | `zee_selection.json` `wp` |
| N | 6,320,097 (60–120 GeV) | Tables/Yields.txt | `zee_fit.json` `data` |
| data/pred. | 0.970 (prefit) | Plots/ee_SR prefit yields | `zee_fit.json` `prefit` |
| counting | (6,320,097 − 44,357)/(0.202 × 16.39 fb⁻¹) = 1896 pb | prefit yields, reference 1954.1 pb | `zee_fit.json` `counting` |
| templates | ±1σ per nuisance parameter (luminosity ±1.2 %, electron ID ±5.8 %, QCD scale, FSR) | Histograms/Zee_fit_histos.root | `zee_fit.json` `templates` |
| pulls | luminosity −0.44 ± 1.01, electron ID +0.59 ± 0.08, QCD scale −1.98 ± 0.39, FSR −0.58 ± 0.50, pile-up −1.58 ± 0.60, L1 prefiring +2.02 ± 0.92, electron reco +0.45 ± 0.98, σ(Z→ττ) −1.81 ± 0.93 | Fits/Zee_fit.txt | `zee_fit.json` `nps` |
| μ_Z | 0.942 ± 0.015 | Fits/Zee_fit.txt, FREEZE.md | `zee_fit.json` `mu` |
| σ(60–120) | 1841 ± 30 pb vs 1954.1 ⁺¹⁵₋₂₁ pb | FREEZE.md; band as in 5-33 | `zee_fit.json` `sigma`, `theory` (`published` is frozen but not drawn) |

## 4. Entry and exit
- Opens on: the last frame of 2-21 `map_three` (v2), not on 3-02 (a cut; 3-01 / 3-02 stay the chapter's openers).
- Ends on: the result frame of the ττ chapter (5-33): σ axis 1750–2250 pb, value above, prediction with band, the ee point in green
  (row of ττ's combined point). No CMS / ATLAS points: the deck owner reveals them only at the end of the talk (17 Sep, 3-16 v2).

## 5. Clips (number = play order)
| # | clip | the one idea |
|---|------|--------------|
| 3-03 | `ee_map_selection` | the ee map grows out of the three maps; into the selection node |
| 3-04 | `ee_sel_event` | a real event: the leading electron fired the trigger; two electrons, opposite charge |
| 3-05 | `ee_sel_pt` | the real electron p_T: the Z electrons sit near m_Z/2; p_T > 20 GeV |
| 3-06 | `ee_sel_eta` | side view: no track beyond |η| = 2.5, an electron there looks like a photon |
| 3-07 | `ee_sel_id` | electron vs jet: shower width, H/E, isolation (real values) |
| 3-08 | `ee_sel_wp` | working points on real m_ee: same-sign pairs drop faster than the peak; purity vs efficiency; medium |
| 3-09 | `ee_sel_mass` | the cut list folds into the funnel; the selected pairs fill m_ee, N = 6,320,097 |
| 3-10 | `ee_sel_mc` | the simulation through the same selection: stack, ratio 0.970 |
| 3-11 | `ee_sel_counting` | N_data = σ·L·εA + N_bkg = … 1896 pb; the histogram flattens: the shape is thrown away |
| 3-12 | `ee_map_fit` | back to the ee map, into the fit node: the shape is back |
| 3-13 | `ee_fit_mu` | ν_i = μ s_i + b_i: μ_Z scales Z→ee only |
| 3-14 | `ee_fit_nuisance` | one nuisance parameter at a time on its real template, stack and ratio move |
| 3-15 | `ee_fit_result` | the likelihood; every parameter to its post-fit value, ratio flat, pulls, μ_Z = 0.942 ± 0.015 |
| 3-16 | `ee_sigma` | σ = 1841 ± 30 pb beside the prediction (v2: no CMS / ATLAS) |

## 6. Look
- Layout: the event view is the inner detector (tracker, ECAL, HCAL) of a 2.4× slice; a cut list on the right with
  cyan ticks through the selection clips; the plot on the left in the selection block, on the right in the fit block
  (μμ defaults) with μ slider and pull rows on the left (below y = 0.22 where they are long).
- Motifs: map zooms of section 2, funnel → histogram rain, stack + data + ratio, slider, pull plot, ττ result frame.
- Colours: ee green for electrons, signal, the funnel; method cyan for ticks, trigger ring, cut lines, the η = 2.5
  lines; backgrounds from `SAMPLE` (Z→ττ red as in μμ, tt̄ slate, VV = `SAMPLE["WW"]`, W+jets = `SAMPLE["Fakes"]`: the
  palette has no W+jets entry, 06 §A5); same-sign pairs grey.
- Words on screen (as in μμ): the working-point names, "efficiency", "purity", nuisance-parameter names.

## 7. Deviations
- Opens on 2-21, not on 3-02 (the map is the chapter's spine).
- W+jets and VV colours borrowed from existing `SAMPLE` entries (no palette change).
- The selection spectra are a 10-file sample of the data (shapes and fractions); totals on screen are the full dataset.

## 8. Process
- One 480p pass, fixes of real failures only, one delivery render; critique by the deck owner on the delivered clips.

## 9. Open questions for the deck owner
- The counting number (1896 pb) sits 3 % above the fit (1841 pb): the fit pulls the electron-ID normalisation
  (+0.59σ of ±5.8 %). Keep the number in 3-11, or show the formula only?
- An uncertainty breakdown clip (luminosity 1.15 % of the 1.53 % total) is not built.
