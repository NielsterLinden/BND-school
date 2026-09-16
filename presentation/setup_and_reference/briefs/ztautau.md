# Chapter brief — Z→ττ (5_ztautau)

Status: approved (2026-09-16, deck owner — the story was given in the request; the remaining
choices below are the director's stated assumptions, to be corrected at the first critique of the
480p preview).
Written by `channel-animation-director` from the request of 16 Sep 2026; the contract for
`channel-story-builder`, `plot-recreator` and the reviewers. Anchors: `../06-chapter-anchors.md`.

The request, verbatim: *"The idea for the Ztautau analysis animation is to highlight the differences
with respect to the other channels, so the neutrinos and how we reconstruct MET as a proxy for them
and how we combine this information with the visible taus to get the corrected mass. Another
highlight are the fakes which are significant compared to the other channels, so add something on
this as well. And also the fact that we trained the BDT in order to better separate the fakes from
the signal. You can also then create a list of other fixes done like the eta reweighting to match
data and MC, etc..."*

## 1. Scope
- Section folder / number: `5_ztautau` (5); chained after the delivered 5-01 `tau_decay`,
  5-02 `tau_jet`, 5-03 `tautau_detector` (`scenes/s5_ztautau.py`, `ZtautauDetector`)
- Time budget: ≈ 4 min in the talk (assumed, as the μμ chapter) → 9 scenes cut into ≈ 24 one-idea
  clips of 5–12 s (the speaker talks over held frames)
- Audience: PhD students, mostly outside CMS
- **The one sentence**: two τ_h and two invisible neutrinos — the MET is folded into the visible
  mass by a per-event likelihood so the Z peak comes back at m_Z; two thirds of the selected events
  are jets faking τ_h, estimated from same-sign data with a fake factor and sorted away from the
  signal by a BDT; the fit gives σ(60–120) = 2082 +222 −194 pb, dominated by the external τ_h ID
  scale factors.

## 2. Story
- Spine nodes this chapter spends time on: **corrections (the mass, the ν/MET step)**,
  **backgrounds (fakes, 64 %)**, **selection-as-categories (the BDT)**, fit, σ. Files, trigger and
  the τ_h ID itself are left to section 2 and to 5-01…5-03.
- Methods explained in depth:
  1. the MET-corrected di-τ mass (`z-tautau/docs/04-ditau-mass.md`): collinear neutrinos, the MET
     as their proxy, the (x₁, x₂) likelihood scan, posterior median, m_vis scale 0.80 → m_ττ 0.99;
  2. the fake factor on the leading τ_h (`docs/05-fake-factors.md`): same-sign determination, the
     era × DM × N_jets × p_T table, the |η(τ₁)| and p_T(τ₂) closure corrections (the "η
     reweighting" of the request), C_OS/SS, the fake template entering the stack;
  3. the k-fold BDT categories (`docs/09-bdt.md`): mass-agnostic inputs, the score, three categories;
  4. the corrections applied to the simulation (`docs/07-corrections-and-systematics.md`) as a
     symbol + value strip (the "list of fixes"; the words go on the PowerPoint slide, §10).
- Left out on purpose: tag-and-probe (section 2), DeepTau internals, the skims, the trigger-rate
  study, the working-point comparison (Medium v2.1), the SigModel discussion, W+jets details.
- Tone: **real data from the first clip of the chain on** (as μμ); the ν/MET geometry and the BDT
  topology are schematic on top of real numbers.

## 3. Numbers on screen
Result version shown: **z-tautau v3, DeepTau Tight on both legs, MC-subtracted fake factor**
(`z-tautau/handoff.md` "Result (v3, nominal)", `output/results.json` `fit.mcsub`,
`fit/results/ztautau_fit_result.json`, 15–16 Sep 2026). Never the v2.1 Medium numbers that
`combination/` quotes.

| symbol | value as printed | source | frozen in |
|--------|------------------|--------|-----------|
| event | Run, Event of one real SR event (OS, both Tight, 1-prong + 3-prong if possible, m_ττ near m_Z) | ntuples_v1 data | `ztautau_events.json` `sr.chosen[0]` |
| p_T(τ₁), p_T(τ₂), p_T^miss | the event's values, 1 decimal, GeV | ntuple | `ztautau_events.json` |
| m_vis, m_ττ | the event's values, 1 decimal | ntuple `m_vis`, `m_tt` | `ztautau_events.json` |
| m_vis / m_ττ scale, resolution | 0.80 → 0.99; 13.1 % → 10.9 %; m_col defined 38 % | docs/04 table | `ztautau_mass.json` `performance` |
| N | 21,160 | docs/02 cutflow, results.json | `ztautau_sr_stack.json` `data.total` |
| data / simulation (no fakes) | 21160 / (20718 − 13592) = 2.97 | RESULTS.md prefit table | `ztautau_sr_stack.json` `totals` |
| data / pred (with fakes) | 1.021 | RESULTS.md prefit table | `ztautau_sr_stack.json` `totals` |
| N_fake | 13592 (64 %) | RESULTS.md, docs/02 | `ztautau_fakes.json` `yields` |
| SS tally | 10 real same-sign pairs (τ₂ Tight): 2 with τ₁ Tight, 8 with τ₁ VVVLoose-not-Tight → f = 2/8 | ntuple | `ztautau_events.json` `ss_ff.chosen` |
| FF map | era G, 0 jets: 4 DM × 5 p_T bins, 2 decimals (of 2 × 4 × 3 × 5 = 120) | fakefactors.json `mcsub.ff.ff` | `ztautau_fakes.json` `ff_map` |
| closure f(\|η(τ₁)\|) | 6 bins: same-sign obs/pred before 0.93, 1.00, 1.11, 1.14, 0.96, 0.84 (v3 Tight; docs/05 quotes the v2 Medium pass 0.92 … 0.93) → after 1.00 ± 0.02–0.03 (recomputed) | fakefactors.json `mcsub.ff.closure`, docs/05 "Closure corrections" | `ztautau_fakes.json` `closure.eta` |
| closure g(p_T(τ₂)) | 5 bins | idem | `ztautau_fakes.json` `closure.pt2` |
| C_OS/SS | per category: 0.96–1.31 (era × N_jets × category), inclusive 1.056 | docs/05 table, fitinputs meta, fakefactors.json | `ztautau_fakes.json` `osss` |
| BDT score | 20 bins: data, fakes, Z→ττ, rest; edges 0.55 / 0.90 | bdt.json `sr_score` | `ztautau_bdt.json` |
| category yields | SR0: 486 / 12351, SR1: 1131 / 1044, SR2: 2141 / 197 (Z→ττ fid / fakes), data 15742 / 2704 / 2714 | RESULTS.md | `ztautau_bdt.json` `categories` |
| AUC | 0.966 | RESULTS.md | `ztautau_bdt.json` |
| SF_ID per DM | 0.902, 0.892, 0.938, 0.812 (DM 0, 1, 10, 11) | docs/07 | `ztautau_corrections.json` |
| TES per DM | 0.993, 0.991, 1.001, 0.997 | docs/07 | idem |
| ⟨w_L1⟩ | 0.990 | docs/07 | idem |
| ⟨w_PU⟩, ⟨SF_trig⟩, ⟨SF_ID⟩ on the fiducial signal | computed by the extractor (3 decimals) | ntuples + `ztautau.analysis.weights` | idem |
| μ_Z | 1.071 +0.114 −0.100 | handoff | `ztautau_fit.json` `poi` |
| impacts | Tau ID 8.4 %, Fakes 6.5 %, Gammas 3.9 %, Tau trigger 3.0 %, Bkg norm 2.6 %, TES 1.5 %, Sig. mod. 1.3 %, MET 0.8 %, Lumi 0.7 %, PU 0.6 %, L1 0.2 %; stat 2.1 % | handoff table | `ztautau_fit.json` `grouped_impact` |
| σ_fid | 4.82 ± 0.09_stat ± 0.47_syst pb (pred 4.50) | handoff | `ztautau_fit.json` `sigma_fid` |
| σ(60–120) | 2082 +222 −194 pb vs 1944.9 pb | handoff, results.json | `ztautau_fit.json` `sigma_60_120` |
| prediction band | 1944.9 +15 −21 pb: the relative uncertainty of the published NNLO+NNLL NNPDF3.1 prediction 1940 +15 −21 pb (stat, PDF, α_s, scales); the repository documents none on the aMC@NLO/FEWZ reference | CMS-SMP-20-004, arXiv:2408.03744, Table 5 | `ztautau_reference.json` `theory` |
| published σ | CMS 1952 ± 49 pb (60–120), ATLAS 1981 ± 57 pb (66–116) | combination/result.md | `ztautau_reference.json` `cms`, `atlas` |

Format: μ_Z three decimals, asymmetric `^{+0.114}_{-0.100}`; σ(60–120) integers, asymmetric
`^{+222}_{-194}`; σ_fid two decimals with `_{\mathrm{stat}}`, `_{\mathrm{syst}}` (no lumi split:
the channel quotes stat/syst); the acceptance uncertainty (± 76 pb) is not printed.

## 4. Entry and exit
- Opens on: the final frame of 5-03 `tautau_detector` (`s5_ztautau.py` `ZtautauDetector`): the
  slice with a 1-prong τ_h at 40° (charge −1, π⁰), a 3-prong τ_h at 218° (charge +1), two dashed
  ν_τ lines at ±6° of them, the p_T^miss arrow on their bisector, **camera zoomed to 0.5 on
  `det.c`**, the outer labels already faded. The first clip is therefore a `MovingCameraScene`
  that starts with `self.camera.frame.scale(0.5).move_to(det.c)` and zooms back out.
- Event display: one real SR event from the data ntuples (chosen by plot-recreator: OS, both Tight,
  p_T > 40, a 1-prong and a 3-prong τ_h if such an event with m_ττ within ±5 GeV of m_Z exists,
  else any DM pair), drawn from φ, charge, DM, p_T; the real p_T^miss direction and size.
- Ends on: σ(60–120) as a point with asymmetric error bars beside the purple prediction line and its
  band, with the published CMS and ATLAS points below, on a **1700–2400 pb** axis (the μμ 1850–2050
  axis cannot hold +222/−194; the tick spacing 100 pb).
- 5-03 ends zoomed (0.5) on (−0.2, −0.25) rather than on the detector centre, so the outer ring stays
  out of the chapter-identifier corner (`s5_ztautau.ZOOM_503`, `A["zoom_503"]`); the zoomed frame still
  fills the title band (accepted before the corner rule; the first 1.4 s of a1 zoom out of it).

## 5. Clips

| scene class | clip (manim section) | the one idea | what moves | numbers used | source |
|-------------|----------------------|--------------|-----------|--------------|--------|
| `TautauEvent` (MovingCameraScene) | `tautau_a1_event` | a real event replaces the schematic one | camera zooms out; schematic tracks / ν / MET fade; the real τ_h tracks (DM → prongs, charge → bend, p_T → curvature), their deposits, the real p_T^miss arrow; p_T values beside each τ_h and the arrow; run/event stamp | event | `ztautau_events.json` |
| | `tautau_a2_visible_mass` | the visible mass alone sits low | the two p_T values converge to `m_vis = … GeV`; a dashed `m_Z` marker shows it is below | m_vis | idem |
| `TautauMass` | `tautau_b1_neutrinos` | the neutrinos are collinear with the τ_h, the MET is their sum | two dashed grey ν arrows grow along each τ_h axis; their vector sum is drawn head-to-tail and lands on the slate p_T^miss arrow | — | docs/04 |
| | `tautau_b2_likelihood` | scan x₁, x₂; weigh by the MET transfer function | to the right an (x₁, x₂) grid (12 × 12, the event's posterior frozen) lights up cell by cell; `m_{\tau\tau} = m_{vis}/\sqrt{x_1 x_2}`; the median cell flashes (SLATE); `m_{\tau\tau} = … GeV` replaces m_vis | posterior grid, m_ττ | `ztautau_events.json` `sr.chosen[0].posterior` |
| | `tautau_b3_shapes` | in simulation: m_vis peaks at 0.80 m_Z, m_ττ at 0.99 m_Z | the grid parks; an m axis 0–200 GeV: **only the Z→ττ signal simulation** (inclusive aMC@NLO sample, one clean peak — the stitched jet-binned shapes had a shoulder), the m_vis shape as a GREY ghost step, the m_ττ shape in red on top, dashed m_Z; `0.80 \to 0.99`, `13\% \to 11\%` | shapes, performance | `ztautau_mass.json` |
| `TautauRain` | `tautau_c1_first_entry` | one event → one entry | the shapes fade; the slice parks left (μμ default); the m_ττ axes (0–350 GeV, 14 fit bins, linear y) appear; the event's m_ττ is the first entry | — | |
| | `tautau_c2_rain` | 21,160 events | clock, rain, the 14 data bins grow to the real counts; `N = 21,160` | N, data counts | `ztautau_sr_stack.json` |
| `TautauStack` | `tautau_d1_simulation` | the simulation explains a third | bars → points; the simulation stack slides in bottom-up (rest / Z→ττ non-fid / Z→ττ fid); key; the stack reaches ~⅓ of the data | mc counts | idem |
| | `tautau_d2_gap` | the missing two thirds | the ratio panel opens at data/pred = 2.97 (y range 0–4) | ratio | idem |
| `TautauFakes` | `tautau_e1_same_sign` | the control region: same-sign pairs are ~99 % jets | plot parks top right; the slice inside a dashed `\tau_h^{\pm}\tau_h^{\pm}` box (DETECTOR_ACCENT); a real SS pair: τ₂ tight (red), τ₁ a jet spray (grey) with a τ_h-like core | pair | `ztautau_events.json` `ss_ff.chosen[0]` |
| | `tautau_e2_tally` | 2 in 10 pass | ten real SS pairs flash through the box, a 10-box tally fills (pass red, fail pale slate); `f = N_T / N_L = 2/8` | 10 pairs | idem |
| | `tautau_e3_ff_map` | the fake factor depends on the decay mode, p_T (and period, N_jets) | the real 4 × 5 map (period G, 0 jets) fades in, rows named **1-prong, 1-prong+π⁰, 3-prong, 3-prong+π⁰** (never DM codes); instead of more tables the dependence is written schematically: `f = f(\mathrm{period},\ \mathrm{DM},\ N_{\mathrm{jets}},\ p_T)` | ff_map | `ztautau_fakes.json` |
| | `tautau_e4_eta_closure` | the leftover η dependence is corrected (the "η reweighting") | **first the disagreement**: a small same-sign plot vs \|η(τ₁)\| — prediction (FF × loose data, pale slate fill) with the data points and a ratio strip that shows 0.93 … 1.14 … 0.84; then the ratio becomes the correction `f(\|\eta_{\tau_1}\|)` and the bars slide to 1; the same for p_T(τ₂) (`g(p_T^{\tau_2})`) | closure (obs, pred, before, after) | idem |
| | `tautau_e5_osss` | measured in SS, applied in OS | no tables: `C_{\mathrm{OS/SS}} = C(\mathrm{period},\ N_{\mathrm{jets}},\ D_{\mathrm{BDT}}) \in [0.96;\ 1.305]`; the box symbol turns `\tau_h^{+}\tau_h^{-}` | osss range | idem |
| `TautauTransfer` | `tautau_f1_apply` | the recipe | `N_{fake} = C_{OS/SS}\, f \times N_{AR}`; the plot returns to the main position | — | |
| | `tautau_f2_template` | the fakes fill the gap | the fake template (pale slate) slides into the bottom of the stack; the ratio drops from 2.97 to 1.02; key gains Fakes; `N_{fake} = 13592` and `64\%` | N_fake, ratio | `ztautau_fakes.json`, `ztautau_sr_stack.json` |
| `TautauCorrections` | `tautau_g1_tau_sf` | the τ_h scale factors | left column, rows per decay mode named 1-prong / 1-prong+π⁰ / 3-prong / 3-prong+π⁰: `SF_{ID} = 0.902, 0.892, 0.938, 0.812`, `TES = 0.993 …`, `\langle SF_{trig}\rangle = …`; the Z→ττ layer of the stack breathes to the corrected height (if the extractor freezes the raw vs corrected signal yield; else the numbers only) | corrections | `ztautau_corrections.json` |
| | `tautau_g2_event_weights` | the event weights | `\langle w_{PU}\rangle`, `\langle w_{L1}\rangle = 0.990`; the MC subtraction of genuine τ from the FF regions (`2.5\%` of the AR) and the fiducial split `38\%` appear as two more rows | idem | idem |
| `TautauBDT` (**its own scene, opens on a blank frame**: everything of g2 fades out first) | `tautau_h1_inputs` | the classifier's inputs | the 16 input symbols appear as a list (by importance: ΔR, p_T(τ₁), p_T(ττ), Δφ, …), the mass variables are visibly absent | features, importance | `ztautau_bdt_inputs.json` |
| | `tautau_h2_shapes` | fakes and signal look different | for the leading inputs (ΔR, Δφ, p_T(τ₁), p_T(ττ), p_T^miss, Δφ(p_T^miss, ττ), N_jets …) a small panel each: the **unit-normalised** distributions of the dominant background (fakes = AR data × FF, pale slate) and the signal (Z→ττ, red), shown one after the other | shapes | idem |
| | `tautau_h3_sketch` | a BDT | schematic: a few decision trees (nodes/branches in DETECTOR_ACCENT) fed by the inputs, summing into one output `D_{\mathrm{BDT}}` | — | docs/09 |
| | `tautau_h4_score` | the output D_BDT | axes `D_{\mathrm{BDT}}` 0–1, **log y with a lower limit high enough that the shapes fill the panel** (start at 10^1.5 ≈ 32, not 10⁰: the smallest bottom-layer bin, the fakes at D_BDT > 0.95, is 55): fakes, Z→ττ, non-fid, rest, data points, **a data/pred ratio panel**; dashed cuts at 0.55 and 0.90; `AUC = 0.966` | sr_score (+ pred) | `ztautau_bdt.json` |
| | `tautau_h5_categories` | three categories, one fit variable | the m_ττ plot splits into three narrow panels SR0 / SR1 / SR2 with their stacks; SR0 bins below 110 GeV greyed (dropped); yields `486/12351`, `1131/1044`, `2141/197` | categories | idem, `ztautau_fit.json` `regions.prefit` |
| `TautauFit` | `tautau_i1_fit` | the profile-likelihood fit | the μ_Z slider (0.8–1.3, ref 1) at the top and below it the pulls of 8 NPs, **labelled by name, not by fit code**: τ_h ID (1-prong / 1-prong+π⁰ / 3-prong / 3-prong+π⁰) = TauID_DM0/1/10/11, C_OS/SS (D_BDT < 0.55) = FakeOSSS_tautau_c0, f closure (D_BDT > 0.90, m_ττ > 110) = FakeClosure_tautau_c2_hi, τ_h trigger (1-prong+π⁰) = TauTrigger_DM1, Luminosity = Lumi; post-fit: pulls move, SR panels go post-fit; `\mu_Z = 1.071^{+0.114}_{-0.100}` | poi, nps, postfit | `ztautau_fit.json` |
| | `tautau_i2_impacts` | why the error is 11 %: external τ_h ID SFs | horizontal bars of the grouped impacts (Tau ID 8.4 … L1 0.2, stat 2.1) growing left to right; group names as the fit spells them **except "Gammas", printed as "Template stat. (γ)"** (the per-bin MC + fake-template statistics) | grouped_impact | idem |
| | `tautau_i3_sigma_fid` | σ_fid | `\sigma_{fid} = 4.82 \pm 0.09_{stat} \pm 0.47_{syst}\ pb` | sigma_fid | idem |
| | `tautau_i4_sigma_total` | the result beside the prediction and the published measurements | a **larger** σ axis (≈ 1700–2400 pb, the label `\sigma_{60-120} = 2082^{+222}_{-194}\ pb` moved **higher**, above the plot): the prediction as a purple dashed line **with its uncertainty band**, labelled `1944.9^{+15}_{-21}` beside the top of the line (`ztautau_reference.json` `theory`, source in §3); the published **CMS** (13 TeV, 60–120 GeV, 1952 ± 49 pb) and **ATLAS** (66–116 GeV, 1981 ± 57 pb) points in slate below the red point, labelled `CMS` / `ATLAS`; this result's red point with its +222/−194 bar | sigma_60_120, reference | idem, `ztautau_reference.json` |

Preview order: 5-03 (rendered `-q l --no-deliver`) → a1 … i4.

## 6. Look
- Layout: the μμ defaults (06 §B3): slice parked at 0.42, (−4.6, −0.7); main plot centred on
  (1.85, −0.25); ratio panel; one-row key with left end (−1.25, −2.92); control-region box with the
  slice at 0.84.
- Motifs reused (06 §B4): event → entry → rain; stack + data + ratio; the data never move; control
  region box → template slides into the stack; tally; value grid; slider; pull plot.
- Colour ideas (06 §C1): chapter ramp in red (`CHANNEL["tautau"]` fills, `CHANNEL_LINE["tautau"]`
  strokes/glyphs, `tint(RED, 0.55)` for the non-fiducial Z/γ*→ττ layer, `tint(RED, 0.80)` boxes);
  invisible things dashed `PARTICLE["nu"]`, p_T^miss the slate dashed arrow; "what was there" as a
  GREY ghost (m_vis shape) vs "what we reconstruct" in red; dominant fakes in `SAMPLE["Fakes"]`
  (pale slate) on purpose; flashes in SLATE (red is the channel).
- Palette requests: **none new** — W+jets, tt̄, tW, VV, Z→ee and the low-mass DY are drawn as one
  "rest" layer in `SAMPLE["TTbar"]` (slate tint) with the striped key swatch the μμ chapter uses;
  the non-fiducial Z/γ*→ττ uses `tint(CHANNEL["tautau"], 0.55)` (a tint of the chapter colour, no
  new entry needed).
- Pace: brisk (5–12 s per clip); the speaker talks over the held frames.

## 7. Deviations from the defaults (06 §B), with the reason
- **Linear y axis** on the m_ττ plots (0–350 GeV, the 14 variable-width fit bins, counts per bin):
  the story is "two thirds of the events are missing from the simulation", which only a linear axis
  shows; the fakes are not a thin wedge here but the bulk. Ratio panel range 0–4 in d, 0.5–1.5 after f.
- **σ axis 1700–2400 pb** (μμ: 1850–2050): the asymmetric ±10 % error does not fit.
- The chain opens on a zoomed camera (5-03 ends with `zoom_in(…, 0.5)`); `TautauEvent` is a
  `MovingCameraScene` and its first move is the zoom-out. The seam is checked at delivery against a
  `-q h --no-deliver` render of 5-03.
- The chapter has **three plots** (m_ττ stack, BDT score, three category panels) instead of one; the
  score plot and the panels replace the main plot at the main position rather than sharing the frame.

## 8. Process
- Critique by: the deck owner (Samuel), on the 480p preview; reviews at delivery: physics / style /
  continuity (yes, once, at delivery).
- Scene file: `scenes/s5_ztautau_story.py`; frozen data `data/ztautau_{events,mass,sr_stack,fakes,bdt,corrections,fit}.json`
  with extractors `data/extract_ztautau_*.py` (LCG env). Untracked channel outputs are read through the
  git-ignored links `z-tautau/output/data` and `z-tautau/fit/results/ztautau/Plots` (→ the main checkout).
- Git: drafts were built on the worktree branch `worktree-ztautau-chapter`; on 16 Sep 2026 the deck
  owner had the chapter copied into `main`, rendered and pushed (clips 5-01…5-03 v3 and the story chain).

## 8b. Critique round 1 (deck owner, 2026-09-16) — applied
1. b3 shows only the Z→ττ signal simulation, one clean peak (inclusive sample, not the stitched jet-binned mix).
2. Decay modes are named (1-prong, 1-prong+π⁰, 3-prong, 3-prong+π⁰), never DM0/DM1/DM10/DM11 — everywhere (FF map, SF rows).
3. e4 first shows the data-vs-prediction disagreement in |η(τ₁)| and p_T(τ₂) (same-sign plot + ratio), then the correction.
4. The dependence of the fake factor and of C on (period, N_jets, D_BDT, …) is written as a schematic function, not as tables.
5. Ranges as `C_{\mathrm{OS/SS}} \in [0.96;\ 1.305]`.
6. The BDT is its own scene starting from a blank frame: input list → normalised fakes-vs-signal shapes per input → BDT sketch →
   `D_{\mathrm{BDT}}` output with a higher log floor and a ratio panel → the categories.
7. Impact group "Gammas" is printed as "Template stat. (γ)"; the σ axis is bigger, the σ label sits higher, and the CMS and ATLAS
   published values and the prediction with its uncertainty band are drawn for comparison.
The general form of these rules is in `../06-chapter-anchors.md` §B5 so every chapter inherits them.

## 8c. Critique round 1, follow-up (deck owner, 2026-09-16) — applied
1. **Chapter-identifier corner**: the slide deck carries a chapter identifier in its top-left 3 cm × 9 cm; like the
   title band nothing is drawn at x < −5.85, y > 0.22 (fixed anchor, `06` §A1, checked with `tools/zonecheck.py`).
   Moved for it: the τ decay diagram of 5-01 (+0.6 in x), the 5-03 zoom centre, the corrections column (g, now
   below y = 0.22), the parked BDT input list (x = −5.75, the shape panels +0.7), N_fake / the recipe (f, below the
   corner), the impact bars (i2), the μ_Z slider (top) and the pulls (below) of i1, the post-fit panels of i3/i4.
2. The prediction is drawn **with** its uncertainty (point 7 above): the band and `1944.9^{+15}_{-21}` (§3 row
   "prediction band"); the σ panel is taller (2.6) with the result label above it.
3. Pull labels are named, not coded (the DM rule of point 2 applied to the fit).
4. `bdt_sketch` used `Circle` without importing it: the h3 clip could not render before this round.

## 9. Open questions (for the first critique)
1. Time budget: ≈ 4 min assumed; cut clips if the chapter must be shorter (candidates: e5, g2, i2).
2. The real event: plot-recreator picks one; say if you prefer another topology (e.g. two 1-prong).
3. The corrections clip animates the signal layer if the extractor can freeze raw vs corrected yields
   in the LCG env; otherwise it prints the values only.
4. Should the v2.1 Medium result appear anywhere (no, by this brief)?

## 10. For the PowerPoint slide beside the corrections clip (words never go into a clip)
The fixes and corrections of the ττ chain (v2 → v3, `REVIEW.md`, `docs/05`, `docs/07`):
- fake-factor closure corrections in |η(τ₁)| (±15 % → ±5 %) and p_T(τ₂) (−7 % → flat): the "η reweighting" of the same-sign data vs prediction;
- fake factor binned per era (G/H trigger paths) and per jet multiplicity (0/1/≥2), not only DM × p_T;
- genuine-τ (MC) subtraction in the fake-factor determination and application regions (6 % of the signal was double counted);
- OS/SS correction C measured per BDT category (the charge correlation depends on the topology);
- W+jets subtracted with uniform weights (one event with weight −63 faked a 75 % non-closure);
- the signal split into fiducial (60 < m < 120 GeV, both vis p_T > 40, |η| < 2.1) and non-fiducial (38 %, a background with its own normalisation);
- SR0 (fake-dominated category) fitted above 110 GeV only (the fake sideband);
- TauPOG UL2016 corrections: τ_h ID SF per decay mode, τ_h energy scale, di-τ trigger leg SFs (genuine legs only), e/μ→τ_h SFs, pileup (69.2 mb), L1 prefiring;
- DeepTau Tight instead of Medium on both legs (fakes 80 % → 64 %, total uncertainty 14.5/12.5 % → 11.4/10.0 %);
- the LO-vs-NLO "SigModel" removed as a nuisance parameter (reported only).
