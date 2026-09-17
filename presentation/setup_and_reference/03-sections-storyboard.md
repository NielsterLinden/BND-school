# 03 — Sections and candidate clips (DRAFT, contents still to be decided)

Six sections, ~30 minutes, ~5 minutes each. Everything below is a proposal of
*what could be animated*; the user decides what goes in. Each candidate is
one clip = one idea. Sources named so the fidelity checker can verify.

## 1 — Theory and relevance (`1_theory`)

17 Sep 2026: "why should we still care about the Z?", told through the **cross section** (not the mass), ~3 min talk,
~70 s of clips, built backwards from the final plot. Brief: `briefs/theory.md`; scene `scenes/s1_theory_story.py`;
numbers `data/theory_reference.json`. Clip number = play order (1-01 is not played in section 1: 1-02 ends on its frame).

| # | clip | what moves | source |
|---|------|-----------|--------|
| 1-01 | `drell_yan` | q q̄ → Z/γ* → ℓ⁺ℓ⁻ built left to right; the base of 3-01 / 4-01 / 5-01 | textbook |
| ✔ 1-02 | `sm_to_process` | the SM table builds; quarks, γ/Z and leptons turn into the Drell–Yan diagram; g_V, g_A at both vertices, g_V = T₃ − 2Q sin²θ_W | textbook |
| ✔ 1-03 | `factorisation` | chained: protons, x₁P / x₂P; σ = Σ∫ f f σ̂ (σ̂ purple = known, f cyan = tested); the α_s series LO 1970 … N³LO 2022; inset x₁x₂ = m_Z²/s, \|y\| < 2.4 lit: x ∈ [6·10⁻⁴; 0.08] | theory_reference `orders`, `x_range` |
| ✔ 1-04 | `prediction` | σ axis; aMC@NLO 1953.9 ⁺⁵⁶·¹₋₈₁.₇ pb, band grows scale → α_s → PDF; NNLO+NNLL by PDF set (CT18 1921, MSHT20 1935, NNPDF3.1 1940, NNPDF4.0 1970), 49 pb spread | combLieke result.json `prediction`; arXiv:2408.03744 Table 5 |
| ✔ 1-05 | `lumi_limited` | chained: CMS 2024 1952 ± 49 pb splits into stat 4 / syst 18 / lumi 45, lumi lit above the PDF spread | combLieke references.json |
| ✔ 1-06 | `z_everywhere` | H→ττ (σ𝓑 / 560, schematic shapes on a log axis) under the Z→ττ tail; small slice: Z→μμ + jet becomes Z→νν + jet | YR4 in theory_reference `higgs` |
| ✔ 1-07 | `lepton_universality` | opens on 1-01's frame; legs split into e / μ / τ; log mass axis ×207 ×17; (g_A, g_V) plane, sin²θ runs to 0.2315, the three leptons stay on one point (g_V^ℓ = −0.037); equal 𝓑 = 3.366 % bars; LEP + SLD Γ_μμ/Γ_ee, Γ_ττ/Γ_ee | PDG, Phys. Rept. 427 |
| ✔ 1-08 | `mass_window` | the generator-level m_ℓℓ spectrum of our aMC@NLO DY sample (log); 60–120 GeV window; σ = 1953.9 pb = 0.965 σ(m > 50) | z-mumu gensums.json `h_lhe_mll` |
| ✔ 1-09 | `final_frame` | chained: the window becomes the prediction line of the final plot's axis (1650–2400 pb), band, empty rows Z→ee, Z→μμ, Z→ττ \| Z→ℓℓ: section 6 fills them | combLieke summary plot |

Retired candidates (15 Sep): `z_lineshape` (mass, not σ), `pdf_x1x2` (→ 1-03), `why_z` (→ 1-05, 1-06).

## 2 — CMS and methods (`2_cms_methods`)

Re-cut of 17 Sep 2026 on the user's notes (plan `notes-to-pipeline-1-dynamic-karp.md`): the one-row pipeline was
too detailed. Section 2 is now **one map with two parallel rows**: theory on top (purple: Drell-Yan → Monte Carlo →
simulated detector → simulated NanoAOD), experiment below (the CMS slice → **Ella's trigger chain** L1 → HLT → RAW →
NanoAOD, `scenes/s2_trigger.py`, drawn at 0.6× as the row itself), the two files meet in a **selection** funnel (added
17 Sep late, v2 of every map clip; section 2 only passes it, the ee chapter zooms into it), generic corrections
(factors from data or simulation, applied to the prediction; the data never move) and a one-parameter fit.
Tag and probe, the fake factor, stacks and pulls are left to the channel chapters (μμ: tag and probe; ττ: fake
factor; ee: the fit in depth). **Clip number = play order.** Scenes: `scenes/s2_map.py` (`MapBuild MapTheory
TriggerInMap MapCorrections MapFit MapThree`), the trigger chain via its `enter`/`leave` hooks. Schematic, no result
numbers; the trigger numbers are Ella's (40 MHz, 100 kHz, 1 kHz, ~440 L1 seeds, O(µs), O(100 ms)).

| # | clip | what moves | source |
|---|------|-----------|--------|
| ✔ 2-01 | `cms_logo_to_slice` | **the detector build** (kept, v3): CMS logo appears, expands to the top right, the quarter layers round up into the full slice, the detail fades in | `CMSLogo`, `logo_rings` |
| ✔ 2-02 | `map_build` | chained on 2-01: the slice parks bottom left; the theory row draws in purple (diagram → MC box → slice in a dashed purple ring → purple file), then the experiment row (L1 → HLT → RAW → NanoAOD), then both files converge into corrections → fit→σ | all channel docs |
| ✔ 2-03 | `theory_prediction` | zoom onto the theory row (full-size row at the bottom, purple panel above): the lineshape dσ/dm_ℓℓ draws, its area = σ_theory | schematic |
| ✔ 2-04 | `theory_generator` | Monte Carlo hit or miss: the first kept point becomes an ℓ⁺ℓ⁻ event; ~90 points thrown, the kept ones drop into a histogram that follows the curve | schematic |
| ✔ 2-05 | `theory_simulation` | the event flies into the simulated slice: purple tracks, a purple file to NanoAOD; a stream of simulated events | z-mumu docs/10 |
| ✔ 2-06 | `trig_a_collisions` | zoom back to the map and onto the experiment row = the trigger spine; bunch crossings, 40 MHz | Ella's outline |
| ✔ 2-07 | `trig_b_l1_primitives` | L1 step 1: coarse calorimeter towers and muon segments light up | Ella's outline |
| ✔ 2-08 | `trig_c_l1_objects` | L1 step 2: candidate objects (μ, e/γ, jet) | Ella's outline |
| ✔ 2-09 | `trig_d_l1_decision` | L1 step 3: menu of ~440 seeds, some fire → accept; O(µs) | Ella's outline |
| ✔ 2-10 | `trig_e_l1_rate` | panel closes; 100 kHz to the HLT | Ella's outline |
| ✔ 2-11 | `trig_f_hlt_paths` | HLT paths of filters, one trigger bit each | Ella's outline |
| ✔ 2-12 | `trig_g_hlt_rate` | O(100 ms), 1 kHz leave the HLT | Ella's outline |
| ✔ 2-13 | `trig_h_raw` | the RAW event: detector data + L1 result + HLT bits + objects | Ella's outline |
| ✔ 2-14 | `trig_i_nanoaod` | the skim: detector data dropped, the rest squeezed into NanoAOD | Ella's outline, docs/10 |
| ✔ 2-15 | `corr_compare` | the panel folds into NanoAOD, zoom out, zoom into corrections: the purple and the grey file drop onto one m_ℓℓ axis: data points, the purple prediction (high, tilted) | schematic |
| ✔ 2-16 | `corr_factor` | dashed control sample: the same quantity in data (grey bar) and in simulation (purple bar); k₁ = grey/purple | docs/11, 12, 13 |
| ✔ 2-17 | `corr_apply` | × k₁, × k₂ (shape), × k₃ fly onto the prediction, which moves toward the data; the data never move | schematic |
| ✔ 2-18 | `fit_model` | zoom out, zoom into the fit: data points with error bars, the purple line μ σ_theory, the μ slider (purple reference at 1) | CONVENTIONS.md §2 |
| ✔ 2-19 | `fit_scan` | μ scans, the line scales through the points, −2Δln L(μ) draws a parabola, settles at μ̂; the Δ = 1 crossings give ±δμ on the slider | schematic |
| ✔ 2-20 | `fit_sigma` | σ = μ̂ σ_theory | CONVENTIONS.md §2 |
| ✔ 2-21 | `map_three` | zoom out to the finished map; it shrinks and repeats in ee green / μμ gold / ττ red with e⁺e⁻, μ⁺μ⁻, τ⁺τ⁻ | palette `CHANNEL` |

Superseded files: `clips/2_cms_methods/00_archive/pre_recut_2026-09-17/` (the one-row pipeline `pipe_a_map …
pipe_h_three`, 2-04…2-24, `scenes/s2_pipeline.py`, and the standalone trigger clips 2-25…2-33), registry rows in
`clips/CLIPLIST_2_cms_methods_pre_recut.tsv`; older still in `00_archive/` (`cms_slice_build`, `cms_signatures`).

## 3 — Z → ee (`3_zee`)

17 Sep 2026 (deck owner's storyline): the chapter is told **on the section-2 map**: it opens on 2-21's three maps, flies
into the **selection** node and later into the **fit** node. Real data throughout; result = the channel's own fit
(z-ee delivery 16 Sep, `FREEZE.md`). Brief `briefs/zee.md`; scene `scenes/s3_zee_story.py` (`ZeeSelection ZeeFit`);
numbers `data/zee_fit.json`, `data/zee_selection.json`. **Clip number = play order.**

| # | clip | what moves | source |
|---|------|-----------|--------|
| ✔ 3-01 | `ee_process` | Drell-Yan diagram in the ee flavour; each electron leg radiates and the photons convert: an EM shower in Feynman style, inside a faint cone | `shower_tree` |
| ✔ 3-02 | `ee_detector` | chained: dissolves to the slice; e⁻/e⁺ tracks bend opposite ways, stop in the ECAL as green clusters, one brem photon, nothing beyond; camera zooms in at the end | `signature("e")` |
| ✔ 3-03 | `ee_map_selection` | opens on 2-21's last frame: μμ and ττ maps fade, the ee map grows to full size, the camera flies into the selection funnel; the inner detector grows out of it | s2_map geometry |
| ✔ 3-04 | `ee_sel_event` | a real selected event (run 279024): two green tracks with deposits, e⁻ 45.5 GeV / e⁺ 41.4 GeV; cyan ring on the trigger electron → cut list "HLT: p_T^e > 27 GeV" ✓; the tracks glow, "e⁺e⁻" ✓ | zee_selection `event` |
| ✔ 3-05 | `ee_sel_pt` | the view parks; the two p_T values drop onto a p_T axis; the real spectrum grows (Z-window pairs green on the non-Z grey); m_Z/2 line; 0–20 GeV veiled, cut line → "p_T > 20 GeV" ✓ | zee_selection `pt` |
| ✔ 3-06 | `ee_sel_eta` | r-z side view (tracker to |η| = 2.5, ECAL barrel + endcaps): the two real electrons with hits; the η = 2.5 lines; a grey electron at η = 2.8 reaches the endcap without hits, "e / γ ?" → "|η| < 2.5" ✓ | zee_selection `cuts`, event η |
| ✔ 3-07 | `ee_sel_id` | η-φ calorimeter panels: a real electron (narrow cluster, one track, empty cone, no H) vs a real failing candidate (wide, tracks in the cone, H/E 0.58); σ_iηiη, H/E, I_rel with the real values; ✓ / ✗ | zee_selection `event`, `fake_candidate` |
| ✔ 3-08 | `ee_sel_wp` | real m_ee of e⁺e⁻ (green) and e^±e^± (grey) pairs steps no ID → veto → loose → medium → tight; purity vs efficiency builds point by point; medium lit → "ID: medium" ✓ | zee_selection `wp` |
| ✔ 3-09 | `ee_sel_mass` | "60 < m_ee < 120 GeV" ✓; the cut list folds into a big green funnel; dots rain from its spout into the m_ee histogram (30 × 2 GeV, log), N counts up to 6,320,097 | zee_fit `data` |
| ✔ 3-10 | `ee_sel_mc` | a purple simulation file drops into the funnel; the prefit stack grows under the data (bars → dots), key, ratio panel, 0.970 | zee_fit `prefit` |
| ✔ 3-11 | `ee_sel_counting` | N_data = σ·L·εA + N_bkg → σ = (6,320,097 − 44,357)/(0.202 × 16.39 fb⁻¹) = 1896 pb; then stack, data and ratio flatten to their averages: the shape is gone | zee_fit `counting` |
| ✔ 3-12 | `ee_map_fit` | everything shrinks back into the funnel node of the ee map; the camera flies into the fit node; the full-shape plot grows out of it | s2_map geometry |
| ✔ 3-13 | `ee_fit_mu` | ν_i = μ s_i + b_i and a μ_Z slider: μ → 1.15 → 0.85 → 1, only the green layer and the ratio move | zee_fit `prefit` |
| ✔ 3-14 | `ee_fit_nuisance` | ν_i = μ s_i(θ) + b_i(θ); pull rows; luminosity, electron ID, QCD scale, FSR each go +1 → −1 → 0 on their real templates, stack and ratio follow | zee_fit `templates` |
| ✔ 3-15 | `ee_fit_result` | L(μ, θ) = Π Pois · Π G; pile-up, L1 prefiring, electron reco, σ(Z→ττ) rows join; all parameters (and the MC-stat γ) move to the post-fit values, the stack lands on TRExFitter's post-fit yields, ratio flat, pull bars; μ_Z = 0.942 ± 0.015 | zee_fit `nps`, `postfit`, `mu` |
| ✔ 3-16 | `ee_sigma` | v2: the 5-33 result frame: σ axis 1750–2250 pb, 1954.1 ⁺¹⁵₋₂₁ pb band, σ_60–120 = 1841 ± 30 pb (green point); no CMS/ATLAS points (revealed only at the end of the talk) | zee_fit `sigma`, `theory` |

## 4 — Z → μμ (`4_zmumu`)

Final pass of 17 Sep 2026 on the frozen result (`FREEZE.md`, tag `zmumu-freeze-2026-09-17`); brief:
`briefs/zmumu.md`. **Clip number = play order** (the section was renumbered once, on the deck owner's request).

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 4-01 | `mumu_process` | Drell-Yan diagram in the μμ flavour; the muon legs simply continue, straight and clean | — |
| ✔ 4-02 | `mumu_into_detector` | chained: the slice fades in behind the diagram, the diagram shrinks into the interaction point, a flash; two outlined gold tracks through everything, MIP dots in both calorimeters, stubs in the four muon stations | `signature("mu")` |
| ✔ 4-03 | `mumu_event_to_mass` | chained (**real data from here on**, `data/zmumu_*.json`): schematic tracks fade, real SR event drawn (φ, charge, p_T), the p_T values converge to m_μμ = 89.7 GeV | data/zmumu_events.json |
| ✔ 4-04 | `mumu_first_entry` | chained: slice parks left, *linear* count axes appear, m_μμ becomes the first entry: a block, not a dot | data/zmumu_events.json |
| ✔ 4-05 | `mumu_fill` | chained: four more real events drop blocks (1,0,1,2,0,0,0,0,1: two share the 89–90 GeV bin); without a pause the axis becomes logarithmic, clock, blocks rain over the whole plot (bins sampled by their visible height), the 60 bins grow to the real counts; N = 10,378,567 | data/zmumu_events.json, zmumu_sr_stack.json |
| ✔ 4-06 | `mumu_prediction` | chained: bars → points; the **uncorrected** prediction slides in (simulation + the data-driven fakes wedge); colour key | data/zmumu_sr_stack.json (`recut` stage `raw`) |
| ✔ 4-07 | `mumu_ratio` | chained: ratio panel opens (0.88–1.12, ticks 0.9/1.0/1.1, data stat error bars — one range for every stage): data/pred = 0.944 | data/zmumu_sr_stack.json |
| ✔ 4-08 | `mumu_tnp_tag_probe` | chained: plot parks top right, the slice grows back; a real Z pair: the tag (solid gold, tight, isolated, cyan trigger ring), the probe (dashed grey); m_μμ counts up to 90.45 GeV and the probe turns gold: a real muon, nothing asked of it | data/zmumu_tnp.json, docs/17 |
| ✔ 4-09 | `mumu_tnp_pass_fail` | chained: the question, tight ID: the probe's muon-station stubs pulse, "3 muon stations ✓", it drops into `pass`; a second real pair, tracker-only probe, "1 muon station ✗", into `fail`; the pass panel is twice as high, own linear scales; all probes of the 40–45 GeV, \|η\| < 0.9 cell rain in; fitted peak on the grey background; ε_data = N_pass/(N_pass+N_fail) = 0.9585 | data/zmumu_tnp.json, docs/12, docs/17 |
| ✔ 4-10 | `mumu_tnp_grid` | chained: the two panels move up = row `data`; the identical procedure on `simulation` below (gold outlines): pass \| fail × data / simulation; ε_data 0.9585, ε_sim 0.9711 above the fail panels | data/zmumu_tnp.json `cells_barrel` |
| ✔ 4-11 | `mumu_tnp_scan` | chained: the grid parks left, ε_ID vs p_T opens; the p_T window steps through the ten barrel cells: the probe's curvature in the parked slice, the lit bin, the four real spectra (at 20–25 GeV the fail panel is mostly background) and both ε change; a data and a simulation point land per cell | data/zmumu_tnp.json |
| ✔ 4-12 | `mumu_tnp_scale_factor` | chained: panel below the ε plot: each data/simulation pair merges into SF = ε_data/ε_sim (0.9585/0.9711 = 0.987) | data/zmumu_tnp.json |
| ✔ 4-13 | `mumu_tnp_table` | chained: the SF points fly into the p_T column of a table; a side view of the detector (`eta_fan`) steps through \|η\| 0.9–1.2, 1.2–2.1, 2.1–2.4, the ε and SF points follow, the table fills to the real 10 × 4 map | data/zmumu_tnp.json |
| ✔ 4-14 | `mumu_tnp_apply` | chained: table beside the m_μμ plot: three real simulated events light their two cells, w = SF·SF (0.974, 0.961, 0.971), a reweighted block lands; the whole prediction steps 0.944 → 0.971; identification, isolation, trigger, reconstruction (same method): → 0.968; row "muon efficiency × 0.976" | data/zmumu_tnp.json `apply_examples`, zmumu_sr_stack.json `recut` |
| ✔ 4-15 | `mumu_corr_pileup` | chained: row "extra collisions (pile-up) × 0.994": → 0.974 | zmumu_sr_stack.json `recut`, docs/11 |
| ✔ 4-16 | `mumu_corr_prefiring` | chained: row "trigger fired early (prefiring) × 0.980": → 0.994 | docs/11 |
| ✔ 4-17 | `mumu_corrected` | chained: the rows collapse to "prediction × 0.950"; the final data vs prediction, 0.994 (= the frozen fit input, = the report plot `z-mumu/summary_deck/figures/sr_mass_log.png`) | zmumu_sr_stack.json `recut.final` |
| ✔ 4-18 | `mumu_rebin` | chained: rebin 60 → 12 (5 GeV), five adjacent bins summed as in `v2_5_fit.py` | data/zmumu_fit.json |
| ✔ 4-19 | `mumu_fit` | chained: pulls + μ slider; post-fit, ratio flat; μ_Z = 0.988 ± 0.014 | data/zmumu_fit.json |
| ✔ 4-20 | `mumu_sigma_fid` | chained: the plot parks left; σ_fid = (N − B)/(C·L) = (10,378,567 − 68,799)/(0.792 × 16.39 fb⁻¹) = 794.5 pb (counting); the fit: σ_fid = 790.2 ± 0.2 ± 6.1 ± 9.6 pb | handoff.md:11-13, 56-58 |
| ✔ 4-21 | `mumu_sigma_total` | chained: side view: pairs with a muon beyond \|η\| = 2.4 or too soft are lost, A = 0.409; σ(60–120) = σ_fid/A = 1931 ± 30 pb on its axis, "this analysis", beside the dashed prediction 1953.9 pb, "prediction, aMC@NLO, NNLO norm." (no CMS/ATLAS points: previous results only in the combination chapter) | handoff.md:13, 99-100, CONVENTIONS.md §6 |

Superseded files: `clips/4_zmumu/00_archive/pre_recut_2026-09-17/` (the numbering of 15–16 Sep: `mumu_a1…h4`,
`mumu_t1…t5`, `mumu_detector`), registry rows in `clips/CLIPLIST_4_zmumu_pre_recut.tsv`; older still in
`00_archive/` (the fake-factor block 4-10…4-16 and `mumu_g3_lepton_sf`, retired 16 Sep: Z→ττ tells the fake factor).
Why tag-and-probe first and in full: the largest single data/MC correction (prediction −2.4 %), measured by us
(20.05 M data pairs, 320 fits per efficiency). The correction order is a choice of the talk; the weights are the
frozen ones and the last stage is the fit input.

## 5 — Z → ττ (`5_ztautau`)

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 5-01 | `tau_decay` | τ⁻ → ν_τ W*⁻, W*⁻ → d ū, hadronisation blob, π⁻ π⁰ (decay mode 1); the only neutrino is the ν_τ | z-tautau/docs/02-selection.md |
| ✔ 5-02 | `tau_jet` | chained: the boosted τ: pions collimate into a narrow cone, π⁰ → γγ | — |
| ✔ 5-03 | `tautau_detector` | chained: dissolves to the slice; a 1-prong τ_h (track, γγ in the ECAL, HCAL) and a 3-prong τ_h, two ν_τ leave unseen, the p_T^miss arrow recoils; camera zooms in (centre (−0.2, −0.25): the chapter-identifier corner stays white) | `signature("tau_h", prongs=...)`, `"met"` |
| ✔ 5-04 | `tautau_a1_event` | chained on 5-03 (**real data from here on**, `data/ztautau_*.json`; brief `briefs/ztautau.md`): camera zooms out, the schematic τ_h / ν / MET fade, a real SR event (τ_h tracks from φ, charge, DM, p_T; the real p_T^miss arrow), p_T values | data/ztautau_events.json |
| ✔ 5-05 | `tautau_a2_visible_mass` | chained: the p_T values converge to m_vis, below the dashed m_Z | data/ztautau_events.json |
| ✔ 5-06 | `tautau_b1_neutrinos` | chained: dashed ν arrows along each τ_h; their vector sum lands on the p_T^miss arrow | docs/04-ditau-mass.md |
| ✔ 5-07 | `tautau_b2_likelihood` | chained: the (x₁, x₂) posterior grid of the event lights up; m_ττ = m_vis/√(x₁x₂), the median cell; m_ττ replaces m_vis | data/ztautau_events.json |
| ✔ 5-08 | `tautau_b3_shapes` | chained: the Z→ττ signal simulation only (inclusive sample, one peak): m_vis (grey ghost, 0.80 m_Z) vs m_ττ (red, 0.99 m_Z); 13 % → 11 % | data/ztautau_mass.json, docs/04 |
| ✔ 5-09 | `tautau_c1_first_entry` | chained: slice parks left, linear m_ττ axes (14 fit bins, 0–350 GeV), first entry | data/ztautau_sr_stack.json |
| ✔ 5-10 | `tautau_c2_rain` | chained: clock, rain, the 14 data bins grow; N = 21,160 | data/ztautau_sr_stack.json |
| ✔ 5-11 | `tautau_d1_simulation` | chained: bars → points; the simulation stack (rest / non-fid Z→ττ / fid Z→ττ) reaches a third of the data; key | data/ztautau_sr_stack.json |
| ✔ 5-12 | `tautau_d2_gap` | chained: ratio panel opens at data/pred = 2.97 | data/ztautau_sr_stack.json |
| ✔ 5-13 | `tautau_e1_same_sign` | chained: plot parks; the slice in the dashed τ_h^±τ_h^± box with a real same-sign pair (τ₁ a jet) | data/ztautau_events.json, docs/05 |
| ✔ 5-14 | `tautau_e2_tally` | chained: ten real same-sign pairs, 2 of 10 pass → f = N_T/N_L = 2/8 | data/ztautau_events.json |
| ✔ 5-15 | `tautau_e3_ff_map` | chained: the real 4 × 5 fake-factor map (period G, 0 jets), rows 1-prong / 1-prong+π⁰ / 3-prong / 3-prong+π⁰; f = f(period, DM, N_jets, p_T) written, no more tables | data/ztautau_fakes.json |
| ✔ 5-16 | `tautau_e4_eta_closure` | chained: first the same-sign data vs FF prediction in \|η(τ₁)\| with its ratio (0.93 … 1.14 … 0.84), then the ratio becomes f(\|η\|) and slides to 1; the same for p_T(τ₂) | data/ztautau_fakes.json, docs/05 |
| ✔ 5-17 | `tautau_e5_osss` | chained: C_OS/SS = C(period, N_jets, D_BDT) ∈ [0.96; 1.305] (no tables); the box turns opposite-sign | data/ztautau_fakes.json |
| ✔ 5-18 | `tautau_f1_apply` | chained: N_fake = C_OS/SS · f × N_AR; the plot returns | docs/05 |
| ✔ 5-19 | `tautau_f2_template` | chained: the fake template fills the gap; ratio 2.97 → 1.02; N_fake = 13592 (64 %) | data/ztautau_fakes.json |
| ✔ 5-20 | `tautau_g1_tau_sf` | chained: SF_ID per DM, TES per DM, ⟨SF_trig⟩; the signal layer breathes to the corrected height | data/ztautau_corrections.json, docs/07 |
| ✔ 5-21 | `tautau_g2_event_weights` | chained: ⟨w_PU⟩, ⟨w_L1⟩ = 0.990, the MC subtraction (2.5 % of the AR), the fiducial split (38 %) | data/ztautau_corrections.json |
| ✔ 5-22 | `tautau_h1_inputs` | own scene from a blank frame: the 16 BDT input symbols by importance (no mass variable) | data/ztautau_bdt_inputs.json |
| ✔ 5-23 | `tautau_h2_shapes` | chained: unit-normalised fakes (pale slate) vs Z→ττ (red) distributions of the leading inputs, one panel after the other | data/ztautau_bdt_inputs.json |
| ✔ 5-24 | `tautau_h3_sketch` | chained: schematic decision trees summing into D_BDT | docs/09-bdt.md |
| ✔ 5-25 | `tautau_h4_score` | chained: the real D_BDT output (log y from 10^1.5, ratio panel): fakes, Z→ττ, non-fid, rest, data; cuts at 0.55 / 0.90; AUC 0.966 | data/ztautau_bdt.json |
| ✔ 5-26 | `tautau_h5_categories` | chained: the plot splits into SR0 / SR1 / SR2 panels; SR0 < 110 GeV greyed; 486/12351, 1131/1044, 2141/197 | data/ztautau_bdt.json, data/ztautau_fit.json |
| ✔ 5-27 | `tautau_i1_fit` | chained: μ_Z slider + pulls labelled by name (τ_h ID (1-prong), C_OS/SS (D_BDT < 0.55), …); post-fit; μ_Z = 1.071 +0.114 −0.100 | data/ztautau_fit.json |
| ✔ 5-28 | `tautau_i2_impacts` | chained: grouped impacts as bars (Tau ID 8.4 % … stat 2.1 %; "Gammas" printed as "Template stat. (γ)") | data/ztautau_fit.json |
| ✔ 5-29 | `tautau_i3_sigma_fid` | chained: σ_fid = 4.82 ± 0.09 ± 0.47 pb | data/ztautau_fit.json |
| ✔ 5-30 | `tautau_i4_sigma_total` | chained: a large σ axis: σ(60–120) = 2082 +222 −194 pb (label above) vs the prediction 1944.9 +15 −21 pb with its band (NNLO+NNLL NNPDF3.1 uncertainty, CMS-SMP-20-004 Table 5) and the published CMS (1952 ± 49) and ATLAS (1981 ± 57, 66–116 GeV) points | data/ztautau_fit.json, data/ztautau_reference.json, combination/result.md |

## 6 — combination and conclusion (`6_combination`)

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 6-01 | `three_to_one` | schematic: three channel points with error bars slide into one slate combined point next to the purple theory line (no numbers yet) | skeleton |
| | `three_to_one_numbers` | the same with the frozen σ values and uncertainties from `data/` | fitting/, handoffs |
| | `correlation_matrix` | correlated nuisance parameters (Lumi, PDF, QCDScale …) light up across the three channels | fitting/CONVENTIONS.md |
| | `chapter_bars` | the six chapter colours as the anchor of the talk, all lit | palette `CHAPTER` |

## 7 — outro (`7_outro`)

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 7-01 | `mic_drop` | the CMS slice; a hand with a microphone comes in from the top, opens, the microphone drops into the centre and explodes (flash, shock rings, sparks, camera shake); the detector cracks into 12 wedges; NL / BE / DE rays (`FLAG`, outro only) come in five staggered bursts at random angles and fly off, the wedges tumble off with them; "Thank you! / Any questions?" centred on the clean slide (the one clip with words) | closing slide |
