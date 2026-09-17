# 03 — Sections and candidate clips (DRAFT, contents still to be decided)

Six sections, ~30 minutes, ~5 minutes each. Everything below is a proposal of
*what could be animated*; the user decides what goes in. Each candidate is
one clip = one idea. Sources named so the fidelity checker can verify.

## 1 — Theory and relevance (`1_theory`)

| # | candidate clip | what moves | source |
|---|----------------|-----------|--------|
| 1-01 | `drell_yan` (done, v1) | q q̄ → Z/γ* → ℓ⁺ℓ⁻ built left to right; ee/μμ/ττ colour variants in the same file | textbook |
| | `z_lineshape` | Breit-Wigner at m_Z, Γ_Z; then the γ* tail and the radiative (FSR) left tail deform it | PDG values |
| | `lepton_universality` | three equal-width bars (BR to ee, μμ, ττ) growing to the same height; the σ×BR per flavour = 6077.22/3 pb | CONVENTIONS.md |
| | `pdf_x1x2` | two protons, a quark and antiquark pulled out with momentum fractions x₁, x₂, ŝ = x₁x₂s | textbook |
| | `why_z` | Z as a standard candle: luminosity / detector calibration / PDF constraint (three icons, animation only) | — |

## 2 — CMS and methods (`2_cms_methods`)

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| archived (2-02) | `cms_slice_build` | inside-out build of the slice, superseded by the logo build (`00_archive/`) | `CMSSlice` |
| archived (2-03) | `cms_signatures` | e, μ, τ_h, jet, γ signature legend on the slice, not used in the talk (`00_archive/`) | `signature()` |
| ✔ 2-01 | `cms_logo_to_slice` | **the detector build**: CMS logo appears, expands to the top right, the quarter layers round up into the full slice, the detail fades in | `CMSLogo`, `logo_rings` |
| ✔ 2-04 | `pipe_a_map` | **the analysis spine**: nodes appear left→right — detector → recorded & simulated files → selection → tag & probe → backgrounds → comparison → fit → σ | all channel docs |
| ✔ 2-05 | `pipe_b1_collisions` | chained: zoom on detector + files; collisions flash, tracks, files fly to a pile | z-mumu docs/10-skims.md |
| ✔ 2-06 | `pipe_b2_skim` | chained: the pile of 8 squeezes to 2 (the skim) | docs/10 |
| ✔ 2-07 | `pipe_b3_simulation` | chained: mini Drell-Yan → the same slice → purple simulated files join the pile | docs/10 |
| ✔ 2-08 | `pipe_c1_trigger` | chained: zoom out, zoom into the funnel; trigger bar | docs/01 |
| ✔ 2-09 | `pipe_c2_lepton_id` | chained: ℓ p_T, ID/iso cut | docs/01 |
| ✔ 2-10 | `pipe_c3_two_leptons` | chained: N_ℓ = 2, then ℓ⁺ℓ⁻ | docs/01 |
| ✔ 2-11 | `pipe_c4_mass_window` | chained: m_ℓℓ window: the surviving bar | docs/01, docs/14 |
| ✔ 2-12 | `pipe_d1_tag_probe` | chained: zoom out, zoom into tag & probe; tag bold, probe dashed | docs/12 |
| ✔ 2-13 | `pipe_d2_efficiency` | chained: probes fall into pass / fail; ε = N_pass/(N_pass+N_fail) | docs/12 |
| ✔ 2-14 | `pipe_d3_scale_factor` | chained: the same in simulation (dashed purple); ε_data/ε_MC slider moves off 1 | docs/12 |
| ✔ 2-15 | `pipe_e1_simulation_stack` | chained: zoom out, zoom into backgrounds: **real Z→μμ SR stack, log y, 2 GeV** builds VV+ττ → tt̄/tW → Z/γ*→ℓℓ (generic legend) | data/zmumu_sr_stack.json (nominal) |
| ✔ 2-16 | `pipe_e2_control_region` | chained: dashed ℓ±ℓ± box, mini slice, same-sign pair with a lepton in a jet | docs/13 |
| ✔ 2-17 | `pipe_e3_fake_factor` | chained: arrow × f into the real Fakes layer at the bottom of the stack | docs/13, data/zmumu_sr_stack.json |
| ✔ 2-18 | `pipe_f1_data` | chained: zoom out, zoom into comparison: real stack + key; real data points | data/zmumu_sr_stack.json |
| ✔ 2-19 | `pipe_f2_ratio` | chained: ratio panel: real data/pred (2 GeV) | data/zmumu_sr_stack.json |
| ✔ 2-20 | `pipe_f3_uncertainty` | chained: the (schematic) systematic band breathes | docs/14 |
| ✔ 2-21 | `pipe_g1_fit_model` | chained: zoom out, zoom into fit: μ_Z slider, θ₁…θ₅ rows, real pre-fit ratio (5 GeV) | data/zmumu_fit.json |
| ✔ 2-22 | `pipe_g2_fit` | chained: the fit: μ moves and tightens, pulls move, ratio goes to the real post-fit (flat) | data/zmumu_fit.json |
| ✔ 2-23 | `pipe_g3_cross_section` | chained: σ = μ_Z σ_pred; σ = (N − B)/(A ε L) lit term by term | CONVENTIONS.md §2, docs/06 |
| ✔ 2-24 | `pipe_h_three` | chained: zoom out to the finished spine; it shrinks and triplicates in ee green / μμ gold / ττ red | palette `CHANNEL` |

(The earlier candidates `open_data_pipeline`, `selection_funnel`, `xsec_formula`, `tag_and_probe`, `fake_factor`,
`profile_likelihood`, `luminosity` are superseded by the `pipe_*` chain, 15 Sep 2026. On 16 Sep the chain was cut
into one-idea clips (manim sections, `tools/deliver_chain.py`; each zoom-out is joined to the next zoom-in) and the
stacks became the real Z→μμ signal region on a log axis, user request.)

## 3 — Z → ee (`3_zee`)

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 3-01 | `ee_process` | Drell-Yan diagram in the ee flavour; each electron leg radiates and the photons convert: an EM shower in Feynman style, inside a faint cone | `shower_tree` |
| ✔ 3-02 | `ee_detector` | chained: dissolves to the slice; e⁻/e⁺ tracks bend opposite ways, stop in the ECAL as green clusters, one brem photon, nothing beyond; camera zooms in at the end | `signature("e")` |
| | `mass_peak_fill_ee` | one event → one bar; rain of events fills the m_ee peak (schematic, seeded) | recipes §4 |
| | `ee_stack` | frozen data/MC stacked histogram (DYee, DYtautau, TTbar, dibosons, fakes) with data points, drawn step by step | data/ee_*.json from z-ee |
| | `ee_result` | the number with stat and syst bars, next to theory | z-ee handoff |

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
