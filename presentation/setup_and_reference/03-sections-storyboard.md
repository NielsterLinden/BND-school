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

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 4-01 | `mumu_process` | Drell-Yan diagram in the μμ flavour; the muon legs simply continue, straight and clean | — |
| ✔ 4-02 | `mumu_detector` | chained: dissolves to the slice; two outlined gold tracks through everything, MIP dots in both calorimeters, stubs in the four muon stations | `signature("mu")` |
| ✔ 4-03 | `mumu_a1_event` | chained: chained on 4-02 (**real data from here on**, `data/zmumu_*.json`): schematic tracks fade, real SR event drawn (φ, charge, p_T), p_T values appear | data/zmumu_events.json |
| ✔ 4-04 | `mumu_a2_mass` | chained: the p_T values converge to m_μμ | data/zmumu_events.json |
| ✔ 4-05 | `mumu_a3_first_entry` | chained: slice parks left, log-y axes appear, m_μμ becomes the first entry | data/zmumu_events.json |
| ✔ 4-06 | `mumu_b1_more_events` | chained: three more real events, each drops an entry | data/zmumu_events.json |
| ✔ 4-07 | `mumu_b2_rain` | chained: clock, rain, the 60 bins grow to the real counts; N = 10,378,567 | data/zmumu_sr_stack.json |
| ✔ 4-08 | `mumu_c1_simulation` | chained: bars → points; the **uncorrected** simulation stack slides in; colour key | data/zmumu_sr_stack.json (stage `raw`) |
| ✔ 4-09 | `mumu_c2_ratio` | chained: ratio panel opens: data/pred = 0.944 | data/zmumu_sr_stack.json |
| ✔ 4-10 | `mumu_d1_control_region` | chained: plot parks top right; slice returns in the dashed μ±μ± box | docs/13-fake-factor.md |
| ✔ 4-11 | `mumu_d2_same_sign_pair` | chained: real same-sign pair: tag, probe inside a jet, charges | data/zmumu_events.json |
| ✔ 4-12 | `mumu_e1_ten_pairs` | chained: ten real same-sign pairs fill a 10-box tally (9 anti-isolated, 1 isolated) | data/zmumu_events.json |
| ✔ 4-13 | `mumu_e2_fake_factor` | chained: f = N_tight/N_anti = 1/9 | docs/13 |
| ✔ 4-14 | `mumu_e3_ff_map` | chained: the real 6×4 fake-factor map, barrel cells highlighted | data/zmumu_fakes.json |
| ✔ 4-15 | `mumu_f1_apply` | chained: back to the SR plot; N_fake = f × N_anti | docs/13 |
| ✔ 4-16 | `mumu_f2_fake_template` | chained: the real fake template enters the stack (key gains Fakes); N_fake = 3870 ± 80 | RESULTS_v2.md, data/zmumu_fakes.json |
| ✔ 4-17 | `mumu_g1_pileup` | chained: ⟨w_PU⟩ = 0.994: the prediction steps down | data/zmumu_sr_stack.json, docs/11 |
| ✔ 4-18 | `mumu_g2_prefiring` | chained: ⟨w_L1⟩ = 0.980 | docs/11 |
| ✔ 4-19 | `mumu_g3_lepton_sf` | chained: SF_ID 0.980, SF_iso 1.006, trigger ε 0.907/0.923 → ratio 0.994; κ_|η| 0.999…1.000 | data/zmumu_corrections.json, docs/12 |
| ✔ 4-20 | `mumu_h1_rebin` | chained: rebin 60 → 12 (5 GeV) | data/zmumu_fit.json |
| ✔ 4-21 | `mumu_h2_fit` | chained: pulls + μ slider; post-fit, ratio flat; μ_Z = 0.988 ± 0.016 | data/zmumu_fit.json |
| ✔ 4-22 | `mumu_h3_sigma_fid` | chained: the plot parks left; μ_Z, σ_fid = 790.1 ± 0.2 ± 8.5 ± 9.6 pb, σ(60–120) = 1931 ± 33 pb | handoff.md |
| ✔ 4-23 | `mumu_h4_sigma_total` | chained: σ(60–120) on its axis beside the prediction 1953.9 pb | CONVENTIONS.md §6 |

(The earlier candidates `tnp_efficiency`, `momentum_calibration`, `mumu_stack`, `mumu_result` are superseded
by the `mumu_*` chain, 15 Sep 2026; cut into one-idea clips on 16 Sep.)

## 5 — Z → ττ (`5_ztautau`)

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 5-01 | `tau_decay` | τ⁻ → ν_τ W*⁻, W*⁻ → d ū, hadronisation blob, π⁻ π⁰ (decay mode 1); the only neutrino is the ν_τ | z-tautau/docs/02-selection.md |
| ✔ 5-02 | `tau_jet` | chained: the boosted τ: pions collimate into a narrow cone, π⁰ → γγ | — |
| ✔ 5-03 | `tautau_detector` | chained: dissolves to the slice; a 1-prong τ_h (track, γγ in the ECAL, HCAL) and a 3-prong τ_h, two ν_τ leave unseen, the p_T^miss arrow recoils; camera zooms in | `signature("tau_h", prongs=...)`, `"met"` |
| | `visible_mass` | why m_vis sits below m_Z: the missing neutrinos, then the likelihood mass | docs/04-ditau-mass.md |
| | `fake_factor_tautau` | same-sign measurement → application region (81 % of the SR are jet→τ_h fakes) | docs/05-fake-factors.md |
| | `tautau_stack` / `tautau_result` | frozen m_ττ stack and σ = 2343 pb result | z-tautau docs/08, handoff |

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
