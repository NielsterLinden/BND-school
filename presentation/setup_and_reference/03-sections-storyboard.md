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
| | `pipe_a_map` | **the analysis spine** (schematic, no numbers): nodes appear left→right — collisions/detector → recorded & simulated events (two rivers merging into one file stack) → selection funnel → corrections (tag & probe) → backgrounds (control region) → comparison (stack + ratio) → fit (μ) → σ | all channel docs |
| | `pipe_b_data` | chained: zoom on nodes 1–2; events flash in the slice and stream into a file stack that shrinks (skim); the simulated river (mini Drell-Yan → the same slice → the same files) merges in | z-mumu docs/10-skims.md |
| | `pipe_c_selection` | chained: the funnel — a bar shrinks through trigger → ℓ ID/isolation → exactly two leptons → opposite sign → mass window (proportions only) | channel handoffs |
| | `pipe_d_tnp` | chained: tag-and-probe — the tag is clean, the probe is tested; pass/fail histograms; ε = pass/(pass+fail); ε_data vs ε_MC → scale-factor slider | z-mumu docs/12 |
| | `pipe_e_backgrounds` | chained: the stack of simulated processes; the one the simulation cannot give (fakes) comes from a control-region box (same-sign) → sliver at the bottom of the stack | docs/13-fake-factor.md, z-tautau docs/05 |
| | `pipe_f_compare` | chained: data points over the stack, the ratio panel opens, the systematic band breathes | docs/14 |
| | `pipe_g_fit` | chained: μ slider, nuisance-parameter pulls, ratio flattens; σ = μ·σ_pred and σ = (N − B)/(A·ε·L) lit term by term | fitting/CONVENTIONS.md §2 |
| | `pipe_h_three` | chained: the spine shrinks to a strip and triplicates in ee green / μμ gold / ττ red → hand-off to section 3 | palette `CHANNEL` |

(The earlier candidates `open_data_pipeline`, `selection_funnel`, `xsec_formula`, `tag_and_probe`, `fake_factor`,
`profile_likelihood`, `luminosity` are superseded by the `pipe_*` chain, 15 Sep 2026.)

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
| | `mumu_a_event` | chained on 4-02 (**real data from here on**, `data/zmumu_*.json`): the schematic tracks fade, a real SR event (φ, charge, p_T from data) is drawn; p_T values peel off and converge to m_μμ; the slice parks left, log-y axes appear right, first entry | data/zmumu_events.json |
| | `mumu_b_rain` | chained: three more real events, a clock spins, dots rain into the plot, the 60 data bins grow to the real counts (10,378,567) | data/zmumu_sr_stack.json |
| | `mumu_c_stack` | chained: the **uncorrected** simulation stack slides in under the data; the ratio panel opens at ≈ 0.94 (raw simulation ~6 % high) | data/zmumu_sr_stack.json (stage `raw`) |
| | `mumu_d_control` | chained: the plot parks; the slice returns with a real same-sign pair (probe inside a jet): the control region | docs/13-fake-factor.md |
| | `mumu_e_tenpairs` | chained: ten real same-sign pairs, a 10-box tally (9 anti-isolated probes, 1 isolated) → f = N_tight/N_anti ≈ 1/9; the real 6×4 fake-factor map | data/zmumu_fakes.json |
| | `mumu_f_transfer` | chained: OS events with one anti-isolated muon × f → the real fake template slides into the stack; 3 870 ± 80 | RESULTS_v2.md |
| | `mumu_g_corrections` | chained: pileup, L1 prefiring 0.980, ID SF 0.980, iso SF 1.006, trigger ε 0.907/0.923, κ — the stack steps through the frozen stages, the ratio ends at the pre-fit 0.965 … 1.01 | data/zmumu_corrections.json, docs/11-12 |
| | `mumu_h_fit` | chained: rebin 60 → 12, pull plot, post-fit; μ_Z = 0.988 ± 0.016 → σ_fid = 790.1 ± 0.2 (stat) ± 8.5 (syst) ± 9.6 (lumi) pb, σ(60–120) = 1931 ± 33 pb beside the prediction 1953.9 pb | data/zmumu_fit.json, handoff.md |

(The earlier candidates `tnp_efficiency`, `momentum_calibration`, `mumu_stack`, `mumu_result` are superseded
by the `mumu_*` chain, 15 Sep 2026.)

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
