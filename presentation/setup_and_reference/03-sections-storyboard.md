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
| | `open_data_pipeline` | Run2016G+H boxes → NanoAOD → skim (259 GB → 31.5 GB) → histograms | z-mumu handoff, docs/10-skims.md |
| | `selection_funnel` | event count shrinking through trigger → object ID → two leptons → OS → mass window | channel handoffs |
| | `xsec_formula` | σ = (N − B) / (A · ε · L) with each term lighting up while its measurement icon appears | docs/06-cross-section.md |
| | `tag_and_probe` | Z → μμ: tag passes, probe tested; pass/fail histograms; ε = pass/(pass+fail) | docs/04, 12 |
| | `fake_factor` | control region → fake factor → transfer to SR | docs/13-fake-factor.md |
| | `profile_likelihood` | the fit: templates, nuisance parameters pulled, μ_Z read off | fitting/, docs/14 |
| | `luminosity` | brilcalc lumi sections summing to 16393.381 pb⁻¹ (G + H) | CLAUDE.md, recid 1059 |

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
| | `tnp_efficiency` | efficiency vs p_T / η map fading in from pass/fail fits | docs/12 |
| | `momentum_calibration` | peak shifts before/after calibration (two step histograms, mean marker slides) | docs on momentum calibration |
| | `mumu_stack` | frozen data/MC stack in mumu_SR | z-mumu fit inputs |
| | `mumu_result` | σ_fid = 776.9 ± 0.2 ± 14.8 pb and the m>50 extrapolation | handoff.md |

## 5 — Z → ττ (`5_ztautau`)

| # | candidate | what moves | source |
|---|-----------|-----------|--------|
| ✔ 5-01 | `tau_decay` | τ⁻ → ν_τ W*⁻, W*⁻ → d ū, hadronisation blob, π⁻ π⁰ (decay mode 1); the only neutrino is the ν_τ | z-tautau/docs/02-selection.md |
| ✔ 5-02 | `tau_jet` | chained: the boosted τ: pions collimate into a narrow cone, π⁰ → γγ | — |
| ✔ 5-03 | `tautau_detector` | chained: dissolves to the slice; a 1-prong τ_h (track, γγ in the ECAL, HCAL) and a 3-prong τ_h, two ν_τ leave unseen, the p_T^miss arrow recoils; camera zooms in | `signature("tau_h", prongs=...)`, `"met"` |
| | `visible_mass` | why m_vis sits below m_Z: the missing neutrinos, then the likelihood mass | docs/04-ditau-mass.md |
| | `fake_factor_tautau` | same-sign measurement → application region (81 % of the SR are jet→τ_h fakes) | docs/05-fake-factors.md |
| | `tautau_stack` / `tautau_result` | frozen m_ττ stack and σ = 2255 pb result | z-tautau docs/08, handoff |

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
