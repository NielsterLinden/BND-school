# FREEZE — BND-school Z cross sections, 17 Sep 2026

The analyses were **frozen on 17 Sep 2026**: there will be no more reruns and no improved analyses. The numbers
below are final. Cite them, and build the talk on them. Git tag: `zmumu-freeze-2026-09-17` (local
annotated tag on the commit that finalised the combination).

Every channel and combination output (`z-*/fit/`, `z-*/output/`, `combination/combLieke/output/`,
`combination/combLieke/checks/`) is read-only from here on. If a number is needed that is not below, read it from
the listed file; do not regenerate it.

> **After the freeze (merge 5b7cc83, 17 Sep):** z-tautau pushed its four-channel v4 measurement (τhτh + μτh + eτh + eμ).
> On `main` it replaces the v3 τhτh inputs: `z-tautau/fit/ztautau.config` and `fit/fitinputs/ztautau.root` are the
> four-channel ones, and `z-tautau/fit/results/` (with `ztautau_fit_result.json`) is gone until z-tautau runs its fits
> (`z-tautau/RESUME.md`). The combination and the ττ numbers below are the **v3** ones; their inputs are at the tag
> (`git show zmumu-freeze-2026-09-17:z-tautau/fit/...`). `combination/combLieke/run.py` no longer reproduces the freeze
> on `main` (the ττ acceptance keys `signal_prediction.A_unc.*` and the v3 fit result are missing): check out the tag to rerun it.

All results: CMS Open Data 2016 G+H, NanoAODv9, √s = 13 TeV, L = 16393.381 pb⁻¹ ± 1.2 %.
Theory reference: aMC@NLO (NNLO-normalised, 6077.22 pb for m > 50 GeV).

## The final numbers

### Combination (Z → ee ⊕ μμ ⊕ τhτh, one TRExFitter v1.8.0 MultiFit)

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1945 ⁺³¹₋₃₀ pb** per lepton flavour, assuming lepton universality
> = 1945.4 ± 0.5 (stat) ⁺³¹·⁰₋₂₉.₇ (syst) pb; μ_Z = 0.9956 ⁺⁰·⁰¹⁵⁹₋₀.₀₁₅₂ × 1953.93 pb

| | value | source |
|---|---|---|
| combined σ(60–120) | **1945.4 ⁺³¹·⁰₋₂₉.₇ pb** | `combination/combLieke/output/result.json` → `combined` |
| aMC@NLO prediction | 1953.9 ⁺⁵⁶·¹₋₈₁.₇ pb | `result.json` → `prediction` |
| channel compatibility | −2 ln(L_common/L_split) = 2.44 for 2 dof, **p = 0.30** | `result.json` → `compatibility` |
| per channel, same likelihood (standalone) | ee 2093.7 ⁺¹²¹·⁹₋₁₁₄.₄, μμ 1931.1 ⁺³⁰·⁸₋₃₀.₁, ττ 2082.5 ⁺²³⁶·⁸₋₂₀₆.₈ pb | `result.json` → `channels.<ch>.standalone` |
| per channel, three-POI joint fit | ee 2097.4, μμ 1940.3, ττ 2082.2 pb | `channels.<ch>.joint` |
| uncertainty groups (covariance decomposition, pb) | luminosity 22.5, acceptance 13.1, L1 prefiring 9.7, muon efficiency 8.4, MC statistics 5.2, signal modelling 3.2, background normalisation 2.5, electron ID 2.3; data statistics 0.47 (stat-only fit) | `combined.grouped_impacts`, `combined.stat_pb` |
| largest ranked parameters (post-fit impact, pb) | Lumi ±22.6, ElectronID ±12.3, L1Prefiring ±9.7, Acc_PDF ±9.4, Acc_PTZ_mumu ±7.0, SigModel_mumu ±7.0, MuonScale ±6.1 | `combined.ranking`, figure `impacts` |
| variations (Δ vs baseline) | without ee 1932.2 (−13.2); without ττ 1943.9 (−1.5); shape/norm split in all channels 1952.8 (+7.4); ee split only shared NPs 1883.9 (−61.5); ee ElectronID ±1.2 % (diagnostic) 1935.8 (−9.6) | `result.json` → `variations` |
| goodness of fit (saturated) | combined 8 × 10⁻³³ (ee peak shape); μμ 0.79, ττ 0.21, ee 4 × 10⁻⁴² | `combined`, `channels.<ch>.standalone` |
| fit health | every fit MIGRAD/HESSE/MINOS status 0, no forced positive-definite covariance; HESSE/MINOS on μ_Z = 1.0015; 75/75 ranking refits converged | `combined.hesse_over_minos` |
| figures | `summary`, `channels_vs_published`, `combined_vs_published`, `breakdown`, `impacts`, `pulls`, `nll_scan`, `variations` (PDF + PNG) | `combination/combLieke/output/plots/` |

Method and discussion: `combination/combLieke/README.md`; systematics vs ATLAS/CMS: `combLieke/docs/systematics.md`;
orthogonality: `combLieke/docs/orthogonality.md`; agent guide: `combination/CLAUDE.md`.

### Z → μμ (`z-mumu/`, v2 + the CMS-SMP-20-004 parity measurements)

> **σ_fid(pp → Z/γ\* → μ⁺μ⁻; dressed, pT > 26/20 GeV, |η| < 2.4, 60 < m < 120 GeV) = 790.2 ± 0.2 (stat) ± 6.1 (syst) ± 9.6 (lumi) pb**
> **σ(Z/γ\* → μμ, 60–120 GeV) = 1931 ± 0.6 (stat) ± 15.0 (syst) ± 13.5 (acceptance) ± 23.4 (lumi) pb = 1931 ± 30 pb**

| | value | source |
|---|---|---|
| μ_Z | 0.9883 ⁺⁰·⁰¹⁴¹₋₀.₀₁₃₈ (0.988 ± 0.014) w.r.t. σ_fid^pred = 799.6 pb | `z-mumu/fit/results/zmumu_fit_result.json` |
| σ(m > 50 GeV) | 2002 ± 32 pb (A = 0.3947) | same, `sigma_m50_*` |
| acceptance A(60–120) | 0.4092 ± 0.70 % (PDF 0.50, boson p_T 0.38, scale 0.29, QED FSR 0.10, MC stat 0.05, generator 0.04, α_s 0.03, PS FSR 0.03 %) | same, `acceptance`; `z-mumu/zmumu/acceptance.py` |
| muon reconstruction SF | 1.0001 ± 0.0013 per muon (measured, tag-and-probe on unskimmed NanoAOD); `MuonReco` ±0.27 % per event | `z-mumu/output/v2/tnp/reco_result.json` |
| goodness of fit | p = 0.79 | fit result `gof` |
| counting cross-check | σ_fid = 794.5 pb (n_obs 10,378,567, n_bkg 68,799) | fit result `meta.counting` |
| binning / model stability | ±0.2 % | `z-mumu/fit/results/stability.json` |
| vs CMS-SMP-20-004 (1952 ± 49 pb) | −21 pb, −0.36σ (budget total 1.59 %, 30.7 pb with the luminosity on the measured value) | `z-mumu/output/v2/cms_parity/parity_zmumu.{json,md}`, `z-mumu/docs/16` |
| report, plots | | `z-mumu/output/v2/RESULTS_v2.md`, `results_v2.json`, `output/v2/plots/` |

Handoff: `z-mumu/handoff.md` (section "Frozen result"). Superseded 15 Sep numbers (reco SF assigned, m > 50 acceptance
uncertainties): σ_fid 790.1 ± 0.2 ± 8.5 ± 9.6 pb, μ_Z 0.988 ± 0.016, σ(60–120) 1931 ± 33 pb, counting 794.7 pb —
kept in `z-mumu/fit/results/zmumu_v2_15sep_fit_result.json`.

### Z → τhτh (`z-tautau/`, v3, unchanged by the freeze)

> **σ(pp → Z/γ\* → ττ, 60–120 GeV) = 2082 ± 41 (stat) ⁺²²²₋₁₉₄ (syst+stat) ± 76 (acc) pb** (prediction 1944.9 pb)
> **σ_fid(τhτh) = 4.82 ± 0.09 (stat) ± 0.47 (syst) pb** (prediction 4.50 pb); μ_Z = 1.071 ⁺⁰·¹¹⁴₋₀.₁₀₀

Source: `z-tautau/fit/results/ztautau_fit_result.json` and `z-tautau/handoff.md` **at the tag** (replaced on `main` by v4, see the top). Acceptance 3.7 %.

### Z → ee (`z-ee/`, delivery of 16 Sep 14:44, unchanged by the freeze)

| | value | source |
|---|---|---|
| the channel's own fit | μ = 0.942 ± 0.015 → σ(60–120) = 1840.8 ± 29.9 pb | `z-ee/Zee_fit.tar.gz` → `Zee_fit/Fits/Zee_fit.txt` |
| in the combination's likelihood (shape and normalisation split) | 2093.7 ⁺¹²¹·⁹₋₁₁₄.₄ pb | `combination/combLieke/output/result.json` |

`z-ee/handoff.md` ("Results so far") still quotes 2051.8 ± 39.5 pb from an earlier fit; the combination uses the
tarball. Known limitations, not fixed before the freeze: the ECAL barrel–endcap gap is not vetoed (so the
`ElectronID` uncertainty is ±5.9 % instead of ~1.2 %), there is no `HLT_Ele27_WPTight_Gsf` scale factor, and there
are no charge-misID or multijet uncertainties. See `combination/combLieke/README.md`, "The ee channel".

## What changed on 17 Sep, before the freeze

1. **Z → μμ** (commit 60602fb). The measurements of the uncertainty-parity pass (`z-mumu/docs/16`) became the
   baseline:
   - the reconstruction SF is measured and applied (it had been assigned);
   - the acceptance block is computed for the 60–120 GeV denominator, with boson-p_T, generator, PS FSR and
     QED FSR rows (`zmumu/acceptance.py`, 0.61 → 0.70 %).
2. **Combination** (the commit carrying this file). It was rerun on those inputs:
   - eight μμ acceptance parameters instead of four: `Acc_PDF`, `Acc_AlphaS`, `Acc_QCDScale`, `Acc_PS_FSR`
     shared with ττ, plus `Acc_PTZ_mumu`, `Acc_Generator_mumu`, `Acc_QEDFSR_mumu`, `AccStat_mumu`;
   - **MINUIT strategy 2** for every fit. With strategy 1 the extra normalisation parameters, all degenerate
     with μ_Z, gave a wrong HESSE covariance at status 0. HESSE was 1.65× below MINOS, the Lumi post-fit error
     0.71, and the luminosity group 1 pb instead of 22.5. MINOS on μ_Z was right either way.
   - Exceptions: the ranking refits use TRExFitter's default strategy escalation (`ranking_fit`; with a
     fixed strategy 2, 27 of 75 failed), and `without_ee` uses strategy 1 (status 1 with strategy 2).
   - `run.py results` now records `combined.hesse_over_minos` and warns if it is off by more than 10 %.

   Result: 1948 ⁺³³₋₃₂ → **1945 ⁺³¹₋₃₀ pb**. The luminosity group went 20.7 → 22.5 pb, acceptance
   10.6 → 13.1, muon efficiency 14.8 → 8.4 pb.

Documentation updated for the freeze:
- the repository `CLAUDE.md` and `README.md`;
- `combination/{CLAUDE,README,handoff}.md`, `combLieke/README.md`, `combLieke/docs/systematics.md`;
- `fitting/CONVENTIONS.md`, `fitting/UNCERTAINTY_PARITY.md`, `.claude/agents/uncertainty-parity-auditor.md`;
- the z-mumu docs;
- the presentation agent definitions in `presentation/setup_and_reference/agents/`.

## For the animation agents: what in `presentation/` still shows pre-freeze numbers

Nothing under `presentation/` was regenerated. Only the agent definitions in
`presentation/setup_and_reference/agents/` now quote the frozen numbers. Against the frozen outputs:

| where | pre-freeze | frozen |
|---|---|---|
| `presentation/data/extract_zmumu_fit.py` asserts (lines ~130–189) | μ_Z 0.988134, 0.988 ± 0.016; σ_fid 790.1, syst 8.5; σ(60–120) 1931 ± 33; C 0.7914 / 0.791415; counting 794.7 pb; stability nominal 0.988134 | μ_Z 0.988286 (+0.014129 −0.013773), 0.988 ± 0.014; σ_fid 790.2, syst 6.1; σ(60–120) 1931 ± 30; C 0.7916 / 0.791567; counting 794.5 pb; stability nominal 0.988286 |
| `presentation/data/extract_zmumu_sr_stack.py:292` | background + fakes 68,786.5 | 68,799 (fit meta `counting.n_bkg`; the simulation carries the reconstruction SF 1.0002 per event) |
| `presentation/data/zmumu_fit.json` | extracted 15 Sep | re-extract: post-fit pulls, `MuonReco` constraint (±0.27 % per event, was ±0.8 %), μ_Z interval |
| clips 4-21 `mumu_h2_fit`, 4-22 `mumu_h3_sigma_fid`, 4-23 `mumu_h4_sigma_total`; pipeline 2-21/2-22 (`pipe_g1_fit_model`, `pipe_g2_fit`) | μ_Z = 0.988 ± 0.016; σ_fid = 790.1 ± 0.2 ± 8.5 ± 9.6 pb; σ(60–120) = 1931 ± 33 pb | μ_Z = 0.988 ± 0.014; σ_fid = 790.2 ± 0.2 ± 6.1 ± 9.6 pb; σ(60–120) = 1931 ± 30 pb |
| `setup_and_reference/briefs/zmumu.md` §1, §3; `03-sections-storyboard.md` rows 4-21, 4-22; `06-chapter-anchors.md` result-format example | as above; "`combination/result.md`" | as above; `combination/combLieke/output/result.json` |
| data/MC ratios of the μμ stack clips (4-08 … 4-28) | 0.944 → … → 0.994 | unchanged at three decimals (the reconstruction SF moves the simulation by 0.02 %) |
| 6-01 `three_to_one` | schematic, no numbers | — |
| any combination clip still to be built | 1948 ⁺³³₋₃₂ pb (16 Sep) if copied from an old doc | 1945 ⁺³¹₋₃₀ pb, groups and variations as above, from `combination/combLieke/output/result.json` |

The git worktrees under `.claude/worktrees/` belong to other sessions and still carry pre-freeze copies of these
docs; they are not the source of truth.

**Status, 17 Sep afternoon — section 4 (Z → μμ) is done** on branch `worktree-tnp-explainer` (not yet merged into
`main`): the three μμ extractors and `zmumu_{fit,sr_stack,tnp}.json` are on the frozen result (all checks pass; the
table above under-reported `extract_zmumu_fit.py`: the GoF, the grouped impacts, the `ratio_prefit` range and the
YAML-vs-`histograms.pkl` comparison had to move too), the chapter was rebuilt as 21 clips numbered in play order
(4-01 … 4-21, `presentation/setup_and_reference/briefs/zmumu.md`), and the brief, the storyboard and the
result-format anchor quote the frozen numbers. The data/MC ladder of those clips is now told tag-and-probe first:
0.944 → 0.971 → 0.968 → 0.974 → 0.994 (`zmumu_sr_stack.json` block `recut`; the old-order stages are kept for
section 2). **Still open from the table:** the pipeline clips 2-21/2-22 and any combination clip.
