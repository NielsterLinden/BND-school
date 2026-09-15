# 01 — Inputs

Everything the combination reads is committed. `comb/inputs.py` is the only module that touches
files, and every `ChannelResult` carries a `provenance` list naming them.

## Z → μ⁺μ⁻

| item | file |
|---|---|
| μ_Z, grouped impacts, σ's, A, C, L, counting yields | `z-mumu/fit/results/zmumu_fit_result.json` |
| grouped impacts of the 1-bin (counting) fit | `z-mumu/review/fitcheck/GroupedImpact_rebin30.txt` |
| μ_Z of the 1-bin fit | `z-mumu/review/fitcheck/zmumu_rebin30.txt` |
| the ±0.7 % lineshape term and why it exists | `z-mumu/REVIEW.md` §4 |

**Which extraction.** The channel's headline v2 number is a 30-bin profile-likelihood fit of
m(μμ), μ_Z = 0.9903 ± 0.0140. Its own review (15 Sep 2026, finding F3) then showed that this fit
is not robust: μ_Z moves between 0.990 and 1.006 depending on the binning and on whether the
`SigModel` template is included, because the data show a 3–4 % deficit relative to aMC@NLO at
62–78 GeV that only that one template can absorb — and that template is itself flawed (F4: the
powheg sample has a generator cut at m_LHE < 120 GeV, depleting its two highest reconstructed
bins by 10 % and 33 %). The review's recommendation 3 is to quote the **counting extraction**,
σ_fid = 794.4 pb (μ_Z = 0.99354), with **±0.7 %** — half the spread of the fit-configuration
table — added as a lineshape-model systematic.

The combination follows that recommendation. Consequences:

- the grouped impacts come from the 1-bin fit, not the 30-bin one. They are *larger* and more
  honest: luminosity 1.204 % instead of 1.164 % (the 30-bin fit "measured" the luminosity from
  the tt̄ sideband shape and constrained it to 0.82σ, which is a fit artefact, F5);
- a `Lineshape model` category, μμ-only, carries the ±0.7 %;
- the data statistical uncertainty is √N_obs/(N_obs − N_bkg) = 0.031 %, unchanged.

Running with the shape fit instead (`--` see `VARIATIONS["mumu_shapefit"]`) shifts the combined
value by −6.5 pb and tightens it by 3 pb. It is reported, not used.

## Z → τhτh

| item | file |
|---|---|
| μ_Z, grouped impacts, ranking, σ's, A and its breakdown, prefit yields | `z-tautau/output/results.json` |

The `nominal` fake-factor variant (no MC subtraction) is the baseline, as the channel published
it; `mcsub` is reported as a variation (−2.0 pb on the combination).

## Two asymmetries between the channels

**1. The two μ_Z do not share a reference.** Both channels define μ_Z as a multiplier on an
aMC@NLO prediction of σ(Z/γ* → ℓℓ, 60 < m < 120 GeV), but they build that prediction differently:

| channel | construction | value |
|---|---|---:|
| μμ | σ_fid^pred / A(60–120) = 799.566 / 0.409209, with σ_fid^pred = 6077.22 pb × Σw_fid/Σw | 1953.9 pb |
| ττ | the LHE-ττ subset of the same sample with 60 < m_LHE < 120 GeV | 1944.9 pb |

They differ by 0.47 %. The difference is a definitional one (the μμ denominator is one third of
the total DY cross section scaled by the LHE-μμ 60–120 fraction; the ττ denominator is the LHE-ττ
60–120 cross section of the sample itself), not a statistical fluctuation — the sample has 72 M
events, so the statistical precision on either is 0.02 %.

**This is why the combination is done on σ, never on μ_Z.** Averaging the two μ_Z would silently
average two slightly different quantities. Each channel's μ_Z is multiplied by *its own*
reference first, and the 0.47 % then shows up where it belongs: as a spread of the theory
reference, drawn as a band in `output/plots/forest.pdf` rather than as a line.

A TRExFitter MultiFit with a single shared `mu_Z` would *not* handle this — it would fit one
number against two different references. At the current precision the effect is negligible
(the ττ channel's weight is 0.5 %, so 0.47 % × 0.005 ≈ 0.002 % on the combination), but the
channels should agree on one construction before a joint fit is quoted.

**2. The acceptance lives outside both fits.** μ_Z is a fiducial signal strength; the conversion
to 60 < m < 120 GeV divides by A, whose uncertainty (0.61 % for μμ, 3.7 % for ττ) is not in
either likelihood. The combination adds it, component by component, with its own correlations
([02](02-correlation-model.md)). The ττ acceptance uncertainty is large because the double 40 GeV
cut on the *visible* τ momentum selects the Z p_T tail, where the QCD-scale dependence is 3.4 %.

## Common to both

| quantity | value | why it matters here |
|---|---|---|
| integrated luminosity | 16 393.381 pb⁻¹ ± 1.2 % | identical normtag value → fully correlated, and the largest single uncertainty |
| DY normalisation | σ(Z/γ*→ℓℓ, m > 50) = 6077.22 pb | same reference cross section |
| signal sample | `DYJetsToLL_M-50` aMC@NLO FxFx (recid 35669) | same events → theory nuisances correlated, MC statistics not |
| fitter | TRExFitter v1.8.0 | same treatment of gammas, smoothing and grouped impacts |
| `Category` strings | `fitting/CONVENTIONS.md` §3 | the keys of the correlation model |
