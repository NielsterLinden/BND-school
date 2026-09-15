# 01 — Inputs

Everything the combination reads is committed. `comb/inputs.py` is the only module that touches
files, and every `ChannelResult` carries a `provenance` list naming them.

## Z → μ⁺μ⁻

| item | file |
|---|---|
| μ_Z, grouped impacts, ranking, σ's, A, C, L, counting yields | `z-mumu/fit/results/zmumu_fit_result.json` |
| μ_Z, MINOS errors, total systematic and GoF of the alternative fit configurations | `z-mumu/fit/results/stability.json` |

**Which extraction.** The channel's v2 measurement was reviewed on 15 Sep 2026 and the review's
findings were fixed the same day (`z-mumu/REVIEW.md` §0). What the review had objected to (F3/F4)
was that the 30 × 2 GeV fit absorbed a 3–4 % data deficit at 62–78 GeV into a `SigModel` template
that itself carried a phase-space artefact — the powheg sample is generated with m_LHE < 120 GeV,
so its two highest reconstructed bins were depleted by 10 % and 33 %. Both are repaired in the
current fit:

- 1 GeV input bins **fitted in 12 × 5 GeV bins**, so the shape nuisance parameters are no longer
  over-constrained; no template smoothing; `UseMinos: all`;
- the `SigModel` template is built with **both generators inside 50 < m_LHE < 120 GeV** and
  mirrored to be two-sided — the last-bin ratio is 1.02 instead of 0.67.

The result is μ_Z = 0.98813 ± 0.0159 with GoF p = 0.79, and the channel's stability table spans
0.9858–0.9901 over the binnings with p > 0.05 — ±0.2 %. **The combination therefore takes the
shape fit** (`mumu="shapefit"`, the default), and the ±0.7 % lineshape term the review had asked
for as a stop-gap is gone: the two-sided `SigModel` carries that uncertainty inside the fit, so
adding it again would double-count. The `Lineshape model` category no longer exists.

The counting extraction (μ_Z = 0.99387, σ_fid = 794.7 pb) is kept as a **variation**, together
with the two alternative binnings the channel accepts:

| variant | configuration | μ_Z | GoF p | shift on the combination |
|---|---|---:|---:|---:|
| `shapefit` *(baseline)* | 12 × 5 GeV, two-sided `SigModel` | 0.98813 | 0.79 | — |
| `bins2gev` | 30 × 2 GeV | 0.99012 | 0.16 | +3.8 pb |
| `bins10gev` | 6 × 10 GeV | 0.98582 | 0.50 | −4.5 pb |
| `counting` | 1 bin | 0.99387 | — | +11.2 pb |

`stab_nosig` and `stab_smooth` are in the channel's table but not offered here: both have
p ≤ 0.01 and the channel rejects them.

The stability table publishes μ_Z, the MINOS errors, the total systematic and the goodness of fit
per configuration, but **not** the grouped impacts per configuration. A variant therefore keeps
the nominal fit's category composition rescaled to its own published total systematic (the largest
rescaling is +5.4 % for the 1-bin fit). It is an approximation, and it is used only for the
variations — never for the baseline, whose impacts are read directly.

## Z → τhτh

| item | file |
|---|---|
| μ_Z, grouped impacts, ranking, σ's, A and its breakdown, prefit yields | `z-tautau/output/results.json` |
| the DeepTau-Tight cross-check, same chain end to end | `z-tautau/variants/tight/output/results.json` |

The baseline is whatever that file calls `nominal_variant`, so the channel stays in charge of its
own nominal. Since v2.1 (per-category OS/SS correction, SR0 used as a fake sideband) that is
**`mcsub`**, the fake factor *with* the genuine-τ MC subtraction — the opposite of v2.0, where the
unsubtracted variant was nominal. Two cross-checks are carried as variations: `nosub`
(μ_Z = 1.166, +0.1 pb on the combination) and `tight` (μ_Z = 1.071, +0.5 pb), the complete re-run
with DeepTau Tight on both legs that the channel recommends as the next iteration's working point
but has not adopted.

One consequence for the correlation model: since v2.1 **z-tautau does not fit a generator nuisance
parameter at all**. Its `SigModel_tautau` (madgraph LO vs aMC@NLO, C_LO/C_NLO = 0.867) is reported
and not used as an uncertainty, because the LO sample is simply the worse model of the visible-τ
p_T spectrum and the number would double-count `QCDScale`/`PS_ISR` (`z-tautau/docs/07`). The
`Signal modelling` category on the ττ side is therefore pure PDF/α_s/scale/PS — all of it
correlated with μμ — and it shrank from 9.8 % to 1.8 % between v2.0 and v2.1.

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
(the ττ channel's weight is 0.04 %, so 0.47 % × 0.0004 ≈ 0.0002 % on the combination), but the
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
