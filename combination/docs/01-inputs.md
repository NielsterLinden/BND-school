# 01 — Inputs

Everything the combination reads is committed. `comb/inputs.py` is the only module that touches
files, and every `ChannelResult` carries a `provenance` list naming them.

| channel | file |
|---|---|
| μμ | `z-mumu/fit/results/zmumu_fit_result.json`, `z-mumu/fit/results/stability.json` |
| ττ | `z-tautau/output/results.json` (the variant its own `nominal_variant` names) |
| ee | `combination/inputs/zee_fit_result.json`, written by `tools/extract_zee.py` from `z-ee/Zee_fit.tar.gz` |

## Z → e⁺e⁻ — and the one thing that has to be fixed

z-ee is the newest input and the only one that does not publish a results file: it commits the
TRExFitter job itself (`z-ee/Zee_fit.tar.gz`, 2.7 MB) and its histogram inputs
(`datasets/z-ee/*.root`). `tools/extract_zee.py` is the adapter — it parses `Fits/Zee_fit.txt`,
the error decomposition, the grouped impacts, the ranking and the prefit yields, computes the
reference cross section and writes `inputs/zee_fit_result.json`. Rerun it whenever z-ee
re-publishes; `run_combination.py --check` asserts every number it produced.

### The problem

`fitting/CONVENTIONS.md` §3 requires the signal theory systematics — `PDF`, `QCDScale`, `PS_ISR`,
`PS_FSR` — to **vary the C factor only, renormalised to a constant fiducial yield**. z-ee builds
`h_mass_PDFUp/Down` and `h_mass_ScaleUp/Down` as raw envelopes over the LHE weight replicas,
normalised by the *nominal* sum of generator weights. The consequence is that `QCDScale` carries
the ±5.87 % uncertainty of the aMC@NLO **total cross section** as a pure normalisation on the
signal template — which is exactly degenerate with the parameter of interest.

What the fit then does with it:

| parameter | pull | post-fit constraint | template normalisation at ±1σ |
|---|---:|---:|---:|
| `QCDScale` | **−1.88** | 0.27 | ±5.87 % |
| `ElectronID` | +0.67 | **0.08** | ±5.83 % |
| `Pileup` | −2.90 | 0.52 | ∓0.79 % |
| `L1Prefiring` | +3.58 | 0.92 | ∓0.47 % |
| `LUMI` | −0.87 | 0.95 | ±1.20 % |

`QCDScale` at −1.88σ rescales the *prediction* by κ = 0.893, and `mu_signal` = 1.0509 is measured
against that rescaled prediction. The honest data-over-prediction ratio is

```
(N_data − N_bkg) / N_sig = (6 320 097 − 44 912) / 6 468 986 = 0.9700
```

so the channel's quoted signal strength is 8 % above what the yields say. The two are consistent:
multiplying μ̂ by the response of *every* nuisance parameter at the best fit gives 0.9706, which
closes with the counting ratio to 0.05 % (`tools/extract_zee.py` does exactly this).

### Which pulls may be absorbed into μ̂ and which may not

Write the expected signal yield as `μ · σ^pred · L · (A·C)`, with each nuisance parameter scaling
one of the factors. A pull is *legitimately* absorbed by μ̂ when it scales the **denominator**:

* `LUMI` — if the fit prefers a lower luminosity, the same observed yield means a higher cross
  section, and σ̂ = μ̂ σ^pred is right.
* `Pileup`, `L1Prefiring`, `ElectronID`, `ElectronRECO` — all modify A·C, the efficiency of
  observing the process. Same argument.
* `PDF`, `QCDScale` **as the other two channels build them**, i.e. varying C only. Same argument.

A pull may **not** be absorbed when it scales σ^pred itself, because then σ̂ = μ̂ σ^pred(nominal)
quotes the measurement against a prediction the fit has already moved. Formally,

```
σ̂ = μ̂ · κ_theory · σ^pred,     κ_theory = Π over the un-renormalised theory NPs
```

and κ_theory ≡ 1 for μμ and ττ by construction. For z-ee, κ_theory = **0.8925**. That is the
`ee_normfix` extraction (`inputs.load_ee("normfix")`): σ(ee) = 1833 pb instead of 2054 pb, and
−99 pb on the combination. κ is computed with the HistFactory `code 4` interpolation — exponential
outside ±1σ — reproduced in `tools/extract_zee.py:kappa()` so the folder needs no ROOT.

### What z-ee should change

1. **Renormalise the `PDF` and `QCDScale` envelopes to a constant yield** before writing them:
   divide each replica by *its own* Σw (the replica-weighted `genEventSumw`), or simply rescale
   each replica histogram to the nominal integral. Only the shape may vary.
2. Accumulate the envelopes **per LHE flavour** rather than over the inclusive DY sample (the
   notebook already flags this as an approximation).
3. Follow §1–§4 of the conventions: POI `mu_Z`, `DY_ee` as `Type: SIGNAL`, region `ee_SR`, sample
   and NP names per §3, one input file `z-ee/fit/fitinputs/zee.root`, and a published
   `z-ee/fit/results/zee_fit_result.json` from `fitting/run_trex.py` — then `tools/extract_zee.py`
   can be deleted.
4. Run the saturated-model **goodness of fit**. z-ee is the only channel without one.
5. Fit in coarser bins than 60 × 1 GeV. The −2.9σ `Pileup` and +3.6σ `L1Prefiring` pulls are a
   data/MC shape mismatch being absorbed by calibration parameters; z-mumu hit the same thing and
   fixed it by going from 30 × 2 GeV to 12 × 5 GeV (`z-mumu/REVIEW.md` F3/F4).
6. `Wjets` has 24 % MC statistics and about thirty bins that were negative before
   `scripts/fix_negative_bins.py` floored them to +10⁻⁶. Flooring keeps TRExFitter happy; it does
   not make the template usable. Smooth it, or merge W+jets into the other backgrounds.

### What z-ee already fixed (16 Sep 2026)

* **The reference cross section.** The notebook now computes
  6077.22 pb × Σw(LHE ee, 60 < m_LHE < 120)/Σw = **1954.1032 pb**, which agrees with the value the
  combination derives independently from `z-mumu/output/v2/gensums.json` (1954.1110 pb) to
  4 × 10⁻⁶. z-ee is now on exactly the same footing as μμ and ττ.
* **Negative bins**, floored by `z-ee/scripts/fix_negative_bins.py`.

### Two structural differences that are *not* problems

**No acceptance term.** μμ and ττ quote a fiducial cross section and divide by A, whose
uncertainty sits outside the likelihood. z-ee fits directly against the 60 < m_LHE < 120 GeV
prediction, so the extrapolation uncertainty is inside the fit as `Signal modelling`.
`ChannelResult.acc` is therefore empty for ee — not missing, *inside the fit*.

**The data statistics are the analytic ones.** TRExFitter's `Stat unc.` row for this fit is 0.92 %
of μ_Z, 23× the Poisson value √N_obs/(N_obs − N_bkg) = 0.040 % that 6.3 M selected events give. It
is the quadrature remainder of the total after the grouped impacts, not a measurement — the
grouped impacts of this fit under-shoot its total MINOS error, which is itself a symptom of the
`QCDScale`/`ElectronID`/μ degeneracy above. The combination therefore carries 0.040 % as
`Data statistics` and the rest as a separate, uncorrelated **`Fit residual`** source, so the
channel's published total (2.02 % on μ_Z) is reproduced exactly while the headline "stat" stays
honest. `run_combination.py --check` asserts the quadrature closure.

## Z → μ⁺μ⁻

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
for as a stop-gap is gone: the two-sided `SigModel` carries that uncertainty inside the fit.

The counting extraction and the two alternative binnings the channel accepts are variations:

| variant | configuration | μ_Z | GoF p | shift on the combination |
|---|---|---:|---:|---:|
| `shapefit` *(baseline)* | 12 × 5 GeV, two-sided `SigModel` | 0.98813 | 0.79 | — |
| `bins2gev` | 30 × 2 GeV | 0.99012 | 0.16 | +0.9 pb |
| `bins10gev` | 6 × 10 GeV | 0.98582 | 0.50 | −2.1 pb |
| `counting` | 1 bin | 0.99387 | — | +11.1 pb |

`stab_nosig` and `stab_smooth` are in the channel's table but not offered here: both have
p ≤ 0.01 and the channel rejects them.

The stability table publishes μ_Z, the MINOS errors, the total systematic and the goodness of fit
per configuration, but **not** the grouped impacts per configuration. A variant therefore keeps
the nominal fit's category composition rescaled to its own published total systematic. It is an
approximation, and it is used only for the variations — never for the baseline.

## Z → τhτh

The baseline is whatever `z-tautau/output/results.json` calls `nominal_variant`, so the channel
stays in charge of its own nominal. Since **v3** that is `mcsub` — the fake factor *with* the
genuine-τ MC subtraction — at the **DeepTau Tight** working point on both legs, which v2.1 had
recommended and v3 adopted: μ_Z = 1.071 ⁺⁰·¹¹⁴₋₀.₁₀₀, σ(60–120) = 2082 ± 253 pb, GoF p = 0.25.
v3 publishes only that variant, so the combination carries no ττ cross-check rows.

One consequence for the correlation model: **z-tautau does not fit a generator nuisance parameter
at all**. Its `SigModel_tautau` (madgraph LO vs aMC@NLO, C_LO/C_NLO = 0.841) is reported and not
used as an uncertainty, because the LO sample is simply the worse model of the visible-τ p_T
spectrum and the number would double-count `QCDScale`/`PS_ISR` (`z-tautau/docs/07`). The
`Signal modelling` category on the ττ side is therefore pure PDF/α_s/scale/PS — all of it
correlated with μμ and ee.

## Three asymmetries between the channels

**1. The three μ_Z do not share a reference.** All three define μ_Z as a multiplier on an aMC@NLO
prediction of σ(Z/γ* → ℓℓ, 60 < m < 120 GeV), and all three now build it the same way — but the
generator's LHE flavour shares are not exactly 1/3:

| channel | Σw(LHE flavour)/Σw | σ^pred(60–120) |
|---|---:|---:|
| ee | 0.33386 | 1954.1 pb |
| μμ | 0.33380 | 1953.9 pb |
| ττ | 0.33234 | 1944.9 pb |

ee and μμ agree to 0.01 %; ττ is 0.44 % lower because the τ mass is in the matrix element. The
0.47 % between the extremes is **physical, not definitional** — the branching fractions of the
generator differ. **This is why the combination is done on σ, never on μ_Z.** Each μ̂ is
multiplied by its own reference first, and the spread shows up where it belongs: as a band on the
theory line in `output/plots/forest.pdf` rather than as a line.

A TRExFitter MultiFit with a single shared `mu_Z` would *not* handle this — it would fit one
number against three references. It is a 0.002 % effect today, but it should be fixed before a
joint fit is quoted.

**2. The acceptance lives outside two of the three fits.** For μμ and ττ, μ_Z is a fiducial signal
strength and the conversion to 60 < m < 120 GeV divides by A, whose uncertainty (0.61 % for μμ,
3.7 % for ττ) is not in either likelihood; the combination adds it component by component with its
own correlations ([02](02-correlation-model.md)). For ee it is inside the fit. The ττ acceptance
uncertainty is large because the double 40 GeV cut on the *visible* τ momentum selects the Z p_T
tail, where the QCD-scale dependence is 3.4 %.

**3. Only μμ fits a generator systematic.** Its `SigModel` is powheg vs aMC@NLO (both NLO);
neither ττ nor ee has one. `ChannelResult.sigmodel` is zero for both, which makes the split in
`comb/model.py` a no-op there and the `sigmodel_correlated` variation exactly null.

## Common to all three

| quantity | value | why it matters here |
|---|---|---|
| integrated luminosity | 16 393.381 pb⁻¹ ± 1.2 % | identical normtag value → fully correlated, and the largest single uncertainty |
| DY normalisation | σ(Z/γ*→ℓℓ, m > 50) = 6077.22 pb | same reference cross section |
| signal sample | `DYJetsToLL_M-50` aMC@NLO FxFx (recid 35669) | same events → theory nuisances correlated, MC statistics not |
| fitter | TRExFitter v1.8.0 | same treatment of gammas, smoothing and grouped impacts |
| `Category` strings | `fitting/CONVENTIONS.md` §3 | the keys of the correlation model |

**Event orthogonality.** μμ reads `SingleMuon`, ττ reads `Tau` and vetoes electrons and muons, ee
reads `Electron`. ττ is disjoint from both others by construction; μμ and ee overlap only through
ZZ → 4ℓ, 61 of the 10.38 M μμ SR events (0.0006 %) and 0.001 % of the ee SR. `Data statistics` is
therefore uncorrelated as a fact, not as a modelling choice. The full check — selections, vetoes,
what would break it — is [05-vs-published.md](05-vs-published.md) §1.
