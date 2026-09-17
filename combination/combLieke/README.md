# Z → ee ⊕ Z → μμ ⊕ Z → τhτh: TRExFitter MultiFit

One profile-likelihood fit (TRExFitter v1.8.0 MultiFit) of the three BND-school channels on CMS Open
Data 2016 G+H (16.4 fb⁻¹, √s = 13 TeV). Shared nuisance parameters are fitted jointly, and the μμ and ττ
acceptance uncertainties are inside the likelihood.

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1948 ⁺³³₋₃₂ pb**
> = 1948 ± 0.5 (stat) ⁺³³₋₃₂ (syst, including 20.7 pb luminosity and 10.6 pb acceptance) pb,
> per lepton flavour, assuming lepton universality.
>
> aMC@NLO (NNLO-normalised): 1954 ⁺⁵⁶₋₈₂ pb. CMS (2024): 1952 ± 49 pb.

All numbers: `output/result.json`. All figures: `output/plots/` (PDF and PNG).

| | figure |
|---|---|
| each channel and the combination against the prediction | `summary` |
| each channel against the published ATLAS/CMS measurement of the same decay | `channels_vs_published` |
| the combination against the published combined measurements | `combined_vs_published` |
| uncertainty by source | `breakdown` |
| the 20 nuisance parameters with the largest impact, with their pulls | `impacts` |
| every nuisance parameter | `pulls` |
| profile likelihood of the combined cross section | `nll_scan` |
| the combination under alternative likelihood models | `variations` |

![summary](output/plots/summary.png)

## Results

| | σ(60–120) [pb] | from |
|---|---|---|
| **combined** | **1947.8 ⁺³³·¹₋₃₂.₄** | common POI, all three channels |
| Z → ee | 2093.8 ⁺¹²¹·⁹₋₁₁₄.₄ | standalone fit, same model |
| Z → μμ | 1930.8 ⁺³³·²₋₃₂.₆ | standalone fit, same model |
| Z → τhτh | 2082.5 ⁺²³⁶·⁸₋₂₀₆.₈ | standalone fit, same model |
| aMC@NLO | 1953.9 ⁺⁵⁶·¹₋₈₁.₇ | scale +2.5/−3.9 %, PDF 0.74 %, α_s 1.3 % |

* **Channel compatibility.** The three-POI fit (one σ per channel, every shared NP profiled together)
  gives ee 2097.5, μμ 1939.8, ττ 2082.2 pb. Against the common POI,
  −2 ln(L_common/L_split) = **2.41 for 2 degrees of freedom, p = 0.30**.
* **Weights.** The result is carried by μμ. Without ee it is 1932.3 ⁺³³·²₋₃₂.₄ pb; without ττ it is
  1945.4 pb (`variations`).
* **Goodness of fit.** Saturated model: μμ p = 0.79, ττ p = 0.21, **ee p = 4 × 10⁻⁴²**. The combined
  p = 8 × 10⁻³³ is entirely the ee peak shape (next section).
* **Fit quality.** Every fit (combined, three-POI, standalone, stat-only, variations) ends with MIGRAD 0,
  HESSE 0, MINOS 0 and no forced positive-definite covariance.
* **The channels' own numbers.** ee 1840.8 ± 29.9 pb (`z-ee/Zee_fit.tar.gz`), μμ 1930.7 pb, ττ 2082.5 pb.
  μμ and ττ are reproduced exactly; their errors are now slightly larger because acceptance is in the fit.
  ee differs, as explained below.

## The likelihood

Every channel enters with **its own TRExFitter config and fit inputs**:

| channel | config | inputs | signal | σ_ref(60–120) |
|---|---|---|---|---|
| ee | `z-ee/fit.config` | `z-ee/Zee_fit.tar.gz` → `Histograms/Zee_fit_histos.root` (bin by bin identical to `datasets/z-ee/zee.root`, checked on every run) | `DYee` | 1954.10 pb (the tarball's `reference_cross_section.txt`) |
| μμ | `z-mumu/fit/zmumu.config` | `z-mumu/fit/fitinputs/zmumu.root` | `DYmumu` | 1953.93 pb |
| ττ | `z-tautau/fit/ztautau.config` | `z-tautau/fit/fitinputs/ztautau.root` | fiducial `DYtautau` | 1944.88 pb |

The only changes, made by `mf/trexcfg.adapt_channel`:

1. **One POI.** Each channel's signal NormFactor is renamed `mu_Z`, and a constant NormFactor
   `xsref_<channel>` = 1953.93 pb / σ_ref(channel) is added on the same samples. So
   σ = μ_Z × 1953.93 pb in every channel, even though the three aMC@NLO references differ (by 0.46 % for ττ).
2. **Acceptance inside the fit.** μμ and ττ quote σ(60–120) = σ_fid / A, with A from aMC@NLO. Their
   δA/A enter as OVERALL parameters on the signal: `Acc_PDF`, `Acc_AlphaS`, `Acc_QCDScale`,
   `Acc_PS_ISR`, `Acc_PS_FSR` (shared μμ/ττ) and `AccStat_<channel>`. The values are read from
   `zmumu.root.meta.json` (total 0.61 %) and `ztautau.root.meta.json` (3.7 %). ee has no acceptance
   term: its template is already normalised to the 60–120 GeV LHE cross section.
3. **Correlations by name,** as `fitting/CONVENTIONS.md` §3 prescribes. Shared are `Lumi`, `Pileup`,
   `L1Prefiring`, `PDF`, `QCDScale`, `PS_FSR`, `PS_ISR`, `XS_TTbar`, `XS_DYtautau`, `XS_WJets`,
   `XS_SingleTop/WW/WZ/ZZ` and the `Acc_*` parameters. Electron, muon, τ, MET, fake-factor and MC-statistics
   parameters are channel-specific. `SigModel` → `SigModel_mumu`.
4. **ee: shape and normalisation as separate parameters.** Every template (HISTO) systematic of ee
   enters twice. Its normalisation keeps the shared name, and its shape gets an ee-only parameter
   `<NP>_eeShape` (TRExFitter `DropShapeIn` / `DropNorm` on the same templates).
5. **ττ: one empty bin dropped.** `tautau_SR2`, 0–40 GeV, has 0 data and 0 predicted events, so its
   MC-statistics γ is unconstrained and sits at 0, which broke HESSE and MINOS. The adapter checks that
   data and every sample are exactly zero before dropping it.

Nothing else is touched: binning, samples, smoothing, symmetrisation, DropBins, MC statistics.

## The ee channel, and why point 4 is needed

With all three channels correlated as delivered, the fit has **no positive-definite minimum**
(TRExFitter retries up to strategy 3 and crashes). This is also what stopped the earlier three-channel
attempt. The cause is in the ee inputs:

* **The ee model does not describe its own peak shape:** saturated-model p = 7 × 10⁻⁴³ with 6.3 M events
  in 30 × 2 GeV bins. The fit absorbs the mismatch by pulling template parameters (`Pileup` −1.6σ,
  `L1Prefiring` +2.0σ, `QCDScale` −2.0σ, `ElectronID` +0.6σ) and constraining them far below their priors.
  Shared with μμ, those pulls move the μμ normalisation. In the three-POI fit without the split,
  μ(μμ) shifts by +2 %.
* **The ee electron-ID uncertainty, ±5.9 %, is an artefact.** The official UL2016postVFP Medium-ID
  scale-factor map gives about 1.2 % per Z → ee event. In the ECAL barrel–endcap gap
  (1.444 < |η_SC| < 1.566) it returns the placeholder sf = 1 ± 1, and z-ee does not veto the gap. The
  shape of this ±100 % variation on ~5 % of events pins the whole 5.9 % normalisation at 0.08σ.
  That is where the channel's quoted ±1.6 % precision comes from (details in `docs/systematics.md`).

Splitting shape from normalisation keeps the correlated normalisations physical. Meanwhile ee's shape
mismodelling stays in ee-only parameters. The price is honest: once the peak shape no longer pins it,
ee's normalisation is limited by the ±5.9 % it was delivered with. Its standalone result becomes
2094 ⁺¹²²₋₁₁₄ pb. The central value sits along an almost flat ElectronID–μ direction, pulled there by the
small backgrounds. ee therefore contributes little to the combined value, and its channel point should
be read with that ±6 % in mind.

How much this matters (`variations`, all converged):

| likelihood | σ [pb] | Δ |
|---|---|---|
| baseline | 1947.8 ⁺³³·¹₋₃₂.₄ | |
| without ee | 1932.3 ⁺³³·²₋₃₂.₄ | −15.6 |
| without ττ | 1945.4 ⁺³³·³₋₃₂.₅ | −2.4 |
| shape/normalisation split in all three channels | 1954.7 ⁺³⁴·⁴₋₃₃.₇ | +6.9 |
| ee: split only the shared parameters (the gap-driven `ElectronID` shape constraint trusted) | 1873.2 ⁺²⁷·³₋₂₆.₈ | −74.6 |
| ee: `ElectronID` normalisation ±1.2 % (official map outside the gap) — *diagnostic* | 1933.5 ⁺³⁰·⁴₋₂₉.₉ | −14.3 |

The last row is not a result: it uses an uncertainty z-ee did not deliver. It shows what to expect
once the gap is vetoed. In that fit ee alone gives 1918 ± 39 pb, in agreement with μμ. The
"shared only" row shows the opposite: if the artefact constraint is trusted, the μμ `MuonReco`
parameter is pulled to +2.5σ.

**For z-ee** (in order of impact): veto the ECAL gap; apply an `HLT_Ele27_WPTight_Gsf` scale factor
with its uncertainty (no trigger correction is applied at all now, which alone could explain the ee
deficit in its own fit); add charge-misidentification and multijet uncertainties; re-check the goodness
of fit, and consider coarser bins if it is still poor.

## Checks

* **Orthogonality** (`docs/orthogonality.md`, measured on data event by event): μμ ∩ ee ≤ 97 events
  (1.5 × 10⁻⁵ of ee), ττ ∩ μμ = 0, ττ ∩ ee ≤ 4 events (2 × 10⁻⁴ of the ττ signal region). The channels
  are statistically independent to this precision.
* **Systematics against ATLAS and CMS** (`docs/systematics.md`): μμ carries the same list as the
  published measurements, with comparable sizes. ττ is complete but lacks CMS's in-situ τh-ID constraint.
  ee lacks trigger, charge-misID and multijet uncertainties, and has the gap artefact above.
* **Published comparisons** (`config/references.json`, each value with paper and table): ATLAS 13 TeV
  ee and μμ (arXiv:1603.09222, moved from 66–116 to 60–120 GeV with the aMC@NLO ratio 1.01425), CMS
  PAS SMP-15-004 ee and μμ, CMS Z → ττ (all five final states, and τhτh alone; arXiv:1801.03535), and the
  combined CMS (arXiv:2408.03744) and ATLAS values. CMS SMP-20-004 publishes no per-channel cross section.

## Uncertainties

The breakdown (`breakdown`) is TRExFitter's covariance decomposition of the combined fit, per Category.
The groups do not add to the total in quadrature, because the post-fit parameters are correlated.
Largest: luminosity 20.7 pb, muon efficiency 14.8, acceptance 10.6, L1 prefiring 9.1, MC statistics 5.3,
electron ID 3.3, signal modelling 3.3 pb. Data statistics, 0.47 pb, come from a separate stat-only fit.
TRExFitter's refit-based grouped impacts (`trex-fitter mi`) fail HESSE in this likelihood and are not
used. The ranking (`impacts`) is refit-based: one `trex-fitter mr Ranking=<NP>` per parameter.

The theory band is the aMC@NLO sample's own uncertainty on σ(60 < m_LHE < 120): 7-point μ_R/μ_F envelope
⊕ NNPDF3.1 Hessian ⊕ α_s ± 0.0015, from `z-mumu/output/v2/gensums.json` (`mf/prediction.py`). The
uncertainty of the NNLO normalisation (6077.22 pb) is not public and not included.

## Run

```bash
source ../../setup.sh
python run.py all                 # ~10 min on stbc-i*; TRExFitter output in work/ (git-ignored)
python tests/test_trexcfg.py      # the adapter
python checks/systematics.py      # -> checks/systematics.json
python checks/orthogonality.py    # -> checks/orthogonality.json (~3 min, reads the data skims)
```

`python run.py <step>` runs one step: `prepare`, `workspaces`, `fits`, `variations`, `impacts`,
`results` or `plots`.

## Files

| path | what |
|---|---|
| `run.py` | the pipeline, step by step |
| `config/channels.json` | the channels, POI, acceptance terms, ee treatment, variations |
| `config/references.json` | published measurements |
| `mf/trexcfg.py` | reads, adapts and writes TRExFitter configs |
| `mf/ee_input.py` | ee histograms from the tarball |
| `mf/prediction.py` | aMC@NLO prediction and uncertainty |
| `mf/results.py`, `mf/plots.py` | `output/result.json`, `output/plots/` |
| `checks/` | orthogonality and systematic-size checks, with their JSON outputs |
| `docs/` | `orthogonality.md`, `systematics.md` |
| `tests/` | adapter tests |
