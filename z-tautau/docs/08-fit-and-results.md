# 08 — The fit and the result

Code: `scripts/step5_fit.py` (writes `fit/ztautau.config` with `fitting/trexconfig.py`, runs TRExFitter,
parses with `fitting/run_trex.py`), `scripts/step6_report.py`. Numbers: `output/RESULTS.md`,
`output/results.json`.

## Set-up

* **TRExFitter v1.8.0**, the build of the git submodule; this worktree uses the main checkout's build,
  the same binary as z-mumu (`config.TREX_FALLBACK_HOME`).
* One region `tautau_SR`, observable m_tt in 14 bins: 0, 40, 60, 70, …, 130, 150, 175, 200, 250, 350 GeV
  (overflow in the last bin). Bins above 130 GeV are ~95 % fakes and fix the fake normalisation; 60–110 GeV
  holds half of the signal.
* Samples: Data; **DYtautau** (signal, `NormFactor mu_Z`, range 0–3); DYee, DYmumu, DYlowmass, WJets, TTbar,
  SingleTop, WW, WZ, ZZ; **Fakes** (data × FF). Histograms follow `fitting/CONVENTIONS.md`
  (`tautau_SR__<sample>[__<syst>Up/Down]`, TH1D with Sumw2, written by uproot).
* Nuisance parameters: `07-corrections-and-systematics.md`. Per-bin γ for the MC and fake-template
  statistics (W+jets excluded: its statistical uncertainty is one normalisation parameter). Shape
  variations of τ ES and MET use TRExFitter smoothing (`Smoothing: 40`).
* Runs: `h w f d p i r` on data, stat-only (`StatOnly=TRUE`), Asimov (`FitBlind=TRUE`), for both
  fake-factor variants.

## Result (nominal: fake factors without MC subtraction, as specified)

```
mu_Z = 1.160 +0.192 -0.163      (data stat. 0.018)           expected: +0.169 -0.140
σ_fid(τhτh; vis pT > 40, |η| < 2.1, 60 < m_ττ < 120)  = 5.22 ± 0.08 (stat) ± 0.80 (syst) pb     pred. 4.50 pb
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV)                   = 2255 ± 36 (stat) +371/−315 (syst) ± 83 (A) pb
                                                                                            pred. 1945 pb (NNLO)
```

Saturated-model goodness of fit: p = 0.18. The measurement agrees with the prediction within 1 σ ((1.16 − 1)/0.163 = 1.0).

Alternative with genuine-τ subtraction in the fake-factor regions: **μ_Z = 1.208 +0.201 −0.169**
(σ = 2348 pb, GoF p = 0.54). The prefit fake yield differs by 4.5 % (1 760 events, 22 % of the signal). The
post-fit μ moves by only 0.048 because the high-mass sideband re-normalises the fakes in both fits
(`FakeOSSS_tautau` is pulled −0.86 σ in the nominal, −0.46 σ in the alternative).

### Uncertainty breakdown (nominal, impact on μ_Z)

| group | impact | largest components |
|---|---:|---|
| **τh ID** | **12.5 %** | TauID_DM0 (+9.7/−8.7 %), DM1 (6.3 %), DM10 (4.1 %), DM11 (3.2 %) |
| signal modelling | 9.8 % | SigModel LO vs NLO (+9.1/−7.9 %), QCD scales (4.3 %), PS FSR (1.7 %) |
| τh trigger | 5.1 % | per DM 1–3.5 % |
| MC statistics (γ) | 4.8 % | signal aMC@NLO in the peak bins (~3 % each) |
| τh energy scale | 3.2 % | DM1 2.6 % |
| luminosity | 1.3 % | |
| pileup | 1.2 % | |
| fakes | 0.8 % | constrained by the sideband |
| MET, bkg. normalisation, prefiring | < 0.7 % each | |
| **data statistics** | **1.8 %** | |
| **total** | **+19.2 / −16.3 %** | |

**This is a systematics-dominated measurement, and the leading systematic is external.** The TauPOG τh ID
scale factors for UL2016 postVFP carry 5–18 % per decay mode, entering squared for two τ. Z→ττ is the
process used to *measure* those scale factors, so a τhτh-only cross section cannot do better than the
scale factors it borrows.

### Pulls and constraints worth knowing

| NP | pull | constraint | interpretation |
|---|---:|---:|---|
| FakeOSSS_tautau | −0.86 | 0.69 | sideband wants 2.6 % fewer fakes — in the direction of the MC-subtracted C |
| FakeClosure_tautau | +0.70 | 0.32 | the same-sign non-closure shape is partly preferred; strongly constrained |
| TauID_DM0 | −0.60 | 0.94 | |
| MCStatNorm_WJets_tautau | −0.59 | 0.80 | W+jets post-fit 339 (prefit 514) |
| PS_FSR | −0.41 | 0.74 | |

Post-fit yields: Z→ττ 8672, fakes 37 668 (−2.8 %), tt̄ 540, W+jets 339, Z→ee 212, dibosons 114;
total 47 616 vs 47 586 observed.

## Robustness checks done along the way

| change | μ_Z |
|---|---:|
| first fit, event-wise theory envelopes (bug), W+jets unsmoothed | 1.17 +0.29/−0.22 |
| theory variations on histograms, W+jets smoothed | 1.16 +0.19/−0.16 |
| + low-mass DY sample (15 events) | 1.158 |
| + fake factors per era | **1.160** |
| same, genuine-τ subtraction in the FF regions | 1.208 |

## For the combination

`fit/ztautau.config` + the workspace in `fit/results/ztautau/RooStats/` enter `fitting/combination_skeleton.config`
directly (job name `ztautau`, POI `mu_Z`). In a combined fit with ee and μμ, μ_Z is fixed by the light
leptons to ~1–2 %. The ττ channel then mostly **constrains the τh ID and trigger nuisance parameters**
rather than μ_Z. That is physically sensible (lepton universality) and a useful by-product.

For a BLUE-style combination (`combination/combination.ipynb`): n_obs = 47 586, n_bkg (prefit) = 40 235,
A·ε = 2.471·10⁻⁴ for σ(60–120), counting estimate 1815 pb. **Do not use the counting numbers**: with
B/S ≈ 5 the result is dominated by the fake normalisation, which only the shape fit constrains; use
σ = 2255 pb, +17/−14 % total (τh ID 12.5 % uncorrelated with ee/μμ; lumi 1.2 % correlated;
acceptance 3.7 % correlated).
