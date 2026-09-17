# Signal extraction, systematics, and the cross section

`zmumu/histograms.py`, `zmumu/momentum.py`, `scripts/v2_4_histograms.py`, `scripts/v2_5_fit.py`,
`fit/zmumu.config`; fit outputs in `fit/results/zmumu/` (TRExFitter v1.8.0).

## Regions and templates

| region | selection | role |
|---|---|---|
| `mumu_SR` | trigger, MET filters, exactly two tight muons (pT > 26 / 20, \|eta\| < 2.4), >= 1 trigger-matched, OS, FSR-recovered 60 < m < 120 GeV; stored in 1 GeV bins, **fitted in 12 bins of 5 GeV** (`--rebin`) | signal region (fit) |
| `mumu_CRemu` | one tight matched muon (pT > 26), no other tight/anti-tight muon, one electron (pT > 20, \|eta\| < 2.5, cutBased medium, not in the gap), OS, 60 < m(e mu) < 120; 12 bins | validation of the flavour-symmetric backgrounds (`--emu-control` promotes it to a control region with `mu_top`) |
| `SS` | as SR, same sign | fake-factor closure, plots |
| `SSemu` | as CRemu, same sign | non-prompt-electron check (plots only) |

Samples: `Data`, `DYmumu` (signal, `NormFactor mu_Z`), `DYtautau`, `DYee`, `TTbar`, `SingleTop`,
`WW`, `WZ`, `ZZ` and the data-driven `Fakes`. In the SR/SS only prompt-prompt MC events enter (the
fake factor replaces the non-prompt muons). In the e-mu regions **all** MC events enter, plus
`WJets` (jet -> e): there is no data-driven estimate of the non-prompt electrons (W+jets, Z -> mumu
+ conversion, Z -> tautau with tau_h -> e, top b -> e), and dropping them was the 9% "excess" of the
first v2 (REVIEW.md F1). The `SSemu` region, dominated by the charge-symmetric part of that
background, is the check that the MC describes it.

## Muon momentum scale and resolution

Rochester corrections for UL2016 are not obtainable here, so the residual data/MC
difference is calibrated on the Z peak: in categories where both muons fall in the same
|eta| bin, the peak is fitted (Breit-Wigner with the PDG Z parameters convolved with a
Gaussian, on a linear background) in data and simulation. MC muon pT is scaled by
`kappa_j = 1 + (mu_data - mu_MC)/m_Z` and smeared by `s_j = sqrt(2(sigma_data^2 - sigma_MC^2))/m_Z`
(`momentum_zpeak_fits.png`, numbers in `RESULTS_v2.md`). `MuonScale` = kappa +- its error
(>= 0.05%), `MuonRes` = smearing off / doubled. Both are re-selections (they move events
across the pT thresholds and the mass window).

## Nuisance parameters

| name | type | applies to | source | Category |
|---|---|---|---|---|
| `Lumi` | norm 1.2% | all MC | CMS-LUM-17-003 | Luminosity |
| `Pileup` | shape | all MC | sigma_minbias +- 4.6% | Pileup |
| `L1Prefiring` | shape | all MC | `L1PreFiringWeight_Up/Dn` | L1 prefiring |
| `MuonID`, `MuonIso`, `MuonTrigger` | shape | all MC | T&P stat (+) syst, coherent | Muon efficiency |
| `MuonReco` | norm 0.27% (T&P, 0.13%/muon, correlated) | all MC | tag-and-probe on the unskimmed NanoAOD (docs/16 §4.1; the SF 1.0001/muon scales all MC); 0.8% assigned until 16 Sep | Muon efficiency |
| `MuonScale`, `MuonRes` | shape (not smoothed) | all MC | Z-peak calibration | Muon momentum |
| `PDF`, `AlphaS`, `QCDScale`, `PS_ISR`, `PS_FSR` | shape | `DYmumu` | LHE/PS weights, renormalised to the fiducial yield (C factor only) | Signal modelling |
| `SigModel` | two-sided shape (SR) | `DYmumu` | powheg / aMC@NLO ratio, both inside the powheg generator window 50 < m_LHE < 120 GeV and normalised to the NLO fiducial prediction in that window; Down = mirrored | Signal modelling |
| `XS_TTbar` 6%, `XS_SingleTop` 10%, `XS_WW/WZ/ZZ` 10%, `XS_DYtautau` 5%, `XS_WJets` 30% (e-mu only) | norm | backgrounds | theory | Background normalisation |
| `FakeStat_mumu`, `FakeMethod_mumu` | shape/norm | `Fakes` (no Sumw2, docs/13) | docs/13 | Fakes |
| `ElectronEff_mumu` 3% | norm | MC in `mumu_CRemu` | no electron T&P here | Electron efficiency |
| MC statistics | gammas per bin | all | `MCstatThreshold: 0` | Gammas |

Correlated across channels (same names): `Lumi`, `Pileup`, `L1Prefiring`, `XS_*`, `PDF`,
`AlphaS`, `QCDScale`, `PS_*`, `SigModel`, and the muon NPs with the mu tau_h channel.

## Fit and extraction

`trex-fitter h w f dp r i` (and `wf "StatOnly=TRUE:Suffix=_statOnly"`); `FitType SPLUSB`,
`FitRegion CRSR`, MINOS on `mu_Z`. Then

```
sigma_fid           = mu_Z x sigma_fid^pred(NLO)
sigma(60 < m < 120) = sigma_fid / A_60_120       sigma(m > 50) = sigma_fid / A_m50
```

with the statistical uncertainty from the stat-only fit, the systematic uncertainty from the
grouped impacts (`FullSyst`; the luminosity is quoted as the external 1.2%, not the profiled
impact), and the acceptance uncertainty added in quadrature for the inclusive cross sections:
0.70% for 60-120 GeV (PDF 0.50%, boson pT 0.38%, scales 0.29%, QED FSR 0.10% estimated, MC stat,
generator, alpha_s, PS FSR < 0.05% each; `zmumu/acceptance.py`, docs/16 §4-5). Until 16 Sep it was
PDF, scale, alpha_s and MC stat for the m > 50 GeV denominator (0.61%, `scripts/mc_acceptance.py`). The counting form `(N_data - N_bkg)/(L C)` is
printed as a cross-check. Results: `fit/results/zmumu_fit_result.json`, `output/v2/RESULTS_v2.md`.

## Why 5 GeV bins, and the stability of the result

The first v2 fit used 30 bins of 2 GeV, smoothed `MuonScale`/`MuonRes` templates and a
one-sided powheg `SigModel` built without the generator window. It was not robust (REVIEW.md
F3-F5): with 10 M events in 2 GeV bins the fit resolves shape differences of 10^-3 per bin,
finer than the physics content and the MC statistics of the variation templates, so
`SigModel`, `MuonScale`, `MuonRes` were constrained to 0.1-0.2 of their priors, `Lumi` to 0.8
and `Pileup` was pulled by -1.5 sigma; mu_Z moved by +-0.8% between binnings and with/without
`SigModel`. The data show a 3-4% deficit against aMC@NLO at 62-78 GeV that only the powheg
template can describe, and the old template also contained the edge migration from
m_LHE > 120 GeV, which the powheg sample lacks.

The fixed configuration keeps a shape fit (the task asks for a binned-likelihood fit of
m(mumu)) but in 5 GeV bins, without smoothing, with a two-sided `SigModel` built inside the
common generator window. `scripts/v2_5_fit_variants.py` repeats the fit with 2 and 10 GeV
bins, one bin (counting), without `SigModel` and with the old smoothing; the table is in
`RESULTS_v2.md` and summarised below.

## Result (frozen 17 Sep 2026: the 15 Sep fit with the measured reconstruction SF and the docs/16 acceptance)

```
mu_Z      = 0.9883 +0.0141 -0.0138        (stat 0.0003, syst 0.0139 of which lumi 0.0120 external / 0.0117 profiled)
sigma_fid = 790.2 +- 0.2 (stat) +- 6.1 (syst) +- 9.6 (lumi) pb       GoF p = 0.79 (12 bins)
sigma(60 < m < 120) = 1931 +- 0.6 (stat) +- 15.0 (syst) +- 13.5 (acc) +- 23.4 (lumi) = 1931 +- 30 pb
sigma(m > 50)       = 2002 +- 32 pb
counting: (10 378 567 - 68 799) / (16393.381 x 0.7916) = 794.5 pb
```

15 Sep (reconstruction SF 1 +- 0.4%/muon assigned, acceptance for m > 50 GeV): mu_Z = 0.9881 +0.0159 -0.0156,
sigma_fid = 790.1 +- 0.2 +- 8.5 +- 9.6 pb, sigma(60-120) = 1931 +- 33 pb (`fit/results/zmumu_v2_15sep_fit_result.json`).

| group | impact on mu_Z |
|---|---:|
| luminosity (external) | 1.20% (profiled 1.18%) |
| L1 prefiring | 0.51% |
| muon efficiency (of which `MuonReco` 0.27%; 0.77% on 15 Sep) | 0.48% |
| signal modelling | 0.46% |
| MC statistics (gammas) | 0.36% |
| muon momentum | 0.35% |
| pileup | 0.15% |
| background normalisation | 0.07% |
| fakes | 0.04% |
| statistical | 0.03% |
| **total systematic** | **1.41%** (0.77% without the luminosity; 1.57% / 1.06% on 15 Sep) |

Stability (`fit/results/stability.json`, MINOS on every parameter):

| configuration | mu_Z | GoF p |
|---|---:|---:|
| **12 x 5 GeV (nominal)** | **0.9883** | 0.79 |
| 30 x 2 GeV | 0.9900 | 0.16 |
| 6 x 10 GeV | 0.9862 | 0.50 |
| 1 bin (= counting) | 0.9937 | -- |
| 12 x 5 GeV without `SigModel` | 1.0066 | 0.001 |
| 12 x 5 GeV, `MuonScale`/`MuonRes` smoothed (old setup) | 1.0011 | 0.014 |
| 30 x 2 GeV without `SigModel` | 1.0048 | 0.0002 |

The binnings that describe the data (p > 0.05) agree within +-0.2%; the counting extraction is
0.55% higher because it does not use the shape information that pulls `SigModel` (+0.66 sigma:
the data prefer a lineshape between aMC@NLO and powheg, `sigmodel_lineshape.png`), `MuonScale`
(-0.48) and `MuonRes` (+0.40). Removing `SigModel` or smoothing the momentum templates gives
fits that do not describe the data (p <= 0.01) and should not be used. Pulls of all other
parameters are within +-0.8 sigma; `Pileup` is now +0.06 (was -1.5 before the N_PV-matched
profile). Constraints: only the three shape parameters (`SigModel` 0.14, `MuonScale` 0.14,
`MuonRes` 0.12) are constrained by the 10 M-event spectrum; `Lumi`, `MuonReco`, `L1Prefiring`
stay at their priors (the first v2 fit reported 0.2-0.8 for them because the numerical Hesse
covariance was ill-conditioned; MINOS fixes that, `UseMinos: all`).

Momentum calibration and the previous v2 result (791.8 +- 6.2 +- 9.3 pb, 30 x 2 GeV, one-sided
`SigModel`, `MuonReco` 0.4%, CSV pileup profile) are in the git history (commit fc9a998) and in
REVIEW.md. Comparison: reviewer 797.2 pb (-0.9%), v1 773.2 pb, aMC@NLO 799.6 pb (ratio 0.988),
madgraph LO 825.4 pb.
