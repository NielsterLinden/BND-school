# Z -> mu+ mu- cross section, v2 (MC-based, TRExFitter fit)

Generated 2026-09-17 from CMS 2016 Open Data (SingleMuon, Run2016G+H, 16.393 fb^-1, normtag).

## Result

**sigma_fid(pp -> Z/gamma* -> mu mu; dressed, pT > 26/20 GeV, |eta| < 2.4, 60 < m < 120 GeV) = 790.2 +- 0.2 (stat) +- 6.1 (syst) +- 9.6 (lumi) pb**

mu_Z = 0.9883 +0.0141 -0.0138; goodness of fit p = 0.792272; signal region in 5 GeV bins

The statistical uncertainty is the data one; the effective statistical limit of the fit is the MC statistics (gammas, see the breakdown), not the data.

| quantity | value |
|---|---:|
| sigma(Z/gamma* -> mu mu, 60 < m < 120 GeV), A = 0.4092 | 1931 +- 0.6 (stat) +- 15.0 (syst) +- 13.5 (acc) +- 23.4 (lumi) = +- 30 pb |
| sigma(Z/gamma* -> mu mu, m > 50 GeV), A = 0.3947 | 2002 +- 0.6 (stat) +- 15.5 (syst) +- 14.5 (acc) +- 24.3 (lumi) = +- 32 pb |
| C factor (reco/fiducial, all corrections) | 0.7916 |
| NLO prediction sigma_fid (6077.22 pb x fiducial fraction) | 799.6 pb |
| counting cross-check (N_obs - N_bkg)/(C L) | 794.5 pb |
| muon reconstruction SF applied to the simulation (T&P, step 7) | 1.0001 per muon, 1.0002 per event |
| v1 (data-only counting) / v1 revised / reviewer | 773.2 / 776.9 / 797.2 pb |

## Uncertainty breakdown (impact on mu_Z, from the grouped-impact fit)

| group | relative |
|---|---:|
| FullSyst | 1.395% |
| Luminosity | 1.166% |
| L1 prefiring | 0.508% |
| Muon efficiency | 0.475% |
| Signal modelling | 0.454% |
| Gammas | 0.353% |
| Muon momentum | 0.350% |
| Pileup | 0.143% |
| Background normalisation | 0.072% |
| Fakes | 0.038% |
| Electron efficiency | 0.000% |
| statistical (stat-only fit) | 0.031% |

Luminosity quoted as the external 1.2% (profiled impact 1.166%).

## Acceptance uncertainties (outside the fit; zmumu/acceptance.py, docs/16)

| source | A(60-120) | A(m>50) |
|---|---:|---:|
| PDF (NNPDF3.1 Hessian) | 0.500% | 0.516% |
| QCD scales (7-point) | 0.290% | 0.315% |
| boson pT: spectrum reweighted to the measured pT(mumu) | 0.375% | 0.375% |
| QED FSR model (estimate, option c) | 0.101% | 0.101% |
| powheg vs aMC@NLO | 0.045% | 0.045% |
| alpha_s | 0.031% | 0.028% |
| PS FSR | 0.027% | 0.028% |
| MC statistics | 0.048% | 0.048% |
| **total** | **0.701%** | **0.723%** |

PS ISR (0.068%) is smaller than the boson-pT row and not added. QED FSR options: a 0.12% (CMS-SMP-20-004) / b 0 / c 0.10%.

## Stability of mu_Z against the fit configuration (scripts/v2_5_fit_variants.py)

| configuration | mu_Z | GoF p | notable pulls (constraint) |
|---|---:|---:|---|
| nominal: 12 x 5 GeV, SigModel two-sided (LHE window), no smoothing | 0.9883 +0.0141 -0.0138 | 0.792272 | MuonRes +0.4 (0.12); MuonScale -0.5 (0.14); SigModel +0.7 (0.14) |
| 30 x 2 GeV | 0.9900 +0.0139 -0.0135 | 0.161256 | MuonRes +0.2 (0.07); MuonScale -0.5 (0.10); PDF +1.4 (0.87); PS_FSR -0.6 (0.52); Pileup -0.4 (0.67) |
| 6 x 10 GeV | 0.9862 +0.0142 -0.0139 | 0.501913 | MuonRes +0.6 (0.53); MuonScale -0.4 (0.29); SigModel +0.7 (0.15) |
| 1 bin (counting) | 0.9937 +0.0149 -0.0146 | 0.0 |  |
| 12 x 5 GeV, no SigModel | 1.0066 +0.0139 -0.0136 | 0.00136029 | MuonIso -2.3 (0.76); MuonRes +0.2 (0.11); MuonScale +0.1 (0.06); MuonTrigger -1.2 (0.94); PS_FSR +1.0 (0.72) |
| 12 x 5 GeV, MuonScale/MuonRes smoothed | 1.0011 +0.0141 -0.0137 | 0.0140223 | MuonIso -1.5 (0.85); MuonRes +1.3 (0.77); MuonScale +0.2 (0.33); PS_FSR +1.1 (0.71); QCDScale -1.0 (0.68) |
| 30 x 2 GeV, no SigModel | 1.0048 +0.0138 -0.0135 | 0.000208255 | MuonIso -3.0 (0.73); MuonRes +0.1 (0.07); MuonScale -0.0 (0.05); MuonTrigger -1.5 (0.94); PDF +1.2 (0.87) |

Spread of mu_Z over the shape fits that describe the data (GoF p > 0.05): 0.9862 - 0.9900 (half-spread 0.19%). The 1-bin fit has no goodness of fit (0 degrees of freedom). Rejected (p <= 0.05): 12 x 5 GeV, no SigModel; 12 x 5 GeV, MuonScale/MuonRes smoothed; 30 x 2 GeV, no SigModel.

## Leading nuisance parameters (ranking)

| NP | pull | constraint | +impact | -impact |
|---|---:|---:|---:|---:|
| Lumi | +0.08 | 1.00 | -1.152% | +1.180% |
| MuonReco | +0.02 | 1.00 | -0.260% | +0.263% |
| Pileup | +0.06 | 0.82 | +0.144% | -0.143% |
| L1Prefiring | -0.06 | 1.00 | +0.514% | -0.505% |
| MuonID | +0.07 | 1.00 | -0.188% | +0.191% |
| MuonIso | +0.01 | 0.90 | -0.294% | +0.298% |
| MuonTrigger | +0.01 | 0.97 | -0.147% | +0.152% |
| MuonScale | -0.48 | 0.14 | +0.336% | -0.333% |
| MuonRes | +0.40 | 0.12 | -0.154% | +0.156% |
| PDF | +0.56 | 0.93 | +0.072% | -0.070% |
| AlphaS | +0.00 | 1.00 | +0.023% | -0.022% |
| QCDScale | +0.76 | 0.79 | -0.096% | +0.101% |

## Yields (mass_fit templates; SR/SS: prompt-prompt MC + data-driven fakes; e-mu regions: all MC)

| region | sample | events |
|---|---|---:|
| SR | Data | 10,378,567.0 |
| SR | DYmumu | 10,375,538.3 |
| SR | TTbar | 31,613.3 |
| SR | DYtautau | 11,025.0 |
| SR | WZ | 9,196.6 |
| SR | ZZ | 6,309.9 |
| SR | Fakes | 3,870.0 |
| SR | WW | 3,848.0 |
| SR | SingleTop | 2,936.2 |
| SS | Data | 2,302.0 |
| SS | WZ | 202.1 |
| SS | DYmumu | 115.7 |
| SS | ZZ | 38.6 |
| SS | TTbar | 0.9 |
| SS | SingleTop | 0.8 |
| SS | WW | 0.2 |
| SS | DYtautau | -1.6 |
| CRemu | Data | 72,357.0 |
| CRemu | TTbar | 43,169.4 |
| CRemu | DYtautau | 13,915.7 |
| CRemu | WW | 5,274.1 |
| CRemu | SingleTop | 4,112.9 |
| CRemu | WJets | 3,878.1 |
| CRemu | DYmumu | 2,220.2 |
| CRemu | WZ | 298.0 |
| CRemu | ZZ | 62.2 |
| CRemu | DYee | 35.4 |
| SSemu | Data | 6,718.0 |
| SSemu | WJets | 2,543.8 |
| SSemu | DYmumu | 1,791.8 |
| SSemu | TTbar | 410.5 |
| SSemu | WZ | 292.6 |
| SSemu | DYtautau | 186.9 |
| SSemu | SingleTop | 63.3 |
| SSemu | ZZ | 52.1 |
| SSemu | WW | 49.3 |
| SSemu | DYee | 38.1 |

Not listed: `DYmumu_powheg` (raw generator weights; SigModel template only); `WJets` (prompt-prompt part is a negative-weight fluctuation; e-mu regions only); `DYother` (LHE flavour not e/mu/tau).

## Fake factor

SR non-prompt estimate 3870 +- 80 (stat); method uncertainty 19% (largest: mcsub_down).
Same-sign closure: predicted 2094 +- 20, observed 1932 +- 79.

## Tag-and-probe

Pairs: 20,050,129 (data). Trigger plateau efficiency per muon: data 0.9071, MC 0.9227.
Isolation SF vs pileup half-spread: 0.0031.

| efficiency | <SF> | <stat> | <syst> |
|---|---:|---:|---:|
| id | 0.9799 | 0.0015 | 0.0031 |
| iso | 1.0059 | 0.0009 | 0.0041 |
| antiiso | 0.8932 | 0.0337 | 0.1349 |

## Momentum calibration (MC pT scale kappa and extra smearing per |eta| bin)

| 0.0-0.9 | 0.9-1.2 | 1.2-2.1 | 2.1-2.4 |
|---:|---:|---:|---:|
| 0.99913 +- 0.00050, 0.00% | 0.99929 +- 0.00050, 0.59% | 0.99876 +- 0.00050, 0.61% | 1.00039 +- 0.00050, 0.78% |
