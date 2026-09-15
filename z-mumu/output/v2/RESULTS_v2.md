# Z -> mu+ mu- cross section, v2 (MC-based, TRExFitter fit)

Generated 2026-09-15 from CMS 2016 Open Data (SingleMuon, Run2016G+H, 16.393 fb^-1, normtag).

## Result

**sigma_fid(pp -> Z/gamma* -> mu mu; dressed, pT > 26/20 GeV, |eta| < 2.4, 60 < m < 120 GeV) = 790.1 +- 0.2 (stat) +- 8.5 (syst) +- 9.6 (lumi) pb**

mu_Z = 0.9881 +0.0159 -0.0156; goodness of fit p = 0.792273; signal region in 5 GeV bins

The statistical uncertainty is the data one; the effective statistical limit of the fit is the MC statistics (gammas, see the breakdown), not the data.

| quantity | value |
|---|---:|
| sigma(Z/gamma* -> mu mu, 60 < m < 120 GeV), A = 0.4092 | 1931 +- 33 pb |
| sigma(Z/gamma* -> mu mu, m > 50 GeV), A = 0.3947 | 2002 +- 34 pb |
| C factor (reco/fiducial, all corrections) | 0.7914 |
| NLO prediction sigma_fid (6077.22 pb x fiducial fraction) | 799.6 pb |
| counting cross-check (N_obs - N_bkg)/(C L) | 794.7 pb |
| v1 (data-only counting) / v1 revised / reviewer | 773.2 / 776.9 / 797.2 pb |

## Uncertainty breakdown (impact on mu_Z, from the grouped-impact fit)

| group | relative |
|---|---:|
| FullSyst | 1.573% |
| Luminosity | 1.163% |
| Muon efficiency | 0.869% |
| L1 prefiring | 0.506% |
| Signal modelling | 0.491% |
| Gammas | 0.382% |
| Muon momentum | 0.367% |
| Pileup | 0.134% |
| Background normalisation | 0.081% |
| Fakes | 0.044% |
| Electron efficiency | 0.000% |
| statistical (stat-only fit) | 0.031% |

Luminosity quoted as the external 1.2% (profiled impact 1.163%).

## Stability of mu_Z against the fit configuration (scripts/v2_5_fit_variants.py)

| configuration | mu_Z | GoF p | notable pulls (constraint) |
|---|---:|---:|---|
| nominal: 12 x 5 GeV, SigModel two-sided (LHE window), no smoothing | 0.9881 +0.0159 -0.0156 | 0.792273 | MuonRes +0.4 (0.12); MuonScale -0.5 (0.14); SigModel +0.7 (0.14) |
| 30 x 2 GeV | 0.9901 +0.0157 -0.0153 | 0.161244 | MuonRes +0.2 (0.07); MuonScale -0.5 (0.10); PDF +1.4 (0.87); PS_FSR -0.6 (0.52); Pileup -0.4 (0.67) |
| 6 x 10 GeV | 0.9858 +0.0160 -0.0156 | 0.502138 | MuonRes +0.6 (0.53); MuonScale -0.4 (0.29); SigModel +0.7 (0.15) |
| 1 bin (counting) | 0.9939 +0.0168 -0.0164 | 0.0 |  |
| 12 x 5 GeV, no SigModel | 1.0073 +0.0158 -0.0154 | 0.00136417 | MuonIso -2.3 (0.76); MuonRes +0.2 (0.11); MuonScale +0.1 (0.06); MuonTrigger -1.2 (0.94); PS_FSR +1.0 (0.72) |
| 12 x 5 GeV, MuonScale/MuonRes smoothed | 1.0012 +0.0160 -0.0156 | 0.014019 | MuonIso -1.5 (0.85); MuonRes +1.3 (0.77); MuonScale +0.2 (0.33); PS_FSR +1.1 (0.71); QCDScale -1.0 (0.68) |
| 30 x 2 GeV, no SigModel | 1.0056 +0.0157 -0.0153 | 0.000209046 | MuonIso -3.0 (0.73); MuonRes +0.1 (0.07); MuonScale -0.0 (0.05); MuonTrigger -1.5 (0.94); PDF +1.2 (0.87) |

Spread of mu_Z over the shape fits that describe the data (GoF p > 0.05): 0.9858 - 0.9901 (half-spread 0.22%). The 1-bin fit has no goodness of fit (0 degrees of freedom). Rejected (p <= 0.05): 12 x 5 GeV, no SigModel; 12 x 5 GeV, MuonScale/MuonRes smoothed; 30 x 2 GeV, no SigModel.

## Leading nuisance parameters (ranking)

| NP | pull | constraint | +impact | -impact |
|---|---:|---:|---:|---:|
| Lumi | +0.08 | 1.00 | -1.151% | +1.175% |
| MuonReco | +0.05 | 1.00 | -0.771% | +0.778% |
| Pileup | +0.06 | 0.82 | +0.129% | -0.137% |
| L1Prefiring | -0.06 | 1.00 | +0.507% | -0.502% |
| MuonID | +0.07 | 1.00 | -0.190% | +0.185% |
| MuonIso | +0.01 | 0.90 | -0.286% | +0.288% |
| MuonTrigger | +0.01 | 0.97 | -0.150% | +0.147% |
| MuonScale | -0.48 | 0.14 | +0.353% | -0.356% |
| MuonRes | +0.40 | 0.12 | -0.153% | +0.151% |
| PDF | +0.56 | 0.93 | +0.085% | -0.088% |
| AlphaS | +0.00 | 1.00 | +0.020% | -0.025% |
| QCDScale | +0.76 | 0.79 | -0.084% | +0.078% |

## Yields (mass_fit templates; SR/SS: prompt-prompt MC + data-driven fakes; e-mu regions: all MC)

| region | sample | events |
|---|---|---:|
| SR | Data | 10,378,567.0 |
| SR | DYmumu | 10,373,541.7 |
| SR | TTbar | 31,607.2 |
| SR | DYtautau | 11,022.9 |
| SR | WZ | 9,194.9 |
| SR | ZZ | 6,308.7 |
| SR | Fakes | 3,870.0 |
| SR | WW | 3,847.3 |
| SR | SingleTop | 2,935.7 |
| SS | Data | 2,302.0 |
| SS | WZ | 202.1 |
| SS | DYmumu | 115.7 |
| SS | ZZ | 38.6 |
| SS | TTbar | 0.9 |
| SS | SingleTop | 0.8 |
| SS | WW | 0.2 |
| SS | DYtautau | -1.6 |
| CRemu | Data | 72,357.0 |
| CRemu | TTbar | 43,165.3 |
| CRemu | DYtautau | 13,914.4 |
| CRemu | WW | 5,273.6 |
| CRemu | SingleTop | 4,112.5 |
| CRemu | WJets | 3,877.7 |
| CRemu | DYmumu | 2,220.0 |
| CRemu | WZ | 297.9 |
| CRemu | ZZ | 62.2 |
| CRemu | DYee | 35.4 |
| SSemu | Data | 6,718.0 |
| SSemu | WJets | 2,543.6 |
| SSemu | DYmumu | 1,791.6 |
| SSemu | TTbar | 410.5 |
| SSemu | WZ | 292.6 |
| SSemu | DYtautau | 186.9 |
| SSemu | SingleTop | 63.2 |
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
