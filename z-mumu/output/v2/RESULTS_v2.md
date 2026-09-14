# Z -> mu+ mu- cross section, v2 (MC-based, TRExFitter fit)

Generated 2026-09-14 from CMS 2016 Open Data (SingleMuon, Run2016G+H, 16.393 fb^-1, normtag).

## Result

**sigma_fid(pp -> Z/gamma* -> mu mu; dressed, pT > 26/20 GeV, |eta| < 2.4, 60 < m < 120 GeV) = 791.8 +- 0.2 (stat) +- 6.2 (syst) +- 9.3 (lumi) pb**

mu_Z = 0.9903 +0.0141 -0.0138; goodness of fit p = 0.32502

| quantity | value |
|---|---:|
| sigma(Z/gamma* -> mu mu, 60 < m < 120 GeV), A = 0.4092 | 1935 +- 30 pb |
| sigma(Z/gamma* -> mu mu, m > 50 GeV), A = 0.3947 | 2006 +- 31 pb |
| C factor (reco/fiducial, all corrections) | 0.7917 |
| NLO prediction sigma_fid (6077.22 pb x fiducial fraction) | 799.6 pb |
| counting cross-check (N_obs - N_bkg)/(C L) | 794.4 pb |
| v1 (data-only counting) / v1 revised / reviewer | 773.2 / 776.9 / 797.2 pb |

## Uncertainty breakdown (impact on mu_Z, from the grouped-impact fit)

| group | relative |
|---|---:|
| FullSyst | 1.397% |
| Luminosity | 1.164% |
| L1 prefiring | 0.506% |
| Muon efficiency | 0.503% |
| Signal modelling | 0.363% |
| Gammas | 0.347% |
| Muon momentum | 0.268% |
| Pileup | 0.158% |
| Background normalisation | 0.093% |
| Fakes | 0.045% |
| Electron efficiency | 0.000% |
| statistical (stat-only fit) | 0.031% |

## Leading nuisance parameters (ranking)

| NP | pull | constraint | +impact | -impact |
|---|---:|---:|---:|---:|
| Lumi | +0.08 | 0.82 | -0.949% | +0.970% |
| MuonReco | +0.03 | 0.97 | -0.372% | +0.376% |
| Pileup | -1.51 | 0.69 | +0.157% | -0.156% |
| L1Prefiring | -0.19 | 1.00 | +0.512% | -0.507% |
| MuonID | +0.20 | 0.99 | -0.177% | +0.180% |
| MuonIso | -1.10 | 0.82 | -0.211% | +0.213% |
| MuonTrigger | -0.67 | 0.94 | -0.107% | +0.106% |
| MuonScale | -0.46 | 0.11 | +0.254% | -0.257% |
| MuonRes | +0.94 | 0.21 | -0.184% | +0.184% |
| PDF | +1.32 | 0.87 | +0.093% | -0.094% |
| AlphaS | -0.07 | 0.99 | +0.029% | -0.028% |
| QCDScale | +0.15 | 0.51 | +0.004% | -0.004% |

## Yields (mass_fit templates, prompt-prompt MC + data-driven fakes)

| region | sample | events |
|---|---|---:|
| SR | DYmumu_powheg | 1,980,364,987.9 |
| SR | Data | 10,378,567.0 |
| SR | DYmumu | 10,377,013.2 |
| SR | TTbar | 31,592.8 |
| SR | DYtautau | 11,112.8 |
| SR | WZ | 9,199.9 |
| SR | ZZ | 6,307.5 |
| SR | WW | 3,844.6 |
| SR | Fakes | 3,824.7 |
| SR | SingleTop | 2,938.3 |
| SR | WJets | -5.0 |
| SS | DYmumu_powheg | 31,752.2 |
| SS | Data | 2,302.0 |
| SS | WZ | 201.9 |
| SS | DYmumu | 114.3 |
| SS | ZZ | 38.6 |
| SS | WJets | 13.5 |
| SS | SingleTop | 0.8 |
| SS | TTbar | 0.8 |
| SS | WW | 0.2 |
| SS | DYtautau | -1.4 |
| CRemu | Data | 72,357.0 |
| CRemu | TTbar | 42,979.8 |
| CRemu | DYmumu_powheg | 26,670.1 |
| CRemu | DYtautau | 13,317.8 |
| CRemu | WW | 5,265.6 |
| CRemu | SingleTop | 4,056.8 |
| CRemu | WZ | 292.8 |
| CRemu | DYmumu | 206.5 |
| CRemu | WJets | 147.2 |
| CRemu | ZZ | 59.9 |
| CRemu | DYee | 22.8 |

## Fake factor

SR non-prompt estimate 3825 +- 80 (stat); method uncertainty 19% (largest: mcsub_down).
Same-sign closure: predicted 2077 +- 20, observed 1922 +- 85.

## Tag-and-probe

Pairs: 20,050,129 (data). Trigger plateau efficiency per muon: data 0.9071, MC 0.9234.
Isolation SF vs pileup half-spread: 0.0024.

| efficiency | <SF> | <stat> | <syst> |
|---|---:|---:|---:|
| id | 0.9797 | 0.0014 | 0.0029 |
| iso | 1.0054 | 0.0009 | 0.0033 |
| antiiso | 0.8996 | 0.0296 | 0.1089 |

## Momentum calibration (MC pT scale kappa and extra smearing per |eta| bin)

| 0.0-0.9 | 0.9-1.2 | 1.2-2.1 | 2.1-2.4 |
|---:|---:|---:|---:|
| 0.99911 +- 0.00050, 0.00% | 0.99923 +- 0.00050, 0.62% | 0.99874 +- 0.00050, 0.63% | 1.00038 +- 0.00050, 0.91% |
