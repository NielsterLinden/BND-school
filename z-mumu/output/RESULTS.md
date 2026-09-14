# Z -> mu+ mu- fiducial cross section -- results

Generated 2026-09-14 from CMS 2016 Open Data (SingleMuon, Run2016G+H).

## Result

**sigma_fid = 773.2 +/- 0.2 (stat) +/- 11.9 (syst) pb**

Total uncertainty 1.54% (0.03% stat, 1.54% syst).

## Fiducial volume

| Requirement | Value |
|---|---|
| muons | exactly 2, opposite sign |
| leading muon pT | > 26.0 GeV |
| subleading muon pT | > 20.0 GeV |
| muon abs(eta) | < 2.4 |
| muon ID | Muon_mediumId |
| muon isolation | pfRelIso04_all < 0.15 |
| dimuon mass | 60.0-120.0 GeV |

Acceptance is 1 by construction: this is the phase space that was measured.

## Cutflow

| Cut | Events | Fraction |
|---|---:|---:|
| skim (>=2 muons in collection) | 80,191,719 | 100.000% |
| golden JSON (certified lumi) | 78,912,289 | 98.405% |
| HLT_IsoMu24 or HLT_IsoTkMu24 | 47,314,794 | 59.002% |
| MET filters + good PV | 47,292,438 | 58.974% |
| exactly 2 selected muons | 11,647,158 | 14.524% |
| leading muon pT > 26 GeV | 11,477,398 | 14.312% |
| opposite sign | 11,472,359 | 14.306% |
| 60 < m(mumu) < 120 GeV | 10,777,373 | 13.440% |

## Efficiencies

| Term | Value | Source |
|---|---:|---|
| reconstruction (per muon) | 0.9960 | external (CMS Muon POG) -- not measurable in NanoAOD |
| medium ID (inclusive) | 0.9775 | tag-and-probe |
| isolation (inclusive) | 0.9495 | tag-and-probe |
| trigger (per event) | 0.9977 | reference-trigger method |
| **total (per event)** | **0.8546** | product, kinematics-weighted |

Tag-and-probe pairs used: 21,510,663.

## Yields

| Quantity | Events |
|---|---:|
| observed (opposite sign) | 10,777,373 |
| same-sign control | 2,375 |
| e-mu control | 23,471 |
| background: non-prompt | 2,375.0 |
| background: flavour-symmetric | 11,735.5 |
| **background: total** | **14,110.5** (0.131%) |
| **signal (counting, nominal)** | **10,763,262** |
| signal (fit, cross-check only) | 10,183,155 |

## Uncertainty breakdown

| Source | Relative |
|---|---:|
| luminosity | 1.200% |
| eff: reconstruction (external) | 0.803% |
| eff: trigger (method bias) | 0.500% |
| muon momentum scale | 0.200% |
| background estimate | 0.056% |
| statistical (data) | 0.031% |
| eff: isolation (T&P stat) | 0.008% |
| eff: ID (tag-and-probe stat) | 0.005% |
| FSR recovery | 0.002% |
| eff: trigger (stat) | 0.001% |
| **total** | **1.542%** |

## Fit validation

- Fitted peak: 90.789 +/- 0.026 GeV (PDG 91.1876)
- Fitted resolution: 1.448 GeV
- chi2/ndf = 104469/115 = 908.42 (poor by construction: a single Voigt profile has no FSR tail)

See `docs/` for the method behind each number.
