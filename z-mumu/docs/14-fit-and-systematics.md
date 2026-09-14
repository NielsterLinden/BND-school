# Signal extraction, systematics, and the cross section

`zmumu/histograms.py`, `zmumu/momentum.py`, `scripts/v2_4_histograms.py`, `scripts/v2_5_fit.py`,
`fit/zmumu.config`; fit outputs in `fit/results/zmumu/` (TRExFitter v1.8.0).

## Regions and templates

| region | selection | role |
|---|---|---|
| `mumu_SR` | trigger, MET filters, exactly two tight muons (pT > 26 / 20, \|eta\| < 2.4), >= 1 trigger-matched, OS, FSR-recovered 60 < m < 120 GeV; 30 bins of 2 GeV | signal region (fit) |
| `mumu_CRemu` | one tight matched muon (pT > 26), no other tight/anti-tight muon, one electron (pT > 20, \|eta\| < 2.5, cutBased medium, not in the gap), OS, 60 < m(e mu) < 120; 12 bins | validation of the flavour-symmetric backgrounds (`--emu-control` promotes it to a control region with `mu_top`) |
| `SS` | as SR, same sign | fake-factor closure, plots |

Samples: `Data`, `DYmumu` (signal, `NormFactor mu_Z`), `DYtautau`, `DYee`, `TTbar`, `SingleTop`,
`WW`, `WZ`, `ZZ` (prompt-prompt events only) and the data-driven `Fakes`.

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
| `MuonReco` | norm 0.4% | all MC | not measurable in NanoAOD | Muon efficiency |
| `MuonScale`, `MuonRes` | shape (smoothed) | all MC | Z-peak calibration | Muon momentum |
| `PDF`, `AlphaS`, `QCDScale`, `PS_ISR`, `PS_FSR` | shape | `DYmumu` | LHE/PS weights, renormalised to the fiducial yield (C factor only) | Signal modelling |
| `SigModel` | one-sided shape | `DYmumu` | powheg template normalised to the NLO fiducial prediction | Signal modelling |
| `XS_TTbar` 6%, `XS_SingleTop` 10%, `XS_WW/WZ/ZZ` 10%, `XS_DYtautau` 5% | norm | backgrounds | theory | Background normalisation |
| `FakeStat_mumu`, `FakeMethod_mumu` | shape/norm | `Fakes` | docs/13 | Fakes |
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
grouped impacts (`FullSyst`, luminosity quoted separately), and the acceptance uncertainty
(PDF 0.52%, scale 0.32%, alpha_s 0.03%, MC stat; from `scripts/mc_acceptance.py`) added in
quadrature for the inclusive cross sections. The counting form `(N_data - N_bkg)/(L C)` is
printed as a cross-check. Results: `fit/results/zmumu_fit_result.json`, `output/v2/RESULTS_v2.md`.

## Result (14 Sep 2026)

```
mu_Z      = 0.9903 +0.0141 -0.0138        (stat 0.0003, syst 0.0140 of which lumi 0.0116)
sigma_fid = 791.8 +- 0.2 (stat) +- 6.2 (syst) +- 9.3 (lumi) pb       GoF p = 0.33
sigma(60 < m < 120) = 1935 +- 30 pb        sigma(m > 50) = 2006 +- 31 pb
counting: (10 378 567 - 68 821) / (16393.381 x 0.7917) = 794.4 pb
```

| group | impact on mu_Z |
|---|---:|
| luminosity | 1.164% |
| L1 prefiring | 0.506% |
| muon efficiency | 0.503% |
| signal modelling | 0.363% |
| MC statistics (gammas) | 0.347% |
| muon momentum | 0.268% |
| pileup | 0.158% |
| background normalisation | 0.093% |
| fakes | 0.045% |
| statistical | 0.031% |
| **total systematic** | **1.397%** |

Momentum calibration: kappa = 0.99911, 0.99923, 0.99874, 1.00038 (+-0.0005) and extra
smearing 0, 0.62, 0.63, 0.91% in the four |eta| bins. Pulls: `Pileup` -1.5 sigma (the data
prefer a lower minimum-bias cross section than 69.2 mb), `MuonIso` -1.1, `PDF` +1.3,
`MuonScale`/`MuonRes` strongly constrained by the 2 GeV mass bins; all others within 1 sigma.
The MINOS uncertainty on `mu_Z` is symmetric within 2%. Comparison: reviewer 797.2 pb (same
volume, -0.7%), v1 773.2 pb, v1 revised 776.9 pb, aMC@NLO 799.6 pb (ratio 0.990), madgraph LO
825.4 pb.
