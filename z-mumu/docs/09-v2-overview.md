# v2 overview -- the MC-based measurement with a TRExFitter fit

`run_v2.py`, `zmumu/`, `scripts/v2_*.py`, `fit/`. Result: `output/v2/RESULTS_v2.md`.

## Why a second version

The data-only measurement (steps 1-6, docs 00-07) was reviewed twice (docs/08 and
`agent_reference/2026.09.14_1400_revisionpoints.md`). The review points that cannot be fixed
inside the data-only design are:

| issue | v1 | v2 |
|---|---|---|
| luminosity without normtag (0.63% low) | 16290.7 pb^-1 | **16393.381 pb^-1** (record 1059, reproduced from the by-LS table) |
| L1 pre-firing (2016: ~2%) not correctable without MC | not applied | `L1PreFiringWeight` on MC (Muon + ECAL, Up/Down) |
| e-mu region from the >= 2-muon skim (x3.2 too few events) | 23 471 events | own skim of the unskimmed parent, 1-muon events kept |
| WZ / ZZ not subtracted | absorbed in a systematic | simulated (WZ -> 3l nu, 2q 2l; ZZ -> 4l, 2l 2nu, 2q 2l) |
| tag-and-probe by counting, tag not trigger-matched | 0.9%/muon bias (tight ID) | simultaneous pass/fail **fits** with MC templates, matched tags |
| trigger efficiency from correlated reference paths | biased high | trigger-object tag-and-probe in fine pT bins |
| "A = 1 by construction" ignores migrations | no C factor | **C factor from DY aMC@NLO** (all corrections applied), A for the combination |
| counting, no likelihood | | **TRExFitter v1.8.0 binned likelihood fit** of m(mu mu) with all systematics as nuisance parameters |
| medium ID | | tight ID (reviewer's convention, CMS standard) |

## The measurement in one line

```
   sigma_fid = mu_Z x sigma_fid^pred        sigma_fid^pred = 6077.22 pb x N_fid^gen / N_all^gen  (aMC@NLO)
```

`mu_Z` is the signal-strength normalisation factor of the DY -> mu mu template in the fit of the
dimuon mass spectrum (60-120 GeV, 2 GeV bins) to the data. Because the template is normalised
to `sigma_th x L x N_sel/N_all`, the product is algebraically `(N_data - N_bkg)/(L x C)` with
`C = N_sel(weighted)/N_fid`: the assumed DY cross section cancels, and the fit only adds the
constraint of the backgrounds and of the nuisance parameters. The counting form is printed
next to the fit as a cross-check.

Fiducial volume (unchanged from v1, but now with a real C factor): exactly two opposite-sign
muons, pT > 26 / 20 GeV, |eta| < 2.4, 60 < m < 120 GeV, **dressed** leptons (GenDressedLepton,
photons within dR < 0.1). The inclusive cross sections follow by dividing by the acceptance
A(60 < m_Born < 120) = 0.4092 (recommended for the combination) or A(m > 50) = 0.3947.

## The chain

| step | script | what it does | doc |
|---|---|---|---|
| 1 | `v2_1_skim.py` | skims of the unskimmed SingleMuon parent and of every MC sample | [10](10-skims.md) |
| 2 | `v2_2_tnp.py` | pileup weights; ID / iso / trigger / anti-iso efficiencies from T&P fits; scale factors | [11](11-mc-weights.md), [12](12-tag-and-probe-fits.md) |
| 3 | `v2_3_control.py` | muon momentum calibration on the Z peak; fake-factor maps and the non-prompt template | [13](13-fake-factor.md), [14](14-fit-and-systematics.md) |
| 4 | `v2_4_histograms.py` | every histogram: samples x regions x variables x systematic variations | [14](14-fit-and-systematics.md) |
| 5 | `v2_5_fit.py` | TRExFitter inputs, config, fit, cross sections | [14](14-fit-and-systematics.md) |
| 6 | `v2_6_report.py` | data/MC plots, `RESULTS_v2.md`, summary | [15](15-combination-inputs.md) |

## Samples

Signal and Z backgrounds: `DYJetsToLL_M-50` aMC@NLO (recid 35669) split by the LHE lepton
flavour into DY -> mu mu / ee / tau tau. Backgrounds: `TTTo2L2Nu`, `ST_tW_top/antitop`
(NoFullyHadronic), `WWTo2L2Nu`, `WZTo3LNu`, `WZTo2Q2L`, `ZZTo4L`, `ZZTo2L2Nu`, `ZZTo2Q2L`.
Alternative signal generator: `ZToMuMu_M-50To120` powheg (recid 75482). `WJetsToLNu` only in
the prompt subtraction of the fake-factor regions. **No QCD simulation**: non-prompt muons are
data-driven. Cross sections and their provenance: `zmumu/samples.py`.

Only *prompt* simulated leptons (`genPartFlav` 1 or 15) enter the signal region, so the
data-driven non-prompt estimate and the simulation never double count.
