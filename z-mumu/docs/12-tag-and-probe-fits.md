# Tag-and-probe with simultaneous pass/fail fits

`zmumu/tnp.py`, `scripts/v2_2_tnp.py`; results `output/v2/tnp/tnp_result.json`,
plots `output/v2/plots/tnp_*.png`.

## Definitions (sequential, Muon POG style, identical in data and simulation)

| efficiency | probe | pass |
|---|---|---|
| ID | loose muon (global or tracker, pT > 20, \|eta\| < 2.4) | `tightId` and \|dxy\| < 0.2, \|dz\| < 0.5 cm |
| isolation | tight ID | `pfRelIso04_all < 0.15` |
| anti-isolation | tight ID | `0.20 < pfRelIso04_all < 1.0` (calibrates the fake-factor prompt subtraction) |
| trigger | tight ID + isolation | dR < 0.1 to an IsoMu24/IsoTkMu24 trigger object (fired events only) |

Tag: tight ID + isolation, pT > 26, **trigger-matched** (v1 did not match the tag; that biased
the probes towards passing by 0.1-0.2%). Both orientations of every opposite-sign pair with
60 < m < 120 GeV; same-sign pairs kept for the counting cross-checks. Bins: pT
[20,25,30,35,40,45,50,60,80,120,200] x |eta| [0,0.9,1.2,2.1,2.4]; trigger in the finer
`TRIG_PT_BINS` across the 24 GeV turn-on.

## The fit

Per (pT, |eta|) cell, a simultaneous extended binned Poisson likelihood of the pass and
fail mass spectra (0.5 GeV bins, 60-120 GeV):

```
pass(m) = N eps      S_p(m) + B_p exp(-lambda_p (m-60))
fail(m) = N (1-eps)  S_f(m) + B_f exp(-lambda_f (m-60))
```

with S_p, S_f the gen-matched DY templates of the same cell (pass and fail separately: failing
probes have a different resolution) convolved with a Gaussian (free shift and width for each).
Ten parameters, `iminuit` (MIGRAD + HESSE); the failing-probe shape is frozen when fewer than
50 failing probes exist. Simulation is fitted the same way (weighted histograms with the
Bohm-Zech scaling), so the fit bias cancels in the scale factor; the gen-matched counting
gives the truth and `|eps_fit - eps_truth|` is a closure systematic.

Systematic variants, applied coherently to data and MC: background `erfc x exp` (CMSShape),
fit range 70-110 GeV, tag pT > 30 with iso < 0.10, templates from the powheg sample. The
scale-factor systematic is the largest deviation, added in quadrature to the closure term.

Trigger: same-sign-subtracted counting in the 70-110 GeV window per cell; the L1-seed veto
variant (objects with 20 <= L1 pT < 22 GeV removed) is added to the data uncertainty.

Isolation vs pileup: the SF is measured in four `PV_npvsGood` bins; the half-spread is
reported (`tnp_iso_vs_npv.png`).

## Application

`tnp.ScaleFactors`: per-muon ID and iso SF, per-event trigger weight with the measured
per-muon efficiencies (a sub-threshold second muon contributes its turn-on efficiency, not
zero). Nuisance parameters `MuonID`, `MuonIso`, `MuonTrigger`: all cells shifted coherently by
+-1 sigma (stat (+) syst) -- conservative for the statistical part, which is far below 0.1%.

## Results

20.0 M data pairs, 14.7 M (weighted) simulated pairs; 320 fits per efficiency.

| | data | simulation | scale factor (mean over cells) | stat | syst |
|---|---:|---:|---:|---:|---:|
| tight ID | 0.957 (20-25 GeV, barrel) ... 0.98 | 0.968 ... 0.99 | **0.980** (0.93-0.99) | 0.0014 | 0.0029 |
| isolation | 0.86 (20-25 GeV) ... 0.99 | | **1.005** (1.00-1.04) | 0.0009 | 0.0033 |
| anti-isolation (0.20-1.0) | | | 0.90 (large spread; only used as a cross-check of the fake-factor prompt subtraction) | | |
| trigger, per muon on the plateau | **0.907** (|eta| bins 0.93 / 0.93 / 0.88 / 0.79) | **0.923** (0.95 / 0.97 / 0.88 / 0.83) | | | L1 veto in the data error |

The data ID efficiency is 1-2% below the simulation everywhere, most in the forward endcap;
the fit-vs-truth closure in simulation is -0.1% (ID) and -0.2% (iso). The isolation scale
factor is flat in pileup (half-spread 0.0024 over four `PV_npvsGood` bins). Plots:
`tnp_eff_{id,iso,antiiso}_vs_pt.png`, `tnp_sf_*_map.png`, `tnp_fits_*_lowpt.png`,
`tnp_trig_eff_vs_pt.png`, `tnp_iso_vs_npv.png`; per-cell numbers in `output/v2/tnp/tnp_result.json`.
