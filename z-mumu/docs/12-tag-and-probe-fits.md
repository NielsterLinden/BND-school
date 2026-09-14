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

Numbers: `RESULTS_v2.md` (mean SFs, plateau efficiencies, per-cell tables in the JSON).
