# Non-prompt muons: the fake-factor method

`zmumu/fakes.py`, `scripts/v2_3_control.py`; results `output/v2/fakes.json`, plots `ff_*.png`.

## Classes

- **tight**: `tightId`, IP cuts, `pfRelIso04_all < 0.15` (the signal muon);
- **anti-tight**: `tightId`, IP cuts, `0.20 < pfRelIso04_all < 1.0` -- same ID, so the fake
  factor is an isolation extrapolation; the 0.15-0.20 gap reduces prompt leakage. Alternative
  window 0.15-0.60 for the systematic.

## Measurement

`FF(pT, |eta|) = [N_tight - N_tight^{MC prompt-prompt}] / [N_anti - N_anti^{MC prompt-prompt}]`

- nominal region: **same-sign** tag + probe pairs (tag tight, pT > 26, trigger-matched; probe
  tight or anti-tight, pT > 20; m > 12 GeV). The probe is not the triggering muon, so its
  isolation is not biased by the online isolation of IsoMu24, and the fake composition
  (W+jets, heavy flavour) is that of the signal region;
- alternative region: single muon + jet (skim category C: exactly one muon, MET < 30, mT < 30,
  >= 1 jet), MC-subtracted, prescale weight 10;
- "prompt-prompt" = both muons `genPartFlav` in {1, 15} in every MC sample (DY, top, diboson,
  W+jets).

## Application

Opposite-sign events with exactly one tight (pT > 26, matched) and one anti-tight muon, no
other tight muon, FSR-recovered 60 < m < 120 GeV:
`N_fake(m) = sum FF(pT, |eta| of the anti-tight muon) x [data - prompt-prompt MC x SF_antiiso]`.
The prompt subtraction is the delicate part: Z -> mu mu with one muon in the anti-isolation
window is far more frequent than genuine fakes, so the MC prompt rate is calibrated with the
anti-isolation scale factor from tag-and-probe. The resulting estimate has a large relative
uncertainty but is a tiny fraction of the signal (v1's same-sign estimate: 2 375 events,
0.02%).

Closure: the fake factor applied to same-sign tight + anti-tight events must reproduce the
same-sign tight-tight yield after prompt subtraction (`ff_closure_and_template.png`).

## Uncertainties (nuisance parameters of the fit)

- `FakeStat_mumu`: statistical uncertainty of the FF map, coherent over cells;
- `FakeMethod_mumu`: largest of {alternative region, alternative anti-isolation window, prompt
  subtraction +-30%, same-sign non-closure}, floored at 30% and capped at 100%.

Numbers and the template: `RESULTS_v2.md`, `fakes.json`.
