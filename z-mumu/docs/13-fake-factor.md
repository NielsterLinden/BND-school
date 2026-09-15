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
other tight muon, FSR-recovered 60 < m < 120 GeV. The FF-weighted mass spectrum of this
region, `D(m) = sum_cells FF x data(m)`, is dominated by Z -> mu mu events with one muon in
the anti-isolation window (~ 90 000 FF-weighted events against a few thousand genuine fakes),
so a plain MC subtraction is unstable (it gave a *negative* yield with the anti-isolation
scale factor from tag-and-probe, and swung by +-25 000 for a +-30% subtraction variation).
The non-prompt yield is therefore obtained from a **two-template fit**
`D(m) = a P(m) + b F(m)` in 2 GeV bins: `P` = FF-weighted prompt-prompt simulation (the Z
peak), `F` = FF-weighted same-sign data minus prompt simulation (the non-prompt shape, which
the closure test validates). The fitted `a` (0.88) is the prompt normalisation -- the
same-sign region and the anti-isolation scale factor (0.945) are its cross-checks -- and
`b F(m)` is the template of the signal region: **3 825 +- 80 events, 0.035% of the data**
(v1's same-sign estimate: 2 375 x R_OS/SS). `b` = 1.84 is the opposite-/same-sign ratio of the
non-prompt background in this region (charge-asymmetric W+jets fakes).

Closure: the fake factor applied to same-sign tight + anti-tight events reproduces the
same-sign tight-tight yield after prompt subtraction: predicted 2 077 +- 20, observed
1 922 +- 85 (`ff_closure_and_template.png`).

The single-muon + jet region gives a fake factor ~ 17 x larger: the IsoMu24 online isolation
depletes its anti-isolated denominator (the probe *is* the triggering muon). It is kept in
`fakes.json` as a documented, biased cross-check and not used as a variant; the non-isolated
prescaled paths (`HLT_Mu8/17`, `HLT_Mu3_PFJet40`, kept in the skim) would be the way to fix it.

## Uncertainties (nuisance parameters of the fit)

- `FakeStat_mumu`: statistical uncertainty of the FF map, coherent over cells;
- `FakeMethod_mumu`: largest of {alternative anti-isolation window (+4%), prompt subtraction
  +-30% in the FF measurement (-19% / +19%), exponential instead of the same-sign non-prompt
  shape (-4%), same-sign non-closure (7%)}, floored at 30% and capped at 100% -> 30%.

The `Fakes` fit template is written without Sumw2: the two-template fit's (sigma_b F)^2 is a
fully correlated normalisation error, already carried by the two nuisance parameters, and would
otherwise be double counted as per-bin gammas (REVIEW.md, section 7).

Numbers and the template: `RESULTS_v2.md`, `fakes.json`.
