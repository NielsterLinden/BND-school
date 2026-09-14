# Systematic uncertainties

Every contribution, where it comes from, and how much it can be trusted.
Numbers are in `output/RESULTS.md`; this page explains them.

| Source | Type | How it is obtained | Trust |
|---|---|---|---|
| luminosity | external | CMS 2016 legacy calibration, 1.2% | solid |
| muon reconstruction | external | Muon POG, 0.4%/muon -> 0.8% | solid but unverifiable here |
| trigger method bias | assigned | flat 0.5%, not measured | **weakest point** |
| background estimate | measured | propagated from control-region stats + transfer factors | solid, and tiny |
| eff(ID), eff(ISO) stat | measured | Clopper-Pearson on tag-and-probe | solid |
| eff(trigger) stat | measured | Clopper-Pearson | solid |
| FSR recovery | measured | bare vs FSR-recovered yield | solid |
| muon momentum scale | assigned | flat 0.2% | approximate |

## The three that dominate

**Luminosity, 1.2%.** The single largest uncertainty, and nothing in this
analysis can reduce it. It is the published CMS 2016 legacy value. Any
cross-section measurement on this dataset is floored here.

**Muon reconstruction efficiency, 0.8%.** An external input, because NanoAOD
does not store `generalTracks` and the denominator for a tracking-efficiency
tag-and-probe therefore does not exist in this format
([04-efficiency.md](04-efficiency.md)). Enters squared, one factor per muon.
Cannot be verified without going back to MiniAOD.

**Trigger method bias, 0.5%.** This is an *assigned* number covering a known
methodological flaw, not a measured one. The reference-trigger method uses
correlated reference paths, because the skim lacks `TrigObj` branches and the
SingleMuon dataset contains no unbiased reference. The 0.5% is a judgement, and
it is the number most likely to be wrong. See the fix in `handoff.md`.

## The ones that turned out not to matter

**Statistical, ~0.03%.** Eleven million signal events. This measurement is
nowhere near statistics-limited; more data would change nothing.

**Background, ~0.06%.** A 0.15% background with a 50% uncertainty on the larger
component gives a negligible contribution. The background work in this analysis
is worth doing because it *demonstrates* the sample is clean, not because the
subtraction moves the answer.

**FSR recovery, ~0.007%.** Surprisingly small, and worth understanding rather
than just accepting: FSR redistributes events *within* the 60-120 GeV window far
more than it moves them *across* the 60 GeV boundary. In a narrow window
(80-100 GeV) this would be one of the leading systematics. See
[02-fsr-photons.md](02-fsr-photons.md).

## What is not included

Honest list of what a publication-grade version would add:

- **Muon momentum scale and resolution (Rochester corrections).** Currently a
  flat 0.2%. The fitted peak sits ~0.4 GeV below PDG, which is partly this.
- **Pileup reweighting.** Irrelevant for a data-only measurement with
  data-driven efficiencies, but it would matter the moment MC enters.
- **Efficiency correlation between the two muons.** Assumed zero.
- **Bin-by-bin efficiency propagation.** The efficiency is applied as a single
  event-level average rather than per event. With a flat efficiency map this is
  exact; with the observed variation it is a sub-0.1% approximation.
- **L1 prefiring.** A real ~1% effect in 2016 for forward objects. Not applied
  -- `L1PreFiringWeight` was not kept in the skim. For muons at `|eta| < 2.4`
  the effect is small, but this is a genuine omission.

## How to re-derive any of them

All systematic inputs are constants in `config.py`:

```python
LUMI_REL_UNC         = 0.012    # luminosity
MU_RECO_EFF_UNC      = 0.0040   # per muon
TRIG_METHOD_REL_UNC  = 0.005    # reference-trigger method bias
MUON_SCALE_REL_UNC   = 0.002    # momentum scale
R_OS_SS_UNC          = 0.50     # same-sign transfer factor
FS_REL_UNC           = 0.50     # e-mu transfer factor
```

Change one, rerun `.venv/bin/python run_all.py --from 3` (steps 1 and 2 do not
need repeating unless the *selection* changed), and the breakdown in
`output/plots/step5_systematics.png` updates.
