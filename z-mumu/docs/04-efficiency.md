# Step 2 -- efficiencies from tag-and-probe

`scripts/step2_efficiency.py`

## Why tag-and-probe

The cross section needs the probability that a Z decay in the fiducial volume
ends up in the selected sample. With no simulation available, that probability
has to come from the data. Tag-and-probe is the standard way: use the Z peak
itself as a source of known-real muons.

Pick one muon that is unambiguously good (the **tag**) and require it to make a
Z-mass pair with a second object (the **probe**) selected *without* the
requirement being measured. The fraction of probes that pass that requirement is
its efficiency. Because the pair mass is constrained to the Z peak, the probe is
a real muon even though nothing was required of it.

## Definitions used here

**Tag** -- a fully selected muon with `pT > 26` GeV (`config.TAG_PT_MIN`), in an
event that fired the analysis trigger.

**Probe** -- deliberately as loose as NanoAOD allows:

```python
(Muon_isTracker | Muon_isGlobal) and pT > 20 GeV and |eta| < 2.4
```

No ID, no isolation, no impact parameter. This is the closest thing to an
unbiased track that the collection offers.

**Pair** -- opposite sign, `70 < m < 110` GeV (tighter than the analysis window,
to keep the pairs pure).

Both orientations of every pair are used, so an event where both muons qualify
as tags contributes two independent measurements. This is handled with
`ak.combinations` followed by a concatenate of the two orderings.

## What is measured

| Efficiency | Denominator | Numerator |
|---|---|---|
| `eff_ID` | every valid probe | probe passes `mediumId` + impact parameters |
| `eff_ISO` | probes that passed the ID | probe also passes `pfRelIso04 < 0.15` |

They are measured **in sequence**, not independently: the isolation denominator
is the ID numerator. That is what makes the product `eff_ID * eff_ISO` the
correct combined efficiency, and it is why the order in `config.py` matters.

Both are stored as 2D maps in (pT, eta) as well as 1D projections. Step 5 uses
the 2D maps, weighted by where the signal muons actually are.

Measured (full dataset): `eff_ID ~ 0.977`, `eff_ISO ~ 0.949`.

## The trigger efficiency, and its caveat

**This is the weakest link in the measurement. Read this before quoting the
number.**

The proper way to measure a trigger efficiency is per-muon, by matching each
offline muon to an HLT trigger object and asking whether *that muon* fired the
path. That requires the `TrigObj_*` branches -- **which were not kept in the
skim.** Nothing in this repository can recover them; it needs a re-skim from the
parent NanoAOD.

The fallback used here is the **reference-trigger method**:

```
    eff_trigger  =  N(reference fired AND IsoMu24/IsoTkMu24 fired)
                    -------------------------------------------
                            N(reference fired)
```

with the reference being `HLT_Mu50`, `HLT_IsoMu27`, `HLT_Mu27`,
`HLT_Mu45_eta2p1` -- single-muon paths that are not the ones being measured, and
that exist in every file.

**Why this is biased.** A proper reference trigger is *independent* of the path
being measured. These are not: they fire on the same muons, in the same event,
through partly the same L1 seeds. The measured number is therefore an
efficiency *relative to* a correlated reference, and it is biased **high** --
events where no muon was triggerable are missing from both numerator and
denominator.

There is no unbiased reference inside the SingleMuon primary dataset, because
the dataset is *defined* by firing a single-muon path. A genuinely independent
measurement needs an orthogonal dataset (`DoubleMuon` is available in the Open
Data release, or a MET/jet-triggered sample).

**What is done about it.** A flat 0.5% systematic
(`config.TRIG_METHOD_REL_UNC`) covers the method bias. The measured value is
~0.998 at the event level, which is plausible -- either of two muons above 26
GeV can fire the trigger, so the event-level efficiency is close to
`1 - (1 - eff_mu)^2` -- but it should be treated as an upper bound until
checked properly. This is the first thing to fix if the measurement is to be
taken further; see `handoff.md`.

## What is *not* measured: reconstruction efficiency

`eff_RECO` is the probability that a muon track becomes a reconstructed `Muon`
object. Measuring it needs a denominator built from tracker tracks
(`generalTracks`), and **NanoAOD does not store them.** There is no way around
this without going back to MiniAOD.

It is therefore taken as an **external input**: `0.996 +/- 0.004` per muon
(`config.MU_RECO_EFF`), the CMS Muon POG value for 2016 UL. Entering squared, it
contributes 0.8% to the total uncertainty -- the second-largest systematic after
the luminosity. It is flagged as external in `RESULTS.md` so nobody mistakes it
for something measured here.

## Outputs

- `output/plots/step2_eff_id_vs_pt.png`, `..._vs_eta.png`
- `output/plots/step2_eff_iso_vs_pt.png`, `..._vs_eta.png`
- `output/plots/step2_eff_trigger_vs_pt.png` -- the turn-on curve
- `output/plots/step2_eff_id_map.png`, `step2_eff_iso_map.png` -- 2D maps
- `output/plots/step2_tagprobe_mass.png` -- passing and failing probe pair mass,
  a purity check on the method itself

Uncertainties are Clopper-Pearson intervals, not `sqrt(eps(1-eps)/N)`, which
would give a zero error in a bin where every probe passes.
