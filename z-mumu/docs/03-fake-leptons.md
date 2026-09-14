# Fake and non-prompt muons

## What "fake" means

A muon in the `Muon` collection did not necessarily come from a W or Z decay.
Three things can put an entry there that the analysis does not want:

| Source | Mechanism | Why it survives cuts |
|---|---|---|
| **Heavy flavour** | b or c hadron decays semileptonically inside a jet | a genuine muon, just not prompt |
| **Decay in flight** | a pion or kaon decays to a muon before the calorimeter | genuine muon track, displaced |
| **Punch-through** | a hadron shower leaks through the calorimeter into the muon chambers | not a muon at all |

The first two are *non-prompt* muons -- real muons from the wrong place. The
third is a true *fake*. The distinction matters for how you suppress them, but
not for how they are estimated here, so this analysis lumps them together and
calls them all fakes.

## Why they are a problem, and why they are a small one here

Fakes live inside jets, so they are suppressed by two of the selection cuts
doing most of the work:

- **Isolation** (`pfRelIso04_all < 0.15`). A muon inside a jet has hadronic
  activity around it. This is the single most powerful cut against fakes, and
  it costs ~5% on real muons (measured: `eff(ISO) = 0.949`).
- **Impact parameters** (`|dxy| < 0.2`, `|dz| < 0.5` cm). Muons from b decays
  come from a displaced vertex.

Together with medium ID, these make the sample very pure. The measured fake
contamination in the Z window is **~0.03%**. This is not a hard measurement to
protect against fakes -- but it is exactly the case where one should *show* that
rather than assume it.

## The estimate: same-sign control region

The method rests on one fact: **fakes carry no charge correlation.** A muon from
a b decay or a punch-through hadron is equally likely to be positive or
negative, independent of the other muon's charge. Real `Z -> mu mu` decays are
almost perfectly opposite-sign (charge misassignment for muons at these momenta
is ~10^-5, negligible).

So, applying the *identical* selection but requiring same-sign muons gives a
sample that is almost pure fakes, and:

```
    N_fake(opposite sign)  =  R_OS/SS  *  N(same sign)
```

`R_OS/SS` is 1.0 if the charge really is uncorrelated. This analysis uses
**1.0 +/- 50%** (`config.R_OS_SS`, `R_OS_SS_UNC`). The 50% is deliberately
generous and is the honest statement that the ratio was not measured in situ --
with a background this small it changes nothing (it contributes 0.06% to the
total uncertainty).

Implemented in `scripts/step1_selection.py` (the `_fill_same_sign` branch, which
shares every cut with the signal region up to the charge requirement) and
combined in `scripts/step3_backgrounds.py`.

## What the same-sign method does *not* catch

This is the part that is easy to get wrong. The same-sign region measures
**charge-uncorrelated** backgrounds. It is blind to any background that produces
genuine **opposite-sign prompt** muon pairs:

- `ttbar -> W+ b W- bbar -> mu+ mu- + X`
- `WW -> mu+ mu-`
- single top (`tW`)
- `Z -> tau tau -> mu mu + neutrinos`

These are all opposite-sign by construction and all have two prompt, isolated
muons. They pass every cut. A same-sign estimate alone would miss them
entirely, and they are the *larger* of the two backgrounds here (0.12% versus
0.03%).

They are estimated separately with the e-mu method -- see
[05-backgrounds.md](05-backgrounds.md).

## Cross-checks worth doing

1. **Same-sign mass shape.** If the same-sign region were pure fakes it should
   be roughly flat, with no Z peak. A visible peak means charge misassignment
   is leaking real Z events in. See `output/plots/step3_control_regions.png`.
2. **Invert isolation.** Selecting non-isolated muons gives a fake-enriched
   sample; check that its same-sign / opposite-sign ratio really is ~1. This is
   how you would *measure* `R_OS/SS` rather than assume it.
3. **ABCD.** Isolation and impact parameter are nearly uncorrelated for fakes,
   so the standard ABCD closure test applies. Not done here (the background is
   too small to be worth it) but it is the natural next step if you tighten the
   selection or move to a lower-pT region where fakes matter more.

## If fakes ever become important

They will if you loosen the isolation cut, lower the pT thresholds, or move to a
mass region away from the Z peak. In that case:

- Measure `R_OS/SS` in a fake-enriched region instead of assuming 1.0.
- Switch to a fake-rate method: measure the probability that a loose muon passes
  the tight selection, in a QCD-dominated control region, binned in pT and eta,
  then apply it to a loose-muon sample.
- Watch for real-muon contamination of the fake-rate measurement region -- it is
  the usual way this method goes wrong.
