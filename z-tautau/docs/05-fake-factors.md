# 05 — Jet → τh fakes: the fake-factor method

Code: `ztautau/fakes.py`, `scripts/step3_fakefactors.py`; plots `output/plots/step3_*.png`.

## Why fakes dominate and why data-driven

QCD multijet production is 10⁶ times larger than Z→ττ, and even at DeepTau Medium ~0.5–1 % of jets pass as
τh. In the signal region **81 % of the events have at least one jet faking τh**. Simulating this needs
jet→τh fake rates at the per-mille level in the tails of the QCD cross section. The QCD MC samples are
statistically far too small for that, and the simulated fake rates are known to be off. So the fakes come
from data, and the QCD simulation is not used.

## The classic fake factor, applied to the leading τ

```
FF(era, DM, N_jets, pT)  =  N(SS, τ1 T, τ2 T)  /  N(SS, τ1 L, τ2 T)                (measured)
N_fakes(SR)              =  Σ_{data in AR: OS, τ1 L, τ2 T}  C_OS/SS(era, N_jets) · FF(era, DM₁, N_jets, pT₁)
```

T = DeepTau VSjet Medium, L = VVVLoose and not Medium.

* **Same-sign (SS) pairs** are the determination region: Z→ττ is always OS, so SS pairs with an isolated
  τ2 are ~99 % QCD.
* **τ2 T in both measurement and application.** In dijet events the isolation of the two fakes is
  correlated, so the FF has to be conditioned on the same τ2 isolation as the SR.
* **Leading τ only.** Every SR event whose leading τ is a jet is covered. SR events with a genuine τ1
  and a jet as τ2 (W+jets, tt̄) are taken from simulation. Simulation keeps only events whose leading τ is
  not a jet, so nothing is counted twice. This is the CMS H→ττ convention for τhτh.

### Binning — driven by closure

| iteration | binning | same-sign closure |
|---|---|---|
| 1 | DM × pT | N_jets: 1.09 / 0.89 / 0.83 / 0.75 for 0/1/2/≥3 jets; m_tt < 110 GeV over-predicted by ~15 % |
| 2 | DM × **N_jets (0, 1, ≥2)** × pT | N_jets flat; m_tt within ±10 % (statistics) |
| 3 (used) | **era** × DM × N_jets × pT | per era 1.000 / 1.000 (was 0.956 G / 1.051 H) |

Quark and gluon jets fake τh at different rates, and the jet multiplicity changes that mixture. The two
eras use different HLT τ isolation. pT bins: 40, 45, 50, 60, 80, ∞ GeV → 2 × 4 × 3 × 5 = 120 bins;
~27 000 SS_T events, median statistical uncertainty per bin 10 %.

Fake factors (era G | H, nominal), N_jets = 0 / 1 / ≥ 2, pT bins 40–45, 45–50, 50–60, 60–80, > 80 GeV:

| | 0 jets | 1 jet | ≥ 2 jets |
|---|---|---|---|
| G DM0 | 0.44 0.39 0.34 0.34 0.33 | 0.42 0.30 0.25 0.26 0.27 | 0.42 0.18 0.21 0.23 0.28 |
| G DM1 | 0.27 0.25 0.24 0.20 0.21 | 0.22 0.23 0.19 0.17 0.17 | 0.23 0.21 0.16 0.15 0.17 |
| G DM10 | 0.18 0.19 0.20 0.23 0.25 | 0.16 0.15 0.15 0.16 0.18 | 0.16 0.15 0.15 0.15 0.13 |
| G DM11 | 0.06 0.07 0.07 0.09 0.08 | 0.05 0.06 0.07 0.07 0.08 | 0.04 0.07 0.06 0.06 0.07 |
| H DM0 | 0.47 0.46 0.41 0.41 0.37 | 0.33 0.29 0.34 0.33 0.37 | 0.53 0.34 0.29 0.22 0.32 |
| H DM1 | 0.29 0.26 0.26 0.24 0.20 | 0.28 0.22 0.20 0.19 0.18 | 0.20 0.16 0.21 0.17 0.15 |
| H DM10 | 0.18 0.21 0.23 0.24 0.25 | 0.12 0.16 0.19 0.20 0.17 | 0.19 0.14 0.15 0.17 0.18 |
| H DM11 | 0.07 0.06 0.08 0.08 0.10 | 0.07 0.06 0.07 0.06 0.08 | 0.11 0.04 0.07 0.08 0.06 |

The DM11 (3-prong + π⁰) FF is 3–5× smaller than for 1-prong: many jets are reconstructed as DM11, and
DeepTau rejects them efficiently.

## OS/SS extrapolation: C_OS/SS

The FF is measured in SS and applied in OS. For quark-initiated jets the charges of the two leading tracks
are correlated, so OS/SS ≠ 1. C is measured in the **τ2 anti-isolated sideband** (τ2 L), which is fake
dominated for both charges: FF measured in SS-AI and applied to OS-AI_L predicts OS-AI_T, and
C = observed / predicted. The table gives C per era and jet multiplicity (statistical uncertainty 0.004–0.011):

| | 0 jets | 1 jet | ≥ 2 jets |
|---|---:|---:|---:|
| Run2016G | 1.068 | 1.095 | 1.107 |
| Run2016H | 1.065 | 1.078 | 1.135 |

As the τ2 sideband is tightened towards Medium, C rises from 1.05 (τ2 VVVLoose-not-VVLoose) to 1.08–1.09
(after subtracting genuine τ). The extrapolation to τ2 = Medium therefore carries a **3 % systematic**,
added in quadrature to the statistical one (`FakeOSSS_tautau`).

## Genuine-τ contamination and the "MC subtraction" question

The determination and application regions contain some genuine τ; simulated yield with a genuine leading τ:

| region | data | genuine τ1 (MC) | fraction |
|---|---:|---:|---:|
| SS_T (FF numerator) | 26,761 | 232 | 0.9 % |
| SS_L (FF denominator) | 140,365 | 108 | 0.1 % |
| AR (application) | 187,973 | 2,386 | 1.3 % |
| OSAI_T (C numerator) | 219,737 | 8,321 | 3.8 % |

As requested for this iteration, the **nominal does not subtract** simulation (`config.FF_SUBTRACT_MC =
False`). Contamination in the AR is mostly Z→ττ with τ1 failing Medium, and it double-counts signal. In the
C numerator it biases C up by ~2.5 %. The **alternative "mcsub"** (subtract simulated genuine τ1 events
everywhere, with the same weights as the SR prediction) is computed by the same code. It gives 36,988
instead of 38,746 prefit fake events (−4.5 %) and C = 1.051 instead of 1.078 (inclusive). Both variants are
fitted; the fit pulls the nominal's `FakeOSSS_tautau` towards the subtracted value, and the final μ differs by
only ~0.05 (`08-fit-and-results.md`). **Recommendation for the next iteration: make the subtracted variant
nominal.**

## Closure and its systematic

* **Same-sign closure** (FF applied to SS_L vs SS_T): exact in pT/DM/N_jets/era by construction. In m_tt,
  bin ratios lie within ±11 % with 4–13 % statistical precision at low mass, and within ±2.5 % above 110 GeV.
  The bin-by-bin ratio in m_tt becomes a shape nuisance parameter (`FakeClosure_tautau`: up = ratio,
  down = 1/ratio). Some residual structure in η(τ1) (±15 %) and pT(τ2) (−7 % at 60–100 GeV) remains.
  It is absorbed by this nuisance and is a candidate for a future closure correction.
* **OS anti-isolated closure** (`step3_OSAI_m_tt_mcsub.png`): after C and genuine-τ subtraction, data/pred
  is within ±5 % over the full m_tt range, so the SS→OS shape extrapolation holds.
* **Statistics of the FF**: one nuisance parameter per decay mode shifts all bins of that DM by their
  statistical uncertainty (`FakeStat_tautau_DM*`). The data statistics of the AR enter the per-bin MC-stat
  γ parameters.
