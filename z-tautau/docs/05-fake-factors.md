# 05 — Jet → τh fakes: the fake-factor method

Code: `ztautau/fakes.py`, `scripts/step3_fakefactors.py`; plots `output/plots/step3_*.png`; numbers
`output/data/fakefactors.json`. The per-category validation in the BDT score is in `09-bdt.md`.

## Why fakes dominate and why data-driven

QCD multijet production is 10⁶ times larger than Z→ττ, and even at DeepTau Tight a few per mille of jets pass as
τh (Medium: 0.5–1 %). In the signal region **64 % of the events have at least one jet faking τh** (80 %
with the Medium working point of v1–v2.1). Simulating this needs
jet→τh fake rates at the per-mille level in the tails of the QCD cross section. The QCD MC samples are
statistically far too small for that, and the simulated fake rates are known to be off. So the fakes come
from data, and the QCD simulation is not used.

## The classic fake factor, applied to the leading τ

```
FF(era, DM, N_jets, pT) = N(SS, τ1 T, τ2 T) − MC  /  N(SS, τ1 L, τ2 T) − MC             (measured)
FF_corr                  = FF · f(|η(τ1)|) · g(pT(τ2))                                 (closure corrections)
N_fakes(SR, cat k)       = Σ_{data in AR, cat k} C_OS/SS(era, N_jets) · FF_corr  −  Σ_{MC genuine τ1 in AR, cat k} w · C · FF_corr
```

T = DeepTau VSjet Tight (v3; Medium in v1–v2.1), L = VVVLoose and not T, "MC" = simulated events whose leading τ is not a
jet (genuine τh, e or μ), weighted like the signal-region prediction.

* **Same-sign (SS) pairs** are the determination region: Z→ττ is always OS, so SS pairs with an isolated
  τ2 are ~99 % QCD.
* **τ2 T in both measurement and application.** In dijet events the isolation of the two fakes is
  correlated, so the FF has to be conditioned on the same τ2 isolation as the SR.
* **Leading τ only.** Every SR event whose leading τ is a jet is covered. SR events with a genuine τ1
  and a jet as τ2 (W+jets, tt̄) are taken from simulation. Simulation keeps only events whose leading τ is
  not a jet, so nothing is counted twice. This is the CMS H→ττ convention for τhτh.
* **MC subtraction (nominal since v2).** The determination and application regions contain genuine τ's
  (table below), mostly Z→ττ with τ1 failing the tight working point in the AR: ~1–2 % of the AR, but ~6 % of the signal,
  which the unsubtracted estimate double counts, and 3.9 % of the C numerator, which biases C up by 2.5 %.
  v1 quoted the unsubtracted variant as nominal; the review (`review/REVIEW_v3.md` 3.2) showed the bias (Δμ = −0.05
  and a 7 % over-prediction of the fakes in the purest fake phase space). Every simulated sample is
  subtracted; W+jets (1114 events with weights 20–200, 16 % of them negative) with **uniform weights**
  (each event carries the sample's mean weight, `analysis.subtraction_weights`): it is a real genuine-τ1
  component of the τ2-anti-isolated sideband (2 % of OSAI_T, worth 1.5 % on C), but with its raw weights a
  single event of weight −63 in the top BDT bin of SS_T faked a 75 % non-closure. The unsubtracted variant
  is only produced on request (`step3_fakefactors.py --with-nosub`, `--ff-variant nosub`).

| region | data | genuine τ1 (MC) | fraction |
|---|---:|---:|---:|
| SS_T (FF numerator) | 9 554 | 164 | 1.7 % |
| SS_L (FF denominator) | 86 929 | 110 | 0.1 % |
| AR (application) | 120 805 | 2 996 | 2.5 % |
| OSAI_T (C numerator) | 144 978 | 8 926 | 6.2 % |

### Binning — driven by closure

| iteration | binning | same-sign closure (Medium working point, v2) |
|---|---|---|
| 1 | DM × pT | N_jets: 1.09 / 0.89 / 0.83 / 0.75 for 0/1/2/≥3 jets; m_tt < 110 GeV over-predicted by ~15 % |
| 2 | DM × **N_jets (0, 1, ≥2)** × pT | N_jets flat; m_tt within ±10 % (statistics) |
| 3 | **era** × DM × N_jets × pT | per era 1.000 / 1.000 (was 0.956 G / 1.051 H) |
| 4 (v2, used) | 3 × **closure corrections f(\|η(τ1)\|) · g(pT(τ2))** | η(τ1) within ±5 % (was ±15 %), pT(τ2) flat (was −7 % at 60–100 GeV) |

Quark and gluon jets fake τh at different rates, and the jet multiplicity changes that mixture. The two
eras use different HLT τ isolation. pT bins: 40, 45, 50, 60, 80, ∞ GeV → 2 × 4 × 3 × 5 = 120 bins;
~9 554 SS_T events, median statistical uncertainty per bin ~10–15 % (all 120 values:
`step3_fakefactors.png`; DM11 (3-prong + π⁰) fake factors are 3–5× smaller than for 1-prong).

### Closure corrections (new in v2)

After the table, the FF still depends on two variables that are not binned (`step3_closure_corrections.png`,
`step3_closure_SS_t1_eta_mcsub_nocorr.png` vs `..._mcsub.png`):

| variable | same-sign obs/pred before correction | why |
|---|---|---|
| \|η(τ1)\| in 6 bins 0–0.4–0.8–1.2–1.5–1.8–2.1 | 0.92, 0.99, 1.08, **1.16**, 0.98, **0.93** | the FF itself varies by ±15 % with \|η\|, with the *same* shape in both eras, all pT and N_jets bins and inside every decay mode (`review/REVIEW_v3.md` section 5): the DeepTau score distribution just below Medium is η-dependent (barrel–endcap transition) and the Medium efficiency for jets drops at the tracker/HLT edge. Not a trigger effect (the eras agree). |
| pT(τ2) in 5 bins 40–45–50–60–80–∞ | 1.02, 0.99, 1.00, **0.94**, **0.95** | the isolation of the two jets in a dijet is correlated beyond what "τ2 Medium" conditions |

Each correction is the same-sign obs/pred ratio (MC subtracted), applied multiplicatively in
`fakes.evaluate`; g is measured after f. The η correction leaves the m_tt template unchanged to < 0.4 %
per bin (m_tt is uncorrelated with η(τ1) in the AR) but it is a prerequisite for the BDT categories, whose
inputs include |η|. The statistical precision of the corrections (2–4 % per bin) is part of the residual
non-closure below.

## OS/SS extrapolation: C_OS/SS

The FF is measured in SS and applied in OS. For quark-initiated jets the charges of the two leading tracks
are correlated, so OS/SS ≠ 1. C is measured in the **τ2 anti-isolated sideband** (τ2 L), which is fake
dominated for both charges: FF measured in SS-AI and applied to OS-AI_L predicts OS-AI_T, and
C = observed / predicted, after subtracting the simulated genuine-τ1 events from both. Per era and jet
multiplicity, inclusive in the BDT score (statistical uncertainty 0.004–0.010):

| | 0 jets | 1 jet | ≥ 2 jets |
|---|---:|---:|---:|
| Run2016G | 1.057 | 1.086 | 1.077 |
| Run2016H | 1.042 | 1.037 | 1.070 |

(inclusive 1.052; without the subtraction 1.078: the difference is the genuine-τ contamination of the
numerator, a third of it W+jets). As the τ2 sideband is tightened towards the tight working point, C rises by a few percent, so the extrapolation
to τ2 = T carries a **3 % systematic**, added in quadrature to the statistical one (`FakeOSSS_tautau`).
`step3_OSAI_m_tt_mcsub.png` shows that after C the SS→OS shape extrapolation holds within ±5 % over the
full m_tt range.

**C depends on the BDT category (since v2.1).** With the inclusive C the same sideband is described to 0.4 % in
the fake-dominated category but under-predicted by 6.5 % in the middle and 3 % in the signal-like
category: the OS/SS charge correlation of the two jets depends on the topology the classifier selects
(more quark-initiated, balanced dijets at high score). In the first v2 fit this showed up as positive
pulls of every fake parameter of SR1 and SR2 (0.6–1.5 σ). Since v2.1, step 4 measures C per
(era, N_jets, BDT category) in the τ2-anti-isolated sideband (`fakes.osss_correction(..., n_cat)`), with
one OS/SS nuisance parameter per category (`FakeOSSS_tautau_c<k>`, statistical ⊕ 3 %). A bin whose statistical
error exceeds 10 % (the 0-jet bins of the signal-like category: 0.58 ± 0.16 in G) takes the N_jets-inclusive
value of its era and category:

| | 0 jets | 1 jet | ≥ 2 jets |
|---|---:|---:|---:|
| category 0, Run2016G | 1.056 | 1.079 | 1.059 |
| category 0, Run2016H | 1.039 | 1.029 | 1.058 |
| category 1, Run2016G | 1.109 | 1.172 | 1.097 |
| category 1, Run2016H | 1.163 | 1.149 | 1.113 |
| category 2, Run2016G | 1.187 | 1.045 | 1.305 |
| category 2, Run2016H | 1.072 | 0.960 | 1.151 |

## The uncertainty model (v2)

| source | v1 | v2 | why |
|---|---|---|---|
| FF statistics | four NPs, one per DM, shifting all 30 bins of a DM coherently | per-event relative error of the FF bin, added to the `Fakes` template variance (γ parameters) | the 120 bin errors are independent; a coherent shift over-states the normalisation part (which the fit measures anyway in the fake-dominated category) and under-states the shape part |
| OS/SS | 3 % ⊕ stat, one NP | C per (era, N_jets, BDT category), one NP per category (`FakeOSSS_tautau_c<k>`) | the charge correlation depends on the topology the BDT selects (3–6 % between categories) |
| non-closure | one NP over all 14 m_tt bins, up = SS ratio, down = 1/ratio | one normalisation-type NP per **BDT category × mass region** (m_tt below / above 110 GeV): `FakeClosure_tautau_c<k>_lo|hi`, size = residual same-sign obs/pred deviation ⊕ its statistical error (the SR0 low-mass bins are not fitted, so `c0_lo` is not used) | v1's single NP tied the precise high-mass bins to the imprecise peak bins and was constrained to 0.32 σ by the former, which silently shrank the peak-bin uncertainty; per category and mass region the constraint comes from the right place |
| η(τ1), pT(τ2) dependence | absorbed by the non-closure NP | corrected | see above |
| W+jets composition of the AR (quark-jet FF ≠ QCD FF) | not covered | not covered, ≲ 1 % of the fakes | noted in `review/REVIEW_v3.md` 3.5 |

The residual same-sign non-closure per category and mass region, and hence the NP sizes, are printed by
step 4 and stored in `output/data/yields.json["meta"]["closure_nps"]` and `output/results.json`.
