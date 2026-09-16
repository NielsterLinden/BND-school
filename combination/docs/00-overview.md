# 00 — What the combination is

All three BND-school channels have now measured a Z cross section on the same CMS Open Data 2016
sample (Run2016G+H, 16 393.381 pb⁻¹): `z-mumu` (Z → μ⁺μ⁻), `z-tautau` (Z → τhτh) and `z-ee`
(Z → e⁺e⁻). This folder turns the three into one number.

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1980 ± 32 pb**, per lepton flavour, assuming lepton
> universality — **with a χ²/ndf of 9.67/2 (p = 0.008)**. See [`../result.md`](../result.md).

The compatibility number is part of the result, not a footnote. The ee and μμ channels are of
comparable precision (1.9 % and 1.8 %) and disagree by 123 pb, and [01](01-inputs.md) shows why:
the z-ee fit lets a ±5.9 % theory *normalisation* of its own signal template float against its
signal strength, which `fitting/CONVENTIONS.md` §3 forbids exactly because the two are degenerate.
Undoing only that (`ee_normfix`) moves the combination by −99 pb and is the largest single number
in the whole cross-check table.

## The chain

```
z-mumu   TRExFitter v1.8.0 fit of m(μμ)   ──►  μ_Z, grouped impacts, A, σ(60–120)
z-tautau TRExFitter v1.8.0 fit of m(ττ)   ──►  ...          │
z-ee     TRExFitter v1.8.0 fit of m(ee)   ──►  ...          │
   │                                                        │
   └─ tools/extract_zee.py (z-ee publishes the raw job)      │
                                                             ├──► comb/inputs.py  read, validate
                                                             │    comb/model.py   assign ρ per category
                                                             │    comb/blue.py    BLUE + χ² + weights
                                                             │    comb/likelihood.py  cross-check + scan
                                                             └──► comb/ratio.py   σ_ττ/σ_μμ, σ_ee/σ_μμ
```

Nothing is refitted and no event is reprocessed. Every input is a number committed to the
repository by the channel that measured it, and `run_combination.py --check` asserts that each one
still equals what that channel published — 26 assertions and 8 closure tests.

## The five documents

| page | what it answers |
|---|---|
| [01-inputs.md](01-inputs.md) | what each channel provides, the three places where they are not symmetric, and the z-ee normalisation problem in full |
| [02-correlation-model.md](02-correlation-model.md) | why each systematic gets ρ = 0 or ρ = 1 |
| [03-method.md](03-method.md) | BLUE, the profile-likelihood cross-check, and what a TRExFitter MultiFit would add |
| [04-results.md](04-results.md) | the numbers, the variations, and what would improve them |
| [05-vs-published.md](05-vs-published.md) | are the channels really orthogonal, which systematics CMS and ATLAS have that we do not, and each channel against the published measurement of the same decay |

## The honest summary

Three things, in order of importance.

**1. The three channels are not compatible at p = 0.008**, and the reason is a convention
violation in one of them, not physics. Everything downstream of that — the central value, the
+3.0σ ee/μμ "universality" ratio — inherits it. The fix is a one-line change in the z-ee notebook
(renormalise the LHE weight envelopes to a constant yield before writing the templates) and a
re-run; it is listed as item 1 of [01](01-inputs.md).

**2. The combination is luminosity-dominated.** 21 pb of the 32 pb total, fully correlated
between the three channels and therefore irreducible by combining. Adding z-ee to μμ ⊕ ττ did
improve the uncertainty — 35.4 → 31.6 pb, 10.7 % — because ee is the first channel of comparable
precision to μμ, and that is exactly what the μμ ⊕ ττ round could not do. But the floor is the
luminosity calibration, and no fourth channel moves it.

**3. What the ττ channel buys is the universality test, not the cross section.** Its weight is
−0.0003. R = σ(ττ)/σ(μμ) = 1.08 ± 0.13, in which the luminosity cancels exactly.
