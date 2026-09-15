# 00 — What the combination is

Two of the three BND-school channels finished a Z cross-section measurement on the same CMS Open
Data 2016 sample (Run2016G+H, 16 393.381 pb⁻¹): `z-mumu` (Z → μ⁺μ⁻) and `z-tautau` (Z → τhτh).
`z-ee` is out of scope. This folder turns the two into one number.

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1931 ± 35 pb**, per lepton flavour, assuming lepton
> universality. See [`../result.md`](../result.md) for the full breakdown.

The inputs are each channel's current baseline as re-published on 15 Sep 2026: z-mumu's fit
rebuilt after its own review (12 × 5 GeV bins, two-sided `SigModel`) and z-tautau v2.1 (MC-
subtracted fake factor as nominal, OS/SS correction per BDT category).

## The chain

```
z-mumu   TRExFitter v1.8.0 fit of m(μμ)  ──►  μ_Z, grouped impacts, A, σ(60–120)
                                                     │
                                                     ├──►  comb/inputs.py   read, validate
z-tautau TRExFitter v1.8.0 fit of m(ττ)  ──►  ...    │     comb/model.py    assign ρ per category
                                                     │     comb/blue.py     BLUE + χ² + weights
                                                     │     comb/likelihood.py  cross-check + scan
                                                     └──►  comb/ratio.py    σ_ττ/σ_μμ
```

Nothing is refitted and no event is reprocessed. Every input is a number committed to the
repository by the channel that measured it, and `run_combination.py --check` asserts that each one
still equals what that channel published.

## The four documents

| page | what it answers |
|---|---|
| [01-inputs.md](01-inputs.md) | what each channel provides, and the two places where the two channels are not symmetric |
| [02-correlation-model.md](02-correlation-model.md) | why each systematic gets ρ = 0 or ρ = 1 |
| [03-method.md](03-method.md) | BLUE, the profile-likelihood cross-check, and what a TRExFitter MultiFit would add |
| [04-results.md](04-results.md) | the numbers, the variations, and how they compare to CMS and ATLAS |

## The honest summary

The combination is **luminosity-dominated** (1.2 %, fully correlated) and the μμ channel is nine
times more precise than the ττ channel, so combining gains 0.00 % on the uncertainty. That is not
a failure of the method: it is the result. What the ττ channel *does* buy is a lepton-universality
test, R = σ(ττ)/σ(μμ) = 1.21 ± 0.17, in which the luminosity cancels exactly.
