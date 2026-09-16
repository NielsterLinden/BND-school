# 03 — Method

## What is combined

σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) **per lepton flavour, assuming lepton universality**.
Each channel contributes σ_c = μ̂_c × σ_c^pred(60–120), and the uncertainties are assembled as
absolute values in pb:

- **in-fit systematics** are impacts on μ, and μ multiplies a fixed prediction, so they scale with
  σ^pred: δσ = s_group × σ^pred. This is the same bookkeeping the channels use (z-mumu:
  `sigma_fid_lumi_pb` = 0.0116267 × 799.566).
- **the acceptance** (μμ and ττ only) enters as σ = σ_fid / A, so its uncertainty is relative to
  the *measured* value: δσ = a × σ.

## BLUE

With x the three measurements, C the covariance and u = (1,1,1)ᵀ (`comb/blue.py`):

```
w = C⁻¹u / (uᵀC⁻¹u)      σ̂ = wᵀx      Var(σ̂) = 1/(uᵀC⁻¹u)
χ² = (x − σ̂u)ᵀ C⁻¹ (x − σ̂u),   ndf = 2
```

C is built source by source, `C = Σ_s C_s` with `C_s = ρ_s · s_i s_j` off-diagonal and `s_i²` on
the diagonal, so the combined uncertainty can be decomposed the same way:
`δ_s = sqrt(wᵀ C_s w)`, and the δ_s add in quadrature back to the total. That is what the
breakdown table in `result.md` is.

**The χ² is now a real test.** With two channels it was one number with one degree of freedom
(1.70/1, p = 0.19). With three it has two, and it gives 9.67/2, **p = 0.008**. Almost all of it is
the ee–μμ pair: they are 123 pb apart with a largely uncorrelated difference. BLUE does not inflate
the uncertainty for a poor χ² (no PDG scale factor is applied here), so the quoted ±32 pb is the
uncertainty *assuming* the inputs are consistent — which at p = 0.008 they are not.
[01](01-inputs.md) identifies the cause; `ee_normfix` takes the χ² to 7.4/2 (p = 0.025), still
poor, because the corrected ee value then sits 98 pb *below* μμ. Either way the two precise
channels do not agree, and that is the finding.

**Iteration.** The acceptance terms are multiplicative, so evaluating them at each channel's own
measured value biases the weights toward whichever channel fluctuated low. `combine` re-evaluates
them at the combined value and re-solves until stable (the iterative BLUE of the LHC top-mass
combinations). It converges in a handful of iterations and moves the result by well under a pb,
because the additive in-fit terms dominate.

**The negative weight.** The ττ weight is −0.0003. This is standard BLUE behaviour whenever the
correlation with a more precise measurement exceeds the ratio of the two uncertainties: the less
precise measurement stops being an average-down and becomes a lever on the shared systematic. It
is mathematically correct, and it is now so small that it moves the central value by 0.06 pb.
μμ and ee both get positive weights (+0.596, +0.404) because ρ(μμ, ee) = 0.44 is well below their
uncertainty ratio, 0.89.

## Profile-likelihood cross-check

`comb/likelihood.py` writes the same combination as an explicit likelihood with one nuisance
parameter per shared source:

```
−2 ln L(σ, θ) = Σ_c (x_c − σ − Σ_s d_{c,s} θ_s)² / u_c²  +  Σ_s θ_s²
```

with u_c the quadrature sum of that channel's uncorrelated sources. For symmetric Gaussian errors
and ρ ∈ {0,1} this is algebraically identical to BLUE, which `run_combination.py --check` asserts
(agreement better than 10⁻³ pb on the value and 2 % on the uncertainty). Two things it buys:

1. the ττ channel's asymmetric uncertainty (+10.7 / −9.3 % on μ_Z) can be carried through with a
   bifurcated response instead of being symmetrised. It changes the combined value by 0.3 pb;
2. a −2Δln L curve to plot (`output/plots/likelihood_scan.pdf`), which shows the parabola and the
   likelihood lying on top of each other — the honest way to say "the asymmetry does not matter
   here".

A source with 0 < ρ < 1 is split into a shared part √ρ·s and an independent part √(1−ρ)·s, which
reproduces its covariance exactly; every ρ in this combination is 0 or 1, so nothing is split.

## Lepton universality

Two ratios, with every ρ = 1 source cancelling:

```
(δR/R)² = Σ_s [ r_X,s² + r_μ,s² − 2 ρ_s r_X,s r_μ,s ]
```

| ratio | value | z from 1 |
|---|---:|---:|
| σ(ττ)/σ(μμ) | 1.079 ± 0.130 | +0.6 |
| σ(ee)/σ(μμ) | 1.064 ± 0.021 | **+3.0** |

Luminosity, pileup, prefiring and the correlated theory terms drop out exactly in both. What
remains in the ττ ratio is that channel's own τh ID, fake factor and MC statistics; what remains
in the ee ratio is the muon efficiency, the ee `Gammas` and `Fit residual`, and the ee electron
efficiency — 2.0 % in total, which is why a 6.4 % difference reads as 3.0σ.

**The ee/μμ ratio is not a measurement of lepton universality.** It is the sharpest available
statement that the two fits disagree, and [01](01-inputs.md) says why. Quote it as a consistency
check between the channels, not as a limit on e–μ universality.

## What a TRExFitter MultiFit would add

The delivered combination is **not** a MultiFit (see `result.md`, "What this is not"). The
configuration is generated and checked by `fit/run_multifit.py`; running it needs the channel
workspaces, which are not all committed. What the joint fit would do that this cannot:

| | covariance combination | MultiFit |
|---|---|---|
| shared NPs | correlated through `Category` totals, using each channel's *own* post-fit impacts | profiled jointly: the μμ and ee data constrain `Lumi`, `Pileup`, `L1Prefiring` and the theory NPs for every channel |
| the z-ee `QCDScale` normalisation | shows up as a −99 pb variation and a p = 0.008 χ² | shows up immediately: one `QCDScale` parameter shared with two channels that renormalise it correctly cannot absorb an ee-only normalisation |
| ττ τh ID | an uncorrelated 12.2 % | could be constrained in situ if eτh/μτh regions were added |
| MC statistics | one `Gammas` number per channel | per-bin γ parameters, correctly uncorrelated bin by bin |
| asymmetry | bifurcated-Gaussian approximation | exact |
| the reference mismatch | handled: each μ̂ multiplied by its own reference | **not** handled: one `mu_Z` against three references |

Expected size of the difference: no longer obviously small. With two channels the μμ fit carried
essentially all the weight and was already profiled; with three, ee carries 40 % and its fit has
a parameter that a joint fit would treat very differently. Running the MultiFit is now the top
item on the "next steps" list, right after z-ee re-publishing.
