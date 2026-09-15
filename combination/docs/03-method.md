# 03 — Method

## What is combined

σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) **per lepton flavour, assuming lepton universality**.
Each channel contributes σ_c = μ̂_c × σ_c^pred(60–120), and the uncertainties are assembled as
absolute values in pb:

- **in-fit systematics** are impacts on μ, and μ multiplies a fixed prediction, so they scale with
  σ^pred: δσ = s_group × σ^pred. This is the same bookkeeping the channels use (z-mumu:
  `sigma_fid_lumi_pb` = 0.0116431 × 799.566).
- **the acceptance** enters as σ = σ_fid / A, so its uncertainty is relative to the *measured*
  value: δσ = a × σ.

## BLUE

With x the two measurements, C the covariance and u = (1,1)ᵀ (`comb/blue.py`):

```
w = C⁻¹u / (uᵀC⁻¹u)      σ̂ = wᵀx      Var(σ̂) = 1/(uᵀC⁻¹u)
χ² = (x − σ̂u)ᵀ C⁻¹ (x − σ̂u),   ndf = 1
```

C is built source by source, `C = Σ_s C_s` with `C_s = ρ_s · s_i s_j` off-diagonal and `s_i²` on
the diagonal, so the combined uncertainty can be decomposed the same way:
`δ_s = sqrt(wᵀ C_s w)`, and the δ_s add in quadrature back to the total. That is what the
breakdown table in `result.md` is.

**Iteration.** The acceptance terms are multiplicative, so evaluating them at each channel's own
measured value biases the weights toward whichever channel fluctuated low. `combine` re-evaluates
them at the combined value and re-solves until stable (the iterative BLUE of the LHC top-mass
combinations). It converges in 4 iterations and moves the result by less than 0.1 pb here, because
the additive in-fit terms dominate.

**The negative weight.** The ττ weight is −0.0004. This is standard BLUE behaviour whenever the
correlation exceeds the ratio of the two uncertainties (0.114 > 0.110 here): the less precise
measurement stops being an average-down and becomes a lever on the shared systematic. It is
mathematically correct and it is now barely there — the two numbers that decide its sign are 4 %
apart, and it moves the central value by 0.2 pb — but it is worth knowing that the combined value
sits slightly *below* the μμ value even though the ττ value is above it. It went from −0.005 in
the previous round to −0.0004 because z-tautau's `Signal modelling` category shrank when the
channel stopped fitting its LO-vs-NLO generator systematic, which took most of what the two
channels had in common out of the ττ uncertainty ([02](02-correlation-model.md)).

## Profile-likelihood cross-check

`comb/likelihood.py` writes the same combination as an explicit likelihood with one nuisance
parameter per shared source:

```
−2 ln L(σ, θ) = Σ_c (x_c − σ − Σ_s d_{c,s} θ_s)² / u_c²  +  Σ_s θ_s²
```

with u_c the quadrature sum of that channel's uncorrelated sources. For symmetric Gaussian errors
and ρ ∈ {0,1} this is algebraically identical to BLUE, which `run_combination.py --check` asserts
(agreement better than 10⁻³ pb on the value and 2 % on the uncertainty). Two things it buys:

1. the ττ channel's asymmetric uncertainty (+12.1 / −10.4 % on μ_Z) can be carried through with a
   bifurcated response instead of being symmetrised. It changes the combined value by 0.7 pb;
2. a −2Δln L curve to plot (`output/plots/likelihood_scan.pdf`), which shows the parabola and the
   likelihood lying on top of each other — the honest way to say "the asymmetry does not matter
   here".

## Lepton universality

R = σ(ττ)/σ(μμ), with every ρ = 1 source cancelling:

```
(δR/R)² = Σ_s [ r_τ,s² + r_μ,s² − 2 ρ_s r_τ,s r_μ,s ]
```

Luminosity, pileup, prefiring and the correlated theory terms drop out exactly. What remains is
dominated by the ττ channel's own τh ID (12.5 %) and generator (8.5 %) uncertainties.

## What a TRExFitter MultiFit would add

The delivered combination is **not** a MultiFit (see `result.md`, "What this is not"). The
configuration is generated and checked by `fit/run_multifit.py`; running it needs the channel
workspaces, which are not committed. What the joint fit would do that this cannot:

| | covariance combination | MultiFit |
|---|---|---|
| shared NPs | correlated through `Category` totals, using each channel's *own* post-fit impacts | profiled jointly: the μμ data can constrain `Lumi`, `Pileup`, `L1Prefiring` and the theory NPs *for the ττ channel too* |
| ττ τh ID | an uncorrelated 12.5 % | could be constrained in situ if eτh/μτh regions were added (the ττ channel's open issue 2) |
| MC statistics | one `Gammas` number per channel | per-bin γ parameters, correctly uncorrelated bin by bin |
| asymmetry | bifurcated-Gaussian approximation | exact |
| the 0.47 % reference mismatch | handled: each μ̂ multiplied by its own reference | **not** handled: one `mu_Z` against two references |

Expected size of the difference: small. The μμ channel carries 99.5 % of the weight and its own
fit is already a profile likelihood, so the only real gain is the ττ channel's nuisance
parameters being constrained by μμ data — and those are the parameters μμ has no sensitivity to
(τh ID, τh trigger, τh energy scale). The claim is not made, which is why this page exists.
