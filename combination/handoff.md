# Handoff – Combination

## Team

- Lieke Gijsen
- Caspar van Eck

## What we did

We combined the two finished channels — `z-mumu` (Z → μ⁺μ⁻) and `z-tautau` (Z → τhτh) — into one
Z cross section. `z-ee` is not included (it has no result yet and was out of scope for this pass).

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1940 ± 35 pb**
> = 1940 ± 0.6 (stat) ± 23 (syst) ± 12 (acc) ± 24 (lumi) pb, per lepton flavour, assuming lepton
> universality.  χ²/ndf = 0.79/1 (p = 0.37), so the two channels are compatible.
>
> **R = σ(Z→ττ)/σ(Z→μμ) = 1.16 ± 0.18** (+0.9σ from 1) — the luminosity cancels exactly here.

Everything: [`result.md`](result.md) (generated), [`docs/`](docs/), [`README.md`](README.md),
[`CLAUDE.md`](CLAUDE.md). Presentation covering the whole analysis: [`../presentation/`](../presentation/).

## Method, in one paragraph

A covariance (BLUE) combination of the two channels' published TRExFitter v1.8.0 results, with
the shared systematics correlated through the `Category` strings of `../fitting/CONVENTIONS.md`
(ρ = 1 for luminosity, pileup, L1 prefiring, background cross sections, PDF/α_s/scale/PS;
ρ = 0 for the object calibrations, MC statistics, fakes and the generator comparison, which is a
different comparison in each channel). Cross-checked against an explicit profile likelihood that
carries the ττ channel's asymmetric uncertainty — same answer. Combination is done on the cross
sections, **not** on `mu_Z`, because the two channels normalise to references that differ by
0.47 % (see `docs/01-inputs.md`).

## Inputs we used

| channel | what we took | value |
|---|---|---|
| μμ | the **counting extraction** + 0.7 % lineshape term, per `z-mumu/REVIEW.md` F3 | μ_Z = 0.9935, σ = 1941 ± 35 pb |
| ττ | the `nominal` fake-factor variant | μ_Z = 1.160 ⁺⁰·¹⁹²₋₀.₁₆₃, σ = 2255 ± 358 pb |

`run_combination.py --check` asserts both against the channels' published files and fails if
either changes. **If a channel re-publishes, run that first** — it will tell you exactly what moved.

## The main message

The combination gains **0.13 %** on the uncertainty over μμ alone. That is not a failure of the
method: luminosity is 24 of the 35 pb and is 100 % correlated between channels, so it cannot be
combined away. The ττ channel is 10× less precise and gets a weight of −0.005 (a negative weight
is normal BLUE behaviour when ρ exceeds the ratio of the uncertainties). What the ττ channel
really contributes is the lepton-universality test.

## Not done

- **No TRExFitter MultiFit.** The channel workspaces are not committed and were not rebuilt, so
  the shared NPs are correlated by hand instead of profiled jointly. The config is generated and
  syntax-checked ([`fit/comb.config`](fit/comb.config)) and
  [`fit/run_multifit.py`](fit/run_multifit.py) has the full recipe — build the submodule, re-run
  `z-mumu --from 4` (~40 min) and `z-tautau --from 3` (~10 min), then one command.
- **No z-ee.** Adding it is a loader in `comb/inputs.py` and a `Fit` block in
  `fit/run_multifit.CHANNELS`, provided it follows `../fitting/CONVENTIONS.md`.
- The caveats we inherited and did not fix are listed at the end of `result.md`
  (`MuonReco` per-muon vs per-event, the pileup profile, the ττ fake-factor double counting).

## Next steps, in order of what they would change

1. Fix the μμ signal extraction (`z-mumu/REVIEW.md` F3/F4): the 0.7 % lineshape term is 13.7 pb,
   the second-largest entry in the whole breakdown.
2. Run the MultiFit and compare. The expected difference is small but it is not measured.
3. Agree one construction of σ^pred(60–120) between the channels, so a shared `mu_Z` means one
   thing (0.47 % today).
4. Add z-ee, and add eτh/μτh to ττ so the τh ID is constrained in situ.
