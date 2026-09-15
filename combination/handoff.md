# Handoff – Combination

## Team

- Lieke Gijsen
- Caspar van Eck

## What we did

We combined the two finished channels — `z-mumu` (Z → μ⁺μ⁻) and `z-tautau` (Z → τhτh) — into one
Z cross section. `z-ee` is not included (it has no result yet and was out of scope for this pass).

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1931 ± 35 pb**
> = 1931 ± 0.6 (stat) ± 24 (syst) ± 12 (acc) ± 23 (lumi) pb, per lepton flavour, assuming lepton
> universality.  χ²/ndf = 1.70/1 (p = 0.19), so the two channels are compatible.
>
> **R = σ(Z→ττ)/σ(Z→μμ) = 1.21 ± 0.17** (+1.3σ from 1) — the luminosity cancels exactly here.

Updated on 15 Sep 2026 after both channels re-published: z-mumu fixed the findings of its own
review and z-tautau moved to v2.1. The previous round gave 1940 ± 35 pb; the −9 pb is almost
entirely the μμ channel switching from the counting extraction to its repaired shape fit.

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
| μμ | the **v2 shape fit** (12 × 5 GeV, two-sided `SigModel`, MINOS), i.e. the channel's own baseline after its review fixes | μ_Z = 0.9881 ± 0.0159, σ = 1931 ± 35 pb |
| ττ | **v2.1**, the variant the channel's `nominal_variant` names (`mcsub`, the MC-subtracted fake factor) | μ_Z = 1.205 ⁺⁰·¹⁴⁵₋₀.₁₂₅, σ = 2343 ± 322 pb |

The ±0.7 % lineshape term the previous round carried is gone: the channel's two-sided `SigModel`
template now covers it inside the fit, and the fit is stable to ±0.2 % between the binnings that
describe the data. The counting extraction and the two alternative binnings survive as variations
(+11.2, +3.8, −4.5 pb).

`run_combination.py --check` asserts both against the channels' published files and fails if
either changes. **If a channel re-publishes, run that first** — it will tell you exactly what moved.

## The main message

The combination gains **0.00 %** on the uncertainty over μμ alone. That is not a failure of the
method: luminosity is 23 of the 35 pb and is 100 % correlated between channels, so it cannot be
combined away. The ττ channel is 9× less precise and gets a weight of −0.0004 (a negative weight
is normal BLUE behaviour when ρ exceeds the ratio of the uncertainties, 0.114 > 0.110 here). What
the ττ channel really contributes is the lepton-universality test.

## Not done

- **No TRExFitter MultiFit.** The channel workspaces are not committed and were not rebuilt, so
  the shared NPs are correlated by hand instead of profiled jointly. The config is generated and
  syntax-checked ([`fit/comb.config`](fit/comb.config)) and
  [`fit/run_multifit.py`](fit/run_multifit.py) has the full recipe — build the submodule, re-run
  `z-mumu --from 4` (~40 min) and `z-tautau --from 3` (~25 min), then one command.
- **No z-ee.** Adding it is a loader in `comb/inputs.py` and a `Fit` block in
  `fit/run_multifit.CHANNELS`, provided it follows `../fitting/CONVENTIONS.md`.
- The caveats we inherited and did not fix are listed at the end of `result.md`. All nine findings
  of `z-mumu/REVIEW.md` were fixed by that channel itself before this round, so what is left is
  each channel's own open-issue list (the μμ 60–80 GeV lineshape, the pileup profile, the in-house
  muon/prefiring calibrations, and z-tautau's recommendation to adopt the Tight working point).

## Next steps, in order of what they would change

1. Settle the μμ 60–80 GeV lineshape (`z-mumu` open issue 3). It is what makes the result move by
   8.3 pb between the binnings the channel accepts and by 11.2 pb if you count instead of fit —
   more than the entire correlation model.
2. Run the MultiFit and compare. The expected difference is small but it is not measured.
3. Agree one construction of σ^pred(60–120) between the channels, so a shared `mu_Z` means one
   thing (0.47 % today).
4. Take z-tautau's Tight working point once that channel adopts it: it moves the combination by
   only +0.5 pb but takes R from 1.21 ± 0.17 to 1.08 ± 0.13.
5. Add z-ee, and add eτh/μτh to ττ so the τh ID is constrained in situ.
