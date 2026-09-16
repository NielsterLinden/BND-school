# 04 — Results

The numbers themselves live in [`../result.md`](../result.md), which is *generated* from
`output/combination_result.json` by `comb/report.py` so it can never drift. This page is the
commentary that does not belong in a generated file.

## The headline

> σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = **1980 ± 32 pb** (1.59 %)
> = 1980 ± 0.5 (stat) ± 22 (syst) ± 7 (acc) ± 21 (lumi) pb
>
> χ²/ndf = **9.67/2, p = 0.008**

Quote the two together. The central value is what BLUE returns from three inputs that are not
mutually consistent, and no scale factor has been applied to the uncertainty to hide that.

## What changed from the μμ ⊕ ττ round

| | μμ ⊕ ττ (15 Sep) | μμ ⊕ ττ ⊕ ee (16 Sep) |
|---|---|---|
| σ | 1931 ± 35 pb | **1980 ± 32 pb** |
| χ²/ndf (p) | 1.70/1 (0.19) | **9.67/2 (0.008)** |
| gain over the best channel | 0.00 % | **10.7 %** |
| ττ input | v2.1, DeepTau Medium, μ_Z = 1.205 | **v3, DeepTau Tight, μ_Z = 1.071** |
| largest variation | `mumu counting`, +11.2 pb | **`ee normfix`, −98.5 pb** |

Two things moved. z-tautau went to v3 and adopted the Tight working point it had recommended,
which brought σ(ττ) from 2343 to 2082 pb and R(ττ/μμ) from 1.21 ± 0.17 to 1.08 ± 0.13 — the
channel's own prediction of what that change would do, confirmed. And z-ee arrived, which is
responsible for the +49 pb and for all of the χ².

## Three things worth saying out loud

**1. Combining finally buys something — and it is not the finding this time.** μμ alone gives
1931 ± 35.4 pb; the three-channel combination gives 1980 ± 31.6 pb, a 10.7 % improvement. That is
what the previous round said would happen as soon as a second channel of comparable precision
appeared, and ee is it (1.9 % against μμ's 1.8 %). The improvement is limited to 10.7 % rather
than the 29 % two independent equal measurements would give, because ρ(μμ, ee) = 0.44 — almost
all of it luminosity, 22.7 × 18.8 pb of the covariance. The floor is unchanged: 21 of the 32 pb is
the luminosity calibration.

**2. The two precise channels do not agree, and the reason is fixable.** σ(ee) = 2054 ± 40 pb
against σ(μμ) = 1931 ± 35 pb, 123 pb apart with a 2.0 % uncorrelated difference: 3.0σ. That is not
a lepton-universality result; it is the z-ee fit's ±5.9 % un-renormalised `QCDScale` normalisation
being pulled −1.88σ, so that `mu_signal` is measured against a prediction the fit itself rescaled
by 0.893 (see [01](01-inputs.md) for the arithmetic and the closure test). Undoing only that gives
σ(ee) = 1833 pb and a combination of 1882 ± 30 pb — and the disagreement flips sign, 2.7σ the
other way, so "just apply the correction" is not the answer either. The answer is a z-ee re-run
with renormalised templates, coarser bins and a goodness of fit. Until then the honest reading of
this measurement is **1931 ± 35 pb from μμ ⊕ ττ, with an ee channel that does not yet agree with
it.**

**3. The measurement is not statistics-limited by four orders of magnitude.** The data statistical
uncertainty on the combined value is 0.5 pb — 0.02 %. The systematic is 44× larger. `z-mumu`'s own
review makes the sharper point (F9): even the *MC* statistical uncertainty (the `Gammas` term,
10.6 pb) is twenty times the data statistics, and on the ee side `Gammas` is 23.8 pb on its own.
More collisions would change nothing; more simulated events would change more.

## Comparison with published results

The per-channel comparison, the orthogonality check and the list of systematics CMS and ATLAS
have that we do not are in [05-vs-published.md](05-vs-published.md). The short version: each of our
two precise channels agrees with the ATLAS measurement of the *same decay* within 1σ, and yet they
disagree with each other by 3.0σ — ATLAS's per-channel uncertainty is ±60 pb against our ±35 and
±40 pb, so a 123 pb gap hides comfortably inside the published error bars.


| measurement | σ·B(Z→ℓℓ) | window | L |
|---|---|---|---|
| this work (μμ + τhτh + ee) | 1980 ± 32 pb | 60–120 GeV | 16.4 fb⁻¹ |
| this work, μμ + τhτh only | 1931 ± 35 pb | 60–120 GeV | 16.4 fb⁻¹ |
| this work, with `ee normfix` | 1882 ± 30 pb | 60–120 GeV | 16.4 fb⁻¹ |
| CMS, [arXiv:2408.03744](https://arxiv.org/abs/2408.03744) (SMP-20-004) | 1952 ± 4 ± 18 ± 45 pb | 60–120 GeV | 206 pb⁻¹ |
| ATLAS, [arXiv:1603.09222](https://arxiv.org/abs/1603.09222) | 1981 ± 7 ± 38 ± 42 pb | 66–116 GeV | 81 pb⁻¹ |
| aMC@NLO reference used here | 1945–1954 pb | 60–120 GeV | — |

The baseline agrees with the CMS measurement of exactly the same quantity within 0.5σ. So does the
μμ ⊕ ττ value (0.4σ) and so, nearly, does the `ee normfix` value (1.2σ) — which is a good
reminder that agreement with a published number at this precision does not discriminate between
the three, and should not be used to pick one. The ATLAS number uses a narrower mass window and is
shown for scale only.

Two remarks. First, the published measurements have a *larger* luminosity uncertainty in relative
terms than this work does (2.3 % and 2.1 % against 1.2 %), because they are early-Run-2
measurements calibrated on small datasets. Second, the agreement is not an independent check of
very much: the same NNLO reference cross section (6077.22 pb) normalises the simulation that all
three channels' C factors come from.

## Variations

Ten variations, in `result.md` and `output/plots/variations.pdf`. The table is now dominated by
two rows that are not modelling choices at all:

| | shift | what it says |
|---|---:|---|
| `ee normfix` | −98.5 pb (−3.12σ) | one channel's signal template normalisation |
| `no ee` | −49.2 pb (−1.56σ) | the size of the ee contribution |
| `mumu counting` | +11.1 pb (+0.35σ) | the μμ signal extraction |
| `rho zero` / `rho one` | +6.2 / −3.7 pb | the *entire* correlation model |
| everything else | < 2.1 pb | — |

The correlation model spans 9.9 pb, a third of the total uncertainty and a tenth of the z-ee
normalisation effect. It is not what limits this measurement.

## What would actually improve this

In order of effect on the combined number:

1. **z-ee re-runs with renormalised `PDF`/`QCDScale` templates.** 99 pb, and it is the difference
   between a χ² of 9.7/2 and something that can be quoted. The full list is [01](01-inputs.md).
   Everything below is second-order until this is done.
2. **Luminosity.** 21 of the 32 pb. Nothing in this repository can improve it; it is the CMS Open
   Data normtag value with its quoted ±1.2 %.
3. **Run the MultiFit.** With ee carrying 40 % of the weight the covariance approximation is no
   longer obviously adequate, and a joint fit would expose the `QCDScale` problem structurally
   rather than as a variation ([03](03-method.md)).
4. **The μμ muon efficiency** (10.1 pb) and the **ee MC statistics** (`Gammas`, 23.8 pb in the
   channel): official Muon-POG scale factors would replace the first (`z-mumu` open issue 4), more
   simulated events the second.
5. **The μμ lineshape** (`z-mumu` open issue 3), worth 11.1 pb between counting and fitting.
6. Everything ττ contributes to the cross section: 0.06 pb. Improving that channel improves the
   *universality test*, not the cross section — and v3 already took R from 1.21 ± 0.17 to
   1.08 ± 0.13 by adopting the Tight working point.
