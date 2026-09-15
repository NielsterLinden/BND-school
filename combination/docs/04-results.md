# 04 — Results

The numbers themselves live in [`../result.md`](../result.md), which is *generated* from
`output/combination_result.json` by `comb/report.py` so it can never drift. This page is the
commentary that does not belong in a generated file.

## The headline

> σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = **1940 ± 35 pb** (1.79 %)
> = 1940 ± 0.6 (stat) ± 23 (syst) ± 12 (acc) ± 24 (lumi) pb

χ²/ndf = 0.79/1, p = 0.37: the two channels are compatible. With only two measurements this is a
weak test — it is one number — but it is the test there is, and it passes.

## Three things worth saying out loud

**1. The combination gains nothing, and that is the finding.** μμ alone gives 1941 ± 35 pb; the
combination gives 1940 ± 35 pb. The uncertainty improves by 0.13 %. The reason is arithmetic: the
luminosity term is 24 pb of the 35 pb total, it is 100 % correlated between the channels, and no
amount of combining reduces a fully correlated uncertainty. A third channel of the same precision
would not change this either. The way to a better Z cross section from this dataset is a better
luminosity calibration, not more channels.

**2. The measurement is not statistics-limited by four orders of magnitude.** The data
statistical uncertainty on the combined value is 0.6 pb — 0.03 %. The systematic is 40× larger.
`z-mumu`'s own review makes the sharper point (F9): even the *MC* statistical uncertainty
(0.35 %, the `Gammas` term) is ten times the data statistics. More collisions would change
nothing; more simulated events would change more.

**3. The ττ channel's value is high but not significantly so.** σ(ττ) = 2255 pb against
σ(μμ) = 1941 pb, i.e. R = 1.16 ± 0.18, +0.9σ from unity. The ττ channel's own open issues point
the same way: its nominal fake factor does not subtract genuine taus (its open issue 1), and
using the `mcsub` variant moves it *further* up, to 2348 pb. The τh ID scale factors (12.5 %) and
the LO-vs-NLO generator systematic (8.5 %) dominate its uncertainty, and the channel's own
recommendation is that both would shrink in a joint fit where μμ fixes μ_Z and ττ measures the
τh ID in situ.

## Comparison with published results

| measurement | σ·B(Z→ℓℓ) | window | L |
|---|---|---|---|
| this work (μμ + τhτh) | 1940 ± 35 pb | 60–120 GeV | 16.4 fb⁻¹ |
| CMS, [arXiv:2408.03744](https://arxiv.org/abs/2408.03744) (SMP-20-004) | 1952 ± 4 ± 18 ± 45 pb | 60–120 GeV | 206 pb⁻¹ |
| ATLAS, [arXiv:1603.09222](https://arxiv.org/abs/1603.09222) | 1981 ± 7 ± 38 ± 42 pb | 66–116 GeV | 81 pb⁻¹ |
| aMC@NLO reference used here | 1945–1954 pb | 60–120 GeV | — |

Agreement with the CMS measurement of exactly the same quantity is within 0.2σ. The ATLAS number
uses a narrower mass window and is shown for scale only, not as a comparison.

Two remarks on that table. First, the published measurements have a *larger* luminosity
uncertainty in relative terms than this work does (2.3 % and 2.1 % against 1.2 %), because they
are early-Run-2 measurements calibrated on small datasets; the CMS Open Data normtag value for
2016 G+H is much better known. Second, the agreement is not an independent check of very much:
the same NNLO reference cross section (6077.22 pb) normalises the simulation that both channels'
C factors come from.

## Variations

All seven variations are in `result.md` and `output/plots/variations.pdf`. The largest is
`mumu_shapefit` (−6.5 pb, −0.19σ): using the μμ channel's 30-bin fit instead of the reviewed
counting extraction. Everything else is below 0.15σ. The two rows that bracket the entire
correlation model — `rho_zero` and `rho_one` — span 6.5 pb, so the modelling choices documented
in [02](02-correlation-model.md) are worth about a fifth of the total uncertainty between them.

## What would actually improve this

In order of effect on the combined number:

1. **Luminosity.** 24 of the 35 pb. Nothing in this repository can improve it; it is the CMS Open
   Data normtag value with its quoted ±1.2 %.
2. **Fix the μμ signal extraction** (`z-mumu/REVIEW.md` F3/F4). The ±0.7 % lineshape term is the
   second-largest entry in the breakdown (13.7 pb) and it exists only because the `SigModel`
   template mixes a generator comparison with a phase-space artefact. Rebuilding that template
   with a common LHE mass window would likely remove most of it.
3. **The μμ muon efficiency** (11.7 pb), and the unresolved per-muon/per-event `MuonReco`
   ambiguity (F8) that could add 0.35 % more.
4. **The acceptance PDF term** (10.0 pb), correlated between the channels and therefore not
   reducible by combining either.
5. Everything ττ contributes: 1.2 pb. Improving the ττ channel improves the *universality test*,
   not the cross section.
