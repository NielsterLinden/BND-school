# 04 — Results

The numbers themselves live in [`../result.md`](../result.md), which is *generated* from
`output/combination_result.json` by `comb/report.py` so it can never drift. This page is the
commentary that does not belong in a generated file.

## The headline

> σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = **1931 ± 35 pb** (1.83 %)
> = 1931 ± 0.6 (stat) ± 24 (syst) ± 12 (acc) ± 23 (lumi) pb

χ²/ndf = 1.70/1, p = 0.19: the two channels are compatible. With only two measurements this is a
weak test — it is one number — but it is the test there is, and it passes.

Both inputs were re-published on 15 Sep 2026 (z-mumu after its review, z-tautau at v2.1). The
combination moved from 1940 ± 35 pb to 1931 ± 35 pb; essentially all of that −9 pb is the μμ
channel switching from the counting extraction to the repaired shape fit, which is exactly the
change the previous round asked for ([01](01-inputs.md)).

## Three things worth saying out loud

**1. The combination gains nothing, and that is the finding.** μμ alone gives 1930.7 ± 35.36 pb;
the combination gives 1930.6 ± 35.36 pb. The uncertainty improves by 0.00 % — less than it did
before, because the ττ channel's `Signal modelling` term shrank and with it the correlation that
gave the ττ measurement its (negative) lever. The reason is arithmetic: the luminosity term is
23 pb of the 35 pb total, it is 100 % correlated between the channels, and no amount of combining
reduces a fully correlated uncertainty. A third channel of the same precision would not change
this either. The way to a better Z cross section from this dataset is a better luminosity
calibration, not more channels.

**2. The measurement is not statistics-limited by four orders of magnitude.** The data
statistical uncertainty on the combined value is 0.6 pb — 0.03 %. The systematic is 40× larger.
`z-mumu`'s own review makes the sharper point (F9): even the *MC* statistical uncertainty
(0.39 %, the `Gammas` term) is ten times the data statistics. More collisions would change
nothing; more simulated events would change more.

**3. The ττ channel's value is high, and slightly more so than before.** σ(ττ) = 2343 pb against
σ(μμ) = 1931 pb, i.e. R = 1.21 ± 0.17, +1.3σ from unity (it was +0.9σ when the channel's nominal
was the unsubtracted fake factor; v2.1 made the MC-subtracted variant nominal, which is the higher
of the two). The τh ID scale factors (11.1 %) and the fake factor (8.6 %) now dominate its
uncertainty — the LO-vs-NLO generator systematic that used to be second largest is gone, dropped
by the channel as a placeholder. The channel's own recommendation points the same way as the
number: its DeepTau **Tight** re-run gives R = 1.08 ± 0.13 (12 % instead of 14 %), +0.6σ from
unity, and in a joint fit where μμ fixes μ_Z the ττ channel would measure `TauID_DM*` in situ.

## Comparison with published results

| measurement | σ·B(Z→ℓℓ) | window | L |
|---|---|---|---|
| this work (μμ + τhτh) | 1931 ± 35 pb | 60–120 GeV | 16.4 fb⁻¹ |
| CMS, [arXiv:2408.03744](https://arxiv.org/abs/2408.03744) (SMP-20-004) | 1952 ± 4 ± 18 ± 45 pb | 60–120 GeV | 206 pb⁻¹ |
| ATLAS, [arXiv:1603.09222](https://arxiv.org/abs/1603.09222) | 1981 ± 7 ± 38 ± 42 pb | 66–116 GeV | 81 pb⁻¹ |
| aMC@NLO reference used here | 1945–1954 pb | 60–120 GeV | — |

Agreement with the CMS measurement of exactly the same quantity is within 0.4σ. The ATLAS number
uses a narrower mass window and is shown for scale only, not as a comparison.

Two remarks on that table. First, the published measurements have a *larger* luminosity
uncertainty in relative terms than this work does (2.3 % and 2.1 % against 1.2 %), because they
are early-Run-2 measurements calibrated on small datasets; the CMS Open Data normtag value for
2016 G+H is much better known. Second, the agreement is not an independent check of very much:
the same NNLO reference cross section (6077.22 pb) normalises the simulation that both channels'
C factors come from.

## Variations

All nine variations are in `result.md` and `output/plots/variations.pdf`. The largest by far is
`mumu_counting` (+11.2 pb, +0.32σ): the 1-bin extraction instead of the 12 × 5 GeV fit. The two
alternative binnings the μμ channel accepts bracket the baseline by +3.8 and −4.5 pb, so the μμ
signal extraction alone spans 8.3 pb — more than the entire correlation model, whose bracketing
rows `rho_zero` and `rho_one` span 8.7 pb but sit symmetrically around the baseline. The two ττ
variations move the result by 0.1 and 0.5 pb: the ττ channel cannot move this number, whatever it
does. `sigmodel_correlated` is exactly null, because z-tautau no longer fits a generator nuisance
parameter for it to be correlated with.

## What would actually improve this

In order of effect on the combined number:

1. **Luminosity.** 23 of the 35 pb. Nothing in this repository can improve it; it is the CMS Open
   Data normtag value with its quoted ±1.2 %.
2. **The μμ muon efficiency** (17.0 pb): ID, isolation, trigger from in-house tag-and-probe plus
   the 0.8 %/event reconstruction term. Official Muon-POG scale factors would replace all of it
   (`z-mumu` open issue 4) but need CERN credentials.
3. **The μμ lineshape** (`z-mumu` open issue 3). It does not show up as a single line in the
   breakdown any more — the `Signal modelling (generator)` term is 7.7 pb — but it is what makes
   the result move by 8.3 pb between accepted binnings and by 11.2 pb if you count instead of fit.
   A third generator, or NLO electroweak corrections, would settle whether the data's 60–80 GeV
   shape is aMC@NLO's or powheg's.
4. **The acceptance PDF term** (10.0 pb), correlated between the channels and therefore not
   reducible by combining either.
5. Everything ττ contributes: 0.2 pb. Improving the ττ channel improves the *universality test*,
   not the cross section — and there the Tight working point already takes R from 1.21 ± 0.17 to
   1.08 ± 0.13.
