# Combination — Z → ee ⊕ Z → μμ ⊕ Z → τhτh

One TRExFitter v1.8.0 MultiFit of the three BND-school channels on CMS Open Data 2016 G+H
(16.4 fb⁻¹, √s = 13 TeV).

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1945 ⁺³¹₋₃₀ pb** per lepton flavour
> (aMC@NLO: 1954 ⁺⁵⁶₋₈₂ pb; channel compatibility p = 0.30)

Frozen on 17 Sep 2026 (git tag `zmumu-freeze-2026-09-17`, [`../FREEZE.md`](../FREEZE.md)); read and
present these numbers, do not rerun.

> **After the freeze: the four-channel Z → ττ.** z-tautau's τhτh + μτh + eτh + eμ measurement has been combined with ee and μμ
> (17 Sep, fits on HTCondor): **1951 ⁺³⁰₋₂₉ pb** (interim), ττ alone 1981 ⁺⁷³₋₇₀ pb, compatibility p = 0.35. Method, decisions,
> checks: [`combLieke/docs/four-channel-tautau.md`](combLieke/docs/four-channel-tautau.md); numbers: `combLieke/interim/result.json`. Everything below, and
> `combLieke/output/`, is still the **frozen** ee ⊕ μμ ⊕ τhτh result until z-tautau's own fit results are in and
> `python run.py results && python run.py plots` has been run.

Everything is in [`combLieke/`](combLieke/):

| | |
|---|---|
| [`combLieke/README.md`](combLieke/README.md) | result, likelihood model, the ee channel, how to run |
| [`combLieke/output/plots/`](combLieke/output/plots/) | figures (PDF for slides, PNG) |
| [`combLieke/output/result.json`](combLieke/output/result.json) | every number |
| [`combLieke/docs/systematics.md`](combLieke/docs/systematics.md) | systematics compared with ATLAS and CMS |
| [`combLieke/docs/orthogonality.md`](combLieke/docs/orthogonality.md) | event-level overlap of the three channels, on data |
| [`handoff.md`](handoff.md), [`CLAUDE.md`](CLAUDE.md) | team handoff, agent guide |

```bash
source ../setup.sh && cd combLieke && python run.py all
```
