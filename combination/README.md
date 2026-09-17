# Combination — Z → ee ⊕ Z → μμ ⊕ Z → τhτh

One TRExFitter v1.8.0 MultiFit of the three BND-school channels on CMS Open Data 2016 G+H
(16.4 fb⁻¹, √s = 13 TeV).

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1948 ⁺³³₋₃₂ pb** per lepton flavour
> (aMC@NLO: 1954 ⁺⁵⁶₋₈₂ pb; channel compatibility p = 0.30)

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
