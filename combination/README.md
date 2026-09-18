# Combination — Z → ee ⊕ Z → μμ ⊕ Z → ττ

One TRExFitter v1.8.0 MultiFit of the three BND-school channels on CMS Open Data 2016 G+H
(16.4 fb⁻¹, √s = 13 TeV).

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1949 ⁺³⁰₋₃₀ pb** per lepton flavour
> (aMC@NLO: 1954 ⁺⁵⁶₋₈₂ pb; channel compatibility p = 0.38)

Z → ττ is z-tautau's four-channel measurement (τhτh + μτh + eτh + eμ); the fits run on HTCondor (17 Sep 2026).
This supersedes the frozen three-channel combination with τhτh only, 1945 ⁺³¹₋₃₀ pb, which is at the git tag
`zmumu-freeze-2026-09-17` ([`../docs/FREEZE.md`](../docs/FREEZE.md)). The ee and μμ inputs are the frozen ones.

Everything is in [`combLieke/`](combLieke/):

| | |
|---|---|
| [`combLieke/README.md`](combLieke/README.md) | result, likelihood model, the ee channel, how to run |
| [`combLieke/output/plots/`](combLieke/output/plots/) | figures (PDF for slides, PNG) |
| [`combLieke/output/result.json`](combLieke/output/result.json) | every number |
| [`combLieke/docs/systematics.md`](combLieke/docs/systematics.md) | systematics compared with ATLAS and CMS |
| [`combLieke/docs/orthogonality.md`](combLieke/docs/orthogonality.md) | event-level overlap of the three channels, on data |
| [`combLieke/docs/four-channel-tautau.md`](combLieke/docs/four-channel-tautau.md) | how the four-channel ττ was brought in |
| [`handoff.md`](handoff.md), [`CLAUDE.md`](CLAUDE.md) | team handoff, agent guide |

```bash
source ../fitting/setup.sh && cd combLieke && python run.py prepare && python run.py condor --submit   # then: results, plots
```
