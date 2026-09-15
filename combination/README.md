# Combination — Z → μμ ⊕ Z → τhτh

One cross section from the two finished BND-school channels, on CMS Open Data 2016
(Run2016G+H, 16 393.381 pb⁻¹, √s = 13 TeV).

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1940 ± 35 pb**
> = 1940 ± 0.6 (stat) ± 23 (syst) ± 12 (acc) ± 24 (lumi) pb, per lepton flavour,
> assuming lepton universality.  χ²/ndf = 0.79/1 (p = 0.37).
>
> Lepton universality: R = σ(Z→ττ)/σ(Z→μμ) = **1.16 ± 0.18** (+0.9σ from 1).

Full numbers and caveats: **[`result.md`](result.md)**. The physics write-up is in
[`docs/`](docs/); start with [`docs/00-overview.md`](docs/00-overview.md).

## Run it

```bash
source ../setup.sh
python run_combination.py            # checks, result.md, output/*.json, output/plots/*
python run_combination.py --check    # input validation and closure tests only (2 s)
python run_combination.py --no-plots # numbers only
python fit/run_multifit.py --dry-run # regenerate and syntax-check the TRExFitter MultiFit config
```

There is nothing to reprocess: every input is a number the channel groups committed. The whole
thing runs in a few seconds.

## What this is — and is not

It **is** a covariance (BLUE) combination of two TRExFitter v1.8.0 profile-likelihood fits, with
the shared systematics correlated through the `Category` strings that `fitting/CONVENTIONS.md`
forces both channels to use, plus a profile-likelihood cross-check and a lepton-universality
ratio.

It is **not** a TRExFitter MultiFit. The channel workspaces are not committed and were not
rebuilt, so the shared nuisance parameters are correlated by hand rather than profiled jointly.
The MultiFit configuration is generated and syntax-checked in
[`fit/comb.config`](fit/comb.config); [`fit/run_multifit.py`](fit/run_multifit.py) carries the
exact recipe to produce the workspaces and run it, and
[`docs/03-method.md`](docs/03-method.md) says what it would add.

## Layout

| path | what |
|---|---|
| `run_combination.py` | orchestrator: checks → BLUE → likelihood → ratio → variations → figures → JSON → `result.md` |
| `comb/inputs.py` | reads and validates the two channels' published results |
| `comb/model.py` | **the correlation model** — the only judgement in the folder |
| `comb/blue.py` | BLUE with iteration, weights, χ², source-by-source breakdown |
| `comb/likelihood.py` | nuisance-parameter likelihood, asymmetric errors, −2Δln L scan |
| `comb/ratio.py` | σ(ττ)/σ(μμ) with the correlated part cancelling |
| `comb/plots.py`, `comb/report.py` | figures, and the generator of `result.md` |
| `fit/comb.config`, `fit/run_multifit.py` | the TRExFitter MultiFit, ready but not run |
| `output/combination_result.json` | every number this folder produces |
| `combination.ipynb` | the same result in a notebook, for reading interactively |

## Headline inputs

| channel | extraction | μ_Z | σ(60–120) |
|---|---|---|---|
| Z → μ⁺μ⁻ | counting (per `z-mumu/REVIEW.md` F3) + 0.7 % lineshape | 0.9935 | 1941 ± 35 pb |
| Z → τhτh | `nominal` fake factor | 1.160 ⁺⁰·¹⁹²₋₀.₁₆₃ | 2255 ± 358 pb |

The presentation covering the whole analysis is in [`../presentation/`](../presentation/).
