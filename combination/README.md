# Combination — Z → μμ ⊕ Z → τhτh ⊕ Z → ee

One cross section from the three BND-school channels, on CMS Open Data 2016
(Run2016G+H, 16 393.381 pb⁻¹, √s = 13 TeV).

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1980 ± 32 pb**
> = 1980 ± 0.5 (stat) ± 22 (syst) ± 7 (acc) ± 21 (lumi) pb, per lepton flavour,
> assuming lepton universality.
>
> **χ²/ndf = 9.67/2 (p = 0.008)** — the three channels are *not* mutually consistent, and that
> is part of the result. σ(ee) = 2054 ± 40 pb sits 3.0σ above σ(μμ) = 1931 ± 35 pb because the
> z-ee fit lets a ±5.9 % theory *normalisation* of its own signal template float against its
> signal strength. Undoing only that moves the combination by −99 pb.
> See [`docs/01-inputs.md`](docs/01-inputs.md) — it is a fixable one-line problem in the z-ee
> notebook, not a physics result.

Full numbers and caveats: **[`result.md`](result.md)**. The physics write-up is in
[`docs/`](docs/); start with [`docs/00-overview.md`](docs/00-overview.md). The orthogonality check,
the systematics CMS and ATLAS carry that we do not, and each channel against the published
measurement of the same decay are in [`docs/05-vs-published.md`](docs/05-vs-published.md). The list
of changes z-ee needs is [`../ASK_Z_EE.md`](../ASK_Z_EE.md).

## Run it

```bash
source ../setup.sh
python tools/extract_zee.py          # only when z-ee re-publishes Zee_fit.tar.gz (needs uproot)
python run_combination.py            # checks, result.md, output/*.json, output/plots/*
python run_combination.py --check    # input validation and closure tests only (2 s)
python run_combination.py --no-plots # numbers only
python fit/run_multifit.py --dry-run # regenerate and syntax-check the TRExFitter MultiFit config
```

There is nothing to reprocess: every input is a number the channel groups committed. The whole
thing runs in a few seconds.

## What this is — and is not

It **is** a covariance (BLUE) combination of three TRExFitter v1.8.0 profile-likelihood fits, with
the shared systematics correlated through the `Category` strings that `fitting/CONVENTIONS.md`
forces every channel to use, plus a profile-likelihood cross-check and two lepton-universality
ratios.

It is **not** a TRExFitter MultiFit. The channel workspaces are not all committed and were not
rebuilt, so the shared nuisance parameters are correlated by hand rather than profiled jointly.
The MultiFit configuration is generated and syntax-checked in
[`fit/comb.config`](fit/comb.config); [`fit/run_multifit.py`](fit/run_multifit.py) carries the
exact recipe to produce the workspaces and run it, and
[`docs/03-method.md`](docs/03-method.md) says what it would add — which, now that z-ee carries
40 % of the weight, is more than it used to be.

## Layout

| path | what |
|---|---|
| `run_combination.py` | orchestrator: checks → BLUE → likelihood → ratios → variations → figures → JSON → `result.md` |
| `tools/extract_zee.py` | adapter: `z-ee/Zee_fit.tar.gz` → `inputs/zee_fit_result.json` (the only thing here that needs uproot) |
| `inputs/zee_fit_result.json` | the z-ee result in the shape the other two channels publish |
| `comb/inputs.py` | reads and validates the three channels' published results |
| `comb/model.py` | **the correlation model** — the only judgement in the folder |
| `comb/blue.py` | BLUE with iteration, weights, χ², source-by-source breakdown |
| `comb/likelihood.py` | nuisance-parameter likelihood, asymmetric errors, −2Δln L scan |
| `comb/ratio.py` | σ_X/σ_μμ with the correlated part cancelling |
| `comb/plots.py`, `comb/report.py` | figures, and the generator of `result.md` |
| `fit/comb.config`, `fit/run_multifit.py` | the TRExFitter MultiFit, ready but not run |
| `output/combination_result.json` | every number this folder produces |
| `combination.ipynb` | the same result in a notebook, for reading interactively |

## Headline inputs

| channel | extraction | μ_Z | σ(60–120) | weight |
|---|---|---|---|---|
| Z → μ⁺μ⁻ | v2 shape fit, 12 × 5 GeV, after the channel's review fixes | 0.9881 ± 0.0159 | 1931 ± 35 pb | +0.596 |
| Z → τhτh | **v3**, DeepTau Tight both legs, MC-subtracted fake factor | 1.071 ⁺⁰·¹¹⁴₋₀.₁₀₀ | 2082 ± 253 pb | −0.0003 |
| Z → e⁺e⁻ | the published TRExFitter job, `Zee_fit.tar.gz` of 16 Sep 2026 | 1.0509 ± 0.0202 | 2054 ± 40 pb | +0.404 |

All three are the channels' own current baselines;
`python run_combination.py --check` asserts every one of those numbers against the files the
channel groups committed and fails loudly if a channel moves again.

The presentation covering the whole analysis is in [`presentation/`](presentation/).
