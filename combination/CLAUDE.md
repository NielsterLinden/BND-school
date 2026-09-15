# CLAUDE.md — combination (Z → μμ ⊕ Z → τhτh)

For agents working in this folder. Read [`README.md`](README.md) for what it is,
[`docs/00-overview.md`](docs/00-overview.md) for the physics, and
[`result.md`](result.md) for the numbers. This file is about not breaking it.

## Environment

```bash
source ../setup.sh          # LCG_110: python 3.13, numpy, scipy, iminuit, matplotlib
python run_combination.py   # everything, a few seconds
```

No ROOT, no uproot, no TRExFitter binary needed: the combination reads JSON and text files that
the channel groups committed. Do not add a dependency that would change that — the point of this
folder is that it reproduces from git alone.

## The one rule

**`result.md` is generated.** `comb/report.py` renders it from
`output/combination_result.json`, which `run_combination.py` writes in the same run. Never edit
`result.md` by hand; edit the template in `comb/report.py`, or the code that produces the number,
and rerun. The same goes for the figures.

## Where each number comes from

| number | source |
|---|---|
| μμ μ_Z, grouped impacts, ranking, A, C, counting yields | `z-mumu/fit/results/zmumu_fit_result.json` |
| μμ alternative fit configurations (binnings, counting) | `z-mumu/fit/results/stability.json` |
| ττ everything | `z-tautau/output/results.json` (the variant its own `nominal_variant` names) |
| ττ Tight-WP cross-check | `z-tautau/variants/tight/output/results.json` |
| correlations | `comb/model.CORRELATION`, `comb/model.ACC_CORRELATION` — *asserted*, justified in `docs/02` |
| published CMS/ATLAS comparisons | `comb/plots.REFERENCES` and `comb/report.REFERENCES`, with arXiv links |

Every `ChannelResult` carries a `provenance` list of the files it was built from; it is written
into the JSON.

## Things that will bite you

1. **Do not combine `mu_Z`.** The two channels normalise to different aMC@NLO references for the
   same quantity — 1953.9 pb (μμ) against 1944.9 pb (ττ), 0.47 % apart, see `docs/01`. Averaging
   the μ's would average two different things. `comb/inputs.ChannelResult.sigma` multiplies each
   μ̂ by *its own* reference; work in pb from there on.
2. **In-fit uncertainties scale with σ^pred, acceptance with σ^measured.** `group_pb()` and
   `acc_pb()` encode this; it is not a detail (they differ by μ = 1.16 for ττ). It matches the
   channels' own bookkeeping.
3. **A new `Category` must be given a ρ.** `model.build` raises `KeyError` rather than defaulting
   to zero. If a channel re-runs and adds a systematic group, the combination fails loudly. That
   is deliberate — do not "fix" it with a `.get(cat, 0.0)`.
4. **`rho_override` does not touch `Data statistics`.** The two channels read disjoint primary
   datasets (`SingleMuon`, `Tau`); its ρ is a fact, not a modelling choice.
5. **The negative ττ weight is correct.** BLUE gives a negative weight whenever ρ exceeds the
   ratio of the two uncertainties (0.114 > 0.110 here — it is marginal now). Do not clip it —
   explain it (`docs/03`).
6. **Each channel decides its own baseline; this folder follows it.** `load_tautau("nominal")`
   reads `nominal_variant` out of the ττ results file rather than hard-coding a variant name —
   that field flipped from `nominal` (unsubtracted) to `mcsub` between v2.0 and v2.1 and the
   combination followed without a code change. On the μμ side the baseline is `"shapefit"`, the
   fit the channel rebuilt after its own review; the counting extraction and the `Lineshape model`
   category that `z-mumu/REVIEW.md` F3 asked for as a stop-gap are **retired** — the two-sided
   `SigModel` carries that uncertainty inside the fit now, so re-adding the term would double
   count. When a channel re-publishes, `run_combination.check()` asserts the published values and
   fails first; fix it there, not by loosening the assertion.
7. **μμ variants rescale, the μμ baseline does not.** `z-mumu/fit/results/stability.json` gives
   μ_Z and the *total* systematic per fit configuration but no per-category breakdown, so
   `load_mumu` gives a variant the nominal category composition scaled to its own total. That is
   fine for `VARIATIONS`; never quote a variant as a result.
8. **Figures are light-background on purpose.** The deck is beamer/metropolis and the channel
   plots are white; a dark figure would be a black rectangle on the slide. (The `z-mumu/review`
   deck is the dark house style — different deck, do not mix the two figure sets.)

## Changing something

* new variation → add to `run_combination.VARIATIONS`; it is one dict entry and it lands in the
  JSON, `result.md` and `output/plots/variations.pdf` automatically;
* changed correlation → `comb/model.py` **and** the justification table in `docs/02`. A ρ with no
  written reason is the one thing this folder should not contain;
* new channel (z-ee) → add a loader in `comb/inputs.py`, a `Fit` block in
  `fit/run_multifit.CHANNELS`, and make sure its `Category` strings are in `model.CORRELATION`.

## Verify after any change

```bash
python run_combination.py --check     # published values + 6 closure tests, all must pass
python fit/run_multifit.py --dry-run  # MultiFit config still generates and validates
```

The closure tests are the safety net: single-channel BLUE reproduces that channel exactly,
identical inputs with ρ = 0 give error/√2 and equal weights, identical inputs with all
systematics correlated leave only the statistical part to average, and the profile likelihood
reproduces BLUE to 10⁻³ pb.
