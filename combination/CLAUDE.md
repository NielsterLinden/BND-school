# CLAUDE.md — combination (Z → μμ ⊕ Z → τhτh ⊕ Z → ee)

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
folder is that it reproduces from git alone. **The one exception is `tools/extract_zee.py`**,
which does need uproot: z-ee publishes a TRExFitter job and histogram files rather than a results
JSON, so that script converts them once into `inputs/zee_fit_result.json` (committed) and
`run_combination.py` reads only the JSON. Rerun it when z-ee re-publishes; delete it the day z-ee
ships a `z-ee/fit/results/zee_fit_result.json` like the other two channels.

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
| ee everything | `combination/inputs/zee_fit_result.json` ← `tools/extract_zee.py` ← `z-ee/Zee_fit.tar.gz` |
| σ^pred(60–120) per flavour | `z-mumu/output/v2/gensums.json`: 6077.22 pb × Σw(LHE flavour, 60–120)/Σw |
| correlations | `comb/model.CORRELATION`, `comb/model.ACC_CORRELATION` — *asserted*, justified in `docs/02` |
| published CMS/ATLAS comparisons | `comb/plots.REFERENCES` and `comb/report.REFERENCES`, with arXiv links |

Every `ChannelResult` carries a `provenance` list of the files it was built from; it is written
into the JSON.

## Things that will bite you

1. **The z-ee input is known to be biased, and the baseline uses it anyway.** Its `PDF`/`QCDScale`
   templates are un-renormalised LHE envelopes, so `QCDScale` is a ±5.9 % *normalisation* on the
   signal that the fit pulls to −1.88σ; `mu_signal` is measured against a prediction the fit
   rescaled by κ = 0.893. `inputs.load_ee("normfix")` divides that out (−99 pb on the
   combination) and is carried as the `ee_normfix` variation. The baseline follows the channel
   (rule 6 below) and says so loudly in `result.md`, `docs/01` and the deck. **If z-ee re-publishes
   with renormalised templates, κ will go to ~1 and `normfix` becomes a no-op — check that before
   assuming the number moved for a physics reason.**
2. **Do not combine `mu_Z`.** The three channels normalise to different aMC@NLO references for the
   same quantity — 1953.9 pb (μμ), 1944.9 pb (ττ), 1954.1 pb (ee) — because the generator's LHE
   flavour shares are not exactly 1/3. Averaging the μ's would average three different things.
   `comb/inputs.ChannelResult.sigma` multiplies each μ̂ by *its own* reference; work in pb from
   there on.
3. **In-fit uncertainties scale with σ^pred, acceptance with σ^measured.** `group_pb()` and
   `acc_pb()` encode this; it is not a detail. It matches the channels' own bookkeeping.
   z-ee has **no** acceptance term at all — it fits against the 60–120 GeV prediction directly, so
   `acc` is empty by design and `--check` asserts that.
4. **A new `Category` must be given a ρ.** `model.build` raises `KeyError` rather than defaulting
   to zero. If a channel re-runs and adds a systematic group, the combination fails loudly. That
   is deliberate — do not "fix" it with a `.get(cat, 0.0)`. The same applies to the z-ee
   `Category`→convention mapping in `tools/extract_zee.CATEGORY_MAP`, which raises on an
   unmapped category.
5. **`rho_override` does not touch `Data statistics` or `Fit residual`.** The three channels read
   disjoint primary datasets (`SingleMuon`, `Tau`, `Electron`; the ZZ→4ℓ overlap is 0.0006 %), and
   `Fit residual` is a property of one fit's own correlation matrix. Both ρ's are facts, not
   modelling choices.
6. **Each channel decides its own baseline; this folder follows it.** `load_tautau("nominal")`
   reads `nominal_variant` out of the ττ results file rather than hard-coding a variant name, and
   `load_ee` reads `nominal_variant` the same way. On the μμ side the baseline is `"shapefit"`,
   the fit the channel rebuilt after its own review. When a channel re-publishes,
   `run_combination.check()` asserts the published values and fails first; fix it there, not by
   loosening the assertion.
7. **μμ variants rescale, the μμ baseline does not.** `z-mumu/fit/results/stability.json` gives
   μ_Z and the *total* systematic per fit configuration but no per-category breakdown, so
   `load_mumu` gives a variant the nominal category composition scaled to its own total. That is
   fine for `VARIATIONS`; never quote a variant as a result. `load_ee("normfix")` makes the same
   approximation for the uncertainty (it rescales every impact by κ rather than refitting).
8. **`Fit residual` is not a fudge factor.** TRExFitter's grouped impacts are quadrature
   differences and do not add up to the total MINOS error. μμ and ττ *over*-shoot (their
   categories sum to more than their published totals, which this folder keeps, conservatively);
   z-ee *under*-shoots by 0.92 % of μ_Z, which its results file carries as `fit_residual` so that
   the published total is reproduced exactly. Do not fold it into `Data statistics`: the analytic
   Poisson number for z-ee is 0.040 %, 23× smaller, and mislabelling it would put ~7 pb of
   "statistics" in the headline of a 6 M-event measurement.
9. **Figures are light-background on purpose.** The deck is beamer/metropolis and the channel
   plots are white; a dark figure would be a black rectangle on the slide. (The `z-mumu/review`
   deck is the dark house style — different deck, do not mix the two figure sets.)
   Channel colours are `comb/plots.COLOUR`: μμ blue, ττ orange, ee purple, combination dark.

## Changing something

* new variation → add to `run_combination.VARIATIONS`; it is one dict entry and it lands in the
  JSON, `result.md` and `output/plots/variations.pdf` automatically. `only=[...]` drops channels;
  `mumu=`/`tautau=`/`ee=` pick a variant; `build=dict(...)` changes the correlation model;
* changed correlation → `comb/model.py` **and** the justification table in `docs/02`. A ρ with no
  written reason is the one thing this folder should not contain;
* a fourth channel → add a loader in `comb/inputs.py`, an entry in `load_channels.loaders`, a
  `Fit` block in `fit/run_multifit.CHANNELS`, a colour in `comb/plots.COLOUR`, names in
  `comb/report.NAMES`/`SHORT`, and make sure its `Category` strings are in `model.CORRELATION`.
  Nothing else in `comb/` is channel-count-specific.

## Verify after any change

```bash
python run_combination.py --check     # published values + 8 closure tests, all must pass
python fit/run_multifit.py --dry-run  # MultiFit config still generates and validates
```

The closure tests are the safety net: single-channel BLUE reproduces that channel exactly,
identical inputs with ρ = 0 give error/√2 and equal weights, identical inputs with all
systematics correlated leave only the statistical part to average, the profile likelihood
reproduces BLUE to 10⁻³ pb, `only=` reproduces the same sub-combination, and the z-ee grouped
impacts plus `Fit residual` plus the analytic statistics reproduce its published MINOS error.
