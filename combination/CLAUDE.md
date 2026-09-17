# CLAUDE.md — combination (TRExFitter MultiFit of Z → ee, Z → μμ, Z → τhτh)

For agents. `combLieke/README.md` is the physics and the result; this file is about changing the
combination without breaking it.

> **Frozen on 17 Sep 2026** (git tag `zmumu-freeze-2026-09-17`, repository `FREEZE.md`). The combination was rerun
> on the frozen channel inputs: Z → μμ with its measured muon reconstruction scale factor and its 60–120 GeV
> acceptance block (eight μμ acceptance parameters instead of four), and MINUIT strategy 2 for every fit.
> `combLieke/output/` is the final result: read and present it, do not regenerate it. Only
> `presentation/` still shows pre-freeze numbers; updating it is the animation agents' job.

The combination is **one TRExFitter v1.8.0 MultiFit** over the three
channel likelihoods. There is no covariance/BLUE combination any more: it and its slide deck were
removed on 16 Sep 2026 (last present in git commit 9501c41).

## Run

```bash
source ../setup.sh               # LCG_110 + trex-fitter v1.8.0 on PATH
cd combLieke
python run.py all                # prepare -> workspaces -> fits -> variations -> impacts -> results -> plots (~10 min)
python run.py results && python run.py plots   # re-parse and redraw only
python tests/test_trexcfg.py     # adapter checks, 2 s
python checks/systematics.py     # input-level systematic sizes -> checks/systematics.json (after prepare)
python checks/orthogonality.py   # event-level overlap on data -> checks/orthogonality.json (~3 min)
```

Everything TRExFitter writes goes to `combLieke/work/` (git-ignored, rebuilt from scratch by `run.py`).
What is kept is `combLieke/output/` (`result.json`, `plots/`) and the `checks/*.json`.

## Where things are decided

| decision | where |
|---|---|
| which channel configs / inputs, references, signal samples | `combLieke/config/channels.json` |
| every change made to a channel config | `mf/trexcfg.adapt_channel` (docstring lists all of them) |
| the ee inputs (histograms from `z-ee/Zee_fit.tar.gz`) | `mf/ee_input.py` |
| acceptance uncertainties of μμ and ττ | `channels.json` → `acceptance`, read by key from the channels' `*.meta.json` (μμ: `z-mumu/zmumu/acceptance.py`) |
| fit options for every channel fit and MultiFit (`FitStrategy: 2`); ranking refits; per-variation overrides | `channels.json` → `fit`, `ranking_fit`, `variations.<name>.fit` |
| correlations | **by nuisance-parameter name only** (TRExFitter). Renames in `rename_systematics`; ee shape parts get `<NP>_eeShape` |
| alternative likelihoods | `channels.json` → `variations` |
| published ATLAS/CMS numbers | `combLieke/config/references.json` (with paper and table) |
| theory prediction and its uncertainty | `mf/prediction.py` |

## Things that will bite you

1. **Do not correlate the full ee templates with the other channels.** Without `split_shape_norm` for ee,
   the fully correlated three-channel fit has no positive-definite minimum (TRExFitter segfaults after
   strategy 2). The reason is physics, not numerics: the ee peak shape pins its ±5.9 % `ElectronID`
   normalisation at 0.08σ and drags `Pileup`/`L1Prefiring`/`QCDScale` into μμ. If z-ee re-delivers
   with the ECAL gap vetoed (`ElectronID` ≈ ±1 %), re-test whether the split is still needed. The
   `ee_split_shared_only` and `ee_electron_id_1p2` variations exist for exactly that.
2. **`NLLOffset: bin` does not work** here: ROOT 6.40 returns NaN for bins with zero observed events,
   which ττ has. The default offset converges once the empty ττ bin is dropped.
3. **`drop_empty_bins` refuses non-empty bins.** Only `tautau_SR2` bin 1 (0 data, 0 prediction) is
   dropped, because its free MC-statistics γ sits at 0 and breaks HESSE/MINOS.
4. **TRExFitter reserves `shape_`, `alpha`, `gamma` in NP names**, so the ee shape parts are `<NP>_eeShape`.
5. **The refit-based grouped impacts (`trex-fitter mi`) return NaN** in this likelihood (HESSE fails
   with groups fixed). The uncertainty groups are TRExFitter's covariance decomposition from `mwf`;
   data statistics come from the separate stat-only fit (`StatOnly=TRUE:Suffix=_statOnly`, using
   `multifit_statonly.config`, which has no likelihood scan so the full scan is not overwritten).
6. **The ranking runs one `trex-fitter mr multifit_ranking.config Ranking=<NP>` per NP in parallel**
   (`--workers`, default 6). TRExFitter matches by substring; `run.py impacts` keeps each NP's own row. A refit
   that fails is recorded as a POI of 0, i.e. a shift of −μ̂ ≈ −0.996: check the ranking for such rows.
7. **σ = μ_Z × 1953.93 pb for every channel.** The constant `xsref_<channel>` NormFactor rescales each
   signal from its own aMC@NLO reference (ee 1954.10, μμ 1953.93, ττ 1944.88 pb) to the POI's. Never
   give two channels one `mu_Z` without it.
8. **The theory band is the generator's own uncertainty** (aMC@NLO weights on σ(60 < m_LHE < 120):
   7-point scale ⊕ NNPDF3.1 Hessian ⊕ α_s). The uncertainty of the NNLO normalisation 6077.22 pb is not
   public and not included.
9. **Strategy 1 HESSE is not reliable in this likelihood.** Every channel normalisation parameter (`Lumi`,
   `MuonReco`, the `Acc_*`) is degenerate with `mu_Z`. With the frozen μμ acceptance block (four more of them),
   MINUIT strategy 1 still returns status 0 and a correct MINOS interval, but its covariance is wrong: HESSE error
   on `mu_Z` 1.65× below MINOS, `Lumi` post-fit error 0.71, `Acc_PDF` 1.19 (above its prior), and the
   uncertainty groups with it (luminosity 1 pb instead of 22). Hence `"FitStrategy": "2"` in `channels.json`
   (8 s per fit). `run.py results` records `combined.hesse_over_minos` and prints a WARNING if it is off by more
   than 10 %. Two exceptions, both in `channels.json`:
   - **the ranking refits (`ranking_fit`, `FitStrategy -1`):** TRExFitter's default starts at strategy 1 and
     retries at 2 and 3. An explicit 2 leaves one retry, and 27 of 75 refits failed (status 3).
     `run.py prepare` writes `work/common/multifit_ranking.config` for them.
   - **`without_ee` (per-variation `fit`, `FitStrategy 1`):** strategy 2 ends with MINUIT and MINOS status 1
     there, strategy 1 with 0. Only its MINOS interval on μ_Z is used.

## Figures

`mf/plots.py`, from `output/result.json` only. Light background and the channel colours μμ `#2B6CB0`,
ττ `#EB811B`, ee `#8E44AD`, combination `#23373B`, prediction `#14B03D`. Header "CMS *Open Data*" and
"16.4 fb⁻¹ (13 TeV)"; keep text inside the axes to a minimum (values in a right-hand column, no titles).

## Verify after any change

```bash
python tests/test_trexcfg.py
python run.py all && grep -E "status|probability" work/common/logs/multifit_mwf.log
```

Every fit in `run.py` must end with MIGRAD status 0, HESSE status 0 and MINOS status 0. `result.json`
records these per fit (`minuit_status`, `hesse_status`, `minos_status`, `pos_def_forced`). Status 0 is not
enough (item 9): `combined.hesse_over_minos` must be within 0.9–1.1, and no constrained parameter may have a
post-fit error above 1.
