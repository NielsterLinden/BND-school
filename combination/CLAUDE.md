# CLAUDE.md — combination (TRExFitter MultiFit of Z → ee, Z → μμ, Z → τhτh)

For agents. `combLieke/README.md` is the physics and the result; this file is about changing the
combination without breaking it. The combination is **one TRExFitter v1.8.0 MultiFit** over the three
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
| acceptance uncertainties of μμ and ττ | `channels.json` → `acceptance`, read from the channels' `*.meta.json` |
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
6. **The ranking runs one `trex-fitter mr ... Ranking=<NP>` per NP in parallel** (`--workers`, default 6).
   TRExFitter matches by substring; `run.py impacts` keeps each NP's own row.
7. **σ = μ_Z × 1953.93 pb for every channel.** The constant `xsref_<channel>` NormFactor rescales each
   signal from its own aMC@NLO reference (ee 1954.10, μμ 1953.93, ττ 1944.88 pb) to the POI's. Never
   give two channels one `mu_Z` without it.
8. **The theory band is the generator's own uncertainty** (aMC@NLO weights on σ(60 < m_LHE < 120):
   7-point scale ⊕ NNPDF3.1 Hessian ⊕ α_s). The uncertainty of the NNLO normalisation 6077.22 pb is not
   public and not included.

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
records these per fit (`minuit_status`, `hesse_status`, `minos_status`, `pos_def_forced`).
