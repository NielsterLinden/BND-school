# CLAUDE.md — combination (TRExFitter MultiFit of Z → ee, Z → μμ, Z → ττ)

For agents. `combLieke/README.md` is the physics and the result; this file is about changing the
combination without breaking it.

> **The result is ee ⊕ μμ ⊕ ττ with the four-channel ττ: 1949 ⁺³⁰₋₃₀ pb** (17 Sep 2026, `combLieke/output/`,
> `combLieke/README.md`). ττ = τhτh + μτh + eτh + eμ of z-tautau (inputs b7ec3d9, results e0b0b1d,
> `z-tautau/docs/11-combination-inputs.md`). The user reopened the combination for it after the freeze; the **frozen**
> three-channel result (ee ⊕ μμ ⊕ τhτh, 1945 ⁺³¹₋₃₀ pb, repository `FREEZE.md`) is at the git tag
> `zmumu-freeze-2026-09-17`, its TRExFitter jobs in `combLieke/work_freeze_2026-09-17/` (git-ignored). The ee and μμ
> inputs are still the frozen ones: do not regenerate those channels. `presentation/` may still quote 1945.

The combination is **one TRExFitter v1.8.0 MultiFit** over the three
channel likelihoods. There is no covariance/BLUE combination any more: it and its slide deck were
removed on 16 Sep 2026 (last present in git commit 9501c41).

## Run

```bash
source ../setup.sh               # LCG_110 + trex-fitter v1.8.0 on PATH
cd combLieke
python run.py prepare            # configs of every likelihood in work/ (seconds)
python run.py condor --submit    # all fits as one HTCondor DAG (mf/condor.py); ends with `results --interim`
python run.py results && python run.py plots   # output/result.json and output/plots/ (see the note above first)
python run.py all                # the same chain on one machine: hours with the four-channel tautau likelihood
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
| which files a result was built from (SHA-256, date, last channel commit) | `mf/provenance.py`, written by `prepare` → `work/status.json` → `result.json` `inputs`. Compare it before calling a rerun "the same inputs" |
| acceptance uncertainties of μμ | `channels.json` → `acceptance`, read by key from `zmumu.root.meta.json` (`z-mumu/zmumu/acceptance.py`). ee and ττ have none outside the fit: their theory templates carry A × ε (`acceptance: null`, `acceptance_note`) |
| the ττ signal samples (15 `DYtautau_tDM*`) | read from `ztautau.root.meta.json` → `signal_samples` (`channels.json` → `tautau.signal`), never hard-coded |
| regions removed from a channel (`mumu_CRemu`) | `channels.json` → `drop_regions`, `mf/trexcfg.drop_regions` |
| ττ lepton parameters decorrelated from μμ / ee; `mu_ttbar` vs `XS_TTbar` | `channels.json` → `tautau.rename_systematics` (+ note), `tautau.ttbar_note`; the alternatives are the variations `tautau_leptons_correlated`, `ttbar_xs_constrained`, `ttbar_cr_for_all_channels` |
| single-channel cross-checks (`mumu_CRemu` as a control region) | `channels.json` → `checks` |
| the HTCondor DAG (nodes, resources, job category) | `mf/condor.py`; the job wrapper is `combLieke/condor/run_step.sh` |
| fit options for every channel fit and MultiFit (`FitStrategy: 3`, item 9); ranking refits; per-variation overrides | `channels.json` → `fit`, `ranking_fit`, `variations.<name>.fit` |
| correlations | **by nuisance-parameter name only** (TRExFitter). Renames in `rename_systematics`; ee shape parts get `<NP>_eeShape` |
| alternative likelihoods | `channels.json` → `variations` (`channels`, `split`, `overall`, `channel_overrides`, `fit`) |
| a modelling uncertainty a channel reports but does not fit | `channels.json` → `<channel>.modelling`: the size is read from two of the channel's published fits (`mf/trexcfg.modelling_nps`), never typed in. Empty in the baseline since z-tautau fits `TauIDpT_tautau` itself; used by `tautau_tauh_only`, whose z-tautau config lacks it |
| sharing a free factor across channels, dropping a prior (the tt̄ variations) | `channel_overrides` with `normfactors`, `drop_systematics`, `normfactor_to_overall` |
| z-tautau's own sub-measurements shown in `tautau_channels` | `channels.json` → `tautau.published_submeasurements` |
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
3. **`drop_empty_bins` refuses non-empty bins, and refuses a region that already has `DropBins`.** ττ v4 drops
   `tautau_SR2` bin 1 (0 data, 0 prediction: its free MC-statistics γ sits at 0 and breaks HESSE/MINOS) in its own
   config, and `tautau_SR0` below 110 GeV (fake sideband), so the manifest lists no bin any more. Re-derive after
   every ττ delivery: a bin needs dropping only if **no sample predicts anything** there; 0 data with a non-zero
   prediction (`tautau_SR1` bin 1, `tautau_SR2` bin 13) is an ordinary Poisson term.
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
9. **Only MINUIT strategy 3 gives a usable covariance; status 0 proves nothing.** Every channel normalisation
   parameter (`Lumi`, `MuonReco`, the `Acc_*`, the ττ `TauIDSF_DM*`, `mu_ttbar`) is degenerate with `mu_Z`, and the
   likelihood has ~320 parameters. With strategy 1 or 2 MINUIT ends with status 0 and a correct MINOS interval on
   `mu_Z`, but HESSE runs on a matrix "forced pos-def by adding to diagonal" and the covariance is wrong:
   HESSE error on `mu_Z` 8× below MINOS (0.0018 vs 0.0152), `Lumi` post-fit error 0.15, `mu_Z` most correlated with
   one ττ MC-statistics γ. The uncertainty groups, the post-fit NP errors and the ±1σ points of the ranking all come
   from that matrix. Strategy 3 (Minuit2 "very high": central differences, no forcing) gives HESSE/MINOS = 1.001 and
   `Lumi` 0.97 at the same minimum. Hence `"FitStrategy": "3"` in `channels.json`.
   - The frozen three-channel fit asked for strategy 2 and was still right: there MIGRAD *failed* at 2 and TRExFitter
     escalated to 3 by itself. The four-channel ττ likelihood "converges" at 2, so nothing escalates. Do not go back.
   - `run.py results` records `combined.hesse_over_minos` and `pos_def_forced` per fit and prints a WARNING if
     HESSE/MINOS is off by more than 10 %. Grep the logs for `forced pos-def` after any change of the model.
     (Until 17 Sep `pos_def_forced` looked for another message and was always false, also in the frozen `result.json`:
     the frozen standalone μμ and ee fits did have a forced matrix; only their MINOS intervals were used.)
   - Strategy 3 needs `MaximumNumberFCNcalls` raised (2 × 10⁶): with TRExFitter's default the ee + μμ likelihood
     exhausts the calls in HESSE (status 300, segfault).
   - **`split_all_channels` (per-variation `fit`, `FitStrategy 2`):** strategy 3 ends in an invalid minimum there;
     strategies 1 and 2 converge and agree on the MINOS interval, the only thing a variation provides.
   - **The ranking refits (`ranking_fit`, `FitStrategy -1`)** keep TRExFitter's default (start at 1, retry at 2 and 3):
     only the POI value of each refit is used, and a fixed high strategy leaves no retry. `run.py prepare` writes
     `work/common/multifit_ranking.config` for them.
10. **The correlation matrix in TRExFitter's `Fits/*.txt` has its rows in reverse parameter order** (row k is
    parameter n−1−k, columns in order). `mf/results._corr` accounts for it; until 17 Sep it did not, so the
    `compatibility.correlations` of the frozen `result.json` are other matrix elements (nothing else used them).
11. **ττ v4 is one `Fit` of the MultiFit with 13 regions** (`z-tautau/fit/ztautau.config`): never add the pT-split
    copies (`*_SRlo/hi_dm*`, same events), never put `mu_Z` on `DYtautau_out`, never add `Acc_*` for ττ (its theory
    templates carry A × ε), and never let `mumu_CRemu` back in (same data as `emu_SR`/`emu_CRtt`). A per-final-state
    POI is not defined in this model: with free `TauIDSF_DM*` a single ℓτh channel cannot separate `mu_Z` from the
    scale factors, so the three-POI fit stays ee / μμ / ττ.
12. **`TauIDpT_tautau` is in the baseline, and since 17 Sep it is z-tautau's own parameter.** Their p_T-split cross-check
    moves μ_Z by 6.3 %; their measurement config carries that as one OVERALL parameter on the ττ signal (size from
    their fits `ztautau_flatsf` and `ztautau_ptsplit`), so `tautau.modelling` is empty: **never add it a second time**
    (`tests/test_trexcfg.py` checks it is there exactly once), and never together with the p_T-split model.
    `tautau_without_tauidpt` drops it (`drop_systematics`), `tautau_tauh_only` adds it to a z-tautau config that lacks it. tt̄: the baseline keeps the free ττ `mu_ttbar` and the
    `XS_TTbar` prior of ee/μμ as two parameters; the two alternatives (±5 pb) and the reasons are in
    `combLieke/README.md`, "tt̄: the eμ control region, and who uses it". Do not share a free factor with μμ or ee
    without a prior: their standalone fits and `without_tautau` would have nothing to determine it.
13. **HTCondor** (`mf/condor.py`, nikhef-condor skill): one DAG, batch name `zcomb`; watch it with
    `/user/sjankovy/.claude/skills/nikhef-condor/condor_watch.sh --prefix zcomb` (it stops on held jobs). `condor_rm`
    of the DAG still runs its FINAL node. The list of parameters to rank is written by the `common` job; the ranking
    node waits for it on the schedd (`condor/wait_for_file.sh`). `run.py condor` removes the `NPRanking_*` files of an
    earlier submission; everything else in `work/` is overwritten job by job, so `rm -rf work` first if the model changed.
    One DAG node per likelihood: DAGMan removes every other job of a cluster when one fails. The FINAL node writes
    `combLieke/interim/result.json` (git-ignored), never `output/`; `python run.py results && python run.py plots` promotes it.

## Figures

`mf/plots.py`, from `output/result.json` only. Light background and the channel colours μμ `#2B6CB0`,
ττ `#EB811B`, ee `#8E44AD`, combination `#23373B`, prediction `#14B03D`. Header "CMS *Open Data*" and
"16.4 fb⁻¹ (13 TeV)"; keep text inside the axes to a minimum (values in a right-hand column, no titles).

## Verify after any change

```bash
python tests/test_trexcfg.py
python run.py prepare && python run.py condor --submit      # wait for the DAG (batch zcomb), then:
grep -E "status|probability|forced pos-def" work/common/logs/multifit_mwf.log && python run.py results
```

Every fit in `run.py` must end with MIGRAD status 0, HESSE status 0 and MINOS status 0. `result.json`
records these per fit (`minuit_status`, `hesse_status`, `minos_status`, `pos_def_forced`). Status 0 is not
enough (item 9): `combined.hesse_over_minos` must be within 0.9–1.1, and no constrained parameter may have a
post-fit error above 1.
