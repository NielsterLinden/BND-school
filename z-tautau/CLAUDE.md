# CLAUDE.md — z-tautau (Z → τhτh)

Guidance for agents working in this folder. Read `README.md` for what the analysis is, `docs/00-overview.md`
for the physics, and `REVIEW.md` for the review that led to the current (v2) set-up; this file is about
working on the code without breaking it.

## Environment

* `source ../setup.sh` (LCG_110: python 3.13, uproot 5.7, awkward 2.9, hist, mplhep 1.1, xgboost 2.1, ROOT 6.40).
  Nothing to pip-install; do not create a venv.
* TRExFitter v1.8.0: `step5_fit.py` finds `trex-fitter` on PATH, else `$TREXFITTER_HOME`, else
  `config.TREX_FALLBACK_HOME`. Every channel must use this one binary. Never build another version.
* Bulk data lives in `$BND_TAUTAU_CACHE` (default `/data/atlas/users/sjankovy/BND-school-cache/ztautau`):
  `skims_v1/` (data and small samples are symlinks to the original skims in
  `/data/atlas/users/nterlind/BND-school-cache/ztautau`, the DY samples were re-skimmed with `LHE_NpNLO`),
  `ntuples_v1/` (the laptop bundle, ~330 MB), `bdt/` (the five fold models), `logs/`.
  `/data/atlas` is ~98 % full and `/project` is over quota: do not duplicate skims, keep plots small.

## Running

```bash
python run_all.py --from 3                  # normal iteration loop (~25 min): FF -> BDT -> histograms -> fits -> report
python scripts/step3b_bdt.py --no-train     # only the BDT validation plots with the existing models
python scripts/step4_histograms.py --no-plots           # fast histogram rebuild (nominal variant)
python scripts/step5_fit.py --skip-ranking              # fast fit (~2 min)
BND_TAUTAU_WP=Medium python run_all.py --from 3        # the same chain with another working point -> variants/medium/
```

* Long jobs (skims): launch detached (`nohup setsid python scripts/step1_skim.py ... > $BND_TAUTAU_CACHE/logs/x.log &`)
  and poll the log. Steps 1 and 2 are resumable/cheap; step 1 skips finished files.
* Every script is standalone and its docstring is its usage. Physics constants are only in
  `ztautau/config.py`; sample definitions only in `ztautau/samples.py`.
* Only the MC-subtracted fake factor (`mcsub`) is produced; `scripts/step3_fakefactors.py --with-nosub` and
  `--ff-variant nosub` in steps 4/5 give the unsubtracted cross-check on request (files with the suffix `_nosub`).
* The fake-dominated category `tautau_SR0` enters the fit only above `config.SIDEBAND_REGION_MTT_MIN`
  (110 GeV, `DropBins` in `step5_fit.py`): it is the fake sideband, not a signal region.
* `BND_TAUTAU_WP` (default `config.NOMINAL_WP` = `Tight` since v3) selects the VSjet working point of both legs and,
  for anything else, redirects `output/`, `fit/` and the BDT models to `variants/<wp>/` (`config.py`). The TauPOG ID
  and trigger scale factors of every working point are in `external/*.json` (step 0).
* The BDT (step 3b) must be retrained after any change to the fake factors, the selection or
  `config.BDT_*`; `bdt.load_models` refuses models trained with other features / folds.

## Conventions you must keep

* Histogram and nuisance-parameter names follow `../fitting/CONVENTIONS.md`: regions `tautau_SR0/1/2`
  (BDT categories), POI `mu_Z` on the **fiducial** signal `DYtautau` only, `DYtautau_nonfid` is a
  background, channel-specific NPs end in `_tautau` (`SigModel_tautau`, `FakeOSSS_tautau_c<k>`,
  `FakeClosure_tautau_c<k>_lo|hi`). The combination (`../fitting/combination_skeleton.config`) expects
  `fit/ztautau.config` and `fit/results/ztautau/`.
* Luminosity 16393.381 pb⁻¹ (normtag), DY cross section 6077.22 pb: identical in all channels.
* Simulation keeps only events whose **leading τ is not a jet** (`analysis.regions(..., is_mc=True)`);
  removing that cut double-counts the fake-factor estimate.
* Nothing that enters the definition of the fake-factor regions (τ₁ isolation, raw DeepTau scores, charge)
  may become a BDT input (`config.BDT_FEATURES`). Every event is scored by the fold model that never saw it
  (`bdt.folds`: event number mod 5). Check `step3b_bdt_closure_SS.png` after any change.
* Commit plots (PNG), `output/RESULTS.md`, `output/results.json`, `fit/*.config`, `fit/bdt_info.json`,
  `external/*.json`, `filelists/`, and — the one exception to the no-ROOT rule, because the combination reads
  them from its own checkout — the fit inputs `fit/fitinputs/ztautau*.root` (+ `.meta.json`), the workspace
  `fit/results/ztautau/RooStats/ztautau_combined_ztautau_model.root` and `fit/results/ztautau_fit_result.json`
  (`git add -f`, they are git-ignored by pattern). Never commit other ROOT files, `output/data/`, the BDT
  models or PDFs of plots.

## Pitfalls that already cost time here (do not reintroduce)

1. **Theory envelopes/Hessians must be combined on histograms, never event by event**
   (`analysis.theory_weights` + `combine_theory`). The event-wise max over scale variations gave ±15 %
   instead of ±4 %.
2. **W+jets (aMC@NLO inclusive) has a handful of events with weights up to ~200.** Bin-by-bin it produced
   "±100 events" of MC stat in empty bins and inflated all γ parameters, and one event with weight −63
   faked a 75 % non-closure of the FF in the top BDT bin. Its SR shape comes from the τ2 ≥ VVVLoose
   selection (`analysis.SMOOTHED_SAMPLES`), its normalisation stat. is `MCStatNorm_WJets_tautau`, and it is
   subtracted from the FF regions with **uniform weights** (`analysis.subtraction_weights`); dropping it instead
   moves C_OS/SS by 1.5 %.
3. **multiprocessing with fork deadlocks** after the parent has used uproot: `step2_ntuples.py` uses the
   `spawn` context. Keep it.
4. **uproot ≥ 5.6 writes `f["name"] = {dict}` as an RNTuple**, which `.arrays(library="np")` cannot read
   back. Write trees with `f.mktree(name, dict_of_arrays)` (this also fills the first chunk).
5. **mplhep 1.1 label helpers (`hep.cms.label/text`) crash with the LCG_110 matplotlib** at draw time; use
   `plotting.label`.
6. **dCache NFS**: stats and reads can hang or return EIO. Use `os.listdir` rather than per-file
   `exists()`, run every file in its own subprocess with a timeout (`batch.run_files`), fall back to EOS.
7. A skim file whose events all fail the preselection has **no `Events` tree** (only `GenSums`); readers
   must check.
8. `pkill -f <pattern>` also kills the shell running it when the pattern is in its own command line.
9. Fake factors: without the **N_jets** binning the same-sign closure is off by 10–25 % per jet bin, without
   the **era** binning by ±5 % per era, and without the **|η(τ₁)| closure correction** by ±15 % in η
   (`fakes.closure_corrections`). Check `step3_closure_SS_*` after any change to the τ selection.
10. **xgboost 2.1 + scikit-learn 1.9**: `XGBClassifier.save_model` raises `_estimator_type undefined`;
    save and load the `Booster` (`bdt.py` does).
11. **The signal template is 38 % non-fiducial** (30 % has m_LHE > 120 GeV and sits in the high-m_tt bins).
    Never merge `DYtautau_nonfid` back into the signal: μ_Z would then scale the high-mass continuum and the
    fake sideband would contain μ-dependent signal (REVIEW.md 3.1).
12. **TRExFitter must really be v1.8.0**: `git submodule status` must show commit `f21d1b36` without a `+`.
    This checkout once had the submodule at a v1.10 commit and a binary built against the StatAnalysis
    ROOT; the symptoms were `Cannot find setting 'UseGammaPulls'` (schema mismatch) and a crash at start-up
    (`TExMap::Add key not unique`, `munmap_chunk`) under LCG. Fix: `git submodule update --init`, then
    `bash ../fitting/build_trexfitter.sh`.
13. **Jet-binned DY stitching** (`analysis.dy_norm`) needs the `sumw_npnlo` GenSums of *every* stitched
    sample, including the inclusive one; a sample skimmed without `LHE_NpNLO` is silently dropped from the
    stitching with a warning. The bin variable is `LHE_NpNLO` (partons of the NLO matrix element), **not**
    `LHE_Njets` (which also counts the real-emission parton: the 0J/1J/2J samples overlap in it and the
    stitched yield came out 7 % high).

## Where numbers come from

| number | source |
|---|---|
| cutflow | `$BND_TAUTAU_CACHE/ntuples_v1/<sample>.meta.json` → `output/results.json["cutflow"]` |
| fake factors, closure corrections, C_OS/SS, contamination | `output/data/fakefactors.json` (step 3) |
| BDT training summary, category yields, closure in the score | `output/data/bdt.json`, `fit/bdt_info.json` (step 3b) |
| prefit yields per category, closure NPs, C_OS/SS per category | `output/data/yields.json` (step 4) |
| μ, impacts, pulls, ranking | `fit/results/ztautau_fit_result.json` (step 5) |
| everything collected | `output/results.json`, `output/RESULTS.md` (step 6) |

## Open issues (good next tasks)

* The τh ID SF uncertainty still dominates: in the combined fit with ee/μμ the ττ channel should be
  used to constrain the τh ID nuisance parameters rather than μ_Z; decay-mode categories would help.
* The DM-binned and pT-binned TauPOG ID prescriptions differ by 7 % per leg (REVIEW.md 3.3); the pT-binned
  result should be quoted as a cross-check.
* Trigger efficiency in situ (μτh tag-and-probe from SingleMuon) instead of the POG turn-on SFs.
* W+jets: replace the inclusive sample by the HT-binned madgraph samples with LHE_HT stitching.
* High-mass DY sample (record 35629) stitched in m_LHE for the non-fiducial template statistics.
* EWK Z→ττ (record 35999) is absent (≈ 0.3 %); decide whether it is signal or background.
