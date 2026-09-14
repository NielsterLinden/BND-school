# CLAUDE.md — z-tautau (Z → τhτh)

Guidance for agents working in this folder. Read `README.md` for what the analysis is and
`docs/00-overview.md` for the physics; this file is about working on the code without breaking it.

## Environment

* `source ../setup.sh` (LCG_110: python 3.13, uproot 5.7, awkward 2.9, hist, mplhep 1.1, ROOT 6.40,
  correctionlib). Nothing to pip-install; do not create a venv.
* TRExFitter v1.8.0: `step5_fit.py` finds `trex-fitter` on PATH, else `$TREXFITTER_HOME`, else
  `config.TREX_FALLBACK_HOME` (the build in the main checkout `/project/atlas/Users/nterlind/BND-school`).
  Every channel must use this one binary. Never build another version.
* Bulk data lives in `$BND_TAUTAU_CACHE` (default `/data/atlas/users/nterlind/BND-school-cache/ztautau`):
  `skims_v1/`, `ntuples_v1/`, `logs/`. `/data/atlas` is ~98 % full: do not duplicate skims.

## Running

```bash
python run_all.py --from 3                  # normal iteration loop (~10 min): FF -> histograms -> fits -> report
python scripts/step4_histograms.py --ff-variant nominal --no-plots   # fast histogram rebuild
python scripts/step5_fit.py --ff-variant nominal --skip-ranking      # fast fit (~40 s)
```

* Long jobs (skims): launch detached (`nohup setsid bash -c "..." > $BND_TAUTAU_CACHE/logs/x.log &`) and poll
  the log. Steps 1 and 2 are resumable/cheap; step 1 skips finished files.
* Every script is standalone and its docstring is its usage. Physics constants are only in
  `ztautau/config.py`; sample definitions only in `ztautau/samples.py`.
* Two fake-factor variants run everywhere: `nominal` (no MC subtraction, as the prompt asked) and `mcsub`.
  Files for the non-nominal variant carry the suffix `_mcsub`. Flip `config.FF_SUBTRACT_MC` to swap roles.

## Conventions you must keep

* Histogram and nuisance-parameter names follow `../fitting/CONVENTIONS.md` (region `tautau_SR`, POI `mu_Z`,
  samples `DYtautau`, `Fakes`, …; channel-specific NPs end in `_tautau`). The combination
  (`../fitting/combination_skeleton.config`) expects `fit/ztautau.config` and `fit/results/ztautau/`.
* Luminosity 16393.381 pb⁻¹ (normtag), DY cross section 6077.22 pb: identical in all channels.
* Simulation keeps only events whose **leading τ is not a jet** (`analysis.regions(..., is_mc=True)`);
  removing that cut double-counts the fake-factor estimate.
* Commit plots (PNG), `output/RESULTS.md`, `output/results.json`, `fit/*.config`, `external/*.json`,
  `filelists/`. Never commit ROOT files, `output/data/` or PDFs of plots.

## Pitfalls that already cost time here (do not reintroduce)

1. **Theory envelopes/Hessians must be combined on histograms, never event by event**
   (`analysis.theory_weights` + `combine_theory`). The event-wise max over scale variations gave ±15 %
   instead of ±4 %.
2. **W+jets (aMC@NLO inclusive) has a handful of events with weights up to ~200.** Bin-by-bin it produced
   "±100 events" of MC stat in empty bins and inflated all γ parameters. Its SR shape comes from the τ2
   ≥ VVVLoose selection (`analysis.SMOOTHED_SAMPLES`); its normalisation stat. is `MCStatNorm_WJets_tautau`.
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
9. Fake factors: without the **N_jets** binning the same-sign closure is off by 10–25 % per jet bin, and
   without the **era** binning by ±5 % per era (different HLT τ isolation in G and H). Check
   `step3_closure_SS_*` after any change to the τ selection.

## Where numbers come from

| number | source |
|---|---|
| cutflow | `$BND_TAUTAU_CACHE/ntuples_v1/<sample>.meta.json` → `output/results.json["cutflow"]` |
| fake factors, C_OS/SS, contamination | `output/data/fakefactors.json` (step 3) |
| prefit yields, prefit variation sizes | `output/data/yields[_mcsub].json` (step 4) |
| μ, impacts, pulls, ranking | `fit/results/ztautau[_mcsub]_fit_result.json` (step 5) |
| everything collected | `output/results.json`, `output/RESULTS.md` (step 6) |

## Open issues (good next tasks)

* Make the genuine-τ-subtracted fake factor nominal (see `docs/05-fake-factors.md`).
* Closure corrections in η(τ1) and pT(τ2) (±15 % / −7 % non-closure, currently only through `FakeClosure_tautau`).
* The τh ID SF uncertainty (12.5 %) dominates: in the combined fit with ee/μμ the ττ channel should be
  used to constrain the τh ID nuisance parameters rather than μ_Z.
* Signal MC statistics (γ: 4.8 %): add the jet-binned DY aMC@NLO samples (records 35577/35595/35613).
* W+jets: replace the inclusive sample by the HT-binned madgraph samples with LHE_HT stitching.
