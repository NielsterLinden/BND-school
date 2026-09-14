# CLAUDE.md — z-mumu (Z → μ⁺μ⁻ cross section)

For agents continuing this channel. Physics write-up: `docs/` (start with `docs/09-v2-overview.md`);
numbers: `output/v2/RESULTS_v2.md` and `output/v2/results_v2.json`; combination inputs: `handoff.md`.

## Two analyses live here

| | v1 (steps 1–6, `run_all.py`) | **v2 (`run_v2.py`, the result)** |
|---|---|---|
| inputs | ≥2-muon skim `/dcache/atlas/kdevries/BND2026/DoubleMuonSkimmed` | own skims of the unskimmed SingleMuon parent + MC (`$BND_SKIM_DIR`) |
| method | data-only counting, T&P by counting, SS + eμ backgrounds | MC-based: T&P **fits** → scale factors, pileup + L1-prefiring weights, fake factor, TRExFitter binned-likelihood fit of m(μμ) |
| muon ID | medium | **tight** |
| result | σ_fid = 773.2 pb (frozen, `output/RESULTS.md`) | `output/v2/RESULTS_v2.md` |

v1 is kept frozen and reproducible (its `.venv` is python 3.9); everything new goes into v2.

## Running v2

```bash
source ../setup.sh                       # LCG_110 + trex-fitter (see ../CLAUDE.md)
python run_v2.py --from 2                # T&P → control → histograms → fit → report (~40 min)
python run_v2.py --from 2 --max-files 2 --workers 4      # smoke test on 2 files per sample
python scripts/v2_1_skim.py              # (re)make the skims: ~1.5 h, resumable, 8 workers
```

| step | script | writes | needs |
|---|---|---|---|
| 0 | `scripts/v2_0_filelists.py` | `filelists/*.sizes.tsv` (committed) | network |
| 1 | `scripts/v2_1_skim.py` | `$BND_SKIM_DIR/<sample>/<file>.root` + `.json`, `manifest.json` | dCache/EOS |
| 2 | `scripts/v2_2_tnp.py` | `output/v2/tnp/tnp_result.json`, `pileup_weights.json`, `tnp_fits.pkl`, plots `tnp_*` | data + DY_NLO skims |
| 3 | `scripts/v2_3_control.py` | `output/v2/momentum.json`, `fakes.json`, `fakes_templates.pkl`, plots `momentum_*`, `ff_*` | all skims, step 2 |
| 4 | `scripts/v2_4_histograms.py` | `output/v2/histograms.pkl` (+ `hist_parts/`), `gensums.json` | steps 2, 3 |
| 5 | `scripts/v2_5_fit.py` | `fit/fitinputs/zmumu.root`, `fit/zmumu.config`, `fit/results/zmumu/`, `fit/results/zmumu_fit_result.json` | step 4, trex-fitter |
| 6 | `scripts/v2_6_report.py` | `output/v2/RESULTS_v2.md`, `results_v2.json`, plots `datamc_*`, `summary_v2.png` | step 5 |

Every per-file step is resumable (parts in `output/v2/*_parts/`; delete a part to redo it) and
has `--summarise-only` / `--merge-only`. Change a cut → rerun from step 2 (T&P) or 4 (histograms).

## Where things are defined (single sources of truth)

- `zmumu/samples.py` — datasets, recids, cross sections (+provenance), dCache/EOS locations.
- `zmumu/regions.py` — muon classes (loose/tight/anti-tight), trigger matching, SR/SS/CRemu/FF regions.
- `zmumu/skim.py` — skim categories A/B/C, branch list, `GenSums`.
- `zmumu/gen.py` — LHE flavour split, dressed/Born fiducial definitions.
- `zmumu/weights.py`, `zmumu/pileup.py` — MC weights; theory variations renormalised to constant fiducial yield.
- `zmumu/tnp.py` — T&P definitions, pass/fail fit model, `ScaleFactors` (event weights incl. trigger).
- `zmumu/fakes.py`, `zmumu/momentum.py` — fake factor; Z-peak momentum calibration.
- `zmumu/histograms.py` — variables, binnings, systematic variations; `zmumu/plotting.py` — data/MC plots.
- `scripts/v2_5_fit.py` — TRExFitter config (NP list, categories) and the σ extraction.
- `../fitting/CONVENTIONS.md` — histogram/NP naming shared with the other channels.

## Traps (all cost time once)

- Awkward: indexing a jagged array with a flat integer array selects *events*; use `objects.take`.
- Skim `TrigObj` is filtered to muons; all other collections are unfiltered (cross-indices valid).
- `GRL.txt` is the full-2016 JSON: the skim restricts to runs 278820–284044.
- MC skims have **no trigger requirement**; data skims do. Category-C events carry `skim_prescale = 10`.
- Only prompt–prompt MC events (`Muon_genPartFlav` 1 or 15) enter SR/SS/CRemu; non-prompt is the fake factor.
- Signal theory variations must keep the fiducial yield fixed (`Weighter.theory_renorm`).
- dCache stalls above ~8 readers; `batch.run_files` retries from EOS. Do not import ROOT in python.
- The powheg `SigModel` template is normalised to the NLO fiducial prediction (its cross section is not used).
- `hist` histograms for TRExFitter need `storage.Weight()` (Sumw2), see `fitting/trexhist.py`.
