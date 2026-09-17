# CLAUDE.md — z-mumu (Z → μ⁺μ⁻ cross section)

For agents continuing this channel. Physics write-up: `docs/` (start with `docs/09-v2-overview.md`);
numbers: `output/v2/RESULTS_v2.md` and `output/v2/results_v2.json`; combination inputs: `handoff.md`.

> **Frozen on 17 Sep 2026 (git tag `zmumu-freeze-2026-09-17`).** The result is v2 with the measured muon
> reconstruction scale factor and the 60–120 GeV acceptance uncertainties of `docs/16`:
> σ_fid = 790.2 ± 0.2 (stat) ± 6.1 (syst) ± 9.6 (lumi) pb, σ(60–120) = 1931 ± 30 pb. No further reruns or
> improvements: read, cite and present these numbers; do not regenerate `fit/`, `output/v2/` or the plots.
> The 15 Sep result (reconstruction SF assigned, m > 50 GeV acceptance uncertainties) is
> `fit/results/zmumu_v2_15sep_fit_result.json`.

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
python run_v2.py --from 2                # T&P → control → histograms → reco T&P → acceptance → fit → report → CMS comparison (~1 h)
python run_v2.py --from 5                # fit, stability variants, report, CMS comparison (~20 min)
python run_v2.py --from 2 --max-files 2 --workers 4      # smoke test on 2 files per sample
python scripts/v2_1_skim.py              # (re)make the skims: ~1.5 h, resumable, 8 workers
```

| step | script | writes | needs |
|---|---|---|---|
| 0 | `scripts/v2_0_filelists.py` | `filelists/*.sizes.tsv` (committed) | network |
| 1 | `scripts/v2_1_skim.py` | `$BND_SKIM_DIR/<sample>/<file>.root` + `.json`, `manifest.json` | dCache/EOS |
| 2a | `scripts/v2_2_pileup.py` | `output/v2/tnp/pileup_weights.json` (N_PV-matched profile), `pileup_npv.pkl`, plots `pu_*` | data + DY_NLO skims |
| 2 | `scripts/v2_2_tnp.py` | `output/v2/tnp/tnp_result.json`, `tnp_fits.pkl`, plots `tnp_*` | step 2a, data + DY_NLO skims |
| 3 | `scripts/v2_3_control.py` | `output/v2/momentum.json`, `fakes.json`, `fakes_templates.pkl`, plots `momentum_*`, `ff_*` | all skims, step 2 |
| 4 | `scripts/v2_4_histograms.py` | `output/v2/histograms.pkl` (+ `hist_parts/`), `gensums.json` | steps 2, 3 |
| 5 | `scripts/v2_5_fit.py` | `fit/fitinputs/zmumu.root`, `fit/zmumu.config`, `fit/results/zmumu/`, `fit/results/zmumu_fit_result.json` with the reconstruction SF applied and the acceptance block of `zmumu/acceptance.py` (`--rebin`, `--tag`, `--no-sigmodel`, `--quick` for variants; `--no-reco-sf` for the 15 Sep setup) | steps 4, 7, 8, trex-fitter |
| 5b | `scripts/v2_5_fit_variants.py` | `fit/results/stability.json` (μ_Z vs binning / SigModel / smoothing) | step 5 |
| 6 | `scripts/v2_6_report.py` | `output/v2/RESULTS_v2.md`, `results_v2.json`, plots `datamc_*`, `summary_v2.png` | step 5 |
| 7 | `scripts/v2_7_reco_tnp.py` | `output/v2/tnp/reco_result.json`, `reco_fits.pkl`, plots `reco_*` (muon reconstruction SF, T&P on the **unskimmed** NanoAOD, ~15 min) | step 2a, dCache parents |
| 8 | `scripts/v2_8_theory_acceptance.py` | `output/v2/cms_parity/theory_acceptance.json`, plot `cms_acc_vs_ptz.png` (A vs pT(Z), bare vs dressed, 25/25 volume) | DY NLO parents from EOS |
| 9 | `scripts/v2_9_cms_parity.py` | `output/v2/cms_parity/parity_zmumu.{json,md}`, `inputs.json`, plots `cms_parity_*` (budget vs CMS-SMP-20-004, three options, comparison) | steps 5, 7, 8 |

`run_v2.py` runs them in the order 0 1 2 3 4 7 8 5 (5b) 6 9: steps 7 and 8 feed the fit.
Every per-file step is resumable (parts in `output/v2/*_parts/`; delete a part to redo it) and
has `--summarise-only` / `--merge-only`. Change a cut → rerun from step 2 (T&P) or 4 (histograms).

## Where things are defined (single sources of truth)

- `zmumu/samples.py` — datasets, recids, cross sections (+provenance), dCache/EOS locations.
- `zmumu/regions.py` — muon classes (loose/tight/anti-tight), trigger matching, SR/SS/CRemu/SSemu/FF regions, lepton-cleaned jets.
- `zmumu/skim.py` — skim categories A/B/C, branch list, `GenSums`.
- `zmumu/gen.py` — LHE flavour split, dressed/Born fiducial definitions.
- `zmumu/weights.py`, `zmumu/pileup.py` — MC weights; theory variations renormalised to constant fiducial yield.
- `zmumu/tnp.py` — T&P definitions, pass/fail fit model, `ScaleFactors` (event weights incl. trigger).
- `zmumu/fakes.py`, `zmumu/momentum.py` — fake factor; Z-peak momentum calibration.
- `zmumu/histograms.py` — variables, binnings, systematic variations; `zmumu/plotting.py` — data/MC plots.
- `scripts/v2_5_fit.py` — TRExFitter config (NP list, categories), the reconstruction SF on the templates, and the σ extraction.
- `zmumu/acceptance.py` — A(60–120), A(m>50) and every acceptance uncertainty (PDF, αs, scales, boson pT, generator, PS FSR, QED FSR estimate, MC stat); the `acceptance` block of the fit result and of `zmumu.root.meta.json`, which the combination reads.
- `zmumu/recoeff.py` — reconstruction-efficiency probes (stand-alone, isolated track) for step 7.
- `../fitting/uncertainty_parity.py`, `../fitting/UNCERTAINTY_PARITY.md` — the comparison with a published budget (docs/16).
- `../fitting/CONVENTIONS.md` — histogram/NP naming shared with the other channels.

## Traps (all cost time once)

- Awkward: indexing a jagged array with a flat integer array selects *events*; use `objects.take`.
- Skim `TrigObj` is filtered to muons; all other collections are unfiltered (cross-indices valid).
- `GRL.txt` is the full-2016 JSON: the skim restricts to runs 278820–284044.
- MC skims have **no trigger requirement**; data skims do. Category-C events carry `skim_prescale = 10`.
- Only prompt–prompt MC events (`Muon_genPartFlav` 1 or 15) enter SR/SS; non-prompt muons are the fake factor. The e-μ regions (CRemu, SSemu) keep *all* MC (W+jets included): nothing replaces the non-prompt electrons there (REVIEW.md F1).
- `mass_fit` is stored in 1 GeV bins; the fit rebins to 5 GeV (`--rebin`). 2 GeV bins over-constrain the shape NPs (REVIEW.md F3/F5).
- Jet multiplicities are lepton-cleaned (`regions.clean_jet_count`).
- Signal theory variations must keep the fiducial yield fixed (`Weighter.theory_renorm`).
- dCache stalls above ~8 readers; `batch.run_files` retries from EOS. Do not import ROOT in python.
- The powheg `SigModel` template is normalised to the NLO fiducial prediction (its cross section is not used).
- `hist` histograms for TRExFitter need `storage.Weight()` (Sumw2), see `fitting/trexhist.py`.
- The reconstruction efficiency *is* measurable in NanoAODv9 (`Muon_isStandalone`, `IsoTrack`), but not in the skims (a failing probe is not a loose muon): step 7 reads the parents.
- `A_60_120` theory uncertainties must use the `lhe_mumu_60_120` denominator; `mc_acceptance.py` (and so the fit JSON until 16 Sep) used m > 50 (docs/16 §4.2). `zmumu/acceptance.py` does it right; the powheg sample only covers 50 < m_LHE < 120, so its generator row is taken from the 60–120 evaluation for σ(m>50) too.
- The combination reads `acceptance.A_*_rel` from `fit/fitinputs/zmumu.root.meta.json` by key (`combination/combLieke/config/channels.json`): renaming a key breaks the MultiFit.
