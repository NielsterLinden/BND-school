# Z → τhτh cross section

CMS 2016 Open Data (Tau dataset, Run2016G+H, 13 TeV, 16.4 fb⁻¹). Both τ decay hadronically. Fakes from a
data-driven fake factor (MC-subtracted, closure-corrected), Z→ττ and the small backgrounds from UL16
simulation, MET-corrected di-τ mass, a **k-fold BDT** that sorts the signal region into three categories,
and a binned profile-likelihood fit of m_ττ in the three categories with **TRExFitter v1.8.0**, using the
shared `fitting/` conventions of all channels. This is the third iteration (v3): the first one was reviewed in [`REVIEW.md`](REVIEW.md), v2 fixed the review
findings, v2.1 fixed the fake-factor pulls, and v3 adopts the DeepTau **Tight** working point on both legs as the
nominal (section "What changed" below).

**Result** (v3, nominal: DeepTau Tight on both legs, MC-subtracted fake factor)

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2082 ± 41 (stat) +222/−194 (syst+stat) ± 76 (acc) pb     NNLO: 1945 pb
σ_fid(τhτh, vis. pT > 40 GeV, |η| < 2.1, 60 < m < 120 GeV) = 4.82 ± 0.09 (stat) ± 0.47 (syst) pb   pred.: 4.50 pb
μ_Z = 1.071 +0.114 -0.100   (v2.1, DeepTau Medium, same chain: 1.205 +0.145 −0.125 — see docs/08 for the comparison)
```

Details: [`output/RESULTS.md`](output/RESULTS.md), [`docs/08-fit-and-results.md`](docs/08-fit-and-results.md),
slides in [`slides/`](slides/).

## Quick start

```bash
source ../setup.sh                  # LCG_110: python 3.13, uproot, awkward, hist, mplhep, xgboost; TRExFitter on PATH
python run_all.py --from 3          # from the existing ntuples: fake factors, BDT, histograms, fits, report (~25 min)
python run_all.py                   # everything, incl. skims from dCache/EOS (~2 h)
python run_all.py --only 5          # just the fits
```

No `pip install` and no venv are needed: LCG_110 has everything. On a laptop, copy the ntuples
(`docs/03-skims.md`), set `BND_TAUTAU_CACHE`, and run `--from 3`.

## Layout

```
ztautau/                  the library -- every physics choice is in config.py
  config.py               selection, fake-factor binning and closure corrections, BDT settings, categories, fit binning, paths
  samples.py              data + simulation registry (incl. jet-binned DY), cross sections, file lists
  batch.py                subprocess-per-file runner: resumable, kills stalled dCache reads, EOS fallback
  io.py  gen.py           golden JSON, trigger OR, MET filters; LHE flavour split, fiducial volume
  skim.py                 NanoAOD -> skim (+ GenSums: generator-weight sums over all events, per LHE_NpNLO bin)
  objects.py              tau / lepton / jet selection, trigger matching, pair choice
  mass.py                 m_vis, collinear mass, MET-likelihood mass m_tt
  corrections.py          TauPOG ID / ES / trigger SFs (genuine legs only), pileup weights
  fakes.py                fake factors (era x DM x N_jets x pT), closure corrections |eta(tau1)|, pT(tau2), OS/SS, per-event FF statistics
  bdt.py                  k-fold XGBoost classifier (features, folds, training, scoring, categories)
  analysis.py             ntuple loading, jet-binned stitching, fiducial split, regions, weights, systematic variations, categories
  plotting.py             CMS-style data / prediction plots
scripts/step0..6_*.py     one standalone script per step (docstring = usage); step3b_bdt.py trains and validates the BDT
run_all.py                runs the steps in order
review/                   the numerical studies behind REVIEW.md (working-point / pT scans, BDT proof of concept, eta non-closure)
external/                 TauPOG inputs converted to JSON (committed, 256 kB)
filelists/                xrootd URL, size, adler32 of every input file
fit/ztautau*.config       generated TRExFitter configs (committed); fit/bdt_info.json (BDT summary); fitinputs/ and results/ are git-ignored
output/plots/             all plots (PNG committed): step3_* fake factors, step3b_* BDT, step4_* control plots, fit_* TRExFitter, step6_* summary
output/RESULTS.md, results.json   the numbers
docs/                     physics documentation, read 00-overview.md first
slides/                   the review deck: make_figures.py (dark vector figures -> figs/) + build_deck.py (PyMuPDF) -> ztautau_slides.pdf
```

## Documentation

| page | topic |
|---|---|
| [00-overview](docs/00-overview.md) | the measurement in one page |
| [01-data-and-samples](docs/01-data-and-samples.md) | data sets, simulation (incl. jet-binned DY stitching, fiducial split), cross sections |
| [02-selection](docs/02-selection.md) | trigger, τh ID, pair choice, vetoes, regions, cutflow |
| [03-skims](docs/03-skims.md) | skims, ntuples, **laptop bundle** |
| [04-ditau-mass](docs/04-ditau-mass.md) | **MET correction of the visible mass**: collinear vs likelihood mass |
| [05-fake-factors](docs/05-fake-factors.md) | **fake factors**, MC subtraction, closure corrections, OS/SS, the uncertainty model |
| [06-cross-section](docs/06-cross-section.md) | fiducial volume, the non-fiducial background, A, C, from μ to σ |
| [07-corrections-and-systematics](docs/07-corrections-and-systematics.md) | TauPOG corrections, every nuisance parameter |
| [08-fit-and-results](docs/08-fit-and-results.md) | fit set-up, result, impacts, pulls, input for the combination |
| [09-bdt](docs/09-bdt.md) | **the k-fold BDT**: inputs, training, validation, categories |
| [REVIEW.md](REVIEW.md) | the review of v1 and what v2 did about it |

## What changed (v1 → v2 → v2.1 → v3)

| review finding | v2 |
|---|---|
| 3.1 signal template 38 % non-fiducial, high-mass tail in the fake sideband | `DYtautau` = fiducial part only (μ_Z), `DYtautau_nonfid` = theory-normalised background |
| 3.2 fake factor without MC subtraction | MC-subtracted FF is the only variant produced (`--ff-variant nosub` on request); W+jets subtracted with uniform weights (pathological raw weights) |
| 3.4 LO-vs-NLO `SigModel` (7 %, correlated with z-mumu) | computed, reported (`SigModel_tautau`, fiducial C), **not** in the fit by default; theory NPs (scales, PS, PDF) + data check of the pT spectra |
| 3.5 fake uncertainty model | closure corrections in \|η(τ1)\| and pT(τ2); FF statistics per event into the template variance; non-closure NPs per category and mass region |
| 3.6 MC statistics | jet-binned DY 0J/1J/2J stitched to the inclusive sample |
| 3.8 trigger SF on jet legs | SF = 1 for `genPartFlav = 0` legs |
| 4.3 BDT | k-fold XGBoost on mass-agnostic kinematics, three categories fitted in m_ττ |
| v2.1: fake NPs pulled 0.6–1.5 σ | C_OS/SS per (era, N_jets, BDT category); the fake-dominated category fitted above 110 GeV only (fake sideband) |
| v3: 4.1 Tight working point | DeepTau Tight on both legs is the nominal (`config.NOMINAL_WP`): 2.7× fewer fakes for 22 % less signal, every uncertainty group smaller, GoF p = 0.25 (`docs/08`, comparison with Medium) |
| 5 η(τ1) non-closure | diagnosed (FF depends on \|η\| by ±15 %), corrected |
| 3.9 second TRExFitter build | rebuilt with `fitting/build_trexfitter.sh` against the LCG ROOT |

Not done (see `CLAUDE.md`, open issues): decay-mode categories, in-situ trigger efficiency, HT-binned W+jets,
the high-mass DY sample for the non-fiducial template, EWK Z.

## Data locations

| what | where | size |
|---|---|---|
| Tau NanoAOD | `/dcache/atlas/sjankovy/BND/collision_data/Tau/Run2016{G__30532,H__30565}` | 160 GB |
| simulation | `/dcache/atlas/sjankovy/BND/mc/` + CERN EOS (tt̄, single top) | |
| skims | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/skims_v1/` (DY re-skimmed; the rest links to the v1 skims of N. ter Linden) | 1.7 GB |
| **ntuples (laptop)** | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/ntuples_v1/` | **330 MB** |
| BDT models | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/bdt/` (Tight nominal; `bdt_medium/` the v2.1 models; recreated by step 3b in ~2 min) | 5 MB |
