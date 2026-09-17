# Z → ττ cross section (τhτh + μτh + eτh + eμ)

CMS 2016 Open Data (Run2016G+H, 13 TeV, 16.4 fb⁻¹). Four ττ final states are measured and fitted together,
following the CMS 13 TeV measurement (arXiv:1801.03535): **τhτh** (Tau dataset, sorted into three BDT
categories), **μτh** (SingleMuon), **eτh** (SingleElectron) and **eμ** (MuonEG) with an eμ tt̄ control
region. The τh identification scale factors (per decay mode) and the τh energy scale are **measured in situ**
by the fit rather than taken from the TauPOG; jet→τh fakes in the lepton channels come from a per-process
fake-factor method (multijet / W+jets / tt̄ fractions), the eμ multijet from same-sign data, tt̄ from the
control region, and lepton efficiencies from the official POG corrections plus in-situ trigger measurements.
There is deliberately **no μμ channel** and every channel vetoes a second muon or electron, so this
measurement is orthogonal to the Z→μμ and Z→ee selections of the other groups (with one exception: our eμ
region overlaps z-mumu's tt̄ control region, see `docs/11-combination-inputs.md`).

Design, selections and the systematics map against the paper's Table 2: [`docs/10-v4-plan.md`](docs/10-v4-plan.md).
Numbers: [`output/RESULTS.md`](output/RESULTS.md). For the combination: [`docs/11-combination-inputs.md`](docs/11-combination-inputs.md).
The measurement was reviewed in [`REVIEW_v4.md`](REVIEW_v4.md); what was changed in answer to it is
[`REVIEW_v4_RESPONSE.md`](REVIEW_v4_RESPONSE.md).

<!-- RESULT:BEGIN -->
**Result** (four channels, τh ID scale factors and energy scale fitted in situ)

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 1981 +73 −70 pb   (stat ±9;  prediction 1945 pb: aMC@NLO acceptance, NNLO normalisation)
μ_Z = 1.019 +0.037 −0.036   (stat ±0.005, syst ±0.036),   expected ±0.036,   GoF p = 0.15,   μ_tt̄ = 1.11 ± 0.04
τh ID SF (Tight):  DM0 0.989 ± 0.040 (POG 0.90 ± 0.13)  DM1 0.963 ± 0.034 (POG 0.89 ± 0.05)  DM10 0.896 ± 0.035 (POG 0.94 ± 0.15)  DM11 0.795 ± 0.050 (POG 0.81 ± 0.15)
τh energy scale:   DM0 -0.6 ± 0.9 %  DM1 -0.2 ± 0.6 %  DM10 +0.6 ± 1.0 %  DM11 +3.2 ± 2.2 %
μ_Z of the two sub-measurements the fit combines:  e mu alone 0.960 +0.042 −0.040;  τ channels alone (SF free) 1.203 +0.091 −0.084   (2.6 σ apart)
per channel alone (POG SFs fixed):  tautau 1.043 +0.079 −0.072  mutau 1.048 +0.055 −0.051  etau 1.016 +0.070 −0.065  emu 0.960 +0.042 −0.040
```

Largest grouped impacts on μ_Z: NormFactors 3.1 %, Emu trigger 2.5 %, Gammas 1.4 %, Tau trigger 1.2 %, Fakes 1.2 %, Electron energy 1.2 %; data statistics 0.5 %. The categories overlap, so their quadrature sum exceeds the MINOS total by a factor 1.46 (`output/RESULTS.md`).

Ranking: mu_ttbar 3.1 %, EmuTrigger 2.6 %, TauIDSF_DM1 2.2 %, TauIDSF_DM10 2.1 %, TauIDSF_DM0 1.9 %, TauIDSF_DM11 1.3 %.
<!-- RESULT:END -->

## Quick start

```bash
source ../setup.sh                  # LCG_110: python 3.13, uproot, awkward, hist, mplhep, xgboost; TRExFitter on PATH
python run_tautau_base.py --from 3  # τhτh fake factors, BDT and base templates (~25 min from the ntuples)
python run_all.py --from 3          # trigger efficiencies, fakes, templates, all fits, report (~4 h with the ranking)
python run_all.py --only 5          # just the fits (~1 h without the ranking: see run_all.py --help)
```

No `pip install` and no venv are needed: LCG_110 has everything. On a laptop, copy the ntuples
(`docs/03-skims.md`), set `BND_TAUTAU_CACHE`, and run from step 3.

## Layout

```
ztautau/                  the library -- every physics choice is in config.py
  config.py               selection (all four channels), fake-factor binning and closure corrections, BDT settings,
                          categories, fit binning, priors, paths
  samples.py              data + simulation registry (incl. jet-binned DY), cross sections, file lists
  batch.py                subprocess-per-file runner: resumable, kills stalled dCache reads, EOS fallback
  io.py  gen.py           golden JSON, trigger OR, MET filters; LHE flavour split, fiducial volume
  skim.py                 NanoAOD -> skim (+ GenSums), v1 (τhτh) and v4 (lepton streams) preselections
  objects.py  leptons.py  τh / lepton / jet selection, trigger matching, pair choice (τhτh; lepton channels)
  mass.py                 m_vis, collinear mass, MET-likelihood mass m_tt
  corrections.py  pog.py  TauPOG ID / ES / trigger SFs, pileup weights; Muon/EGM/BTV/JME POG json readers
  fakes.py  fakes_v4.py   τhτh fake factors; per-process fake factors of the lepton channels
  bdt.py                  k-fold XGBoost classifier (features, folds, training, scoring, categories)
  analysis.py             τhτh ntuples, stitching, regions, weights, systematic variations, categories
  analysis_v4.py          the four-channel layer: lepton-channel ntuples, regions, weights, in-situ trigger SFs,
                          theory variations renormalised to σ(60-120), systematic registry
  plotting.py             CMS-style data / prediction plots
scripts/step*.py          one standalone script per step (docstring = usage)
run_tautau_base.py        the τhτh base chain (steps 0-4a) -> fit/fitinputs/tautau_base.root
run_all.py                the measurement (steps 1-6), including every cross-check fit
review/                   the numerical studies behind REVIEW.md
external/                 TauPOG / POG inputs converted to JSON, and the in-situ trigger tables (committed)
filelists/                xrootd URL, size, adler32 of every input file
fit/                      ztautau*.config (generated, committed), comb.config (MultiFit), bdt_info.json,
                          fitinputs/*.root (committed: the combination reads them), results/ (fit outputs)
output/plots/             all plots (PNG committed): step3_* fake factors, step3b_* BDT, step4_* control plots,
                          fit_* TRExFitter, step6_* summary
output/RESULTS.md, results.json   the numbers
docs/                     physics documentation, read 00-overview.md first
slides/                   the review deck: make_figures.py (dark vector figures -> figs/) + build_deck.py -> ztautau_slides.pdf
```

## Documentation

| page | topic |
|---|---|
| [00-overview](docs/00-overview.md) | the measurement in one page |
| [01-data-and-samples](docs/01-data-and-samples.md) | data sets, simulation (incl. jet-binned DY stitching, fiducial split), cross sections |
| [02-selection](docs/02-selection.md) | τhτh trigger, τh ID, pair choice, vetoes, regions, cutflow |
| [03-skims](docs/03-skims.md) | skims, ntuples, **laptop bundle** |
| [04-ditau-mass](docs/04-ditau-mass.md) | **MET correction of the visible mass**: collinear vs likelihood mass |
| [05-fake-factors](docs/05-fake-factors.md) | τhτh **fake factors**, MC subtraction, closure corrections, OS/SS, the uncertainty model |
| [06-cross-section](docs/06-cross-section.md) | fiducial volume, the non-fiducial background, A, C, from μ to σ |
| [07-corrections-and-systematics](docs/07-corrections-and-systematics.md) | TauPOG corrections, every τhτh nuisance parameter |
| [08-fit-and-results](docs/08-fit-and-results.md) | fit set-up and the τhτh-only result (superseded as *the* result by 10) |
| [09-bdt](docs/09-bdt.md) | **the k-fold BDT**: inputs, training, validation, categories |
| [10-v4-plan](docs/10-v4-plan.md) | **the four-channel measurement**: channels, objects, orthogonality, fakes per channel, in-situ τh ID / ES, systematics vs the CMS paper, the fit model, the result |
| [11-combination-inputs](docs/11-combination-inputs.md) | **what the combination gets**: files, region and NP names, correlations, overlaps, the two routes, the caveats |
| [REVIEW.md](REVIEW.md), [REVIEW_v4.md](REVIEW_v4.md), [REVIEW_v4_RESPONSE.md](REVIEW_v4_RESPONSE.md) | the reviews and what was done about them |

---

## The τhτh part in detail (the backbone of the measurement)

CMS 2016 Open Data (Tau dataset, Run2016G+H, 13 TeV, 16.4 fb⁻¹). Both τ decay hadronically. Fakes from a
data-driven fake factor (MC-subtracted, closure-corrected), Z→ττ and the small backgrounds from UL16
simulation, MET-corrected di-τ mass, a **k-fold BDT** that sorts the signal region into three categories,
and a binned profile-likelihood fit of m_ττ in the three categories with **TRExFitter v1.8.0**, using the
shared `fitting/` conventions of all channels. This is the third iteration (v3): the first one was reviewed in [`REVIEW.md`](REVIEW.md), v2 fixed the review
findings, v2.1 fixed the fake-factor pulls, and v3 adopts the DeepTau **Tight** working point on both legs as the
nominal (section "What changed" below).

**Historical** — the τhτh channel fitted on its own, with the TauPOG scale factors as priors and the τhτh
visible fiducial volume as the signal (v3, superseded by the four-channel fit above):

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2082 ± 41 (stat) +222/−194 (syst+stat) ± 76 (acc) pb     NNLO: 1945 pb
σ_fid(τhτh, vis. pT > 40 GeV, |η| < 2.1, 60 < m < 120 GeV) = 4.82 ± 0.09 (stat) ± 0.47 (syst) pb   pred.: 4.50 pb
μ_Z = 1.071 +0.114 -0.100   (DeepTau Medium, same chain: 1.205 +0.145 −0.125 — see docs/08 for the comparison)
```

Method details: [`docs/08-fit-and-results.md`](docs/08-fit-and-results.md), slides in [`slides/`](slides/).
The τhτh-only fit is no longer produced; `run_tautau_base.py` stops at the templates.

## What changed (v1 → v2 → v2.1 → v3 → v4)

| review finding | what was done |
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

| v4: τh ID scale factors dominate (11 %) | decay-mode regions in μτh / eτh and free `TauIDSF_DM*` NormFactors: the scale factors are measured in situ at 4–5 % |
| v4: trigger efficiency from the POG turn-on | in-situ Ele27 and eμ cross-trigger efficiencies from the other stream's events (`step3c_trigger.py`) |
| `REVIEW_v4.md` (1, 5, 9) | Asimov MINOS, the order of the TRExFitter steps, the over-shoot of the grouped impacts, report hygiene |
| `REVIEW_v4.md` (2, 4) | the eμ and τ-channel sub-measurements are quoted separately; the pT dependence of the τh ID scale factor is measured (`ztautau_ptsplit`) |
| `REVIEW_v4.md` (3, 6, 7) | the eμ trigger, ℓτh fake-composition and Ele27 turn-on priors corrected (`REVIEW_v4_RESPONSE.md`) |

Still not done (see `CLAUDE.md`, open issues): HT-binned W+jets, the high-mass DY sample for the
non-fiducial template, EWK Z→ττ, SM H→ττ.

## Data locations

| what | where | size |
|---|---|---|
| Tau NanoAOD | `/dcache/atlas/sjankovy/BND/collision_data/Tau/Run2016{G__30532,H__30565}` | 160 GB |
| simulation | `/dcache/atlas/sjankovy/BND/mc/` + CERN EOS (tt̄, single top) | |
| skims | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/skims_v1/` (DY re-skimmed; the rest links to the v1 skims of N. ter Linden) | 1.7 GB |
| **ntuples (laptop)** | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/ntuples_v1/` | **330 MB** |
| BDT models | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/bdt/` (Tight nominal; `bdt_medium/` the v2.1 models; recreated by step 3b in ~2 min) | 5 MB |
| v4 skims | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/skims_v4/` (SingleMuon, SingleElectron, MuonEG G+H, all simulation; `skim_cat` bitmask) | 17 GB |
| **v4 ntuples** | `/data/atlas/users/sjankovy/BND-school-cache/ztautau/ntuples_v4/` (`<sample>_<channel>.root`, channels mutau / etau / emu / emu_mu / emu_el) | 2.9 GB |
| SingleElectron NanoAOD | CERN EOS only (`root://eospublic.cern.ch//eos/opendata/cms/Run2016{G,H}/SingleElectron/...`, 151 files, 257 GB; no dCache copy) | |
