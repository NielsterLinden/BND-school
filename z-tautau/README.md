# Z → τhτh cross section

CMS 2016 Open Data (Tau dataset, Run2016G+H, 13 TeV, 16.4 fb⁻¹). Both τ decay hadronically. Fakes from a
data-driven fake factor, backgrounds and signal from UL16 simulation. MET-corrected di-τ mass, binned
profile-likelihood fit with **TRExFitter v1.8.0**, using the shared `fitting/` conventions of all channels.

**Result** (nominal: fake factors without MC subtraction)

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2255 ± 36 (stat) +371/−315 (syst) ± 83 (acceptance) pb     NNLO: 1945 pb
σ_fid(τhτh, vis. pT > 40 GeV, |η| < 2.1, 60 < m < 120 GeV) = 5.22 ± 0.08 (stat) ± 0.80 (syst) pb   pred.: 4.50 pb
μ_Z = 1.160 +0.192 −0.163
```

With genuine-τ subtraction in the fake-factor regions: μ_Z = 1.208 +0.201 −0.169. The uncertainty is
dominated by the TauPOG τh ID scale factors (12.5 %), signal modelling (9.8 %) and trigger scale factors
(5.1 %); data statistics are 1.8 %. Details: [`output/RESULTS.md`](output/RESULTS.md),
[`docs/08-fit-and-results.md`](docs/08-fit-and-results.md), slides in [`slides/`](slides/).

## Quick start

```bash
source ../setup.sh                  # LCG_110: python 3.13, uproot, awkward, hist, mplhep, ROOT; TRExFitter on PATH
python run_all.py --from 3          # from the existing ntuples: fake factors, histograms, fits, report (~10 min)
python run_all.py                   # everything, incl. skims from dCache/EOS (~1 h)
python run_all.py --only 5          # just the fits
```

No `pip install` and no venv are needed: LCG_110 has everything. On a laptop, copy the ntuples
(`docs/03-skims.md`), set `BND_TAUTAU_CACHE`, and run `--from 3`.

## Layout

```
ztautau/                  the library -- every physics choice is in config.py
  config.py               selection, fake-factor binning, fit binning, paths
  samples.py              data + simulation registry, cross sections, file lists (dCache manifest / Open Data API)
  batch.py                subprocess-per-file runner: resumable, kills stalled dCache reads, EOS fallback
  io.py  gen.py           golden JSON, trigger OR, MET filters; LHE flavour split, fiducial volume
  skim.py                 NanoAOD -> skim (+ GenSums: generator-weight sums over all events)
  objects.py              tau / lepton / jet selection, trigger matching, pair choice
  mass.py                 m_vis, collinear mass, MET-likelihood mass m_tt
  corrections.py          TauPOG ID / ES / trigger SFs, pileup weights
  fakes.py                fake factors (era x DM x N_jets x pT), OS/SS correction
  analysis.py             ntuple loading, regions, MC weights, systematic variations
  plotting.py             CMS-style data / prediction plots
scripts/step0..6_*.py     one standalone script per step (docstring = usage)
run_all.py                runs the steps in order
external/                 TauPOG inputs converted to JSON (committed, 256 kB)
filelists/                xrootd URL, size, adler32 of every input file
fit/ztautau*.config       generated TRExFitter configs (committed); fitinputs/ and results/ are git-ignored
output/plots/             all plots (PNG committed): step3_* fake factors, step4_* control plots, fit_* TRExFitter, step6_* summary
output/RESULTS.md, results.json   the numbers
docs/                     physics documentation, read 00-overview.md first
slides/                   the summary slide deck (PDF + LaTeX source)
```

## Documentation

| page | topic |
|---|---|
| [00-overview](docs/00-overview.md) | the measurement in one page |
| [01-data-and-samples](docs/01-data-and-samples.md) | data sets, simulation, cross sections |
| [02-selection](docs/02-selection.md) | trigger, τh ID, pair choice, vetoes, regions, cutflow |
| [03-skims](docs/03-skims.md) | skims, ntuples, **laptop bundle** |
| [04-ditau-mass](docs/04-ditau-mass.md) | **MET correction of the visible mass**: collinear vs likelihood mass |
| [05-fake-factors](docs/05-fake-factors.md) | **fake factors**, OS/SS, closure, the MC-subtraction question |
| [06-cross-section](docs/06-cross-section.md) | fiducial volume, A, C, from μ to σ |
| [07-corrections-and-systematics](docs/07-corrections-and-systematics.md) | TauPOG corrections, every nuisance parameter |
| [08-fit-and-results](docs/08-fit-and-results.md) | fit set-up, result, impacts, pulls, input for the combination |

## Data locations

| what | where | size |
|---|---|---|
| Tau NanoAOD | `/dcache/atlas/sjankovy/BND/collision_data/Tau/Run2016{G__30532,H__30565}` | 160 GB |
| simulation | `/dcache/atlas/sjankovy/BND/mc/` + CERN EOS (tt̄, single top) | |
| skims | `/data/atlas/users/nterlind/BND-school-cache/ztautau/skims_v1/` | 1.6 GB |
| **ntuples (laptop)** | `/data/atlas/users/nterlind/BND-school-cache/ztautau/ntuples_v1/` | **300 MB** |
