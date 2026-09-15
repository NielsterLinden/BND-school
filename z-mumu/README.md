# Z -> mu+ mu- cross section

CMS 2016 Open Data (SingleMuon, Run2016G+H, 13 TeV, 16.39 fb^-1). Two measurements live in
this folder; **v2 is the result**, v1 is kept as the documented first iteration.

## v2 result (MC-based, TRExFitter binned-likelihood fit) -- `run_v2.py`

```
sigma_fid(pp -> Z/gamma* -> mu+ mu-) = 791.8 +/- 0.2 (stat) +/- 6.2 (syst) +/- 9.3 (lumi) pb
```

in the fiducial volume: two opposite-sign muons, `pT > 26/20 GeV`, `|eta| < 2.4`,
`60 < m(mu mu) < 120 GeV`, dressed leptons (dR < 0.1). Signal strength `mu_Z = 0.990 +/- 0.014`
against the aMC@NLO prediction of 799.6 pb; goodness of fit p = 0.33; counting cross-check
`(N_obs - N_bkg)/(C L)` = 794.4 pb. Inclusive: `sigma(60 < m < 120 GeV) = 1935 +/- 30 pb`
(A = 0.4092), `sigma(m > 50 GeV) = 2006 +/- 31 pb` (A = 0.3947). Reviewer's independent
number for the same volume: 797.2 pb; v1: 773.2 pb.

What v2 does that v1 could not (details in [`docs/09-v2-overview.md`](docs/09-v2-overview.md)):
own skims of the unskimmed NanoAOD (trigger objects, prefiring weights, 1-muon e-mu events),
simulation for signal and every prompt background (DY, ttbar, tW, WW, WZ, ZZ), pileup and
L1-prefiring weights, tight-ID / isolation / trigger scale factors from tag-and-probe **fits**,
a data-driven **fake-factor** estimate of non-prompt muons (no QCD simulation), a Z-peak muon
momentum calibration, and a **TRExFitter v1.8.0** fit of the mass spectrum with 24 nuisance
parameters. Full numbers: [`output/v2/RESULTS_v2.md`](output/v2/RESULTS_v2.md); machine-readable:
`output/v2/results_v2.json`, `fit/results/zmumu_fit_result.json`; plots: `output/v2/plots/`;
fit outputs: `fit/results/zmumu/`. **Review (15 Sep 2026): [`REVIEW.md`](REVIEW.md), slides `review/deck/zmumu_review.pdf` -- the fit's central value depends on the shape model at the 1% level; see there before quoting the number.** Twelve-slide summary deck: <https://claude.ai/code/artifact/8072b961-01d5-47ed-9804-ee21f70ac985>. Agents: read [`CLAUDE.md`](CLAUDE.md).

Uncertainty budget (impact on the cross section): luminosity 1.16%, L1 prefiring 0.51%,
muon efficiencies 0.50%, signal modelling (PDF, scales, parton shower, generator) 0.36%,
MC statistics 0.35%, muon momentum 0.27%, pileup 0.16%, backgrounds 0.09%, fakes 0.05%,
statistics 0.03%. Total 1.4% (0.8% without the luminosity).

## v1 result (data-only counting) -- `run_all.py`

```
sigma_fid(pp -> Z -> mu+ mu-) = 773.2 +/- 0.2 (stat) +/- 11.9 (syst) pb      (frozen)
```

Same volume, no simulation, medium ID, luminosity without normtag (0.63% low). Full numbers in
[`output/RESULTS.md`](output/RESULTS.md), plots in `output/plots/`. The review
([`docs/08-review-and-crosschecks.md`](docs/08-review-and-crosschecks.md)) found the trigger
efficiency biased high, background in the tag-and-probe, the normtag luminosity, a 1.1%
migration effect, and computed the acceptance; the external reviewer added the L1-prefiring
(2%) and e-mu skim issues (`agent_reference/2026.09.14_1400_revisionpoints.md`). All of these
are addressed in v2.

## Quick start

```bash
# v2 (LCG_110 environment + TRExFitter; see ../CLAUDE.md for the one-time build)
source ../setup.sh
python run_v2.py --from 2                 # T&P -> control -> histograms -> fit -> report, ~40 min
python run_v2.py --from 2 --max-files 2   # smoke test
python scripts/v2_1_skim.py               # (re)make the skims from the NanoAOD parents, ~1.5 h

# v1 (its own python 3.9 venv, frozen)
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python run_all.py               # ~80 s on 12 cores
```

Do **not** `pip install` into the system or user environment.

## Layout

```
zmumu/                 the library (v1 modules + the v2 modules listed in CLAUDE.md)
  config.py            selection cuts, luminosity, binning, systematics
  io.py                file discovery, chunked reading, golden-JSON lumi mask
  objects.py           muon selection, FSR recovery, dimuon kinematics
  hists.py             histogram booking, persistence, CMS-style plotting
  stats.py             Clopper-Pearson intervals, error propagation
  batch.py             per-file subprocess runner (resumable) + CERN catalogue helpers

scripts/               one script per component, each runnable on its own
  step1_selection.py   OS / SS / e-mu regions, cutflow, mass spectra
  step2_efficiency.py  tag-and-probe ID, isolation, trigger efficiencies
  step3_backgrounds.py same-sign and e-mu data-driven backgrounds
  step4_signal.py      background-subtracted yield + fit validation
  step5_crosssection.py  combine and propagate uncertainties
  step6_report.py      results.json, RESULTS.md, summary plot
  xcheck_efficiency.py review: trigger-object T&P, T&P background, lumi (parent NanoAOD)
  mc_acceptance.py     review: acceptance A and method closure in DY simulation
  v2_0_filelists.py .. v2_6_report.py   the v2 chain (skim, T&P fits, control regions,
                       histograms, TRExFitter fit, report); v2_bundle.py laptop skims

run_all.py             runs v1 steps 1-6; run_v2.py runs the v2 chain
fit/                   v2 TRExFitter inputs (fitinputs/), config (zmumu.config) and results/
filelists/             cernopendata-client listings (xrootd URL, size, adler32) for the review scripts
docs/                  the physics documentation (read 00-overview.md first)
output/                plots + results (committed); output/data/ pickles are git-ignored
```

Each script has a full docstring explaining what it does and why; each has a
matching page in `docs/`.

## Documentation

| Page | Topic |
|---|---|
| [00-overview.md](docs/00-overview.md) | the measurement in one page |
| [01-selection.md](docs/01-selection.md) | cuts, the three regions, the trigger-menu wrinkle |
| [02-fsr-photons.md](docs/02-fsr-photons.md) | **FSR photons**: what they do and how they are recovered |
| [03-fake-leptons.md](docs/03-fake-leptons.md) | **fake / non-prompt muons**: sources and the same-sign estimate |
| [04-efficiency.md](docs/04-efficiency.md) | tag-and-probe, and what cannot be measured from NanoAOD |
| [05-backgrounds.md](docs/05-backgrounds.md) | both background classes and why you need both |
| [06-cross-section.md](docs/06-cross-section.md) | fiducial vs total, signal extraction, the formula |
| [07-systematics.md](docs/07-systematics.md) | every uncertainty and how far to trust it |
| [08-review-and-crosschecks.md](docs/08-review-and-crosschecks.md) | review: what was tested, what changes, inputs for the combination |
| [09-v2-overview.md](docs/09-v2-overview.md) | **v2**: why and how, the chain, the samples |
| [10-skims.md](docs/10-skims.md) | v2 skims: categories, branches, `GenSums`, laptop bundle |
| [11-mc-weights.md](docs/11-mc-weights.md) | normalisation, pileup from the luminosity record, L1 prefiring, scale factors |
| [12-tag-and-probe-fits.md](docs/12-tag-and-probe-fits.md) | pass/fail fits, systematics, trigger efficiency, application |
| [13-fake-factor.md](docs/13-fake-factor.md) | non-prompt muons: fake factor, template fit, closure |
| [14-fit-and-systematics.md](docs/14-fit-and-systematics.md) | regions, momentum calibration, nuisance parameters, extraction |
| [15-combination-inputs.md](docs/15-combination-inputs.md) | what the combination gets |

## Input data

v2 reads its own skims of the unskimmed NanoAOD parents (data: `/dcache/atlas/sjankovy/BND/collision_data/SingleMuon`,
MC: `/dcache/atlas/sjankovy/BND/mc/` and, for ttbar/tW, CERN EOS over xrootd), stored at
`$BND_SKIM_DIR` = `/data/atlas/users/nterlind/BND-school-cache/skims_v2` (24 GB, one file per
parent file, `manifest.json`) with a per-sample laptop bundle at
`/project/atlas/users/nterlind/BND-school-skims-lite/` (see `handoff.md`). v1 reads the
>= 2-muon skim `/dcache/atlas/kdevries/BND2026/DoubleMuonSkimmed/` (152 files, 31.5 GB).

## What limits the result

Not statistics (0.03%). The luminosity (1.2%) dominates; the next terms -- L1 prefiring,
muon efficiencies, signal modelling, MC statistics -- are each ~0.3-0.5% and could be reduced
with the official CMS correction files (prefiring maps, POG scale factors), more simulated
events, and Rochester corrections, none of which are reachable from this cluster. See
`docs/14-fit-and-systematics.md` and the open issues in [handoff.md](handoff.md).
