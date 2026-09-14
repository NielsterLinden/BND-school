# Z -> mu+ mu- fiducial cross section

CMS 2016 Open Data (SingleMuon, Run2016G+H, 13 TeV, 16.3 fb^-1).
A complete, data-only measurement: no simulation is used anywhere.

**Result**

```
sigma_fid(pp -> Z -> mu+ mu-) = 773.2 +/- 0.2 (stat) +/- 11.9 (syst) pb
```

in the fiducial volume: exactly two opposite-sign muons, `pT > 26/20 GeV`,
`|eta| < 2.4`, `60 < m(mu mu) < 120 GeV`, dressed lepton level.

Full numbers in [`output/RESULTS.md`](output/RESULTS.md), plots in
`output/plots/`, machine-readable in `output/results.json`.

**Review (Niels, 14 Sep 2026).** Steps 1-6 reproduce exactly. Cross-checks with
the unskimmed NanoAOD (trigger objects) and DY simulation find a biased trigger
efficiency, background in the ID/iso tag-and-probe, a luminosity without
normtag and a 1.1% migration effect. With those corrections, and the acceptance
from DY NLO for the combination:

```
sigma_fid = 776.9 +/- 0.2 (stat) +/- 14.8 (syst) pb        (DY NLO prediction: 799.6 pb)
sigma(pp -> Z/gamma* -> mu mu, m > 50 GeV) = 1968 pb       (A = 0.3947)
```

Not yet adopted in steps 1-6 -- see [`docs/08-review-and-crosschecks.md`](docs/08-review-and-crosschecks.md).

## Quick start

```bash
# one-time setup
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# the whole measurement (~80 s on 12 cores)
.venv/bin/python run_all.py

# a fast smoke test on 4 files
.venv/bin/python run_all.py --max-files 4

# rerun only the cheap steps after changing a systematic in config.py
.venv/bin/python run_all.py --from 3
```

Do **not** `pip install` into the system or user environment -- everything runs
from `.venv`.

## Layout

```
config.py is inside zmumu/ -- start there, every physics choice lives in it

zmumu/                 the library
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

run_all.py             runs steps 1-6 in order (the review scripts run separately)
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

## Input data

Skimmed NanoAOD at

```
/dcache/atlas/kdevries/BND2026/DoubleMuonSkimmed/
    Run2016G__30530/   70 files
    Run2016H__30563/   82 files
```

152 files, 31.5 GB, 80,191,719 events -- every event with at least two entries
in the `Muon` collection, from the 259 GB / 324 M event parent dataset at
`/dcache/atlas/sjankovy/BND/collision_data/SingleMuon`. The skim keeps ~300 of
the original 1363 branches.

The review scripts read the parent NanoAOD and the DY samples (recids 35669,
35671) from CERN EOS over xrootd (`xrootd` and `fsspec-xrootd` are in
`requirements.txt`). dCache NFS reads of the parent stalled under parallel
access on 14 Sep 2026, although the files match the catalogue checksums.

Note that 39 trigger branches are not present in every file (the 2016 HLT menu
changed during data-taking); `zmumu/io.py` handles this. The analysis triggers
`HLT_IsoMu24` and `HLT_IsoTkMu24` are present in all 152 files.

## What limits the result

Not statistics. With 10.8 million signal events the statistical uncertainty is
0.03%. The 1.5% total is dominated by the luminosity (1.2%), the external muon
reconstruction efficiency (0.8%), and an assigned trigger-method systematic
(0.5%). See [07-systematics.md](docs/07-systematics.md) and the open issues in
[handoff.md](handoff.md).
