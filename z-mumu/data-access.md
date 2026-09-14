# Z → μ⁺μ⁻: cluster environment and data access

Status on 2026-09-14. Setup and access tests only. No analysis has been run yet.

## 1. Cluster environment

| | |
|---|---|
| Node | `lvp036login.nikhef.nl`: 16 cores, 7.5 GiB RAM, shared with other users |
| OS | AlmaLinux 9.8 (Olive Jaguar), kernel `5.14.0-687.31.1.el9_8.x86_64` |
| Per-user limits on this node | memory cgroup `memory.max` = **2.00 GiB** (hard), no CPU quota |
| CVMFS | `/cvmfs/sft.cern.ch` mounted |
| Batch | **no HTCondor client on this node** (`condor_q`/`condor_submit` not found) |

**This is not a stoomboot node.** It is a login node, not `stbc-i*`. From here, `stbc-i1`, `stbc-i2` and `stbc-i3` accept connections on ssh port 22 (`stbc-i4` does not). One commit in this repo was made from `stbc-i3`. I did not check HTCondor on those nodes.

Environment (run this in every new shell):

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc14-opt/setup.sh
source /project/atlas/users/nterlind/venvs/bnd-opendata/bin/activate   # adds cernopendata-client
```

| Package | Version |
|---|---|
| Python | 3.13.11 |
| uproot | 5.7.1 |
| awkward | 2.9.0 |
| numpy | 2.4.4 |
| matplotlib | 3.11.0 |
| ROOT | 6.40.02 |
| xrdcp | v6.0.3 (the XRootD Python bindings also import) |
| cernopendata-client | 1.0.2 (in the venv) |

`LCG_110` is the newest numbered release on CVMFS. I only tested the `x86_64-el9-gcc14-opt` flavour. The venv sits on top of the LCG view, so uproot, awkward and ROOT still come from LCG. Only `cernopendata-client` and its dependencies are installed in the venv. To recreate it:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc14-opt/setup.sh
python3 -m venv --system-site-packages /project/atlas/users/nterlind/venvs/bnd-opendata
source /project/atlas/users/nterlind/venvs/bnd-opendata/bin/activate
pip install cernopendata-client
```

The venv is on `/project/atlas`, not in home, because the home quota is 1.2 GB (see §2). `pip install --user` would write into home.

Fallbacks (untested): the system `python3` is 3.9.25 and has no xrdcp. There is a personal miniconda at `/data/atlas/users/nterlind/venvs/miniconda3` (its `condabin` is on `PATH`). There is no mamba.

## 2. Storage

| Area | Free on filesystem | My quota (soft / hard) | Used by me | Use for |
|---|---|---|---|---|
| `/user/nterlind` (home) | 974 G | **1.2 G / 1.2 G** | 793 M | config only |
| `/project/atlas/users/nterlind` | 4.1 T of 12 T | 150 G / 1024 G | 11.4 G | repo, venvs, plots, small outputs |
| `/data/atlas/users/nterlind` | 1.4 T of 42 T (**97 % full**) | 2048 G / 3072 G | 167 G | **data cache** |
| `/tmp` | 9.3 G (root filesystem) | – | – | nothing large |
| `/dcache` | not mounted on this node | | | |

**Recommended cache directory: `/data/atlas/users/nterlind/BND-school-cache`** (already created). My quota there has about 1.9 TB left, which is enough for all 448 GB. The risk is that the filesystem is shared and 97 % full, so its 1.4 TB of free space can shrink.

## 3. Data inventory

| Dataset | recid | Events | Files | Size | File lists (in `z-mumu/filelists/`) |
|---|---|---:|---:|---:|---|
| SingleMuon Run2016G | 30530 | 149,916,849 | 70 | 119.60 GB | `SingleMuon_Run2016G_30530.{xrootd,http}.txt` |
| SingleMuon Run2016H | 30563 | 174,035,164 | 82 | 140.06 GB | `SingleMuon_Run2016H_30563.{xrootd,http}.txt` |
| DYJetsToLL_M-50 amcatnloFXFX (NLO) | 35669 | 71,839,442 | 41 | 91.18 GB | `DYJetsToLL_M-50_NLO_amcatnloFXFX_35669.{xrootd,http}.txt` |
| DYJetsToLL_M-50 madgraphMLM (LO) | 35671 | 82,448,537 | 61 | 97.27 GB | `DYJetsToLL_M-50_LO_madgraphMLM_35671.{xrootd,http}.txt` |
| **Total** | | **478,239,992** | **254** | **448.11 GB** | |

Each `*.txt` file has one URL per line. Each `*.sizes.tsv` file has three columns: xrootd URL, size in bytes, and adler32 checksum. Sizes are summed from the `.sizes.tsv` files, and event counts come from `datasets/datasets.tsv`. File counts and sizes match `datasets/SIZES.md` exactly. To regenerate a list:

```bash
cernopendata-client get-file-locations --recid 30530 --protocol xrootd            # or --protocol http
cernopendata-client get-file-locations --recid 30530 --protocol xrootd --verbose  # + size, checksum
```

## 4. Access test results

**xrootd works from this node.** TCP port 1094 to `eospublic.cern.ch` is open, and uproot opens `root://` URLs directly. The opendata.cern.ch HTTP API also answers (HTTP 200 in 0.19 s), but I did not test reading files with uproot over HTTP because xrootd was not blocked.

Test files:

- Data: `root://eospublic.cern.ch//eos/opendata/cms/Run2016G/SingleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/130000/0A4230E2-0C75-604D-890F-A4CE5E5C164E.root` (2.35 GB)
- MC: `root://eospublic.cern.ch//eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/30000/0082C29D-E74C-024A-BE9B-97B29EE7A4A2.root`

| | SingleMuon Run2016G file | DY NLO file |
|---|---:|---:|
| Events | **2,939,781** | **1,933,726** |
| Branches in `Events` | 1380 | 1504 |
| Trees | `tag, Events, LuminosityBlocks, Runs, MetaData, ParameterSets` | same |
| `HLT_*Mu*` branches | 178 | 191 |
| `uproot.open` time | 0.45 s | 0.39 s |

### Branch checks

All requested branches are present.

- Both files have `run`, `luminosityBlock`, `Muon_pt`, `Muon_eta`, `Muon_phi`, `Muon_charge`, `Muon_mass`, `Muon_tightId`, `Muon_mediumId`, `HLT_IsoMu24` and `HLT_IsoTkMu24`. The task listed `oMu24`, which I read as `HLT_IsoMu24`.
- The DY NLO file also has `genWeight`, `L1PreFiringWeight_Nom`, `LHEPart_pdgId`, `LHEPart_status` and `Pileup_nTrueInt`.
- The DY NLO `Runs` tree has `genEventSumw`, along with `genEventCount`, `genEventSumw2`, `LHEPdfSumw`, `LHEScaleSumw`, `nLHEPdfSumw`, `nLHEScaleSumw` and `run`. That file has 1 `Runs` entry with `genEventSumw` = 3.2887e10.
- The data `Runs` tree contains only `run`, as expected for data.

### Timings (Run2016G test file, 2.94 M events)

"Analysis columns" means these 11 branches: `run, luminosityBlock, HLT_IsoMu24, HLT_IsoTkMu24, nMuon, Muon_pt, Muon_eta, Muon_phi, Muon_charge, Muon_mass, Muon_tightId`.

| Test | Time | Notes |
|---|---:|---|
| Read `Muon_pt` over xrootd | 2.7 – 9.1 s (4 runs) | 15.5 MB compressed |
| Read the analysis columns over xrootd | 16.5 s, 12.2 s | 46.4 MB compressed, **~2 % of the file**; 180–240 kHz; peak memory 0.67 GB |
| Same columns via `uproot.iterate(step_size=500_000)` over xrootd | 18.3 s | peak memory **0.28 GB** |
| `xrdcp` the whole file to the cache | **29.9 s** | 2.35 GB at **78.8 MB/s**; adler32 matches the catalogue |
| Read the analysis columns from the cached copy | 9.1 s | 321 kHz. Best case: the file was probably still in page cache right after copying |
| DY NLO: analysis columns + weights + LHE (16 columns) over xrootd | 15.7 s | 29.3 MB compressed, 123 kHz, peak memory 0.69 GB. `Muon_pt` alone took 14.9 s (1 run) |

Parallel streaming test: several threads in one process, reading the analysis columns from Run2016G files over xrootd.

| Workers = files | Wall time | Events | Throughput |
|---:|---:|---:|---|
| 1 | 12–17 s | 2.94 M | 180–240 kHz |
| 2 | 12.9 s | 4.03 M | 312 kHz |
| 4 | 25.0 s | 10.0 M | 400 kHz |
| 8 | no result after 6 min 20 s | – | killed by the out-of-memory killer (2 GiB per-user limit) |

**Cached file:** `/data/atlas/users/nterlind/BND-school-cache/SingleMuon_Run2016G/0A4230E2-0C75-604D-890F-A4CE5E5C164E.root` (2,351,733,902 bytes)

## 5. Recommendation

**Stream over xrootd. Don't download the full 448 GB.**

- The Z→μμ selection only needs 11–16 columns, about 2 % of each file's bytes. Streaming the analysis columns of one data file (12–17 s) is faster than copying the file (30 s) and then reading it locally (9 s).
- Estimated time for one full pass on this node with 4 workers: about 14 min for the data (324 M events at 400 kHz), plus about 10 min for the MC (154 M events). The MC number is not measured: I assumed ~250 kHz, and a single MC file read at 123 kHz. A single-worker pass over the data alone takes about 27 min.
- On this login node, use `uproot.iterate` (step size ≤ 500k), with at most 4 workers or 2 whole-file readers. Otherwise you will hit the 2 GiB memory limit. Other processes you run on the node (VS Code server etc.) count against the same limit.
- For batch jobs, move to `stbc-i1/2/3` and check HTCondor there. One job per file, each streaming its own file, avoids the login node's limits.
- Suggestion: skim once (trigger + ≥ 2 muons, keeping only the needed columns) into the cache directory, then iterate on the skim. Also add retries on xrootd errors: I saw none, but I only read a handful of files.

**Download estimate for everything:** 448 GB at the measured 78.8 MB/s (single `xrdcp` stream) ≈ **95 min**. The data alone (260 GB) ≈ 55 min. Parallel `xrdcp` was not tested.

## 6. What failed or was surprising

1. **Wrong node type.** The hostname is `lvp036login.nikhef.nl`, not an stbc interactive node, and there is no HTCondor here.
2. **2 GiB memory limit per user.** The 8-worker test printed nothing. The kernel log shows one out-of-memory kill of a `python3` process (1.65 GB resident memory) in my user slice. The kill was triggered by another of my processes (`wandb-core`), which shares the same limit.
3. **Home quota is 1.2 GB** and 793 MB is already used. `~/.local/bin` is on `PATH`, so `pip install --user` would fill it quickly.
4. **`/project/atlas` has a 150 GB soft quota** even though the filesystem has 4.1 TB free, so it is not the place for bulk data. `/data/atlas` has the large quota, but that filesystem is 97 % full.
5. `/dcache` is not mounted on this node.
6. **Timings vary a lot.** `Muon_pt` alone on the same file took anywhere from 2.7 to 9.1 s, and 14.9 s on the MC file. Each number comes from one or a few runs on a shared node (load average 7–9 at the time), so treat them as order of magnitude.
7. The `HLT_*Mu*` branch lists differ: 14 branches exist only in the MC file, and 1 only in data (`HLT_L2Mu35_NoVertex_3Sta_NoBPTX3BX`). `HLT_IsoMu24` and `HLT_IsoTkMu24` are in both. Full lists: `filelists/hlt_mu_branches_SingleMuon_Run2016G.txt` and `filelists/hlt_mu_branches_DY_NLO.txt`.
8. `datasets/datasets.tsv` lists the four DY records three times, once per group. Deduplicate on recid when reading it.
9. `cernopendata-client` has no `--version` flag. Use `pip show cernopendata-client`.
10. Git had no `user.name`/`user.email` and no push credentials on this cluster.

## Appendix: all `HLT_*Mu*` branches in the SingleMuon Run2016G file (178)

<details>
<summary>show list</summary>

```
```

</details>
