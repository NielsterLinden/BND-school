# CLAUDE.md — BND-school Z cross-section project (repository level)

Read this first, then the `CLAUDE.md` of the folder you work in (`z-mumu/`, `z-tautau/`, `combination/`,
`presentation/`). Three channel measurements (Z→ee, Z→μμ, Z→ττ) plus their combination, on CMS 2016 Open Data
(NanoAODv9, Run2016G+H). Human entry point: the root `README.md`.

> **The project is finished (18 Sep 2026).** The analyses were frozen on 17 Sep 2026, the talk was given on 18 Sep
> with `2026.09.18 -- Presentation BND -- OpenData.pptx` (repository root, Git LFS). No reruns and no improved
> analyses: channel and combination outputs are final; cite them, never regenerate them. Every final number with its
> source file: `docs/FREEZE.md`.

**Final result:** σ(Z/γ*→ℓℓ, 60–120 GeV) = **1949 ⁺³⁰₋₃₀ pb** per flavour (compatibility p = 0.38), one TRExFitter
v1.8.0 MultiFit of ee ⊕ μμ ⊕ ττ in `combination/combLieke/` (`README.md`, `output/result.json`, `output/plots/`).
ττ is z-tautau's four-channel measurement (τhτh + μτh + eτh + eμ), 1981 ⁺¹⁵³₋₁₃₆ pb, which carries its τh ID p_T
dependence (`TauIDpT_tautau`; ⁺⁷³₋₇₀ without). μμ: 1931 ± 30 pb. `mumu_CRemu` is dropped in the combination (same data
as the ττ eμ channel). The earlier three-channel freeze with τhτh only, 1945 ⁺³¹₋₃₀ pb, is the git tag
`zmumu-freeze-2026-09-17`; `docs/FREEZE.md` is its record. The combination fits ran on HTCondor
(`python run.py prepare && python run.py condor --submit`, then `results`, `plots`).

> **Known limitation of the ee channel:** z-ee (delivery of 16 Sep 14:44) does not veto the ECAL barrel–endcap gap.
> The official electron-ID SF map is a placeholder (sf = 1 ± 1) there, which makes the ee `ElectronID` uncertainty
> ±5.9 % instead of ~1.2 %; its shape then over-constrains the channel (saturated GoF p = 4e-42). No
> `HLT_Ele27_WPTight_Gsf` scale factor is applied either. The combination splits ee shape from normalisation
> parameters to stay stable: "The ee channel" in `combination/combLieke/README.md`.

## Where things are

| | |
|---|---|
| root | only `README.md`, the talk (`.pptx`) and folders. Keep it that way. |
| `z-ee/`, `z-mumu/`, `z-tautau/` | one channel each: `README.md`, `CLAUDE.md`, `handoff.md` (the hand-in to the combination), `docs/NN-*.md`, code, `fit/`, `output/`, `review/` |
| `combination/combLieke/` | the MultiFit (`run.py`, `mf/`, `config/`, `output/`, `docs/`, `checks/`); agent guide `combination/CLAUDE.md` |
| `fitting/` | shared python package (`trexhist`, `trexconfig`, `run_trex`, `uncertainty_parity`), `setup.sh`, `build_trexfitter.sh`, `selftest.py`. Imported as `from fitting import …` with the repository root on `sys.path`: do not move it. |
| `docs/` | shared documents: `CONVENTIONS.md` (the fit-input contract), `FREEZE.md`, `UNCERTAINTY_PARITY.md`, `OpenDataIntro.pdf` |
| `datasets/` | `datasets.tsv`, `SIZES.md`, `GRL/GRL.txt`, `corrections/`, the z-ee histograms |
| `presentation/` | the Manim clips of the talk: `README.md`, `CLAUDE.md`, `docs/01…07`, `docs/briefs/`, `clips/`, `scenes/`, `style/`, `data/`, `tools/` |
| `.claude/` | this guide, `agents/` (the animation agents and `uncertainty-parity-auditor`), `prompts/` (the prompts the analyses and the chapter animations started from) |

Documentation pattern, everywhere: a folder's `README.md` is the way in, details are numbered pages in its `docs/`,
a `CLAUDE.md` is the agent guide. New prose goes into a `docs/`, not next to the code and not into the root.

## Environment (one line, every shell)

```bash
source fitting/setup.sh    # LCG_110: ROOT 6.40, python 3.13, uproot 5.7, awkward 2.9, hist, mplhep, iminuit
                           # + trex-fitter on PATH (from $TREXFITTER_HOME) + $BND_ROOT + $BND_SKIM_DIR
```

Never `pip install` into the system or user environment. Analysis code is python + uproot; do **not** import ROOT in
analysis processes (uproot writes the TH1D inputs; TRExFitter is a binary). The Manim environment of `presentation/`
is separate (`presentation/CLAUDE.md`); never mix the two in one shell.

## TRExFitter (mandatory version 1.8.0)

TRExFitter is **not in this repository** (it was a submodule until 18 Sep 2026). `bash fitting/build_trexfitter.sh`
clones tag v1.8.0 into `$TREXFITTER_HOME` (default `../TRExFitter-v1.8.0`, next to the checkout) and builds it with
C++20 to match the LCG ROOT; it patches the single C++20 incompatibility of v1.8.0 (`u8string()` in `Root/Common.cc`)
idempotently. Fallback route: `setupATLAS; asetup StatAnalysis,0.7.2`. Self-test of the whole chain:
`python fitting/selftest.py`. Every channel must use v1.8.0: the combination relies on it.

## Conventions shared by the three channels (`docs/CONVENTIONS.md` is the contract)

- Fit inputs: one ROOT file per channel `<channel>/fit/fitinputs/<channel>.root` with `TH1D` (Sumw2 stored) named
  `<region>__<sample>` and `<region>__<sample>__<syst>Up|Down`; regions channel-prefixed (`mumu_SR`, `ee_SR`,
  `tautau_SR`); samples `Data, DYmumu, DYee, DYtautau, TTbar, SingleTop, WW, WZ, ZZ, Fakes`; POI `mu_Z`; correlated
  NP names identical across channels (`Lumi`, `Pileup`, `L1Prefiring`, `XS_*`, `PDF`, `QCDScale`, `PS_*`, `SigModel`,
  `Muon*`, `Electron*`, `Tau*`), uncorrelated ones suffixed by channel (`FakeStat_mumu`). The combination's
  settled conventions are `docs/CONVENTIONS.md` §6 (its references to `combination/docs/` point to the removed
  covariance/BLUE combination, which is in the git history).
- **Luminosity: 16393.381 pb⁻¹** (normtag, CMS Open Data record 1059) for every channel; the earlier
  16290.713 pb⁻¹ (brilcalc without `--normtag`) is 0.63 % low.
- Signal normalisation: σ(Z/γ*→ℓℓ, m > 50) = 6077.22/3 pb per flavour; the combination quotes the
  60 < m < 120 GeV Born-level cross section (NLO acceptance).

## Data

- Catalogue: `datasets/datasets.tsv`, `datasets/SIZES.md`; golden JSON `datasets/GRL/GRL.txt` (full 2016 — restrict
  to runs 278820–284044).
- dCache (`/dcache/atlas/sjankovy/BND/`): all six primary datasets (G+H) and 120 MC samples; **no ttbar / single
  top** there (stream from EOS: `root://eospublic.cern.ch//eos/opendata/...`, file lists via the Open Data API — see
  `z-mumu/zmumu/samples.py`).
- v2 skims (all channels may use them): `$BND_SKIM_DIR` = `/data/atlas/users/nterlind/BND-school-cache/skims_v2`
  (NanoAOD branch names, `z-mumu/docs/10-skims.md`).
- Cluster: run on `stbc-i*` (dCache mounted); ≤ 8 parallel dCache readers; every per-file job falls back to EOS.

## Git rules

Work on `main`, `git pull --rebase` before pushing. Never commit data files; the committed ROOT files are the
deliberate exceptions (the channels' fit inputs and workspaces the combination hashes in `output/result.json`, the
z-ee histograms). `*.pptx` goes through Git LFS (`.gitattributes`). Commit messages carry the attribution lines used
so far. History anchors: tags `zmumu-freeze-2026-09-17` (frozen three-channel state) and `pre-cleanup-2026-09-18`
(everything the final cleanup removed: TRExFitter submodule, superseded and unused clips, archived cuts, old scenes).
