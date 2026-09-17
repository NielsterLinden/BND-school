# CLAUDE.md — BND-school Z cross-section project (repository level)

Read this first, then the channel `CLAUDE.md` you are working on. Shared workspace of three
subgroups (Z→ee, Z→μμ, Z→ττ) plus a combination, on CMS 2016 Open Data (NanoAODv9, Run2016G+H).

> **Frozen on 17 Sep 2026 — read `FREEZE.md` first** (git tag `zmumu-freeze-2026-09-17`). No more reruns or
> improved analyses: the channel and combination outputs are final; cite them. The talk (`presentation/`) is
> built on them by the animation agents; `FREEZE.md` lists what there still shows pre-freeze numbers.

**State (17 Sep 2026):** all three channels have a result. `combination/combLieke/` combines them in
one **TRExFitter MultiFit** (σ(Z/γ*→ℓℓ, 60–120 GeV) = 1945 ⁺³¹₋₃₀ pb, channel compatibility p = 0.30 —
`combination/combLieke/README.md`, `output/result.json`, `output/plots/`). `combination/CLAUDE.md` is
the agent guide. The earlier covariance (BLUE) combination and its slide deck were removed; they are in
the git history. The conventions the combination settled are `fitting/CONVENTIONS.md` §6. Its
references to `combination/docs/` point to that removed version.

> **Open, limiting the ee channel:** z-ee (re-delivered 16 Sep 14:44, renormalised PDF/QCDScale) does not
> veto the ECAL barrel–endcap gap. The official electron-ID SF map is a placeholder (sf = 1 ± 1) there,
> which makes the ee `ElectronID` uncertainty ±5.9 % instead of ~1.2 %. Its shape then over-constrains
> the channel (saturated GoF p = 7e-43). No `HLT_Ele27_WPTight_Gsf` scale factor is applied either. The
> combination splits ee shape from normalisation parameters to stay stable; the list for z-ee is at the end
> of "The ee channel" in `combination/combLieke/README.md`. `ASK_Z_EE.md` is about the previous delivery.

## Environment (one line, every shell)

```bash
source setup.sh            # LCG_110: ROOT 6.40, python 3.13, uproot 5.7, awkward 2.9, hist, mplhep, iminuit
                           # + the TRExFitter v1.8.0 build on PATH (trex-fitter) + $BND_SKIM_DIR
```

Never `pip install` into the system or user environment. Analysis code is python + uproot; do
**not** import ROOT in analysis processes (uproot writes the TH1D inputs; TRExFitter is a binary).

## TRExFitter (mandatory version 1.8.0)

- `TRExFitter-v1.8.0/` is a **git submodule pinned at tag v1.8.0** (HTTPS URL). After cloning
  or pulling: `git submodule update --init --recursive`.
- Build once: `bash fitting/build_trexfitter.sh` (C++20 to match the LCG ROOT; it patches the
  single C++20 incompatibility of v1.8.0, `u8string()` in `Root/Common.cc`, idempotently — the
  submodule then shows that file as modified, which is expected). Fallback route documented in
  `TRExFitter-v1.8.0/setup.sh`: `setupATLAS; asetup StatAnalysis,0.7.2`.
- Self-test of the whole chain: `python fitting/selftest.py`.

## Conventions shared by the three channels (`fitting/CONVENTIONS.md` is the contract)

- Fit inputs: one ROOT file per channel `<channel>/fit/fitinputs/<channel>.root` with `TH1D`
  (Sumw2 stored) named `<region>__<sample>` and `<region>__<sample>__<syst>Up|Down`; regions
  channel-prefixed (`mumu_SR`, `ee_SR`, `tautau_SR`); samples `Data, DYmumu, DYee, DYtautau,
  TTbar, SingleTop, WW, WZ, ZZ, Fakes`; POI `mu_Z`; correlated NP names identical across
  channels (`Lumi`, `Pileup`, `L1Prefiring`, `XS_*`, `PDF`, `QCDScale`, `PS_*`, `SigModel`,
  `Muon*`, `Electron*`, `Tau*`), uncorrelated ones suffixed by channel (`FakeStat_mumu`).
- Helpers: `fitting/trexhist.py` (write/check inputs), `fitting/trexconfig.py` (config from
  python), `fitting/run_trex.py` (run + parse), `fitting/combination_skeleton.config` (MultiFit).
- **Luminosity: 16393.381 pb⁻¹** (normtag, CMS Open Data record 1059) for every channel; the
  earlier 16290.713 pb⁻¹ (brilcalc without `--normtag`) is 0.63% low.
- Signal normalisation: σ(Z/γ*→ℓℓ, m>50) = 6077.22/3 pb per flavour; acceptance convention
  for the combination to be agreed (recommended: 60 < m < 120 GeV Born-level denominator, NLO).

## Data

- Catalogue: `datasets/datasets.tsv`, `datasets/SIZES.md`; golden JSON `datasets/GRL/GRL.txt`
  (full 2016 — restrict to runs 278820–284044).
- dCache (`/dcache/atlas/sjankovy/BND/`): all six primary datasets (G+H) and 120 MC samples;
  **no ttbar / single top** there (stream from EOS: `root://eospublic.cern.ch//eos/opendata/...`,
  file lists via the Open Data API — see `z-mumu/zmumu/samples.py`).
- v2 skims (all channels may use them): `$BND_SKIM_DIR` = `/data/atlas/users/nterlind/BND-school-cache/skims_v2`
  (NanoAOD branch names, see `z-mumu/docs/10-skims.md`).
- Cluster: run on `stbc-i*` (dCache mounted); ≤ 8 parallel dCache readers; every per-file job
  falls back to EOS.

## Git rules (from README.md)

Work on `main`, `git pull --rebase` before pushing, never commit ROOT/data files, keep each
channel's `handoff.md` current. Commit messages carry the attribution lines used so far.
