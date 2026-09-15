# Handoff – Z → e⁺e⁻

## Team

- Luke
- 

## What we worked on

- Built a NanoAOD-based Z → e⁺e⁻ analysis using CMS 2016 (Legacy/UL) Open Data, from raw file download through a stacked, luminosity-normalized invariant-mass plot with a data/MC ratio panel.
- Added LHE-truth-level decay channel splitting (ee / μμ / ττ) for the inclusive DY MC, so the ττ feed-down into the ee selection is visible separately from the ee signal.
- Added diboson (WW/WZ/ZZ), ttbar, and W+jets backgrounds to the MC stack.
- Exported per-sample weighted histograms (`boost_histogram` → ROOT via `uproot`) and set up a TRExFitter config to fit the DY→ee signal-strength normalization (`mu_signal`) against data in a signal region spanning the full 60–120 GeV mass window.

## Data and simulation used

| Dataset | Type (data / MC) | Record / DOI | Notes |
|---|---|---|---|
| Run2016G | Data | recid 30529 | Electron primary dataset, UL2016 NanoAODv9 |
| Run2016H | Data | recid 30562 | Electron primary dataset, UL2016 NanoAODv9 |
| DY_inclusive_NLO | MC | recid 35669 | `/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/.../NANOAODSIM` — used as the DY sample in the current plot/fit |
| DY_inclusive_LO | MC | recid 35671 | `/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/.../NANOAODSIM` — downloaded, not currently used in the plot or fit |
| WW_2L2Nu | MC | recid 72676 | `WWTo2L2Nu` |
| WZ_3LNu | MC | recid 72752 | `WZTo3LNu` |
| WZ_2Q2L | MC | recid 72742 | `WZTo2Q2L` |
| ZZ_2L2Nu | MC | recid 75567 | `ZZTo2L2Nu` |
| ZZ_2Q2L | MC | recid 75573 | `ZZTo2Q2L` |
| ZZ_4L | MC | recid 75589 | `ZZTo4L` |
| TTTo2L2Nu | MC | recid 67801 | powheg ttbar, dilepton |
| TTToSemiLeptonic | MC | recid 67993 | powheg ttbar, semileptonic |
| TTToHadronic | MC | recid 67841 | powheg ttbar, all-hadronic |
| WJetsToLNu | MC | recid 69745 | amcatnloFXFX inclusive |

Certification JSON: `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt` (from [opendata.cern.ch/record/14220](https://opendata.cern.ch/record/14220)), restricted to the Run2016G+H run range (278820–284044) for the luminosity calculation.

## Selection and method

- **Trigger:** `HLT_Ele27_WPTight_Gsf`
- **Object selection:** electrons with `pt > 20 GeV`, `|eta| < 2.5`, `Electron_cutBased >= 3` (Medium working point); exactly two such electrons per event, opposite sign.
- **Event selection / mass window:** invariant mass histogrammed over 60–120 GeV (121 bins, 0.5 GeV wide); data events additionally filtered against the certified-runs JSON via a `(run, luminosityBlock)` lumi mask.
- **Corrections / weights:** `genWeight × L1PreFiringWeight_Nom × (xsec × lumi / sum_genweight)`. `sum_genweight` is the full per-sample sum of `genEventSumw` from the `Runs` tree (not restricted to selected events). Cross-section used for **both** DY_inclusive_NLO and DY_inclusive_LO: **6077.22 pb** (official NNLO DY M-50 value, applied identically to both samples so LO/NLO generator differences are absorbed via each sample's own `sum_genweight` rather than via separate k-factors). Diboson/ttbar/W+jets cross-sections are commonly-cited 13 TeV literature values (see notebook `XSEC_PB` dict) — **not yet individually verified against each record's own quoted cross-section**, same caveat as the DY value. Integrated luminosity: **16290.713420 pb⁻¹** (Run2016G+H, from `brilcalc lumi -c web --begin 278820 --end 284044`).
- **Channel splitting:** DY MC events are labeled `ee` / `mumu` / `tautau` / `other` using LHE-level truth (`LHEPart_pdgId`, `LHEPart_status`, outgoing leptons only), so the plotted DY stack separates the ee signal from ττ feed-down (μμ is essentially zero given the electron selection, as expected). Non-DY backgrounds (diboson, ttbar, W+jets) are each labeled by process group rather than LHE-classified, since the LHE-channel split is DY-specific.

## Histogram export

- Each MC component (DY→ee, DY→ττ, Diboson, ttbar, W+jets) and the data are written as weighted `boost_histogram` objects to individual ROOT files under `../datasets/z-ee/` (one file per sample, histogram named `h_mass`), for use as TRExFitter inputs.
- Data histogram is filled unweighted (Poisson-error convention); MC histograms use the full per-event weight so their `Weight()` storage carries the correct statistical uncertainty (`sumw2`) for the fit.

## Fit setup

- **Tool:** TRExFitter, config `fit.config`.
- **Fit type:** SPLUSB (signal+background), fit region `SR` = full 60–120 GeV window (rebinned by 2, so ~0.6M events in the peak bin after rebinning).
- **POI:** `mu_signal`, a floating normalization factor applied **only to the `DY_ee` sample**.
- **Backgrounds (fixed normalization, not floated):** `DY_tautau`, `Diboson`, `ttbar`, `Wjets`.
- **No systematics defined yet** — the `% - SYSTEMATICS - %` block in the config is empty.

## How to run

1. Open and run the analysis notebook top-to-bottom.
Downloads files via cernopendata-client, computes sum_genweight,
processes NanoAOD in chunks with uproot, produces the mass histogram,
the data/MC ratio plot, and exports per-sample ROOT histograms.

jupyter notebook z-ee.ipynb

2. Run the TRExFitter fit (histogram-only, no ntuple production step needed
since ReadFrom: HIST reads directly from the exported ROOT files).

trex-fitter hwdfp fit.config

## Results so far

- Full Run2016G+H dataset: ~6.56M raw selected data events (3,108,273 + 3,452,445) passing trigger + two-medium-electron + opposite-sign selection, before any mass-window cut.
- Z-peak invariant mass distribution built and plotted (log-scale y-axis), with MC stacked by process (ttbar, W+jets, Diboson, DY→ττ, DY→ee) and overlaid with data points including Poisson error bars.
- Added a data/MC ratio panel with a hatched gray band showing MC statistical uncertainty (from `sum(weight_i²)` per bin, since MC events are weighted) and Poisson error bars on the ratio points.
- Peak bin (~91 GeV) reaches roughly 500,000 events per 0.5 GeV bin — consistent in order of magnitude with the ~6.5M total selected data events concentrated in a narrow window around the Z mass.
- ttbar and W+jets yields in the selected sample are visibly very small relative to DY and diboson (e.g. `TTToHadronic` and `WJetsToLNu` each returned well under 300 selected events out of tens of millions processed) — expected given the tight two-medium-electron opposite-sign requirement, but worth a sanity check that these aren't being over- or under-selected due to a selection/branch issue specific to those samples.
- TRExFitter config written. Current best fit gives cross-section of 1848.6 pb.

## Open issues / next steps

- **Cross-section verification** — diboson/ttbar/W+jets cross-sections need to be checked against each record's own quoted value (same open caveat as the shared DY NNLO xsec).
- **No pileup reweighting applied yet** (only `L1PreFiringWeight_Nom` is used).
- **No electron ID/reco/trigger scale factors applied** — a second real normalization gap alongside pileup reweighting.
- **No systematics in the fit config** — statistical-only fit for now; will need at minimum a luminosity uncertainty and MC-statistical uncertainty treatment before any final number is quotable.
- **σ(Z→ee) extraction**: once the fit is run, remember the assumed DY cross-section (6077.22 pb) is the **inclusive DY→ℓℓ** value (summed over e/μ/τ), so `μ_signal × 6077.22` gives the inclusive measured cross-section, not σ(Z→ee) — divide by 3 to get the ee-channel cross-section.
- Consider adding a similar validation step (compare raw LHE-level channel fractions pre-selection) to confirm the ee/μμ/ττ classifier is correct, rather than relying only on post-selection sanity checks.