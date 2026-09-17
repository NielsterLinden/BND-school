# Handoff – Z → e⁺e⁻

## Team

- Luke
- Noémie
- Mathias

## What we worked on

- Built a NanoAOD-based Z → e⁺e⁻ analysis using CMS 2016 (Legacy/UL) Open Data, from raw file download through a stacked, luminosity-normalized invariant-mass plot with a data/MC ratio panel.
- Added LHE-truth-level decay channel splitting (ee / μμ / ττ) for the inclusive DY MC, so the ττ feed-down into the ee selection is visible separately from the ee signal.
- Added diboson (WW/WZ/ZZ), ttbar, and W+jets backgrounds to the MC stack.
- Applied pileup reweighting, L1 pre-firing weights, and electron Reco/ID scale factors from the CMS jsonpog-integration corrections (correctionlib).
- Added per-event weight systematics (PU, L1Prefire, EleReco, EleID) and PDF/scale envelope histograms accumulated in-stream.
- Exported per-sample weighted histograms (`boost_histogram` → ROOT via `uproot`), including nominal and all systematic variations, for use as TRExFitter inputs.
- Set up a TRExFitter config to fit the DY→ee signal-strength normalization (`mu_signal`) against data in a signal region spanning the full 60–120 GeV mass window.

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

Corrections: CMS jsonpog-integration-2016post (`POG/LUM/2016postVFP_UL/puWeights.json.gz` and `POG/EGM/2016postVFP_UL/electron.json.gz`), loaded via `correctionlib`.

## Selection and method

- **Trigger:** `HLT_Ele27_WPTight_Gsf`
- **Object selection:** electrons with `pt > 20 GeV`, `|η| < 2.5`, `Electron_cutBased >= 3` (Medium working point); exactly two such electrons per event, opposite sign.
- **Event selection / mass window:** invariant mass histogrammed over 60–120 GeV (121 bins, 0.5 GeV wide); data events additionally filtered against the certified-runs JSON via a `(run, luminosityBlock)` lumi mask.
- **Corrections / weights:** `genWeight × L1PreFiringWeight_Nom × PU_weight × RecoSF(e₁) × RecoSF(e₂) × IDSF(e₁) × IDSF(e₂) × (xsec × lumi / sum_genweight)`. `sum_genweight` is the full per-sample sum of `genEventSumw` from the `Runs` tree (not restricted to selected events). Cross-section used for DY_inclusive_NLO: **6077.22 pb** (official NNLO DY M-50 value). Diboson/ttbar/W+jets cross-sections are commonly-cited 13 TeV literature values (see notebook `XSEC_PB` dict) — **not yet individually verified against each record's own quoted cross-section**. Integrated luminosity: **16393.381 pb⁻¹** (Run2016G+H, computed with `brilcalc lumi -c web --normtag` using the `normtag_PHYSICS` value; the earlier value of 16290.713420 pb⁻¹ was derived without `--normtag` and was ~0.63% low).
- **Scale factors:** PU weights from `Collisions16_UltraLegacy_goldenJSON` (central, up, down); electron reco and ID SFs from `UL-Electron-ID-SF` for year tag `2016postVFP`, using working points `RecoAbove20` (reco) and `Medium` (ID). Per-event SF is the product of the leading and subleading electron SFs.
- **Channel splitting:** DY MC events are labeled `ee` / `mumu` / `tautau` / `other` using LHE-level truth (`LHEPart_pdgId`, `LHEPart_status`, outgoing leptons only), so the plotted DY stack separates the ee signal from ττ feed-down (μμ is essentially zero given the electron selection, as expected). Non-DY backgrounds are each labeled by process group rather than LHE-classified.

## Histogram export

Each MC component (DY→ee, DY→ττ, Diboson, ttbar, W+jets) and the data are written as weighted `boost_histogram` objects to individual ROOT files under `../datasets/z-ee/` (one file per component, histogram named `h_mass`), for use as TRExFitter inputs.

Each MC file contains:
- `h_mass` — nominal weighted histogram
- `h_mass_{syst}` for each of: `PUUp`, `PUDown`, `L1PrefireUp`, `L1PrefireDown`, `EleRecoUp`, `EleRecoDown`, `EleIDUp`, `EleIDDown`

The `DY_ee.root` file additionally contains `h_mass_PDFUp`, `h_mass_PDFDown`, `h_mass_ScaleUp`, `h_mass_ScaleDown` (computed as max-deviation envelopes over the LHE PDF and scale weight replicas, accumulated in-stream). Note: these envelopes are currently computed over the full inclusive DY sample before the ee/μμ/ττ channel split, which is dominated by the ee channel post-selection and is a good approximation; an exact per-channel accumulation can be added if needed.

Data histogram is filled unweighted (Poisson-error convention); MC histograms use the full per-event weight so their `Weight()` storage carries the correct statistical uncertainty (`sumw2`) for the fit.

## Fit setup

- **Tool:** TRExFitter, config `fit.config`.
- **Fit type:** SPLUSB (signal+background), fit region `SR` = full 60–120 GeV window (rebinned by 2).
- **POI:** `mu_signal`, a floating normalization factor applied **only to the `DY_ee` sample**.
- **Backgrounds (fixed normalization, not floated):** `DY_tautau`, `Diboson`, `ttbar`, `Wjets`.
- **Systematics:** histogram inputs include weight systematics and PDF/scale envelopes; **the fit config has not yet been updated to include them** — it is currently a statistical-only fit.

## How to run

1. Open and run the analysis notebook top-to-bottom. Downloads files via `cernopendata-client`, computes `sum_genweight`, processes NanoAOD in chunks with `uproot`, produces the mass histogram, the data/MC ratio plot, and exports per-sample ROOT histograms with systematics.

```
jupyter notebook z-ee.ipynb
```

2. Run the TRExFitter fit (histogram-only, no ntuple production step needed since `ReadFrom: HIST` reads directly from the exported ROOT files).

```
trex-fitter hwdfp fit.config
```

## Results so far

- Full Run2016G+H dataset: ~6.56M raw selected data events (3,108,273 + 3,452,445) passing trigger + two-medium-electron + opposite-sign selection, before any mass-window cut.
- Z-peak invariant mass distribution built and plotted (log-scale y-axis), with MC stacked by process (ttbar, W+jets, Diboson, DY→ττ, DY→ee) and overlaid with data points including Poisson error bars.
- Data/MC ratio panel with hatched gray band for MC statistical uncertainty and Poisson error bars on the ratio.
- Peak bin (~91 GeV) reaches roughly 500,000 events per 0.5 GeV bin — consistent with ~6.5M total selected data events concentrated in a narrow window around the Z mass.
- ttbar and W+jets yields in the selected sample are very small relative to DY and diboson (e.g. `TTToHadronic` returned 19 selected events out of ~107M processed, `WJetsToLNu` returned 216) — expected given the tight two-medium-electron opposite-sign requirement, but worth a sanity check that these are not affected by a selection or branch issue.
- TRExFitter config written. Current best fit (with systematics) gives a cross-section of **1840.8 pb ± 29.9**.

## Open issues / next steps

