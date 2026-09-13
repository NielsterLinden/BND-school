# Handoff – Z → e⁺e⁻

## Team

- Luke
- 

## What we worked on

- Built a NanoAOD-based Z → e⁺e⁻ analysis using CMS 2016 (Legacy/UL) Open Data, from raw file download through a stacked, luminosity-normalized invariant-mass plot with a data/MC ratio panel.
- Added LHE-truth-level decay channel splitting (ee / μμ / ττ) for the inclusive DY MC, so the ττ feed-down into the ee selection is visible separately from the ee signal.

## Data and simulation used

| Dataset | Type (data / MC) | Record / DOI | Notes |
|---|---|---|---|
| Run2016G | Data | recid 30529 | Electron primary dataset, UL2016 NanoAODv9 |
| Run2016H | Data | recid 30562 | Electron primary dataset, UL2016 NanoAODv9 |
| DY_inclusive_NLO | MC | recid 35669 | `/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| DY_inclusive_LO | MC | recid 35671 | `/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

Certification JSON: `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt` (from [opendata.cern.ch/record/14220](https://opendata.cern.ch/record/14220)), restricted to the Run2016G+H run range (278820–284044) for the luminosity calculation.

## Selection and method

- **Trigger:** `HLT_Ele27_WPTight_Gsf`
- **Object selection:** electrons with `pt > 20 GeV`, `|eta| < 2.5`, `Electron_cutBased >= 3` (Medium working point); exactly two such electrons per event, opposite sign.
- **Event selection / mass window:** invariant mass histogrammed over 60–120 GeV (121 bins, 0.5 GeV wide); data events additionally filtered against the certified-runs JSON via a `(run, luminosityBlock)` lumi mask.
- **Corrections / weights:** `genWeight × L1PreFiringWeight_Nom × (xsec × lumi / sum_genweight)`. `sum_genweight` is the full per-sample sum of `genEventSumw` from the `Runs` tree (not restricted to selected events). Cross-section used for **both** DY_inclusive_NLO and DY_inclusive_LO: **6077.22 pb** (official NNLO DY M-50 value, applied identically to both samples so LO/NLO generator differences are absorbed via each sample's own `sum_genweight` rather than via separate k-factors). Integrated luminosity: **16290.713420 pb⁻¹** (Run2016G+H, from `brilcalc lumi -c web --begin 278820 --end 284044`).
- **Channel splitting:** DY MC events are labeled `ee` / `mumu` / `tautau` / `other` using LHE-level truth (`LHEPart_pdgId`, `LHEPart_status`, outgoing leptons only), so the plotted DY stack separates the ee signal from ττ feed-down (μμ is essentially zero given the electron selection, as expected).

## How to run

```
# 1. Open and run the analysis notebook top-to-bottom.
#    Downloads files via cernopendata-client, computes sum_genweight,
#    processes NanoAOD in chunks with uproot, produces the mass histogram
#    and data/MC ratio plot.
jupyter notebook z-ee.ipynb
```

## Results so far

- Full Run2016G+H dataset: ~6.56M raw selected data events (3,108,273 + 3,452,445) passing trigger + two-medium-electron + opposite-sign selection, before any mass-window cut.
- Z-peak invariant mass distribution built and plotted (log-scale y-axis), with MC stacked by true decay channel (ee dominant, small ττ tail, μμ ≈ 0) and overlaid with data points including Poisson error bars.
- Added a data/MC ratio panel with a hatched gray band showing MC statistical uncertainty (from `sum(weight_i²)` per bin, since MC events are weighted) and Poisson error bars on the ratio points.
- Peak bin (~91 GeV) reaches roughly 500,000 events per 0.5 GeV bin — consistent in order of magnitude with the ~6.5M total selected data events concentrated in a narrow window around the Z mass.

## Open issues / next steps

- No non-DY backgrounds included yet (e.g. tt̄, diboson, single top) — current MC stack is DY-only.
- The shared NNLO cross-section (6077.22 pb) for both LO and NLO DY samples should be double-checked against the current version of the CMS `StandardModelCrossSectionsat13TeV` TWiki before it's used in any final result; it was sourced from a downstream analysis repo rather than the TWiki directly.
- No pileup reweighting applied yet (only `L1PreFiringWeight_Nom` is used).
- Consider adding a similar validation step (compare raw LHE-level channel fractions pre-selection) to confirm the ee/μμ/ττ classifier is correct, rather than relying only on post-selection sanity checks.