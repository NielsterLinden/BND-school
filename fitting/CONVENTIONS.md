# Fit-input conventions shared by the three channels

Every channel (z-ee, z-mumu, z-tautau) produces its fit inputs in the same form, so the
combination is a TRExFitter *MultiFit* over three independent workspaces and nothing has to be
re-edited. TRExFitter **v1.8.0** (the git submodule `TRExFitter-v1.8.0/`, build with
`fitting/build_trexfitter.sh`) is used everywhere.

## 1. Files and histogram names (`ReadFrom: HIST`)

One ROOT file per channel, written with uproot (no ROOT import in analysis code):

```
<channel-dir>/fit/fitinputs/<channel>.root        e.g. z-mumu/fit/fitinputs/zmumu.root
    <region>__<sample>                            nominal, e.g. mumu_SR__DYmumu
    <region>__Data                                data
    <region>__<sample>__<syst>Up                  +1 sigma variation (shape and/or normalisation)
    <region>__<sample>__<syst>Down                -1 sigma variation
```

TRExFitter builds the histogram path as
`<HistoPath>/<HistoFile>.root/<HistoName><HistoNameSuff(region)><HistoNameSuff(sample)><HistoNameSufUp/Down(syst)>`,
which the config in `fitting/trexconfig.py` maps onto exactly these names. Rules:

- every histogram is a `TH1D` with the sum of weights squared stored (`hist.storage.Weight()`
  before writing; `fitting/trexhist.py` refuses histograms without it);
- all samples of a region share one binning; no dots (`.`) in any name; no negative bins in the
  nominal templates (`trexhist.write` clips and reports);
- region names are channel-prefixed and globally unique: `mumu_SR`, `mumu_CRemu`, `ee_SR`,
  `tautau_SR`, ... (MultiFit does not merge regions with the same name, it warns);
  ττ additionally uses the sample `DYlowmass` (Z/γ*→ℓℓ, 10 < m < 50 GeV);
- sample names: `Data`, `DYmumu`, `DYee`, `DYtautau`, `TTbar`, `SingleTop`, `WW`, `WZ`, `ZZ`,
  `Fakes` (data-driven; carries its own systematics and its own Sumw2);
- one-sided variations provide only `__<syst>Up` (declare `Symmetrisation: ONESIDED`);
- the histograms are already in events: **never set `Lumi:` in the Job block**.

## 2. Parameter of interest and normalisation

- POI name `mu_Z` in every channel: a `NormFactor` (nominal 1) on the channel's Drell-Yan
  samples. Signal templates are normalised to the same prediction in every channel,
  σ(Z/γ*→ℓℓ, m > 50 GeV) = 6077.22 / 3 pb per flavour, times L = **16393.381 pb⁻¹**
  (normtag; CMS Open Data record 1059), divided by the sample's sum of generator weights.
- The channel cross section is σ_fid = μ̂ × σ_fid^pred, and σ(60 < m_Born < 120) = σ_fid / A
  with A from the same DY sample (dressed leptons, ΔR < 0.1). Quote both; the combination must
  agree on the denominator (recommended: 60 < m < 120 GeV) and on NLO (aMC@NLO) for A.

## 3. Systematic names (identical strings = correlated nuisance parameters)

| name | type | correlated across channels | Category |
|---|---|---|---|
| `Lumi` | OVERALL ±1.2 % on all MC | yes | Luminosity |
| `Pileup` | HISTO | yes | Pileup |
| `L1Prefiring` | HISTO | yes | L1 prefiring |
| `XS_TTbar` 6 %, `XS_SingleTop` 10 %, `XS_WW` 10 %, `XS_WZ` 10 %, `XS_ZZ` 10 %, `XS_DYtautau` 5 % | OVERALL | yes | Background normalisation |
| `PDF`, `QCDScale`, `PS_ISR`, `PS_FSR`, `SigModel` | HISTO on the DY samples (vary the C factor only: renormalised to a constant fiducial yield) | yes | Signal modelling |
| `MuonID`, `MuonIso`, `MuonTrigger`, `MuonReco` (OVERALL 0.4 %), `MuonScale`, `MuonRes` | HISTO/OVERALL | yes between channels that use muons (μμ, μτ_h) | Muon efficiency / Muon momentum |
| `ElectronID`, `ElectronReco`, `ElectronTrigger`, `ElectronScale` | | yes between channels that use electrons | Electron efficiency / Electron energy |
| `TauID_DM0/1/10/11`, `TauTrigger_DM*`, `TauES_DM*`, `TauFakeEle`, `TauFakeMu` | HISTO | ττ only (would be shared with eτh/μτh) | Tau (SubCategory Tau ID / Tau trigger / Tau energy scale) |
| `MET_Unclustered` | HISTO (shape) | channels that use MET | MET |
| `XS_DYll` 5 % (Z→ee/μμ as a background in ττ), `XS_DYlowmass` 10 %, `XS_WJets` 10 % | OVERALL | yes | Background normalisation |
| `MCStatNorm_<sample>_<channel>` | OVERALL (statistical normalisation of a smoothed template) | no | Background normalisation |
| `FakeStat_<channel>`, `FakeMethod_<channel>` | HISTO on `Fakes` | no (channel suffix) | Fakes |
| ττ: `FakeStat_tautau_DM*`, `FakeOSSS_tautau`, `FakeClosure_tautau` | HISTO on `Fakes` | no | Fakes |
| MC statistics | per-bin gammas (`MCstatThreshold: 0`) | no | Gammas |

Use the same `Category` strings, so the grouped-impact tables of the individual fits and of the
combination have the same rows.

## 4. Outputs every channel provides

- `<channel>/fit/fitinputs/<channel>.root` (not committed; path in the channel's handoff.md),
- `<channel>/fit/<channel>.config` (committed),
- `<channel>/fit/results/<job>/RooStats/<job>_combined_<job>_model.root` (workspace, not committed),
- `<channel>/fit/results/<channel>_fit_result.json` written by `fitting/run_trex.py`
  (μ, uncertainties, grouped impacts, σ_fid, σ(60–120), n_obs, n_bkg, A, C, L),
- the same numbers in the channel's `handoff.md`.

## 5. Combination

`fitting/combination_skeleton.config` is a MultiFit with `Combine: TRUE` and one `Fit:` block per
channel pointing at the channel's config and results directory. Run from `combination/`:

```
trex-fitter mwf comb.config      # combined workspace + fit
trex-fitter mdp comb.config      # plots
trex-fitter mr  comb.config      # ranking
trex-fitter mi  comb.config      # grouped impacts
```

Stat-only: run every channel with `"StatOnly=TRUE:Suffix=_statOnly"` first, then the same option
string on the MultiFit.
