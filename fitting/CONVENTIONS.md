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
  `tautau_SR0`, `tautau_SR1`, `tautau_SR2` (the ττ BDT categories), ... (MultiFit does not merge
  regions with the same name, it warns);
  ττ additionally uses the samples `DYlowmass` (Z/γ*→ℓℓ, 10 < m < 50 GeV) and `DYtautau_nonfid`
  (Z/γ*→ττ outside the fiducial volume, mostly m > 120 GeV: a theory-normalised background, *not*
  scaled by `mu_Z`);
- sample names: `Data`, `DYmumu`, `DYee`, `DYtautau`, `TTbar`, `SingleTop`, `WW`, `WZ`, `ZZ`,
  `Fakes` (data-driven; carries its own systematics and its own Sumw2);
- one-sided variations provide only `__<syst>Up` (declare `Symmetrisation: ONESIDED`);
- the histograms are already in events: **never set `Lumi:` in the Job block**.

## 2. Parameter of interest and normalisation

- POI name `mu_Z` in every channel: a `NormFactor` (nominal 1) on the channel's Drell-Yan
  samples (in ττ on the fiducial part `DYtautau` only). Signal templates are normalised to the same prediction in every channel,
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
| `SigModel_tautau` | HISTO (normalisation only, one-sided) on `DYtautau`: fiducial C with LO madgraph instead of aMC@NLO | **no** (the μμ `SigModel` is a powheg-vs-aMC@NLO shape: a different quantity, must not be correlated) | Signal modelling |
| `MuonID`, `MuonIso`, `MuonTrigger`, `MuonReco` (OVERALL 0.4 % per muon, correlated: 0.8 % per event in μμ), `MuonScale`, `MuonRes` | HISTO/OVERALL | yes between channels that use muons (μμ, μτ_h) | Muon efficiency / Muon momentum |
| `ElectronID`, `ElectronReco`, `ElectronTrigger`, `ElectronScale` | | yes between channels that use electrons | Electron efficiency / Electron energy |
| `TauID_DM0/1/10/11`, `TauTrigger_DM*`, `TauES_DM*`, `TauFakeEle`, `TauFakeMu` | HISTO | ττ only (would be shared with eτh/μτh) | Tau (SubCategory Tau ID / Tau trigger / Tau energy scale) |
| `MET_Unclustered` | HISTO (shape) | channels that use MET | MET |
| `XS_DYll` 5 % (Z→ee/μμ as a background in ττ), `XS_DYlowmass` 10 %, `XS_WJets` 10 %, `XS_DYtautau_nonfid` 5 % | OVERALL | yes | Background normalisation |
| `MCStatNorm_<sample>_<channel>` | OVERALL (statistical normalisation of a smoothed template) | no | Background normalisation |
| `FakeStat_<channel>`, `FakeMethod_<channel>` | HISTO on `Fakes` | no (channel suffix) | Fakes |
| ττ: `FakeOSSS_tautau_c<k>`, `FakeClosure_tautau_c<k>_lo`, `FakeClosure_tautau_c<k>_hi` (one per BDT category k, the closure ones also per m_tt region below / above 110 GeV) | HISTO on `Fakes` (the FF statistics are in the `Fakes` Sumw2, i.e. in the γ parameters) | no | Fakes |
| MC statistics | per-bin gammas (`MCstatThreshold: 0`) | no | Gammas |

Use the same `Category` strings, so the grouped-impact tables of the individual fits and of the
combination have the same rows.

## 4. Outputs every channel provides

- `<channel>/fit/fitinputs/<channel>.root` (path in the channel's handoff.md; ττ commits it, 0.8 MB, so the
  combination works from any checkout — `git add -f`),
- `<channel>/fit/<channel>.config` (committed),
- `<channel>/fit/results/<job>/RooStats/<job>_combined_<job>_model.root` (workspace; ττ commits it. Note that
  TRExFitter names it `<job>_allBinsFitRegions_combined_<job>_model.root` when a region uses `DropBins`; ττ
  keeps a copy under the conventional name),
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

## 6. Decisions taken by the combination (15 Sep 2026)

§2 and §5 left two conventions "to be agreed" and the channels have since diverged in one place.
The combination (`combination/`, μμ ⊕ ττ; z-ee not yet available) settled them as follows.

- **Acceptance denominator: 60 < m < 120 GeV, NLO (aMC@NLO).** The recommendation of §2, and what
  both finished channels already quote. The μμ fiducial volume is defined with *dressed* leptons
  (ΔR < 0.1), the ττ one with LHE (Born-like) mass and `GenVisTau` momenta; the two definitions
  are not mixed because each channel's A is used only to convert its own σ_fid.
- **σ^pred(60–120) is not the same number in the two channels.** μμ builds it as
  σ_fid^pred / A(60–120) = 799.566 / 0.409209 = **1953.9 pb**; ττ uses the LHE-ττ subset of the
  same sample with 60 < m_LHE < 120 GeV, **1944.9 pb**. The two differ by **0.47 %** — a
  definitional difference, not a statistical one. Until the channels agree on one construction:
  - the combination multiplies each μ̂ by *its own* reference and combines the **cross sections**,
    never `mu_Z`;
  - a TRExFitter MultiFit with one shared `mu_Z` (§5) would fit one parameter against two
    references. At the current precision the bias is ~0.002 % on the combination, but a channel
    adding a third measurement of comparable weight should fix this first.
- **`SigModel` must be decorrelated between channels.** Both channels use the name, but μμ
  compares powheg with aMC@NLO (0.2 % on C) and ττ compares madgraph LO with aMC@NLO (7.3 %).
  §3 correlates by name, so a MultiFit needs `DecorrSysts: "SigModel"` + `DecorrSuff: "_<channel>"`
  in each Job block before the workspace is built. Every other name in §3 stays correlated.
  *Since z-tautau v2 (15 Sep 2026) this is moot: the ττ parameter is `SigModel_tautau` and is not in the
  workspace by default; the ττ signal is the fiducial `DYtautau` only, the fit has three regions
  `tautau_SR0/1/2` (SR0 above 110 GeV only), the ττ working point is DeepTau Tight (v3), and the ττ reference
  σ^pred(60–120) = 1944.9 pb is unchanged.*
- **Data statistics are uncorrelated between channels** — `SingleMuon` and `Tau` are disjoint
  primary datasets selected by orthogonal triggers. Worth stating because a future eτh/μτh channel
  would *not* be orthogonal to μμ and would need an overlap treatment.
- **Z/γ*→ττ as a background in the μμ signal region** carries its own `XS_DYtautau` (5 %) rather
  than scaling with `mu_Z`, although under lepton universality it is the same process. 11.1k of
  10.4M events: 0.005 % on the combination, documented rather than fixed.

The rationale for each is in `combination/docs/02-correlation-model.md`; the resulting numbers in
`combination/result.md`.
