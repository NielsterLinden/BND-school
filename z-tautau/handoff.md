# Handoff – Z → τ⁺τ⁻: v4 (τhτh + μτh + eτh + eμ) on top of v3 (τhτh)

## v4 in short (17 September 2026)

- Four channels fitted together (`run_v4.py`, `docs/10-v4-plan.md`): the v3 τhτh categories plus μτh (SingleMuon,
  IsoMu24), eτh (SingleElectron, Ele27_WPTight; EOS only) and eμ (MuonEG cross triggers) with an eμ tt̄ control region.
  **No μμ channel**; second-muon / second-electron vetoes make every channel orthogonal to the z-mumu and z-ee
  selections (the eμ channel shares events with z-mumu's tt̄ *control* region `mumu_CRemu`, which a joint fit must drop).
- Signal: Z/γ*→ττ with 60 < m_LHE < 120 GeV in every decay (one μ_Z); the rest of the DY ττ simulation is a background.
- τh ID scale factors free per decay mode (TRExFitter NormFactors, products on the τhτh templates), τh energy scale with a
  3 % prior: both measured in situ as in CMS arXiv:1801.03535. Every source of the paper's Table 2 is implemented from an
  official correction, measured in situ, or replaced by a stated estimate (`docs/10-v4-plan.md` section 9).
- Lepton-channel fakes: per-process fake factors (multijet / W+jets / tt̄) with AR fractions, OS/SS and m_T corrections,
  same-sign validation (0.97 ± 0.01 μτh, 1.03 ± 0.02 eτh). Two things learned: the isolated-lepton same-sign region is ~50 %
  W+jets, and the W+jets fake factor is charge-correlated (OS 0.08, SS 0.04).
- Outputs: `output_v4/RESULTS.md`, `output_v4/results.json`, `output_v4/plots/`, `fit_v4/ztautau_v4.config`,
  `fit_v4/fitinputs/ztautau_v4.root`, `fit_v4/results/ztautau_v4_fit_result.json`, `external/trigger_insitu_v4.json`.
- Caches: `skims_v4/` (17 GB), `ntuples_v4/` (2.9 GB) under `/data/atlas/users/sjankovy/BND-school-cache/ztautau/`.

**Result (v4, four channels, τh ID scale factors and energy scale fitted in situ)**

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2057 +90 −86 pb   (stat ±9;  prediction 1945 pb: aMC@NLO acceptance, NNLO normalisation)
μ_Z = 1.058 +0.046 −0.044   (stat ±0.005, syst ±0.045),   goodness of fit p = 0.09,   μ_tt̄ = 1.16 ± 0.05
τh ID SF (Tight):  DM0 0.959 ± 0.043 (POG 0.90 ± 0.13)  DM1 0.934 ± 0.038 (POG 0.89 ± 0.05)  DM10 0.866 ± 0.038 (POG 0.94 ± 0.15)  DM11 0.768 ± 0.051 (POG 0.81 ± 0.15)
τh energy scale:   DM0 -0.7 ± 0.9 %  DM1 -0.2 ± 0.6 %  DM10 +0.5 ± 1.0 %  DM11 +3.4 ± 2.2 %
per channel alone (POG SFs fixed):  tautau 1.043 +0.079 −0.072  mutau 1.049 +0.056 −0.051  etau 1.017 +0.068 −0.063  emu 0.978 +0.056 −0.053
v3 (τhτh alone):  μ_Z = 1.071 +0.114 −0.100,  σ = 2082 +222 −194 pb
```

Largest impacts on μ_Z: NormFactors 4.1 %, Electron efficiency 3.9 %, Tau trigger 2.0 %, Gammas 1.8 %, Fakes 1.7 %, MET 1.5 %; data statistics 0.5 %.

---

# Handoff – Z → τ⁺τ⁻ (τhτh), v3

## Team

- Niels ter Linden (v1, with Claude Code agents); v2, v2.1 and v3 by Claude Code (Fable 5.1) for Samuel Jankovych,
  after the review in `REVIEW.md`.

## What we worked on

- Complete **Z/γ* → τhτh cross-section measurement** on the CMS 2016 Tau dataset (Run2016G+H):
  - jet→τh fakes from a data-driven fake factor (MC subtracted, closure corrected, no QCD simulation);
  - UL16 simulation for the **fiducial** Z→ττ signal (inclusive + jet-binned aMC@NLO, stitched), the
    non-fiducial Z/γ*→ττ as a theory-normalised background, and the small backgrounds, with TauPOG corrections;
  - a **MET-corrected di-τ mass**;
  - a **k-fold BDT** (XGBoost, mass-agnostic kinematics) sorting the signal region into three categories;
  - a **TRExFitter v1.8.0** profile-likelihood fit of m_ττ in the three categories, in the shared `fitting/` conventions.
- Skims (1.7 GB), flat ntuples (330 MB, laptop bundle) of all data and simulation, the five BDT fold models.
- Documentation: `README.md`, `CLAUDE.md`, `REVIEW.md` (the review of v1 and the numerical studies behind it in
  `review/`), `docs/00–09`, slides `slides/ztautau_slides.pdf`.

**Result (v3, nominal: DeepTau Tight on both legs, MC-subtracted fake factor):**

> **σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2082 ± 41 (stat) +222/−194 (syst+stat) ± 76 (acc) pb** (NNLO 1945 pb)
> **σ_fid(τhτh) = 4.82 ± 0.09 (stat) ± 0.47 (syst) pb** (prediction 4.50 pb) — μ_Z = 1.071 +0.114 -0.100

Cross-check without the MC subtraction in the fake-factor regions: μ_Z = 1.166 +0.142 -0.124.
The same chain with DeepTau Medium (v2.1) gave μ_Z = 1.205 +0.145 −0.125, σ(60–120) = 2343 pb (`docs/08`, working-point comparison).

## Data and simulation used

| Dataset | Type (data / MC) | Record / DOI | Notes |
|---|---|---|---|
| Tau Run2016G | data | 30532 | 45 files, 79.6 M events, 7653.3 pb⁻¹ |
| Tau Run2016H | data | 30565 | 55 files, 76.8 M events, 8740.1 pb⁻¹ |
| DYJetsToLL_M-50 amcatnloFXFX | MC | 35669 | signal (fiducial LHE ττ) + non-fiducial ττ + Z→ee/μμ; σ = 6077.22 pb; normalisation and acceptance |
| DYJetsToLL_0J / 1J / 2J amcatnloFXFX | MC | 35577 / 35595 / 35613 | signal statistics, stitched per LHE_NpNLO bin (no extra cross sections) |
| DYJetsToLL_M-50 madgraphMLM | MC | 35671 | generator cross-check of the fiducial C (24/61 files) |
| DYJetsToLL_M-10to50 amcatnloFXFX | MC | 35631 | Z/γ*→ℓℓ with m < 50 GeV (16 events in the SR) |
| WJetsToLNu amcatnloFXFX | MC | 69745 | genuine-τ1 part only; subtracted from the FF regions with uniform weights |
| TTTo2L2Nu / TTToSemiLeptonic powheg | MC | 67801 / 67993 | subsets, from EOS |
| ST_tW top / antitop | MC | 64895 / 64839 | from EOS |
| WW / WZ / ZZ inclusive pythia8 | MC | 72696 / 72754 / 75593 | |
| TauPOG SFs, UL2016postVFP, DeepTau2017v2p1 | external | GitHub cms-tau-pog/TauIDSFs, TauTriggerSFs | `external/tau_pog_UL2016postVFP.json` |

- Skims: `/data/atlas/users/sjankovy/BND-school-cache/ztautau/skims_v1/` (DY samples re-skimmed with `LHE_NpNLO`;
  the rest are links to `/data/atlas/users/nterlind/BND-school-cache/ztautau/skims_v1/`)
- **Ntuples (use these):** `/data/atlas/users/sjankovy/BND-school-cache/ztautau/ntuples_v1/`
- Luminosity: **16393.381 pb⁻¹ ± 1.2 %** (normtag, as in z-mumu v2 and `fitting/CONVENTIONS.md`).

## Selection and method

- **Trigger:** `HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg` (G) / `HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg`
  (H). Both legs are matched to HLT τ objects (filterBits 2, pT > 35, ΔR < 0.5).
- **Object selection:** τh pT > 40 GeV, |η| < 2.1, |dz| < 0.2, DM 0/1/10/11; DeepTau2017v2p1 VSe ≥ VVLoose,
  VSmu ≥ VLoose; VSjet Tight (SR) or VVVLoose-not-Tight (fake-factor regions); most isolated pair with
  ΔR > 0.5; extra e/μ veto (orthogonal to eτ, μτ, ee, μμ).
- **Event selection:** OS, both Tight → 21 160 events (64 % fakes). No mass cut; m_ττ from 0 to 350 GeV is
  fitted in three BDT categories.
- **Signal definition:** `DYtautau` = LHE ττ with 60 < m_LHE < 120 GeV, both τ hadronic, both visible τ pT > 40,
  |η| < 2.1 (62 % of the selected Z/γ*→ττ). The rest (`DYtautau_nonfid`, 30 % of it m_LHE > 120 GeV) is a
  background with a 5 % normalisation NP and the theory variations on the non-fid/fid ratio.
- **Fakes:** FF(era, DM, N_jets, pT) = [N(SS, τ1 Tight) − MC] / [N(SS, τ1 VVVLoose!Tight) − MC], τ2 Tight,
  times closure corrections f(|η(τ1)|) (±15 %) and g(pT(τ2)) (−6 %); applied to OS events with τ1 failing Tight
  (MC with a genuine τ1 subtracted, W+jets with uniform weights); multiplied by C_OS/SS(era, N_jets, BDT category) = 1.04–1.15 from the
  τ2 anti-isolated sideband. FF statistics per event into the template variance; residual non-closure per
  category × mass region as nuisance parameters.
- **Mass:** m_ττ is a per-event likelihood over the τ visible-energy fractions using the MET covariance
  (posterior median, 1/m² prior). Scale 0.99 and resolution 11 %, against 0.80 and 13 % for m_vis.
- **BDT:** XGBoost on 16 mass-agnostic inputs (ΔR, Δφ, pT balance of ττ + MET, MET significance, DMs, …),
  signal = fiducial Z→ττ MC, background = AR data × FF; 5 folds by event number; held-out AUC 0.970.
  Categories: score < 0.55 / 0.55–0.90 / > 0.90 (S/B ≈ 5.3). Same-sign closure of the FF in the score verified.
- **Corrections / weights:** TauPOG τh ID (per DM), energy scale, trigger leg SFs (genuine legs only), e/μ→τh SFs,
  pileup (Open Data lumi table, 69.2 mb), L1 prefiring.
- **Fit:** TRExFitter v1.8.0, m_ττ in 3 × 14 bins, POI `mu_Z` on `DYtautau`, 39 NPs + 42 γ.

## How to run

```bash
cd z-tautau
source ../setup.sh
python run_all.py --from 3        # from the ntuples: FF, BDT, histograms, fits, report (~25 min)
python run_all.py                 # from NanoAOD: ~2 h
rm -f slides/figs/*; python slides/make_figures.py      # slide figures (dark, vector) from the outputs above
env -u PYTHONPATH -u LD_LIBRARY_PATH -u PYTHONHOME <betterplottingtool venv>/bin/python slides/build_deck.py   # the deck
```

Every plot (`output/plots/`, `slides/figs/`) carries the stamp `v3: DeepTau Tight τh` (`config.PLOT_TAG`) and every
slide a `v3 | DeepTau Tight` footer, so figures of different versions or working points cannot be mixed up.

## Results so far

| quantity | v3 (Tight, nominal) | v2.1 (Medium, same chain) |
|---|---:|---:|
| μ_Z | 1.071 +0.114 -0.100 | 1.166 +0.142 -0.124 |
| expected (Asimov) | +0.106 −0.093 | +11.8 / −10.2 % |
| σ(60–120) [pb] | 2082 +222 −194 | 2357 |
| σ_fid [pb] | 4.82 ± 0.09 ± 0.47 | 4.82 ± 0.09 ± 0.47 |
| GoF p | 0.22 | 0.16 |

Uncertainty on μ_Z (nominal, grouped impacts):

| source | impact |
|---|---:|
| Tau ID | 8.4 % |
| Fakes | 6.5 % |
| Gammas | 3.9 % |
| Tau trigger | 3.0 % |
| Background normalisation | 2.6 % |
| Tau energy scale | 1.5 % |
| Signal modelling | 1.3 % |
| MET | 0.8 % |
| Luminosity | 0.7 % |
| Pileup | 0.6 % |
| L1 prefiring | 0.2 % |
| **data statistics** | **2.1 %** |

Prediction: σ(Z/γ*→ττ, 60–120) = 1944.9 pb; σ_fid = 4.501 pb; A = 0.002314 ± 3.7 %
(scale 3.4 %, PDF 0.5 %, α_s 0.6 %, ISR 0.9 %, FSR 0.2 %); C = 0.0510.

## For the combination

**Profile-likelihood combination (recommended).**

- Inputs, **all committed** (exception to the no-ROOT rule, they are 0.8–1.2 MB and the combination needs them
  from any checkout): the histograms `fit/fitinputs/ztautau.root` (+ the `.meta.json` sidecar), the config `fit/ztautau.config`, the workspace
  `fit/results/ztautau/RooStats/ztautau_combined_ztautau_model.root` (a copy of TRExFitter's
  `ztautau_allBinsFitRegions_combined_ztautau_model.root`, the name it uses when a region has `DropBins`)
  and `fit/results/ztautau_fit_result.json`. Job `ztautau`, POI `mu_Z`, regions `tautau_SR0`, `tautau_SR1`,
  `tautau_SR2` (SR0 fitted above 110 GeV only). The entry in `fitting/combination_skeleton.config` already
  points there. Regenerate with `python run_all.py --from 4 --to 6` (~15 min from the ntuples).
- Correlated NPs (shared names): `Lumi`, `Pileup`, `L1Prefiring`, `QCDScale`, `PDF`, `PS_ISR`, `PS_FSR`,
  `XS_TTbar`, `XS_SingleTop`, `XS_WW`, `XS_WZ`, `XS_ZZ`.
- ττ-only NPs: `TauID_DM*`, `TauTrigger_DM*`, `TauES_DM*`, `TauFakeEle`, `TauFakeMu`, `MET_Unclustered`,
  `XS_DYll`, `XS_DYlowmass`, `XS_DYtautau_nonfid`, `XS_WJets`, `MCStatNorm_WJets_tautau`, `FakeOSSS_tautau_c*`,
  `FakeClosure_tautau_c*_lo|hi`. **`SigModel_tautau` is reported, not fitted, and must never be correlated
  with z-mumu's `SigModel`** (different quantities).
- In the combination μ_Z will be fixed by ee/μμ, and this channel will mainly constrain the τh ID parameters;
  present the fitted `TauID_DM*` as a result. The 60–120 GeV denominator of σ is Born level here (m_LHE);
  the acceptance uncertainty (3.7 %) is correlated.

**BLUE / counting (`combination.ipynb`).**

| variable | value |
|---|---|
| `n_obs` | 21 160 |
| `n_bkg` (prefit, incl. non-fiducial DY) | 41459 |
| `acc_eff` | 0.0001179 (for σ(60–120), L = 16393.381, fiducial signal) |

> ⚠️ Counting is useless here (B/S ≈ 8 over the whole signal region): use the fitted σ = 2082 +222 −194 pb with the
> uncertainties above. Correlated: lumi 1.2 %, acceptance 3.7 %. Uncorrelated: τh ID, trigger, fakes, MC stat.

## Open issues / next steps

1. τh ID scale factors dominate: decay-mode categories and/or the combined fit to measure `TauID_DM*` in situ;
   quote the pT-binned POG prescription (14 % different on the yield) as a cross-check.
2. Trigger efficiency in situ (μτh tag-and-probe from SingleMuon) instead of the POG turn-on SFs.
3. W+jets: HT-binned madgraph samples with LHE_HT stitching (the inclusive sample has pathological weights).
4. High-mass DY sample (record 35629) stitched in m_LHE for the non-fiducial template; EWK Z→ττ.
5. Add eτh and μτh (the skims would need looser lepton content).
