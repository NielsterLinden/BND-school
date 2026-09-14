# Handoff – Z → τ⁺τ⁻ (τhτh)

## Team

- Niels ter Linden (with Claude Code agents)

## What we worked on

- Complete **Z/γ* → τhτh cross-section measurement** on the CMS 2016 Tau dataset (Run2016G+H):
  - jet→τh fakes from a data-driven classic fake factor (no QCD simulation);
  - UL16 simulation for Z→ττ and the small backgrounds, with TauPOG corrections;
  - a **MET-corrected di-τ mass**;
  - a **TRExFitter v1.8.0** profile-likelihood fit in the shared `fitting/` conventions.
- Skims (1.6 GB) and flat ntuples (300 MB, laptop bundle) of all data and simulation.
- Documentation: `README.md`, `CLAUDE.md`, `docs/00–08`, slides `slides/ztautau_slides.pdf`.

**Result (nominal = fake factors without MC subtraction):**

> **σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2255 ± 36 (stat) +371/−315 (syst) ± 83 (acc) pb** (NNLO 1945 pb)
> **σ_fid(τhτh) = 5.22 ± 0.08 (stat) ± 0.80 (syst) pb** (prediction 4.50 pb) — μ_Z = 1.160 +0.192 −0.163

Alternative (genuine-τ subtraction in the fake-factor regions): μ_Z = 1.208 +0.201 −0.169, σ = 2348 pb.

## Data and simulation used

| Dataset | Type (data / MC) | Record / DOI | Notes |
|---|---|---|---|
| Tau Run2016G | data | 30532 | 45 files, 79.6 M events, 7653.3 pb⁻¹ |
| Tau Run2016H | data | 30565 | 55 files, 76.8 M events, 8740.1 pb⁻¹ |
| DYJetsToLL_M-50 amcatnloFXFX | MC | 35669 | signal (LHE ττ) + Z→ee/μμ; σ = 6077.22 pb |
| DYJetsToLL_M-50 madgraphMLM | MC | 35671 | generator systematic (24/61 files) |
| DYJetsToLL_M-10to50 amcatnloFXFX | MC | 35631 | 15 events in the SR |
| WJetsToLNu amcatnloFXFX | MC | 69745 | genuine-τ1 part only |
| TTTo2L2Nu / TTToSemiLeptonic powheg | MC | 67801 / 67993 | subsets, from EOS |
| ST_tW top / antitop | MC | 64895 / 64839 | from EOS |
| WW / WZ / ZZ inclusive pythia8 | MC | 72696 / 72754 / 75593 | |
| TauPOG SFs, UL2016postVFP, DeepTau2017v2p1 | external | GitHub cms-tau-pog/TauIDSFs, TauTriggerSFs | `external/tau_pog_UL2016postVFP.json` |

- Skims: `/data/atlas/users/nterlind/BND-school-cache/ztautau/skims_v1/`
- **Ntuples (use these):** `/data/atlas/users/nterlind/BND-school-cache/ztautau/ntuples_v1/`
- Luminosity: **16393.381 pb⁻¹ ± 1.2 %** (normtag, as in z-mumu v2 and `fitting/CONVENTIONS.md`).

## Selection and method

- **Trigger:** `HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg` (G) / `HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg`
  (H). Both legs are matched to HLT τ objects (filterBits 2, pT > 35, ΔR < 0.5).
- **Object selection:**
  - τh pT > 40 GeV, |η| < 2.1, |dz| < 0.2, DM 0/1/10/11;
  - DeepTau2017v2p1 VSe ≥ VVLoose, VSmu ≥ VLoose;
  - VSjet Medium (SR) or VVVLoose-not-Medium (fake-factor regions);
  - most isolated pair with ΔR > 0.5;
  - extra e/μ veto (orthogonal to eτ, μτ, ee, μμ).
- **Event selection / mass window:** OS, both Medium → 47 586 events (81 % fakes, 16 % Z→ττ). No mass cut;
  m_tt from 0 to 350 GeV is fitted.
- **Fakes:**
  - FF(era, DM, N_jets, pT) = N(SS, τ1 Medium) / N(SS, τ1 VVVLoose!Medium), with τ2 Medium;
  - applied to OS events with τ1 failing Medium;
  - multiplied by C_OS/SS(era, N_jets) = 1.07–1.13, from the τ2 anti-isolated sideband;
  - no MC subtraction in the nominal (as requested).
- **Mass:** m_tt is a per-event likelihood over the τ visible-energy fractions, using the MET covariance
  (posterior median, 1/m² prior). Scale 0.99 and resolution 11 %, against 0.80 and 13 % for m_vis.
- **Corrections / weights:** TauPOG τh ID (per DM), energy scale, trigger leg SFs, e/μ→τh SFs, pileup
  (Open Data lumi table, 69.2 mb), L1 prefiring.
- **Fit:** TRExFitter v1.8.0, m_tt in the SR (14 bins), POI `mu_Z` on `DYtautau`, 38 NPs + 14 γ.

## How to run

```bash
cd z-tautau
source ../setup.sh
python run_all.py --from 3        # from the ntuples: ~10 min
python run_all.py                 # from NanoAOD: ~1 h
```

## Results so far

| quantity | nominal | MC-subtracted FF |
|---|---:|---:|
| μ_Z | 1.160 +0.192 −0.163 | 1.208 +0.201 −0.169 |
| expected (Asimov) | +0.169 −0.140 | +0.169 −0.140 |
| σ(60–120) [pb] | 2255 | 2348 |
| σ_fid [pb] | 5.22 | 5.43 |
| GoF p | 0.18 | 0.54 |

Uncertainty on μ_Z (nominal):

| source | impact |
|---|---:|
| τh ID | 12.5 % |
| signal modelling (LO vs NLO 8.5 %, scales 4.3 %) | 9.8 % |
| τh trigger | 5.1 % |
| MC statistics | 4.8 % |
| τh energy scale | 3.2 % |
| luminosity | 1.3 % |
| pileup | 1.2 % |
| fakes | 0.8 % |
| **data statistics** | **1.8 %** |

Prediction: σ(Z/γ*→ττ, 60–120) = 1944.9 pb; σ_fid = 4.501 pb; A = 0.002314 ± 3.7 %
(scale 3.4 %, PDF 0.5 %, α_s 0.6 %, ISR 0.9 %, FSR 0.2 %); C = 0.107.

## For the combination

**Profile-likelihood combination (recommended).**

- Inputs: `fit/ztautau.config` and the workspace `fit/results/ztautau/RooStats/ztautau_combined_ztautau_model.root`.
  Job `ztautau`, POI `mu_Z`, region `tautau_SR`. The entry in `fitting/combination_skeleton.config` already
  points there.
- Correlated NPs (shared names): `Lumi`, `Pileup`, `L1Prefiring`, `QCDScale`, `PDF`, `PS_ISR`, `PS_FSR`, `SigModel`.
- ττ-only NPs: `TauID_DM*`, `TauTrigger_DM*`, `TauES_DM*`, `TauFakeEle`, `TauFakeMu`, `MET_Unclustered`,
  `XS_DYll`, `XS_DYlowmass`, `MCStatNorm_WJets_tautau`, `Fake*_tautau`.
- In the combination μ_Z will be fixed by ee/μμ, and this channel will mainly constrain the τh ID parameters.

**BLUE / counting (`combination.ipynb`).**

| variable | value |
|---|---|
| `n_obs` | 47 586 |
| `n_bkg` | 40 235 (prefit) |
| `acc_eff` | 2.471e-4 (for σ(60–120), L = 16393.381) |

> ⚠️ Counting gives 1815 pb: B/S ≈ 5, and only the shape fit fixes the fake normalisation (−2.8 %).
> Use the fitted σ = 2255 pb, +16.6/−14.1 %. Correlated: lumi 1.2 %, acceptance 3.7 %. Uncorrelated:
> τh ID 12.5 %, trigger 5.1 %, fakes, MC stat.

## Open issues / next steps

1. Make the MC-subtracted fake factor nominal. It removes a known double counting and raises μ by 0.05.
2. τh ID scale factors dominate. In a combined fit, let the ττ channel constrain `TauID_DM*` (lepton
   universality) instead of μ_Z.
3. The LO-vs-NLO generator systematic (8.5 %) is conservative. An NLO-only prescription (scales + PS)
   would give ~5 %.
4. Fake-factor closure corrections in η(τ1) and pT(τ2).
5. More signal MC statistics (jet-binned DY) and HT-binned W+jets.
6. Add eτh and μτh (the skims would need looser lepton content). They measure the τh ID SF in situ.
