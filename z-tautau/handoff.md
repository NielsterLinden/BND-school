# Handoff — Z → τ⁺τ⁻ (τhτh + μτh + eτh + eμ)

## Team

- Niels ter Linden (v1, with Claude Code agents); v2–v4 by Claude Code (Fable 5.1) for Samuel Jankovych,
  after the reviews in `REVIEW.md` (the τhτh iteration) and `REVIEW_v4.md` (the four-channel measurement).

## State (17 September 2026)

- **One measurement, one set of results.** Four channels are fitted together (`run_all.py`,
  `docs/10-v4-plan.md`): the τhτh BDT categories plus μτh (SingleMuon, IsoMu24), eτh (SingleElectron,
  Ele27_WPTight; EOS only) and eμ (MuonEG cross triggers) with an eμ tt̄ control region. The τhτh-only fit
  of v3 no longer exists; what remains of that chain is the τhτh Data and fake templates
  (`fit/fitinputs/tautau_base.root`, `run_tautau_base.py`).
- **No μμ channel**; second-muon / second-electron vetoes make every channel orthogonal to the z-mumu and
  z-ee selections (the eμ channel shares events with z-mumu's tt̄ *control* region `mumu_CRemu`, which a
  joint fit must drop).
- Signal: Z/γ*→ττ with 60 < m_LHE < 120 GeV in every decay (one μ_Z); the rest of the DY ττ simulation is a
  background. σ^pred(60–120) = 1944.9 pb, ours (z-mumu and z-ee have their own).
- τh ID scale factors free per decay mode (TRExFitter NormFactors, products on the τhτh templates), τh energy
  scale with a 3 % prior: both measured in situ as in CMS arXiv:1801.03535. Every source of the paper's
  Table 2 is implemented from an official correction, measured in situ, or replaced by a stated estimate
  (`docs/10-v4-plan.md` §9).
- Lepton-channel fakes: per-process fake factors (multijet / W+jets / tt̄) with AR fractions, OS/SS and m_T
  corrections, same-sign validation (0.97 ± 0.01 μτh, 1.03 ± 0.02 eτh). Two things learned: the isolated-lepton
  same-sign region is ~50 % W+jets, and the W+jets fake factor is charge-correlated (OS 0.08, SS 0.04).
- Outputs: `output/RESULTS.md`, `output/results.json`, `output/plots/`, `fit/ztautau.config`,
  `fit/fitinputs/ztautau.root`, `fit/results/ztautau_fit_result.json`, `external/trigger_insitu_v4.json`.
- Caches: `skims_v1/`, `ntuples_v1/` (τhτh), `skims_v4/` (17 GB), `ntuples_v4/` (2.9 GB) under
  `/data/atlas/users/sjankovy/BND-school-cache/ztautau/`.

<!-- RESULT:BEGIN -->
**Result** (four channels, τh ID scale factors and energy scale fitted in situ)

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 1981 +73 −70 pb   (stat ±9;  prediction 1945 pb: aMC@NLO acceptance, NNLO normalisation)
μ_Z = 1.019 +0.037 −0.036   (stat ±0.005, syst ±0.036),   expected ±0.036,   GoF p = 0.15,   μ_tt̄ = 1.11 ± 0.04
τh ID SF (Tight):  DM0 0.989 ± 0.040 (POG 0.90 ± 0.13)  DM1 0.963 ± 0.034 (POG 0.89 ± 0.05)  DM10 0.896 ± 0.035 (POG 0.94 ± 0.15)  DM11 0.795 ± 0.050 (POG 0.81 ± 0.15)
τh energy scale:   DM0 -0.6 ± 0.9 %  DM1 -0.2 ± 0.6 %  DM10 +0.6 ± 1.0 %  DM11 +3.2 ± 2.2 %
μ_Z of the two sub-measurements the fit combines:  e mu alone 0.960 +0.042 −0.040;  τ channels alone (SF free) 1.203 +0.091 −0.084   (2.6 σ apart)
per channel alone (POG SFs fixed):  tautau 1.043 +0.079 −0.072  mutau 1.048 +0.055 −0.051  etau 1.016 +0.070 −0.065  emu 0.960 +0.042 −0.040
```

Largest grouped impacts on μ_Z: NormFactors 3.1 %, Emu trigger 2.5 %, Gammas 1.4 %, Tau trigger 1.2 %, Fakes 1.2 %, Electron energy 1.2 %; data statistics 0.5 %. The categories overlap, so their quadrature sum exceeds the MINOS total by a factor 1.46 (`output/RESULTS.md`).

Ranking: mu_ttbar 3.1 %, EmuTrigger 2.6 %, TauIDSF_DM1 2.2 %, TauIDSF_DM10 2.1 %, TauIDSF_DM0 1.9 %, TauIDSF_DM11 1.3 %.
<!-- RESULT:END -->

## What the review changed

`REVIEW_v4.md` reviewed the four-channel measurement; `REVIEW_v4_RESPONSE.md` is the point-by-point answer.
In short:

* The eμ trigger prior applied the paper's 2 % **once per leg**, so it was 5.4 % instead of 2.7 %. Applying
  it once per event moved the result from μ_Z = 1.055 ± 0.045 to **1.019 ± 0.037**; the cross-check fit
  `ztautau_emutrig2x` restores the old treatment and gives 1.055 back, so the whole shift is that prior.
* The expected uncertainty was the HESSE error of a failed Asimov MINOS (±0.009). The cause was one bin
  with no data and no prediction, whose unconstrained MC-statistics γ made the per-bin offset log(0); such
  bins are now dropped and the Asimov fit converges at **±0.036**, next to the observed ±0.037.
* The eτh Ele27 turn-on and the ℓτh fake-composition priors were too tight and were corrected.
* The eμ and τ-channel sub-measurements are quoted separately, because they do not agree: 0.960 against
  1.203, 2.6 σ apart in the uncorrelated limit.
* **The τh ID scale factor is not flat in pT.** `ztautau_ptsplit` gives the ℓτh regions below 40 GeV their
  own scale factors and finds them 6 / 0 / 6 / 14 % higher; μ_Z then moves to 0.954, 1.8 times the total
  uncertainty, and the τ channels come down onto the eμ value. That 6 % is the size of the
  single-scale-factor assumption and it is larger than any experimental systematic we quote.
* The grouped impacts are reported with their (real) over-shoot instead of as if they added up to the total.

## For the combination

Everything is in **`docs/11-combination-inputs.md`** — files, region and nuisance-parameter names, which
NPs may be correlated with z-mumu / z-ee and which must not, the `mumu_CRemu` overlap, the two combination
routes, and the caveats that must be understood before our number is used. The short version:

| what | path |
|---|---|
| fit inputs (all channels / per channel) | `fit/fitinputs/ztautau.root`, `fit/fitinputs/ztautau_<ch>.root` (+ `.meta.json`) |
| configs | `fit/ztautau.config`, `fit/ztautau_<ch>.config`, MultiFit `fit/comb.config` |
| workspaces | `fit/results/ztautau[_<ch>]/RooStats/*_combined_*_model.root` |
| results | `fit/results/*_fit_result.json`, `output/results.json` → `for_combination` |

`output/results.json["for_combination"]["channel_result"]` has exactly the fields of
`combination/comb/inputs.ChannelResult`, including `sigma_pred`, `acc = {}` (the acceptance is profiled
inside the fit) and `groups_rescaled` (the grouped impacts scaled so that their quadrature sum equals the
MINOS total). The combination's `load_tautau` still expects the v3 layout and has to be updated once, and
its `--check` assertions on the ττ line will fail until then; §6 of `docs/11` gives the replacement.

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

- v4 skims / ntuples also cover SingleMuon, SingleElectron and MuonEG (Run2016G+H) and every simulation
  sample in the μτh / eτh / eμ selections; `docs/10-v4-plan.md` §10.

## Selection and method (τhτh; the lepton channels are `docs/10-v4-plan.md` §3–§5)

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
python run_tautau_base.py --from 3   # τhτh fake factors, BDT, base templates (~25 min from the ntuples)
python run_all.py --from 3           # trigger efficiencies, fakes, templates, all fits, report (~4 h)
python run_all.py --only 4           # templates only (~50 min)
python scripts/step5_fit.py --skip-ranking --skip-asimov     # the four-channel fit alone (~5 min)
```

Every plot carries the stamp `config.PLOT_TAG` (`v4: DeepTau Tight τh`), so figures of different versions
or working points cannot be mixed up. The slide deck is rebuilt with

```bash
python slides/make_figures.py                                  # the 25 dark vector figures (LCG)
env -u PYTHONPATH -u LD_LIBRARY_PATH -u PYTHONHOME \
    /project/atlas/users/sjankovy/boostHHbbtautau/HHARD_workfolder/betterplottingtool/venv/bin/python \
    slides/build_deck.py                                        # slides/ztautau_slides.pdf (needs PyMuPDF)
```

Every number in it is read from `output/results.json`, so the deck cannot drift from the measurement.

## Open issues / next steps

1. The eμ and τ-channel sub-measurements do not agree well (`RESULTS.md`, `docs/11` §7). The candidates are
   the eμ trigger efficiency, the pT dependence of the τh ID scale factor (measured now, see
   `ztautau_ptsplit`) and the τhτh trigger and fake modelling. This is the first thing to look at.
2. The eμ channel alone has the worst goodness of fit of the four; the 130–150 GeV m_tt bin of the lepton
   channels is where the leptonic likelihood mass sits ~10 % above m_LHE.
3. W+jets: HT-binned madgraph samples with LHE_HT stitching (the inclusive sample has pathological weights).
4. High-mass DY sample (record 35629) stitched in m_LHE for the `DYtautau_out` template; EWK Z→ττ; SM H→ττ.
5. τh trigger efficiency in situ (μτh tag-and-probe from SingleMuon) instead of the POG turn-on curves —
   the single-lepton and cross triggers are done, the di-τ one is not.
6. A τh trigger scale factor measured in situ would close the last external correction in the τhτh channel.
