# 00 — The Z → ττ measurement in one page

**Goal.** Measure σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) in CMS 2016 Open Data (Run2016G+H, 16.4 fb⁻¹) with a
binned profile-likelihood fit (TRExFitter v1.8.0) in the same framework as Z→μμ and Z→ee, so the three
channels can be fitted together. Four ττ final states are fitted simultaneously: **τhτh** (Tau dataset, three
BDT categories), **μτh** (SingleMuon), **eτh** (SingleElectron) and **eμ** (MuonEG, with a tt̄ control
region). The τh identification scale factors and the τh energy scale are measured *in situ* by that fit
instead of being taken from the TauPOG. History: v1 was reviewed (`REVIEW.md`) and rebuilt, v3 was the
τhτh-only measurement, v4 added the three lepton channels and was reviewed in turn (`REVIEW_v4.md`,
answered in `REVIEW_v4_RESPONSE.md`). Sections 1–9 of these docs describe the τhτh part, which is still the
backbone; `10-v4-plan.md` describes the four-channel measurement and `11-combination-inputs.md` what the
combination gets.

**Chain.** `run_tautau_base.py` builds the τhτh base, `run_all.py` the measurement:

| step | script | what | doc |
|---|---|---|---|
| 0 | `step0_external.py` | file lists; TauPOG SFs → `external/*.json` | 07 |
| 1 | `step1_skim.py` (`--v4`) | 156 M τhτh + 4 lepton-stream data/MC events → NanoAOD-format skims | 03, 10 |
| 2 | `step2_ntuples_tautau.py`, `step2_ntuples_lepton.py` | pair selection, trigger matching, vetoes, **di-τ masses** → flat ntuples | 02, 04, 10 |
| 3 | `step3_fakefactors.py` | τhτh **fake factors** (era × DM × N_jets × pT, MC subtracted), closure, C_OS/SS | 05 |
| 3b | `step3b_bdt.py` | **k-fold BDT** training and validation (closure of the FF in the score) | 09 |
| 3c | `step3c_trigger.py` | **in-situ** single-electron and eμ cross-trigger efficiencies | 10 |
| 3d | `step3d_fakes_lepton.py` | per-process fake factors of μτh / eτh and the eμ multijet estimate | 10 |
| 4a | `step4a_tautau_base.py` | τhτh Data, fake estimate and BDT categories → `fit/fitinputs/tautau_base.root` | 05, 09 |
| 4 | `step4_histograms.py` | templates of all four channels + ~90 systematic variations → `fit/fitinputs/ztautau.root` | 06, 07, 10 |
| 4b | `step4b_export_channels.py` | one fit-input file per channel (for the per-channel workspaces) | 11 |
| 5 | `step5_fit.py` | TRExFitter config + fit, stat-only, Asimov, impacts, ranking, the cross-check fits | 10 |
| 5b | `step5b_multifit.py` | MultiFit of the four per-channel workspaces | 11 |
| 6 | `step6_report.py` | `output/results.json`, `RESULTS.md`, summary plots | 06, 10, 11 |

**Physics in eight lines.**
1. **τhτh** (Tau dataset): two τh with pT > 40 GeV, |η| < 2.1, DeepTau Tight, opposite sign, di-τ trigger
   matched, lepton vetoes — 21 160 events, 64 % of them jet→τh fakes.
2. **μτh / eτh** (SingleMuon, SingleElectron): one isolated lepton on its single-lepton trigger plateau, one
   τh with pT > 30 GeV, opposite sign, m_T(ℓ, MET) < 40 GeV, b-jet-free; split into one region per τh decay
   mode, because that is what makes the per-decay-mode τh identification scale factors measurable.
3. **eμ** (MuonEG cross triggers): one muon and one electron, opposite sign, D_ζ > −20 GeV, b veto — no τh
   at all, so this channel measures the cross section without any τh identification scale factor; a second
   region with D_ζ < −40 GeV and MET > 80 GeV normalises tt̄.
4. Jet→τh fakes come from data: same-sign fake factors in τhτh, per-process fake factors (multijet /
   W+jets / tt̄, with the fractions of the application region) in μτh and eτh, and same-sign data scaled by
   an OS/SS factor for the eμ multijet. No QCD simulation anywhere.
5. The MET is folded into the visible mass with a per-event likelihood (MET covariance + τ decay phase
   space). m_tt peaks at m_Z with 11 % resolution; fakes sit at 130–300 GeV.
6. The signal is Z/γ*→ττ with 60 < m_LHE < 120 GeV in **all** decays (aMC@NLO inclusive + jet-binned);
   Z/γ*→ττ outside that window is a theory-normalised background, like tt̄, W+jets, Z→ee/μμ and dibosons.
   The theory variations are renormalised to a fixed σ(60–120), so they vary the acceptance only.
7. A k-fold BDT on mass-agnostic kinematics sorts the τhτh signal region into three categories
   (S/B 0.04, 1.1, ≈ 10.9); the lowest one is fitted above 110 GeV only, as a fake sideband.
8. One simultaneous fit of m_tt in all 13 regions gives μ_Z, the four τh identification scale factors and
   the four τh energy scales; σ = μ_Z × σ^pred(60–120) = μ_Z × 1944.9 pb.

<!-- RESULT:BEGIN -->
**Result** (four channels, τh ID scale factors and energy scale fitted in situ)

```
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 1981 +153 −136 pb   (stat ±9;  prediction 1945 pb: aMC@NLO acceptance, NNLO normalisation)
μ_Z = 1.019 +0.079 −0.070   (stat ±0.005, syst ±0.074),   expected ±0.073,   GoF p = 0.15,   μ_tt̄ = 1.11 ± 0.04
includes the τh ID scale-factor pT dependence (TauIDpT_tautau, ±6.3 % on the signal: μ_Z 1.019 flat, 0.954 pT-split);   without it 1981 +73 −70 pb
τh ID SF (Tight):  DM0 0.989 ± 0.040 (POG 0.90 ± 0.13)  DM1 0.963 ± 0.034 (POG 0.89 ± 0.05)  DM10 0.896 ± 0.035 (POG 0.94 ± 0.15)  DM11 0.795 ± 0.050 (POG 0.81 ± 0.15)
τh energy scale:   DM0 -0.6 ± 0.9 %  DM1 -0.2 ± 0.6 %  DM10 +0.6 ± 1.0 %  DM11 +3.2 ± 2.2 %
μ_Z of the two sub-measurements the fit combines:  e mu alone 0.960 +0.042 −0.040;  τ channels alone (SF free) 1.203 +0.091 −0.084   (2.6 σ apart)
per channel alone (POG SFs fixed):  tautau 1.043 +0.079 −0.072  mutau 1.048 +0.055 −0.051  etau 1.016 +0.070 −0.065  emu 0.960 +0.042 −0.040
```

Largest grouped impacts on μ_Z: Tau ID pT dependence 6.5 %, NormFactors 3.1 %, Emu trigger 2.5 %, Gammas 1.4 %, Tau trigger 1.2 %, Fakes 1.2 %; data statistics 0.5 %. The categories overlap, so their quadrature sum exceeds the MINOS total by a factor 1.13 (`output/RESULTS.md`).

Ranking: TauIDpT_tautau 6.8 %, mu_ttbar 3.1 %, EmuTrigger 2.6 %, TauIDSF_DM1 2.2 %, TauIDSF_DM10 2.1 %, TauIDSF_DM0 1.9 %.
<!-- RESULT:END -->

**Read next:** `10-v4-plan.md` for the four-channel measurement (the result), `02-selection.md` for the τhτh
selection, `05-fake-factors.md` for the background method, `09-bdt.md` for the categories, and
`11-combination-inputs.md` for what the combination gets and what it must know before using it.
