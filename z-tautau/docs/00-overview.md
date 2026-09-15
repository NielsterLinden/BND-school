# 00 — The Z → τhτh measurement in one page

**Goal.** Measure σ(pp → Z/γ* → ττ) with both τ decaying hadronically, in CMS 2016 Open Data (Run2016G+H,
Tau dataset, 16.4 fb⁻¹), with a binned profile-likelihood fit (TRExFitter v1.8.0) in the same framework as
Z→μμ and Z→ee, so the three channels can be fitted together. This is v2: v1 was reviewed (`REVIEW.md`) and
rebuilt.

**Chain.** `run_all.py` runs eight standalone scripts:

| step | script | what | doc |
|---|---|---|---|
| 0 | `step0_external.py` | file lists; TauPOG SFs → `external/*.json` | 07 |
| 1 | `step1_skim.py` | 156 M data + MC events (incl. jet-binned DY) → 1.7 GB NanoAOD-format skims | 03 |
| 2 | `step2_ntuples.py` | pair selection, trigger matching, vetoes, **di-τ masses** → 330 MB flat ntuples | 02, 04 |
| 3 | `step3_fakefactors.py` | **fake factors** (era × DM × N_jets × pT, MC subtracted), closure corrections, C_OS/SS, closure | 05 |
| 3b | `step3b_bdt.py` | **k-fold BDT** training and validation (closure of the FF in the score) | 09 |
| 4 | `step4_histograms.py` | templates per BDT category + ~70 systematic variations, control plots | 06, 07 |
| 5 | `step5_fit.py` | TRExFitter config + fit, stat-only, Asimov, impacts, ranking | 08 |
| 6 | `step6_report.py` | `output/results.json`, `RESULTS.md`, summary plots | 06, 08 |

**Physics in six lines.**
1. Two τh with pT > 40 GeV, |η| < 2.1, DeepTau Medium, opposite sign, di-τ trigger matched, lepton vetoes:
   47 586 events.
2. 80 % are jet→τh fakes: estimated from data with fake factors measured in same-sign events (simulated
   genuine τ subtracted, closure corrected in |η(τ1)| and pT(τ2)) and applied to events whose leading τ fails
   Medium (no QCD simulation).
3. The MET is folded into the visible mass with a per-event likelihood (MET covariance + τ decay phase
   space). m_tt peaks at m_Z with 11 % resolution; fakes sit at 130–300 GeV.
4. The signal is the **fiducial** Z→ττ (60 < m < 120 GeV, both visible τ pT > 40, |η| < 2.1; aMC@NLO
   inclusive + jet-binned, TauPOG corrections); the non-fiducial Z/γ*→ττ (38 % of the selected DY, mostly
   m > 120 GeV) is a theory-normalised background, like tt̄, W+jets, Z→ee and dibosons.
5. A k-fold BDT on mass-agnostic kinematics sorts the signal region into three categories (S/B 0.02, 0.6, ≈ 5.3).
6. μ_Z from the m_tt fit in the three categories; σ = μ_Z × prediction.

**Result.**

> **σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2358 ± 40 (stat) +282/−242 (syst+stat) ± 86 (acc) pb**
> (NNLO 1945 pb; μ_Z = 1.212 +0.145 -0.125) — **σ_fid(τhτh) = 5.46 ± 0.09 (stat) ± 0.60 (syst) pb** (prediction 4.50 pb)

Limited by the external τh identification scale factors (11.2 %); data statistics 2.0 % on μ_Z. The goodness of fit is poor (p = 0.01): data/prediction rises with the visible-τ pT in the signal-dominated category (`08-fit-and-results.md`).

**Cross-check with DeepTau Tight** on both legs (`variants/tight/`): μ_Z = 1.099 +0.116 −0.102, σ(60–120) = 2137 +226 −198 pb, GoF p = 0.32: more precise and a better fit (`08-fit-and-results.md`).

**Read next:** `02-selection.md` for the selection, `05-fake-factors.md` for the background method,
`09-bdt.md` for the categories, and `08-fit-and-results.md` for the result and what the combination should use.
