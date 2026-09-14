# 00 — The Z → τhτh measurement in one page

**Goal.** Measure σ(pp → Z/γ* → ττ) with both τ decaying hadronically, in CMS 2016 Open Data (Run2016G+H,
Tau dataset, 16.4 fb⁻¹), with a binned profile-likelihood fit (TRExFitter v1.8.0) in the same framework as
Z→μμ and Z→ee, so the three channels can be fitted together.

**Chain.** `run_all.py` runs seven standalone scripts:

| step | script | what | doc |
|---|---|---|---|
| 0 | `step0_external.py` | file lists; TauPOG SFs → `external/*.json` | 07 |
| 1 | `step1_skim.py` | 156 M data + MC events → 1.6 GB NanoAOD-format skims | 03 |
| 2 | `step2_ntuples.py` | pair selection, trigger matching, vetoes, **di-τ masses** → 300 MB flat ntuples | 02, 04 |
| 3 | `step3_fakefactors.py` | **fake factors** (era × DM × N_jets × pT), C_OS/SS, closure | 05 |
| 4 | `step4_histograms.py` | templates + ~70 systematic variations, control plots | 07 |
| 5 | `step5_fit.py` | TRExFitter config + fit, stat-only, Asimov, impacts, ranking | 08 |
| 6 | `step6_report.py` | `output/results.json`, `RESULTS.md`, summary plots | 06, 08 |

**Physics in five lines.**
1. Two τh with pT > 40 GeV, |η| < 2.1, DeepTau Medium, opposite sign, di-τ trigger matched, lepton vetoes:
   47 586 events.
2. 81 % are jet→τh fakes: estimated from data with fake factors measured in same-sign events and applied to
   events whose leading τ fails Medium (no QCD simulation).
3. The MET is folded into the visible mass with a per-event likelihood (MET covariance + τ decay phase
   space). m_tt peaks at m_Z with 11 % resolution; fakes sit at 130–300 GeV.
4. Simulated Z→ττ (aMC@NLO) with TauPOG ID, energy-scale and trigger corrections; small backgrounds from
   simulation (tt̄, W+jets, Z→ee, dibosons).
5. μ_Z from the m_tt fit; σ = μ_Z × prediction.

**Result.**

> **σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2255 ± 36 (stat) +371/−315 (syst) ± 83 (acc) pb**
> (NNLO 1945 pb; μ_Z = 1.16 +0.19/−0.16) — **σ_fid(τhτh) = 5.22 ± 0.08 (stat) ± 0.80 (syst) pb**

Limited by the external τh identification scale factors (12.5 %) and signal modelling (9.8 %); data
statistics 1.8 %.

**Read next:** `02-selection.md` for the selection, `05-fake-factors.md` for the background method, and
`08-fit-and-results.md` for the result and what the combination should use.
