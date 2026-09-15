# 08 — The fit and the result

Code: `scripts/step5_fit.py` (writes `fit/ztautau.config` with `fitting/trexconfig.py`, runs TRExFitter,
parses with `fitting/run_trex.py`), `scripts/step6_report.py`. Numbers: `output/RESULTS.md`,
`output/results.json`.

## Set-up

* **TRExFitter v1.8.0**, the build of the git submodule made with `fitting/build_trexfitter.sh` against the
  LCG_110 ROOT (a build against another ROOT crashes at start-up, `CLAUDE.md`).
* **Three regions** `tautau_SR0`, `tautau_SR1`, `tautau_SR2`: the BDT categories (score < 0.55, 0.55–0.90,
  > 0.90; `09-bdt.md`), each with the observable m_tt in 14 bins: 0, 40, 60, 70, …, 130, 150, 175, 200, 250,
  350 GeV (overflow in the last bin). SR0 is 90 % fakes and fixes their normalisation and shape; SR2 is a
  nearly background-free Z peak (S/B ≈ 5.3).
* Samples: Data; **DYtautau** (fiducial signal, `NormFactor mu_Z`, range 0–3); DYtautau_nonfid, DYee, DYmumu,
  DYlowmass, WJets, TTbar, SingleTop, WW, WZ, ZZ; **Fakes** (data × FF, MC subtracted). Histograms follow
  `fitting/CONVENTIONS.md` (`tautau_SR<k>__<sample>[__<syst>Up/Down]`, TH1D with Sumw2, written by uproot).
* Nuisance parameters: `07-corrections-and-systematics.md` (39 NPs). Per-bin γ for the MC and fake-template
  statistics (the `Fakes` variance includes the FF statistics; W+jets has no per-bin γ: its statistical
  uncertainty is one normalisation parameter). Shape variations of τ ES and MET use TRExFitter smoothing (`Smoothing: 40`).
* Runs: `h w f d p i r` on data, stat-only (`StatOnly=TRUE`), Asimov (`FitBlind=TRUE`), for both
  fake-factor variants.

## Result (nominal: MC-subtracted, closure-corrected fake factor)

```
mu_Z = 1.212 +0.145 -0.125      (data stat. 0.020)           expected: +0.117 −0.100
σ_fid(τhτh; vis pT > 40, |η| < 2.1, 60 < m_ττ < 120)  = 5.46 ± 0.09 (stat) ± 0.60 (syst) pb     pred. 4.50 pb
σ(pp → Z/γ* → ττ, 60 < m < 120 GeV)                   = 2358 ± 40 (stat) +282/−242 (syst+stat) ± 86 (acc) pb        pred. 1945 pb (NNLO)
```

Saturated-model goodness of fit: p = 0.01. The tension is a shape effect in the signal-dominated category:
data/prediction rises with the visible-τ pT from 1.0 at the 40 GeV threshold to ≈ 1.2 above 70 GeV
(`step4_tautau_SR2_t1_pt.png`, `step4_tautau_SR2_t2_pt.png`), which the m_tt fit sees as a distorted peak
(data/pred 0.85 at 50–60 GeV, 1.15 at 120–130 GeV in `fit_SR2_postfit.png`). It is also where part of μ > 1
comes from: at the threshold μ would be ≈ 1.05. The candidates are the trigger turn-on modelling (both legs
sit on it, the POG SF uncertainty is pT-flat per DM) and the NLO Z pT spectrum; an in-situ trigger
efficiency measurement is the first thing to do next (`CLAUDE.md`, open issues). Cross-check without the MC subtraction in the fake-factor
regions: μ_Z = 1.212 +0.144 -0.125 (GoF p = 0.01); the difference is the known double counting of the
unsubtracted estimate (`05-fake-factors.md`), not an uncertainty.

v1 (one region, unsubtracted FF, whole Z/γ*→ττ template as signal, LO/NLO NP) gave μ_Z = 1.160 +0.192 −0.163.

### Uncertainty breakdown (nominal, grouped impacts on μ_Z)

| source | impact |
|---|---:|
| Tau ID | 11.2 % |
| Fakes | 8.8 % |
| Gammas | 4.7 % |
| Tau trigger | 4.2 % |
| Tau energy scale | 2.5 % |
| Background normalisation | 1.8 % |
| Signal modelling | 1.8 % |
| MET | 1.2 % |
| Luminosity | 1.0 % |
| Pileup | 0.9 % |
| L1 prefiring | 0.3 % |
| **data statistics** | **2.0 %** |

**This is a systematics-dominated measurement, and the leading systematic is external.** The TauPOG τh ID
scale factors for UL2016 postVFP carry 5–16 % per decay mode, entering squared for two τ, and the two POG
prescriptions differ by 14 % on the yield. Z→ττ is the process used to *measure* those scale factors, so
a τhτh-only cross section cannot do better than the scale factors it borrows. The fake-related uncertainty
is now what the same-sign closure per category says it is, not what a single over-constrained parameter
allowed.

### Pulls and constraints worth knowing

The pulls (`fit_pulls.png`) and the ranking (`fit_ranking.png`): `FakeOSSS_tautau` is pulled +0.64 σ, i.e. the
fake-dominated category measures the fake normalisation itself; the closure parameters of SR0 stay at their
prior. The post-fit m_tt distributions per category are `fit_SR0_postfit.png`, `fit_SR1_postfit.png`,
`fit_SR2_postfit.png`; the summary of all three regions `fit_summary_postfit.png`.


## Working-point cross-check: DeepTau Tight on both legs

The review (section 4.1) noted that Tight removes 63 % of the fakes for 24 % of the signal at no cost in
statistical precision. The whole chain was rerun with `BND_TAUTAU_WP=Tight` (`variants/tight/`): the
fake factors, closure corrections and C are re-derived with T = Tight, L = VVVLoose & !Tight, the BDT is
retrained, the TauPOG Tight ID and trigger scale factors are used (`external/*.json`).

| quantity | Medium (nominal) | Tight |
|---|---:|---:|
| μ_Z | 1.212 +0.145 −0.125 | 1.099 +0.116 −0.102 |
| total uncertainty on μ_Z | +14.5 / −12.5 % | +11.6 / −10.2 % |
| expected (Asimov) | +11.7 / −10.0 % | +10.4 / −9.1 % |
| data statistics | 2.0 % | 2.1 % |
| σ(60–120) [pb] | 2358 +282 −242 | 2137 +226 −198 |
| goodness of fit p | 0.012 | 0.317 |
| impact: Tau ID | 11.2 % | 8.6 % |
| impact: Fakes | 8.8 % | 7.0 % |
| impact: Gammas | 4.7 % | 3.9 % |
| impact: Tau trigger | 4.2 % | 3.1 % |
| impact: Tau energy scale | 2.5 % | 1.7 % |
| impact: Signal modelling | 1.8 % | 1.4 % |
| impact: Background normalisation | 1.8 % | 2.5 % |
| prefit fiducial signal / fakes (all categories) | 4827 / 36931 | 3759 / 13569 |
| SR2: signal / fakes | 2790 / 529 | 2141 / 186 |
| BDT held-out AUC | 0.970 | 0.966 |
| C_OS/SS (inclusive) | 1.052 | 1.056 |
| closure NPs (c0/lo … c2/hi) | 15 %, 1 %, 5 %, 7 %, 12 %, 55 % | 14 %, 1 %, 7 %, 21 %, 15 %, 44 % |
| τh ID SF per DM (0/1/10/11) | 0.92±0.14, 0.88±0.05, 0.87±0.09, 0.90±0.16 | 0.90±0.13, 0.89±0.05, 0.94±0.15, 0.81±0.15 |

**Reading.** The Tight working point gives a *tighter* result (total +11.6/−10.2 % instead of
+14.5/−12.5 %): every group shrinks, the τh ID because the Tight SF uncertainties are smaller for the
decay modes that carry most of the signal (DM1, and DM0 slightly), the fakes because there are 2.7× fewer
of them, the MC statistics because the BDT categories are purer. It also has an acceptable goodness of fit
(p = 0.32 instead of 0.01) and a central value 1.4 σ closer to NNLO and to Z→μμ; the visible-pT slope in
the signal-dominated category is somewhat weaker with Tight (data/pred 1.0 → 1.15 instead of 1.0 → 1.25,
`variants/tight/output/plots/step4_tautau_SR2_t1_pt.png`), so part of it is fake- or WP-related and part
remains (trigger turn-on or Z pT modelling). The two results share
the data and most systematics, so their difference (Δμ = 0.11, i.e. 1 σ of either) is not an additional
uncertainty; it is a cross-check that the Medium result passes at the 1 σ level and that identifies the
better operating point. **Recommendation for the next iteration and the combination: adopt Tight as the
nominal working point** (this commit keeps Medium nominal so that every number, plot and page refers to one
configuration; switching is one environment variable and 25 minutes).

## For the combination

`fit/ztautau.config` + the workspace in `fit/results/ztautau/RooStats/` enter `fitting/combination_skeleton.config`
directly (job name `ztautau`, POI `mu_Z`, regions `tautau_SR0/1/2`). In a combined fit with ee and μμ,
μ_Z is fixed by the light leptons to ~1–2 %. The ττ channel then mostly **constrains the τh ID and trigger
nuisance parameters** rather than μ_Z. That is physically sensible (lepton universality) and a useful
by-product. `SigModel_tautau` is not in the workspace by default and must not be correlated with the μμ
`SigModel` if it is switched on.

For a BLUE-style combination (`combination/combination.ipynb`): n_obs = 47 586, n_bkg (prefit, including the
non-fiducial DY) = 41459, A·ε = 0.0001514 for σ(60–120). **Do not use the counting numbers**: with B/S ≈ 8
the result is dominated by the fake normalisation, which only the shape fit in the categories constrains;
use σ = 2358 +282 −242 pb (τh ID and fakes uncorrelated with ee/μμ; lumi 1.2 % correlated; acceptance 3.7 %
correlated).
