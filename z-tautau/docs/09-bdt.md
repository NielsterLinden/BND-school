# 09 — The k-fold BDT and the category fit

Code: `ztautau/bdt.py`, `scripts/step3b_bdt.py`; the categories are used by `step4_histograms.py` and
`step5_fit.py`. Numbers: `fit/bdt_info.json`, `output/data/bdt.json`; plots `output/plots/step3b_*.png`.

## Why

With DeepTau Medium on both legs and no kinematic selection, the signal region is 80 % jet → τh fakes
(S/B = 0.2). The fakes come from the application region of the data, so every fake-related uncertainty
(closure, OS/SS extrapolation, FF statistics, W+jets composition) acts on ~38 000 events and scales with B/S.
The review (`REVIEW.md` section 4) quantified three ways out; the group chose a classifier, because it
keeps the full signal and lets the fit see the fakes where they dominate and the signal where it dominates.

A classifier on **mass-agnostic** kinematics separates Z → τhτh from QCD dijets at AUC ≈ 0.97: a Z decays
back to back with the neutrinos along the visible τ's (small ΔR for a 40 GeV threshold, balanced ττ + MET
system), a dijet has the MET pointing elsewhere or nowhere and a wider opening angle. The score is used as
a **category** variable; m_tt stays the fit variable in every category, so the fake normalisation is still
measured in the high-mass sideband of every category and the Z peak is what is fitted.

## Inputs (`config.BDT_FEATURES`)

| input | why |
|---|---|
| pT(τ1), pT(τ2), pT(τ2)/pT(τ1) | signal τ's are soft and balanced, fakes are a steeply falling dijet spectrum |
| \|η(τ1)\|, \|η(τ2)\| | central Z vs forward-ish dijets (after the η closure correction of the FF, docs/05) |
| ΔR(τ1, τ2), Δφ(τ1, τ2) | the most powerful input: a 90 GeV Z with two 40 GeV τ's is close to back to back, dijets are back to back with ISR |
| MET, MET significance, Δφ(MET, ττ), pT(ττ + MET), pT(vis. ττ) | neutrinos along the τ's vs MET from mismeasurement |
| N_jets, leading-jet pT | ISR activity |
| DM(τ1), DM(τ2) | fake rates and ID efficiencies are decay-mode dependent |

**Forbidden inputs**: anything that enters the definition of the fake-factor regions. The τ1 raw DeepTau score,
the isolation sums, `leadTkPtOverTauPt`, the charge (OS vs SS), and any mass variable (m_vis, m_tt, m_col,
m_T^tot): the first group would let the classifier learn the difference between the application region and
the signal region (breaking the FF method), the charge would break the same-sign closure, and a mass input
would make the score a mass variable and the categories a mass selection.

## Training (`bdt.train`)

* **Signal**: simulated Z → ττ in the signal region, *fiducial part only* (60 < m_LHE < 120 GeV, both
  visible τ pT > 40, |η| < 2.1; `config.BDT_TRAIN_ON_FIDUCIAL`), from every available Drell-Yan sample
  (inclusive + jet-binned), positive weights. The non-fiducial high-mass tail is not a target: it goes to
  the lower categories, where it is a background with its own normalisation parameter.
* **Background**: the fake estimate itself — application-region data weighted by the nominal fake factor —
  so the classifier learns the fakes exactly as the analysis models them. (Same-sign data would be the
  alternative; it is used as the closure test instead.)
* Both classes carry the same total weight. XGBoost, depth 4, 300 trees, learning rate 0.05, subsampling
  0.8, single thread with a fixed seed (reproducible), `config.BDT_PARAMS`.
* **k = 5 folds by event number** (`event mod 5`). Model k is trained on the other four folds and applied
  to fold k of *everything*: application-region data (the fake template), same-sign data (closure),
  signal-region data, every simulated sample and every kinematic variation (τ ES, MET: the features are
  recomputed, `analysis.categories(..., kin, syst, direction)`). No event is ever scored by a model that
  saw it, which is what makes it legitimate to train on the very data events that form the fake template.
* Models: `$BND_TAUTAU_CACHE/bdt/fold<k>.ubj` + `bdt.json` (features, parameters, AUC per fold,
  importance). `bdt.load_models` refuses models trained with other settings. Saved through the `Booster`
  API (the sklearn wrapper's `save_model` fails with scikit-learn 1.9).

## Validation (`step3b_bdt.py`)

1. **Held-out AUC per fold** and train/test agreement (over-training check).
2. **Same-sign closure of the FF in the score** (`step3b_bdt_closure_SS.png`): SS_T observed vs FF × SS_L
   predicted, MC subtracted. This is *the* test that the fake factor does not depend on the score, i.e.
   that the categories do not bias the fake estimate. The residual per category and mass region becomes a
   nuisance parameter in the fit (`FakeClosure_tautau_c<k>_lo|hi`, docs/05). The OS/SS correction is
   measured per category as well (the charge correlation depends on the topology, docs/05). A first version of this plot
   showed a 75 % non-closure in the top bin caused by a single W+jets event with weight −63 in the MC
   subtraction; W+jets is therefore subtracted with uniform weights (`analysis.subtraction_weights`).
3. **Prefit signal-region score** (`step3b_bdt_SR_score[_log].png`): data vs fakes + simulation. The
   fake-dominated bins close at the percent level; the excess in the signal-like bins is the signal strength.
4. **Category yields and the stat-only sensitivity** of an m_tt fit inclusive vs in categories.

## Categories (`config.BDT_CATEGORY_EDGES`, regions `tautau_SR0/1/2`)

| region | score | role |
|---|---|---|
| `tautau_SR0` | < 0.55 | fake dominated: fixes the fake normalisation and shape (only its bins above 110 GeV are fitted, docs/08), holds most of the non-fiducial DY |
| `tautau_SR1` | 0.55–0.90 | mixed |
| `tautau_SR2` | > 0.90 | signal dominated (S/B ≈ 5): the Z peak with little background, also a clean τ energy-scale control |

The boundaries were chosen on the signal score distribution (roughly its upper tercile and the fake-dominated
bulk); they are not tuned on data. All three regions share the 14 m_tt bins of the inclusive fit.

## What it does and does not buy

* Statistically the gain is modest (the measurement is systematics limited): the stat-only δμ improves by
  ~10 %.
* The fake-related uncertainties now act on a few hundred events in the signal-dominated category instead
  of on 38 000, and the fake normalisation is measured in a region with essentially no signal.
* The τh ID, trigger and signal-modelling uncertainties are multiplicative on the signal and are *not*
  reduced by a classifier; those need the combination with ee/μμ (τh ID constrained in situ) and the
  fixes of `REVIEW.md` 3.3–3.4.
