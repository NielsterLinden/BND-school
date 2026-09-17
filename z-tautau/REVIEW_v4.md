# Review of the v4 four-channel Z → ττ measurement (τhτh + μτh + eτh + eμ)

Reviewer: Claude (Fable 5.1), 17 September 2026, at the request of Samuel Jankovych.
Reviewed state: commits `2a989a0` … `cb3753a` (v4 measurement, ranking and report, per-channel exports and MultiFit, deck).

What was reviewed: `docs/10-v4-plan.md`, `output_v4/RESULTS.md`, `handoff.md`, `ztautau/analysis_v4.py`,
`scripts/step4_histograms_v4.py`, `scripts/step5_fit_v4.py`, `fit_v4/*.config`, the committed fit inputs
`fit_v4/fitinputs/ztautau_v4.root`, the fit results and logs in `fit_v4/results/`, the in-situ trigger tables
`external/trigger_insitu_v4.json`, the post-fit plots, and the shared `fitting/CONVENTIONS.md` and
`combination/` (`docs/01-inputs.md` is the recipe for the uncertainty check of section 2). Five cross-check fits
of the committed model were run with TRExFitter v1.8.0 (section 3) and the three-channel combination was re-run
with the v4 input (section 5). Scripts and printouts: `review/v4_uncertainty_checks.{py,txt}`,
`review/v4_combination.{py,txt}`, `review/v4_fit_variants.py`. Nothing in `output_v4/` or `fit_v4/` was modified.

---

## 1. Verdict in short

The v4 chain is a large and careful piece of work: three new channels with their own data streams, a per-process
fake-factor method that found two real effects (the same-sign region is half W+jets; the W+jets fake factor is
charge-correlated), in-situ trigger efficiencies, per-decay-mode regions so that the τh ID scale factors and
energy scales are measured by the data, and per-channel exports whose MultiFit reproduces the single-file fit.
The fit inputs follow `fitting/CONVENTIONS.md` §3 exactly where z-ee does not: every theory variation is
renormalised to the generator yield in 60 < m_LHE < 120 GeV, so no nuisance parameter scales the prediction
that μ_Z multiplies (section 2). The fit converges, MINOS succeeds, the goodness of fit is p = 0.085.

The headline number, however, is not what the report says it is. Three findings change what the result means:

| # | severity | finding | effect |
|---|---|---|---|
| 1 | **critical** | **The quoted expected uncertainty (±0.009) is wrong.** In the Asimov fit MINOS failed on `mu_Z` ("Invalid lower error", Hessian "forced pos-def by adding to diagonal 0.185"); the ±0.009 is the HESSE error of a non-positive-definite matrix. The MINOS side that did converge gave +0.043. The four per-channel Asimov fits are fine (MINOS status 0, errors equal to the observed ones). | `RESULTS.md`, `handoff.md`, `docs/10` §12 quote "expected ±0.009" against an observed ±0.045; the real expected uncertainty is ≈ ±0.044. `results/ztautau_v4/CorrelationMatrix.yaml`, `CorrMatrix.png` and the `*_asimov.png` files are those of that failed fit (the Asimov `wf` ran last and overwrote them). |
| 2 | **critical** | **The four-channel μ_Z = 1.055 ± 0.045 is a compromise between two measurements that disagree by 2.2σ**, and the fit pays for it with the eμ trigger scale factor. The eμ channel alone gives μ_Z = 0.978 +0.056 −0.053 (no τh, no ID scale factor). The three τ channels alone, with the ID scale factors free, give μ_Z = **1.210 +0.091 −0.084** (section 3): with free scale factors the τhτh regions measure SF² μ and the ℓτh regions SF μ, so μ comes from the ratio τhτh / (ℓτh)². In the combined fit the eμ prediction is lowered by 9.7 % through nuisance parameters (`EmuTrigger` −1.62σ = −8.7 %, constrained to 0.71), the ℓτh predictions by 3–8 % (`MET_Unclustered` +1.05σ constrained to 0.32, `Pileup`, `ElectronTrigger`), the τhτh prediction raised by 2–3 % (`QCDScale` +0.5σ). `mu_ttbar` = 1.15 is a by-product: the control region fixes μ_tt̄ × ε_trig, the signal region μ_Z × ε_trig, so the eμ channel really measures μ_Z/μ_tt̄ and the τ channels push ε_trig down. | The statement "the four channels agree with each other" (`docs/10` §12) holds only with the POG scale factors fixed (1.04 / 1.05 / 1.02 / 0.98). With the scale factors free, which is the point of v4, the τ channels and eμ do not agree. The 2.4× gain in precision over v3 is real only if the tension is understood. |
| 3 | **major** | **The eμ trigger prior is 5.4 % and it sets the answer.** `analysis_v4._insitu_sf` multiplies each leg by (SF ± stat) × (1 ± 2 %) with the same sign on both legs, so the 2 % of the paper (per channel, Table 2) enters twice and linearly: the template moves by ±5.4 % (plus `ElectronID` 1.2 %, `ElectronReco` 1.4 %, `MuonID/Iso` 1 % ⇒ ≈ 5.9 % on the eμ normalisation against 2.8 % in the paper). With the variation scaled by 0.5 (2.7 %, i.e. the 2 % applied once) the fit gives **μ_Z = 1.008 +0.034 −0.033** (σ = 1961 pb), `EmuTrigger` −1.31σ, μ_tt̄ = 1.10, p = 0.06. The central value moves by one full standard deviation with a change of the prior that is arguably the correct one. | The result is prior-dominated; the quoted uncertainty does not contain this. |
| 4 | major | **The τhτh / (ℓτh)² lever assumes the same τh ID scale factor at pT > 30 GeV (ℓτh) and pT > 40 GeV (τhτh).** Most of the ℓτh signal sits at 30–45 GeV. The TauPOG DM-binned and pT-binned prescriptions already differ by 7 % per leg (`REVIEW.md` 3.3, `docs/07`); a 5 % lower scale factor at 30–40 GeV would make the lever's μ_Z ≈ 10 % too high, which is the size of the effect in finding 2. `docs/10` §11 names this as a limitation but nothing tests it. | the single model assumption that decides whether finding 2 is physics or an artefact |
| 5 | major | **The grouped impacts over-shoot the total by a factor 1.6** (quadrature sum of the categories + stat = 7.2 % vs MINOS 4.5 %): `NormFactors` 4.0 % and `Electron efficiency` 3.9 % are the same degeneracy (μ_tt̄ ↔ `EmuTrigger` ↔ μ_Z) counted twice. Fixing μ_tt̄ at ±1σ moves μ_Z by ∓0.04 and halves its error, so the ranking entry (4.1 %) is right, but "NormFactors" is not the τh ID scale factors, it is the eμ normalisation. | in the covariance combination (rule: keep the over-shoot) v4 enters with ±141 pb instead of ±88 pb; either way its weight is ≤ 7 % and the three-channel result moves by < 5 pb (section 5) |
| 6 | minor | `FakeFrac_etau` is pulled +2.6σ (constraint 0.62) and `FakeWmT_etau` +1.6σ: the eτh fake composition wants to move by 2–3 times its prior (−1.3 … −6.5 %), most in DM10/DM11. The eμ channel alone has p(GoF) = 0.036. | the eτh fake prior is too tight or the composition is off; the eμ shape has a mild problem |
| 7 | minor | eτh electrons at pT > 29 GeV sit on the `Ele27_WPTight` turn-on: data efficiency in the 29–35 GeV bin is 0.62 (barrel) and 0.42 / 0.28 (endcap), scale factors 0.91 / 0.97 / 0.78 / 0.57, covered by one ±3 % normalisation NP (`ElectronTrigger`, pulled −0.8σ, ρ(μ_Z) ≈ +0.3). | a pT-dependent modelling error is treated as a flat 3 % |
| 8 | minor | `MET_Unclustered` changes the ℓτh signal *normalisation* by −3.8 % / +2.9 % (m_T < 40 GeV cut and the likelihood mass), is pulled +1.05σ and constrained to 0.32: a shape systematic used as a normalisation knob. | 1.5 % post-fit impact, and part of the mechanism of finding 2 |
| 9 | minor | Report hygiene: `RESULTS.md` prints `mu_ttbar` as a raw list and `stat 0.00456806`; the "expected" row (finding 1); `output_v4/plots/fit_corrmatrix.png` should be checked to be the observed-data matrix, not the Asimov one. | cosmetic / traceability |

For the combination (section 5): v4 does **not** change the three-channel result today (1980 → 1981–1985 pb,
χ²/ndf unchanged at ≈ 10/2), whichever v4 variant is used, because its weight is small and the μμ–ee discrepancy
dominates the χ². What v4 changes is the ττ line: 2053 +90 −86 pb (four channels), 1902 +108 −102 pb (eμ alone)
or 2354 +178 −164 pb (τ channels alone), against μμ 1931 ± 35 and ee 2054 ± 40. The eμ-alone number sides with
μμ, the τ-lever number with nothing.

---

## 2. The uncertainty check (`combination/docs/01-inputs.md`, `CONVENTIONS.md` §3)

The z-ee problem was a theory variation that changed the *total* signal yield by the cross-section uncertainty and was
therefore degenerate with the POI: the fit pulled it, and μ̂ was measured against a prediction the fit had moved. The
recipe: (1) the prefit normalisation effect of every systematic on the signal template, (2) the pulls and constraints,
(3) which pulls μ̂ absorbed, and whether they scale the *denominator* (luminosity, efficiency, acceptance — legitimate)
or the *prediction* (not legitimate). Full tables: `review/v4_uncertainty_checks.txt`.

### 2.1 Template normalisation at ±1σ

`analysis_v4.theory_weights` divides every scale / PDF / PS member by its own Σw over the generated 60 < m_LHE < 120 GeV
Z→ττ events (`gensums["scale_lhe_tautau_60_120"]` etc.), so each member has the same σ(60–120) and only A × ε per
region varies; `combine_theory` builds the envelope / Hessian sum **on histograms**. The templates confirm it:

| systematic (signal `DYtautau`) | τhτh SR0/1/2 | μτh dm0/1/10/11 | eτh dm0/1/10/11 | eμ SR |
|---|---|---|---|---|
| `QCDScale` | +5.5/−5.7, +5.8/−5.2, +4.9/−4.9 | 1.0–1.5, dm11 2.8 | 1.8–2.8 | +0.7/−0.8 |
| `PDF` | ±1.5, ±1.4, ±1.3 | ±0.5–0.7 | ±0.6–0.8 | ±0.5 |
| `PS_ISR` | ≤ 1.5 | ≤ 1.0 | ≤ 1.1 | 0.4 |
| `PS_FSR` | +0.7/−1.3, +1.0/−2.9, +1.3/−1.7 | ≤ 1.9 | ≤ 2.3 | 0.5 |

These are acceptance-sized (the τhτh double 40 GeV cut selects the Z pT tail: 5–6 %, as v3's 3.4 % on A), not the
5.9 % total-cross-section envelope of z-ee; in eμ, where the acceptance is almost inclusive, `QCDScale` is 0.7 %.
κ_theory ≡ 1 by construction. **The v4 inputs pass the z-ee test.** The variations on `DYtautau_out` (15–21 % in
τhτh SR1/2) are the scale dependence of the out/in ratio, as in v3.

Experimental normalisation effects worth knowing (signal, per region): `EmuTrigger` ±5.4 % (eμ, finding 3);
`ElectronTrigger` ±3.0 % (eτh); `MET_Unclustered` −3.8/+2.9 % (ℓτh, finding 8); `Pileup` −1.8…−3.0 / +1.7…+2.7 %
(ℓτh, asymmetric); `TauES_DM1` +33/−21 % on the signal in τhτh SR0 (the fake sideband above 110 GeV: migration of the
signal tail); `TauFakeMu` ±21 % on Z→μμ in μτh; `TauFakeEle` ±15 % on Z→ee in eτh. `BTag` ±4.7 % on tt̄ in the eμ SR.

### 2.2 Pulls and constraints

|pull| > 1 or constraint < 0.5 (all others are within ±0.8σ and above 0.6):

| NP | pull ± constraint | what it moves |
|---|---|---|
| `FakeFrac_etau` | +2.60 ± 0.62 | eτh fakes −1.3…−6.5 % per σ |
| `EmuTrigger` | −1.62 ± 0.71 | everything in eμ, −8.7 % |
| `FakeWmT_etau` | +1.57 ± 0.62 | eτh fakes +6…+10 % per σ |
| `FakeFrac_mutau` | +1.48 ± 0.68 | μτh fakes +2…+7 % per σ |
| `TopPt` | +1.17 ± 0.30 | tt̄ shape |
| `TauES_DM11` | +1.09 ± 0.73 | 3 % prior → +3.3 % |
| `FakeClosure_tautau_c2_hi` | +1.08 ± 0.87 | τhτh SR2 fakes above 110 GeV |
| `MET_Unclustered` | +1.05 ± 0.32 | ℓτh signal −3.8 % |
| `TauFakeMu` | +1.02 ± 0.49 | Z→μμ in μτh +21 % |
| `Pileup` | +0.25 ± 0.47 | ℓτh signal −2 % per σ |
| `PS_FSR` | +0.39 ± 0.39 | τhτh shape |
| `QCDScale` | +0.50 ± 0.58 | τhτh acceptance +2.8 %, ℓτh +0.7 % |
| `mu_ttbar` | 1.153 ± 0.048 | tt̄ everywhere |
| `TauIDSF_DM0/1/10/11` | 0.962 / 0.935 / 0.867 / 0.771 (± 0.04–0.05) | POG 0.902 / 0.892 / 0.938 / 0.812 |

None of the pulled parameters scales σ^pred, so σ̂ = μ̂ σ^pred is the right formula (unlike z-ee). But `QCDScale`,
`PS_FSR`, `Pileup`, `MET_Unclustered` are constrained to 0.3–0.6 of their priors by a fit that has no region where these
are the leading effect: they are being used as normalisation knobs between the channels (2.3).

### 2.3 Which pulls μ̂ absorbed (post-fit / pre-fit signal yield per region)

| region | post/pre | /μ_Z | /μ_Z/(SF/POG) | reading |
|---|---|---|---|---|
| τhτh SR0 / SR1 / SR2 | 1.09 / 1.08 / 1.08 | 1.03 / 1.02 / 1.02 | (SF products) | prediction raised 2–3 % beyond μ: `QCDScale`, SF² |
| μτh dm0 / dm1 / dm10 / dm11 | 1.07 / 1.08 / 0.97 / 1.06 | 1.02 / 1.02 / 0.92 / 1.00 | 0.95 / 0.97 / 1.00 / 1.05 | NPs lower the prediction 0–5 % |
| eτh dm0 / dm1 / dm10 / dm11 | 1.04 / 1.03 / 0.93 / 1.01 | 0.98 / 0.97 / 0.88 / 0.96 | 0.92 / 0.93 / 0.95 / 1.01 | NPs lower the prediction 5–8 % |
| eμ SR | 0.953 | **0.903** | — | NPs lower the prediction 9.7 %; data are 4.7 % *below* the prefit prediction |
| eμ CRtt | 1.20 | 1.14 | — | μ_tt̄ |

Every one of these moves is of the "denominator" kind (efficiency, acceptance), so the rule of `01-inputs.md` says μ̂
may absorb them. The pattern, however, is one-directional: all three lepton channels are pulled *down* by their
efficiency parameters, the τhτh channel *up*, which is what a fit does when the τhτh / (ℓτh)² lever wants a higher μ
than the eμ normalisation allows. That is finding 2, quantified in section 3.

---

## 3. Cross-check fits of the committed model (`review/v4_fit_variants.py`, TRExFitter v1.8.0, MINOS status 0 in all)

| fit | μ_Z | σ(60–120) [pb] | SF DM0 / DM1 / DM10 / DM11 | `EmuTrigger` | μ_tt̄ | GoF p |
|---|---|---|---|---|---|---|
| **nominal, four channels** | **1.055 +0.046 −0.044** | 2053 +90 −86 | 0.96 / 0.94 / 0.87 / 0.77 | −1.62 ± 0.71 | 1.153 ± 0.048 | 0.085 |
| eμ alone (`ztautau_v4_emu`) | 0.978 +0.056 −0.053 | 1902 +108 −102 | — | | | 0.036 |
| τhτh + μτh + eτh, SFs free (no eμ) | **1.210 +0.091 −0.084** | 2354 +178 −164 | 0.88 / 0.86 / 0.79 / 0.70 | — | | |
| four channels, μ_tt̄ fixed at 1.2016 (+1σ) | 1.097 ± 0.021 | 2133 | | −2.20 | fixed | |
| four channels, μ_tt̄ fixed at 1.1052 (−1σ) | 1.015 ± 0.020 | 1975 | | −1.03 | fixed | |
| four channels, `EmuTrigger` variation × 0.5 (2.7 %) | **1.008 +0.034 −0.033** | 1961 +67 −65 | 1.00 / 0.97 / 0.90 / 0.80 | −1.31 ± 0.89 (= −3.5 %) | 1.101 ± 0.035 | 0.062 |
| per channel, POG SFs fixed (`*_fixedid`, from `RESULTS.md`) | ττ 1.043, μτ 1.049, eτ 1.017, eμ 0.978 | | | | | |
| v3 (τhτh alone, POG priors) | 1.071 +0.114 −0.100 | 2082 | | | | 0.22 |

Reading the table:

* **eμ vs τ channels: (1.210 − 0.978) / √(0.088² + 0.055²) = 2.2σ.** The nominal 1.055 is their weighted compromise, and
  its ±0.045 is smaller than either input because the fit lets the eμ trigger efficiency and the ℓτh acceptance
  parameters move to reconcile them.
* **μ_tt̄ ↔ `EmuTrigger` ↔ μ_Z.** Both eμ regions scale with ε_trig; the control region fixes μ_tt̄ · ε_trig and the
  signal region μ_Z · ε_trig, so the eμ channel measures the *ratio* μ_Z/μ_tt̄ = 0.915 to ±2 % and leaves the absolute
  normalisation to the 5.4 % prior. Fixing μ_tt̄ at +1σ pushes `EmuTrigger` to −2.2σ and μ_Z to 1.10; the CMS tt̄
  cross section is not 15 % off, the trigger parameter is. (The ±0.021 with μ_tt̄ fixed also says that "NormFactors
  4.0 %" is this chain, not the τh ID scale factors.)
* **The prior decides.** Halving the trigger variation moves the central value by −0.047, one full σ, and the τh ID
  scale factors up by 0.03–0.04 (DM0 to 1.00). Whether 5.4 % or 2.7 % is right is a choice the analysis has to defend;
  the paper's Table 2 has 2 % for "electron ID and trigger" and 2 % for "muon ID and trigger", i.e. about 2.8 %
  for the eμ efficiency in total, against ≈ 5.9 % here.

---

## 4. Why the τ-lever may be biased: the pT range of the scale factors (finding 4)

The lever works because the same `TauIDSF_DM*` multiplies the ℓτh templates once and the τhτh templates twice. The
ℓτh τh has pT > 30 GeV and its spectrum peaks at 30–45 GeV; the τhτh legs have pT > 40 GeV on the trigger turn-on.
If SF(30–40) is lower than SF(> 40) by δ, the ℓτh regions want SF·μ lower by δ, the τhτh regions do not, and the
fit returns μ higher by ≈ 2δ / (1 + ℓτh weight …), of order δ–2δ. The TauPOG pT-binned prescription for pT > 40 GeV is
already 7 % below the DM-binned one (`docs/07`), and the POG's own pT-binned SFs are not flat between 20 and 40 GeV.
The fitted DM11 scale factor (0.77, and 0.70 with the τ channels alone) and the τES pull of DM11 (+1.1σ, anticorrelated
with its SF at −0.53) point the same way: the parameters are absorbing a pT dependence.

The test is cheap relative to what was built: split the ℓτh regions (or only the μτh ones) at pT(τh) = 40 GeV, and
either fit one SF per decay mode from the pT > 40 GeV regions of all channels (letting the 30–40 GeV regions float with
their own scale factors) or simply raise the ℓτh τh threshold to 40 GeV for the scale-factor constraint. If the τ-only μ_Z
then drops towards 1.0, finding 2 is an artefact of the assumption; if it stays at 1.2, the eμ channel (or the τhτh
trigger and fake modelling) has to be looked at instead.

---

## 5. In the light of the full combination (`review/v4_combination.py`, the combination's own BLUE and ρ model)

The combination currently reads v3 (`z-tautau/output/results.json`, weight −0.0003). The v4 results file has a different
layout (no `nominal_variant`, no `prediction.A_unc`: the acceptance is inside the fit, as for z-ee), so a v4 loader
would be needed; the numbers below construct the `ChannelResult` by hand with `acc = {}` and ρ = 0 for the categories the
model does not know (`NormFactors`, `Electron energy`, `Jets`, `Background modelling`, `b tagging`).

| ττ input | ττ σ [pb] | combination [pb] | χ²/ndf, p | weights μμ / ττ / ee |
|---|---|---|---|---|
| v3 τhτh (baseline of `result.md`) | 2082 ± 253 | 1980.3 ± 31.6 | 9.67/2, 0.008 | 0.596 / −0.000 / 0.404 |
| v4 four channels, sized by the combination's rule (groups in quadrature, ±141) | 2053 ± 141 | 1981.5 ± 31.5 | 9.66/2, 0.008 | 0.587 / 0.017 / 0.396 |
| v4 four channels, sized by its MINOS error (±88) | 2053 ± 88 | 1985.4 ± 30.9 | 10.02/2, 0.007 | 0.554 / 0.074 / 0.372 |
| v4 eμ alone (MINOS ±105) | 1902 ± 105 | 1976.2 ± 31.2 | 10.10/2, 0.006 | 0.571 / 0.048 / 0.382 |
| v4 τ channels alone (MINOS ±170) | 2354 ± 171 | 1982.8 ± 31.6 | 14.10/2, 0.001 | 0.592 / 0.007 / 0.401 |
| v4 four channels @ MINOS, ee normfix | 2053 ± 88 | 1892.8 ± 29.4 | 10.61/2, 0.005 | 0.469 / 0.064 / 0.467 |
| v4 eμ alone @ MINOS, ee normfix | 1902 ± 105 | 1882.8 ± 29.6 | 6.91/2, 0.032 | 0.482 / 0.041 / 0.478 |
| μμ + v4 four channels (no ee) | | 1934.6 ± 35.1 | 0.75/1, 0.39 | |
| μμ + v4 eμ alone (no ee) | | 1930.0 ± 35.1 | 0.04/1, 0.85 | |

* **v4 does not move the combination**: +1 to +5 pb on 1980, with χ²/ndf ≈ 10/2 in every case, because the μμ–ee
  discrepancy (the z-ee `QCDScale` normalisation, `ASK_Z_EE.md`) dominates and the ττ weight stays ≤ 7 %.
* **v4 does change the universality picture.** The eμ-alone value (1902 ± 105) agrees with μμ (1931 ± 35) to 0.3σ;
  the four-channel value (2053 ± 88) sits with ee; the τ-lever value (2354 ± 171) agrees with neither (2.4σ from μμ).
  Which ττ number is quoted should wait for the resolution of findings 2–4.
* **A TRExFitter MultiFit of v4 with z-mumu / z-ee** (`comb_v4.config` + `Fit:` blocks) needs, beyond what `docs/10` §13
  says: `DecorrSysts` for `MuonID/Iso/Trigger/Scale` and `ElectronID/Reco/Scale` unless the groups agree that the
  POG-json (v4) and in-house tag-and-probe (z-mumu) / z-ee measurements are the same nuisance parameter — they are
  not; `mumu_CRemu` dropped (shares events with `emu_SR`); a decision on tt̄ (`XS_TTbar` 6 % in μμ vs the free `mu_ttbar`
  here); and the awareness that ττ would constrain the shared `QCDScale` to 0.58 of its prior through the τhτh/ℓτh
  acceptance ratio, which then tightens the μμ and ee acceptance uncertainties — legitimate only if that constraint is
  believed (finding 2).
* The grouped-impact over-shoot (finding 5) means `comb/inputs.ChannelResult` would carry v4 at ±141 pb by the
  combination's conservative rule; the MINOS ±88 pb is the honest total. Both are in the table.

---

## 6. Recommendations, in order

1. **Fix the expected uncertainty** (finding 1): rerun the Asimov fit with MINOS (strategy 2 or `FitBlind` after `f`)
   and quote it; if MINOS fails, quote nothing rather than the HESSE number; regenerate `RESULTS.md`, `handoff.md`,
   `docs/10` §12 and the deck. Re-run `f` (or `d`) last so the correlation matrix and plots are the observed-data ones.
2. **Report the tension** (finding 2): quote eμ alone (0.978 +0.056 −0.053) and the τ channels alone with free scale
   factors (1.210 +0.091 −0.084) next to the combined number, say that they differ by 2.2σ, and replace "the four channels
   agree" by that statement. The combined number is not wrong as a likelihood result, but it is not a measurement
   of μ_Z independent of the eμ trigger prior.
3. **Apply the 2 % trigger uncertainty once per event** (finding 3), on top of the statistical error of the in-situ
   tables (≈ 2.3–2.7 % in total), or justify 5.4 %; quote both fits until then (1.008 vs 1.055).
4. **Test the scale-factor pT assumption** (finding 4, section 4) before the τ-lever number is used for anything.
5. Revisit the eτh fake composition prior (`FakeFrac_etau` +2.6σ) and the electron trigger turn-on (raise the eτh
   electron pT to ≥ 32–33 GeV, or bin `ElectronTrigger` in pT); look at the eμ shape (p = 0.036 alone; the 130–150 GeV
   bin noted in `docs/10` §12).
6. Split `MET_Unclustered` into shape and normalisation, or understand why the unclustered-energy variation changes
   the ℓτh acceptance by 3.8 % asymmetrically (m_T cut).
7. For the combination: no action needed now; when v4 goes in, add a v4 loader with `acc = {}` and ρ entries for the
   new categories, and decide the Muon*/Electron* correlation question above.
8. Hygiene: `mu_ttbar` formatting in `RESULTS.md`; state in `docs/10` §12 that "NormFactors 4.0 %" is the
   μ_tt̄–trigger chain, not the τh ID scale factors; update the `CLAUDE.md` open-issues list (the first two items are
   done by v4, the τ trigger in situ is not).

## 7. What is fine and should not be touched

The theory-template renormalisation (2.1); the per-decay-mode regions and the `Expression` products; the τES
prior and its in-situ constraint (0.6–1.0 %, DM11 2.2 %); the per-channel exports and the MultiFit closure; the
fake-factor findings of `docs/10` §5.1 and the same-sign closures (0.97 ± 0.01, 1.03 ± 0.02); the veto-based
orthogonality to z-mumu / z-ee and the `mumu_CRemu` caveat; the GoF of the combined fit (0.085); MINOS on the observed
data (status 0). The luminosity, cross-section and sample conventions match `fitting/CONVENTIONS.md`.
