# Review of the Z → τhτh cross-section measurement

Reviewer: Claude (Fable 5.1), 15 September 2026, at the request of Samuel Jankovych.
Reviewed state: commit `6a45615` ("z-tautau: complete Z -> tau_h tau_h cross-section measurement").

What was reviewed: every file in `ztautau/`, `scripts/`, `docs/`, `fit/*.config`, `output/results.json`,
`output/RESULTS.md`, `handoff.md`, the committed plots, the TauPOG inputs in `external/`, and the shared
`fitting/CONVENTIONS.md` and z-mumu conventions this channel has to match. The ntuples in
`$BND_TAUTAU_CACHE/ntuples_v1` were used to re-derive the nominal fake factors (reproduced exactly:
38 758 fake events) and to run three short numerical studies on the fake background. The study scripts,
their printouts and plots are in `review/` (`fake_studies.py`, `eta_study.py`, `*.txt`, `review_bdt.png`,
`review_eta.png`); nothing in `output/` or `fit/` was modified.

---

## 1. Verdict in short

The analysis is complete, well engineered and unusually well documented. The chain reproduces, the
fake-factor method is implemented correctly for what it claims to do, the MET-likelihood mass is sound,
and the fit is set up in the shared conventions. The result,

    σ(pp → Z/γ* → ττ, 60 < m < 120 GeV) = 2255 ± 36 (stat) +371/−315 (syst) ± 83 (acc) pb,  μ_Z = 1.16 +0.19/−0.16,

is compatible with the z-mumu result (1935 ± 30 pb, μ_Z = 0.990) within one standard deviation. That is
the only place where the two channels can be compared, and they agree.

However, several *choices* are wrong or at least not defensible as they stand, and two of them move the
central value or change what it means:

| # | severity | finding | effect |
|---|---|---|---|
| 1 | **critical** | 38 % of the selected Z→ττ is outside the fiducial volume, 30 % has m_LHE > 120 GeV, and 31 % of the signal template sits in the m_tt > 130 GeV "fake sideband". μ_Z scales all of it. | the quoted σ(60–120) is a measurement of the whole DY spectrum weighted by the SM shape; the sideband that fixes the fakes is 7.5 % signal |
| 2 | **critical** | the nominal fake factor does not subtract genuine τ from the determination/application regions, although the team knows this double-counts signal | central value biased low by Δμ = −0.05 (a third of the systematic uncertainty); fakes over-predicted by 7 % in the purest fake phase space |
| 3 | major | τh ID scale factors: the dominant uncertainty (12.5 %), and the two TauPOG prescriptions (DM-binned 0.87–0.92, pT-binned 0.82) differ by 7 % per leg = 14 % on the yield, more than any nuisance parameter | the choice of prescription is a hidden systematic; only the combination can constrain it |
| 4 | major | `SigModel` = madgraph LO vs aMC@NLO, 7.3 % one-sided normalisation, is (a) not a credible alternative model, (b) double counts `QCDScale`/`PS_*`, (c) shares its name (and hence is correlated) with z-mumu's powheg-vs-aMC@NLO NP that means something else and is 0.2 % | second-largest uncertainty (9.8 %) is an artefact of the choice of alternative generator |
| 5 | major | fake uncertainty model: one non-closure NP links the ±2.5 % high-mass bins to the ±11 % peak bins and is constrained to 0.32 σ; FF statistics as four fully-coherent NPs; the ±15 % η(τ1) non-closure is a real ±15 % |η| dependence of the FF (diagnosed in section 5; no effect on μ, but every η plot is wrong) and the −7 % pT(τ2) non-closure is only "absorbed" | the quoted fake uncertainty (0.8 %) is not credible; realistic is a few % |
| 6 | major | MC statistics 4.8 %, third-largest, while the jet-binned (35577/35595/35613), high-mass (35629) and pT_Z-binned DY samples are already on dCache | avoidable |
| 7 | minor | both τ sit on the trigger turn-on (data efficiency 0.3–0.6 at 40 GeV); raising the threshold does not reduce the SF uncertainty (checked), so the 5.1 % is external | in-situ constraint needed (combination, or μτh T&P) |
| 8 | minor | trigger SFs applied to jet legs of W+jets/tt̄, W+jets fake composition not in the FF, no PU-jet ID in N_jets, trigger filter bit to be verified, EWK Z→ττ absent | each ≲ 1 % |
| 9 | minor | an untracked second TRExFitter build (`z-tautau_build_trex.sh`, `asetup StatAnalysis,0.8.0`) exists at the repository root, against CLAUDE.md's "never build another version" | reproducibility across checkouts |

The large fake fraction (81 %) is a consequence of the loosest reasonable working point (DeepTau Medium
on both legs) and of no kinematic selection at all. Section 4 quantifies three ways to reduce it: a tighter
working point (Tight: −63 % fakes, −24 % signal, unchanged statistical sensitivity), a k-fold BDT on
mass-agnostic kinematics (AUC 0.89; top quartile of the signal at S/B ≈ 10), and decay-mode categories.
An ML classifier is worthwhile but, honestly, it addresses the fake-related uncertainties and not the two
that dominate (τh ID, signal modelling); those need the fixes in findings 3, 4 and 6 and the combination.

---

## 2. What is good (keep it)

* **Structure and documentation.** One physics choice per constant in `config.py`, one standalone script per
  step, `docs/00–08` that explain *why*, a CLAUDE.md pitfall list that saves the next person days, and a
  300 MB laptop bundle. I could reproduce the fake-factor yield to the event in 6 s from the ntuples.
* **Statistical independence and orthogonality.** Tau primary dataset (independent of SingleMuon/Electron),
  explicit e/μ vetoes so that eτh/μτh/ee/μμ never overlap. Required for the combination and often forgotten.
* **The MET-likelihood mass.** Correct scale (0.99), 11 % resolution, always defined, per-event MET
  covariance, sensible prior, and an honest statement that it separates signal from fakes no better than
  m_vis here. Choosing the physical variable anyway is right for a Z measurement.
* **Closure-driven FF binning** (N_jets, era) with the same-sign closure plots that justify it, and the τ2-T
  conditioning of the FF, which is the CMS τhτh convention.
* **Theory variations combined on histograms**, fiducial-normalised so they move only C; envelope / Hessian
  done right after having been done wrong once — and the wrong way documented.
* **W+jets large-weight handling** (shape from the relaxed region, one normalisation NP) instead of letting
  ±100-event γ parameters pollute the fit.
* **Both FF variants fitted and reported side by side**, with the difference tracked through the fit.

---

## 3. Findings in detail

### 3.1 [critical] The signal template is 38 % non-fiducial and the fake sideband contains signal

Numbers from the DY_NLO ntuple (SR, LHE ττ, all corrections):

| component of the selected Z→ττ | fraction |
|---|---:|
| inside the fiducial volume (60 < m_LHE < 120, both vis. τ pT > 40, \|η\| < 2.1) | 62.4 % |
| m_LHE > 120 GeV | **30.0 %** |
| 60–120 GeV, both genuine τh, visible pT/η outside the fiducial cuts (median gen pT(τ2) = 38.4 GeV) | 6.0 % |
| 60–120 GeV, a leg that is not a hadronic τ (τ→e/μ, flav 1–4) | 1.5 % |
| m_LHE < 60 GeV | 0.2 % |

and per fit bin (signal template vs data, from `review/`):

| m_tt bin [GeV] | signal | of which m_LHE > 120 | data |
|---|---:|---:|---:|
| 60–110 | 3992 | 0 % | 6941 |
| 110–130 | 1394 | 15 % | 7674 |
| **130–350** | **2440** | **88 %** | 32 205 |

Why this matters:

* The two 40 GeV visible-pT cuts have an acceptance of 0.23 % at the Z peak but tens of percent at
  m > 200 GeV, so the γ* continuum is enriched by two orders of magnitude relative to the peak. `docs/06`
  states "μ_Z multiplies the whole ττ template, including events outside any fiducial volume" as if this
  were a detail. It is not: 30 % of what is being fitted is *not the process whose cross section is quoted*,
  and its normalisation relative to the peak is taken from theory (NNLO/aMC@NLO) with no uncertainty
  assigned to that ratio.
* The bins above 130 GeV are described as "~95 % fakes and fix the fake normalisation". They are
  92.5 % fakes and 7.5 % high-mass signal *scaled by μ_Z*. The fake normalisation NP (`FakeOSSS_tautau`,
  pulled −0.86 σ) and μ_Z are therefore partially degenerate through the sideband. The correlation matrix
  shows only −2 % between them, because the peak bins dominate μ, but the sideband no longer does what
  the documentation says.
* The out-of-fiducial fraction is also why the theory NPs are large here and tiny in z-mumu:
  `QCDScale` moves C by 4.1 % (0.3 % in μμ) and the LO/NLO C difference is 7.3 % (0.2 % in μμ). Those
  are not uncertainties on the efficiency of a fiducial measurement, they are the theory uncertainty on
  the high-mass tail and on the migration across a cut placed on a steeply falling spectrum, dressed up as
  "signal modelling".
* The MC statistics of the high-mass tail (2718 raw events for 30 % of the template, 971 with
  m_LHE > 200) is part of the 4.8 % γ impact.

What to do:

1. Split `DYtautau` into `DYtautau_fid` (the fiducial part, `NormFactor mu_Z`) and `DYtautau_nonfid`
   (m_LHE > 120, m_LHE < 60, out-of-acceptance) treated as a background with its own theory
   normalisation uncertainty (high-mass DY is known to ~3–5 %; take `QCDScale`+`PDF` on the ratio
   non-fid/fid, computed from the GenSums which already exist). This is the standard ATLAS/CMS
   treatment and makes σ_fid = μ̂ · σ_fid^pred literally true. The non-fiducial part could alternatively be
   tied to μ_Z with a free ratio NP; either is honest, the current setup is not.
2. Use the high-mass sample DYJetsToLL_M-100to200 (record 35629, on dCache) stitched to the inclusive one
   for the tail, which also fixes a good part of finding 3.6.
3. Reconsider the fiducial visible-τ pT cut: with the reco cut at 40 GeV and a 10 % τ energy resolution, a
   gen cut at 40 GeV puts 6 % of the in-window signal just below it (median 38.4 GeV). Common practice is
   to keep them equal, so this is acceptable, but then the 6 % must be counted as the "migration" part of C
   and the τES NPs must be allowed to move it (they are; fine).
4. The CONVENTIONS file defines A with dressed leptons for ee/μμ and this channel uses LHE (Born) mass.
   For the total σ(60 < m_Born < 120) that is what the combination wants, so the ττ choice is right and
   the *other* channels have to make sure their 60–120 denominator is Born-level too. Flag it in the
   combination handoff explicitly; it is a 1–2 % effect for ee/μμ.

### 3.2 [critical] The nominal fake factor knowingly double-counts signal

`config.FF_SUBTRACT_MC = False` (config.py:126). The genuine-τ1 contamination measured by the team
itself: 1.3 % of the application region (2386 weighted MC events, almost all Z→ττ with τ1 failing Medium,
i.e. under the peak), 3.8 % of the C_OS/SS numerator, 0.9 % of the FF numerator. Consequences the team
documented: the prefit fake yield is 4.5 % (1760 events) too high, C is 1.078 instead of 1.051, μ is
1.160 instead of 1.208.

Why it is bad and not merely "an alternative":

* Subtracting simulated genuine-τ events from a fake-factor sideband is not an option in the method, it is
  part of its definition; the "old-time normal fake factor" of the CMS/ATLAS notes always includes it.
  Without it the estimate is by construction (data-driven fakes) + (a fraction of the signal), and the fit
  then removes signal from μ. The written prompt (`prompts/z_tautau.md`) only asks for "the old time normal
  Fake Factor"; it does not ask for no subtraction. If a verbal instruction did, it should be recorded, and
  the physically correct variant should still be what is quoted.
* Δμ = −0.048 is 30 % of the total systematic uncertainty and is a *bias*, not a noise term. Quoting the
  biased number as central and the correct one as "alternative" inverts the logic. Nor is the difference
  between the two a systematic uncertainty (it is not treated as one either, which is consistent, but then
  the bias is simply uncorrected).
* Independent evidence from the review BDT (section 4.3): in the purest fake phase space (BDT score < 0.1,
  11 900 fakes and 180 signal events) the nominal prediction is 7 % above the data
  (data/pred = 0.93 ± 0.01). The MC-subtracted C alone (−2.5 %) recovers a third of that; the rest is
  non-closure (3.3 below). The fit sees the same thing: `FakeOSSS_tautau` is pulled −0.86 σ, "in the
  direction of the MC-subtracted C".

What to do: flip the nominal (the code is ready; `handoff.md` item 1 already says so). Then re-run the
closure plots: the SS closure in m_tt with subtraction (`step3_closure_SS_m_tt_mcsub.png`) is the one
that should define `FakeClosure_tautau`.

### 3.3 [major] The τh ID scale factors: dominant, external, and internally inconsistent

The per-DM Medium SFs used (TauPOG `TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root`, read
correctly, bin errors as uncertainties): DM0 0.923 ± 0.142, DM1 0.880 ± 0.054, DM10 0.868 ± 0.085,
DM11 0.898 ± 0.163. Applied to both legs, fully correlated between legs, uncorrelated between DMs.
Combined impact 12.5 %, DM0 alone 9.7 %.

Observations:

* The pT-binned prescription from the same POG release gives 0.823 +0.086/−0.036 for pT > 40 GeV
  (recorded in `external/*.json` as a cross-check). That is 7 % lower per leg than the DM-binned average,
  i.e. ~14 % on the τhτh yield. Had the pT-binned SFs been used, μ_Z would be ≈ 1.33, 2 σ from the SM
  and from z-mumu. The docs choose DM-binned because "TauPOG recommends it for τhτh" (true: the
  DM-binned SFs are meant for selections with DM-dependent trigger/ID, such as the di-τ trigger); the
  choice is defensible, but the 14 % swing between two official prescriptions is larger than any nuisance
  parameter in the fit and is nowhere discussed. It has to be, in `docs/07` and in the handoff.
* These SFs are measured by the POG in Z→ττ (μτh) events. A τhτh Z cross section that takes them as
  external input is not a measurement of σ, it is a consistency check of the SF at the 12 % level, which the
  documentation admits. The interesting deliverable of this channel is therefore the *combination*, where
  ee/μμ fix μ_Z at 1–2 % and the ττ channel measures the four `TauID_DM*` parameters in situ.
  `handoff.md` says this; the combination config should make it explicit (e.g. decorrelate `SigModel`
  first, see 3.4, otherwise the ττ SFs would be measured against a 7.3 % "generator" prior).
* Uncorrelated-between-DM but fully-correlated-between-legs is the right correlation model. The DM0 and
  DM11 uncertainties (14 %, 16 %) are so large that a fit in decay-mode categories (section 4.4) would
  constrain them from the relative category yields even without the combination.

What to do: (i) document the DM- vs pT-binned discrepancy and quote the pT-binned result as a
cross-check; (ii) categorise the SR by decay mode so that the SF NPs are constrained in situ; (iii) in the
combination, present the fitted `TauID_DM*` values as a result.

### 3.4 [major] `SigModel` (LO vs NLO) is not a credible uncertainty and is wrongly correlated with z-mumu

`step4_histograms.py:211` builds the `SigModel` template from madgraph MLM (LO) normalised to the NLO
fiducial cross section; `step5_fit.py:107` makes it a one-sided, symmetrised, normalisation-only NP of
7.3 %. It is the second-largest uncertainty (9.8 % after the fit, correlated +48 % with μ_Z).

Why it is bad:

* LO madgraph is a *worse* model of the Z pT spectrum than aMC@NLO FxFx, and the entire 7.3 % arises
  from the migration across the 40 GeV visible-pT cut, i.e. from the Z pT / τ pT slope (the docs say so).
  An uncertainty built as "nominal minus a model known to be worse" is a conservative placeholder, not a
  measurement of ignorance. z-mumu handled this correctly: NLO nominal, powheg (NLO) alternative, 0.2 %.
* It double counts `QCDScale` (4.1 % on C, same physics: Z pT slope) and `PS_ISR`.
* It carries the name `SigModel`, which `fitting/CONVENTIONS.md` correlates across channels. In the
  combination the ττ 7.3 % normalisation NP and the μμ 0.2 % powheg-shape NP would be *the same
  parameter*. They are physically unrelated; the ττ one must be renamed (`SigModel_tautau`) or, better,
  replaced.

What to do (in order of preference):

1. Replace it by a data-driven Z pT uncertainty: the Z→μμ channel measures the Z pT spectrum in the same
   data set with per-mille statistics. Reweight the ττ aMC@NLO sample to the measured μμ Z pT spectrum
   (a standard CMS procedure for exactly this sample) and take the reweighting as the shape/normalisation
   uncertainty. This is what "the three channels share the workspace" should buy you.
2. Or use an NLO-only prescription (scales + PS + PDF), which the handoff estimates at ~5 %, and add
   the pT_Z-binned aMC@NLO samples on dCache (35615–35625) for statistics in the tail.
3. At minimum rename the NP so the combination does not correlate it with z-mumu.

### 3.5 [major] The fake-uncertainty model under-estimates the fake uncertainty

The fit reports "fakes 0.8 %". The ingredients:

* **`FakeClosure_tautau`** (`step4_histograms.py:91–114`): up = same-sign closure ratio per m_tt bin
  (clipped to 0.5–1.5), down = 1/ratio, one NP for all 14 bins. The ratio is ±11 % with 4–13 % statistical
  precision in the peak bins and ±2.5 % with sub-percent precision above 110 GeV. Post-fit: pull +0.70,
  constraint 0.32.

  Reasoning: a single parameter ties the (statistically precise, physically tiny) non-closure in the
  high-mass bins to the (imprecise, 10 %) non-closure under the Z peak. The data constrain the NP to a
  third of its prior through the high-mass bins, and that constraint is then applied to the peak bins where
  the prior was mostly statistical noise. The peak-bin fake uncertainty that survives is ~3 % instead of
  ~10 %. The 36 % correlation of this NP with `γ_bin6` (100–110 GeV) and −25 % with `γ_bin10` shows it
  is fitting bin fluctuations. The correlation with μ_Z is small (−2 %), so the *central value* is not
  endangered, but the 0.8 % is.

  Fix: either decorrelate the closure NP in two or three mass regions (< 110, 110–150, > 150 GeV), or
  replace the bin-by-bin ratio by a smooth (linear) fit to the SS non-closure with its two parameters as
  NPs, or — cleanest — derive a closure *correction* and keep only its statistical uncertainty.

* **`FakeStat_tautau_DM*`**: one NP per DM shifts all 30 (era × N_jets × pT) bins of that DM coherently
  by their statistical error (`fakes.py:101`). The statistical errors of 120 independent bins are not
  correlated; treating them as fully correlated over-states the coherent normalisation part and
  under-states the shape part. Since the fake normalisation is anyway fixed by the sideband, the net effect
  is an under-estimate of the shape uncertainty. Fix: propagate the FF statistical uncertainty bin-by-bin
  into the Fakes template's Sumw2 (it is then handled by the γ parameters, which already carry the AR
  statistics), or use a few eigen-variations.

* **η(τ1) non-closure ±15 %** (`step3_closure_SS_t1_eta_nominal.png`, and identically in the OS
  τ2-anti-isolated region and the SR): diagnosed in section 5. The fake factor itself varies by ±15 % with
  |η(τ1)| and the binning has no η. It integrates out of the m_tt template (< 0.4 % per bin, −0.03 % on the
  yield), so it does not affect μ_Z, but it makes every η(τ1) plot wrong and would bias any categorisation
  (ML, section 4) that uses η. Fix with a factorised closure correction in |η(τ1)| (the shape is universal).

* **pT(τ2) non-closure −7 % at 60–100 GeV**: the FF depends on the second leg because the isolation of the
  two jets in a dijet is correlated. The same physics makes the BDT closure (4.3) trend by ±4 %. A
  correction in pT(τ2) or ΔR is the usual fix.

* **`FakeOSSS_tautau`**: 3 % flat, well motivated by the C vs τ2-WP trend; pulled −0.86 σ in the nominal
  because of 3.2. Fine once 3.2 is fixed. One caveat: C is measured in the τ2-anti-isolated sideband where
  W+jets (charge-asymmetric) is a larger fraction than in the SR, and the quark/gluon mixture differs. The
  3 % is a reasonable envelope; keep.

* **W+jets in the application region**: SR events with a jet as τ1 and a genuine τ2 from W→τν are
  estimated with the QCD fake factor. The W+jets jet→τh FF is typically 1.5–2× the QCD one at the same pT
  (quark jets). At the few-percent level of the AR this is ≲ 1 % of the fakes and is not covered by any NP.
  Mention it in `docs/07`; a W-enriched check (m_T(τ2, MET) > 70 GeV in the AR) would bound it.

Realistic fake uncertainty after these fixes: a few percent on μ, still far below the τh ID.

### 3.6 [major] MC statistics (4.8 %) is avoidable

Third-largest impact, from the aMC@NLO signal in the peak bins (negative weights halve the effective
statistics; 10 437 raw SR events) and from the high-mass tail. On dCache and unused: DYJetsToLL_0J/1J/2J
(35577/35595/35613), M-100to200 (35629), the pT_Z-binned set (35615–35625). Stitching 0J/1J/2J with the
inclusive sample is standard and would multiply the SR statistics by ~4. The handoff lists it as a next
step; it should have been done before quoting a 4.8 % uncertainty from it.

### 3.7 [minor] Trigger: both legs on the turn-on

TauPOG data efficiencies at pT = 40 GeV: DM0 0.49, DM1 0.63, DM10 0.28, DM11 0.31 (SF uncertainties
3.5–10 % per leg). The claim in `docs/02` that 40 GeV is "above the 35 GeV HLT threshold where trigger SFs
are reliable" is optimistic; the plateau starts at ~60 GeV. I checked whether raising the threshold helps
(section 4.2): the signal-weighted mean SF uncertainty is flat at 6.3–6.5 % from 40 to 60 GeV because the
POG uncertainties do not shrink on the plateau for this year. So the 5.1 % is external like the ID, and the
cure is the same: in situ (combination) or an independent trigger-efficiency measurement with μτh events
from SingleMuon (not in the skims). Not worth a threshold change on its own.

### 3.8 [minor] Small implementation points

* `analysis.py:193`: the di-τ trigger SF (a genuine-τ efficiency ratio) is applied to legs with
  `genPartFlav = 0` (jets faking τ2 in W+jets/tt̄) and to e/μ→τh legs. Harmless at 1.2 % of the SR, but
  the W+jets normalisation NP is then the only cover. Apply SF = 1 to jet legs, or note it.
* `objects.py:127–129`: N_jets uses tight jet ID, no pileup-jet ID, pT > 30, |η| < 4.7. As a binning
  variable of a data-driven FF this is self-consistent; for the MC subtraction (`mcsub`) the pileup-jet
  content of MC and data differ slightly. Negligible now, worth a PU-ID if N_jets ever becomes a fit
  category.
* `config.py:57`: trigger matching uses `filterBits & 2` ("MediumChargedIso"). In NanoAODv9 the τ
  filter-bit assignment differs between 2016 and 2017/18 and the di-τ leg bit is not always filled for 2016;
  97 % match rate suggests it works, but a one-line check that bit 2 corresponds to the
  `*MediumIsoPFTau35*` legs in the 2016 UL TrigObj documentation belongs in `docs/02`.
* EWK Z→ττ (EWKZ2Jets, record 35999, on dCache) is absent; ≈ 0.3 % of DY inclusively, more in the
  ≥ 2-jet FF bin. Small, but decide whether it is signal or background and write it down.
* `TauFakeMu` SF of 3.55 ± 0.29 in 1.7 < |η| < 2.3 is the real POG value; fine, and it affects nothing.
* `fakes.py:87`: the FF statistical error uses the *unsubtracted* counts; correct for the nominal, slightly
  optimistic for `mcsub` (the MC subtraction's own statistics is ignored). Fine at 1 %.
* Repository root: `z-tautau_build_trex.sh/.log` (untracked) build TRExFitter with
  `asetup StatAnalysis,0.8.0`, a different ROOT than LCG_110 and a different version than the documented
  fallback (0.7.2). CLAUDE.md forbids a second build. Delete the scripts or move the recipe into
  `fitting/`, and record which binary produced the committed `results.json`.

### 3.9 Documentation nits

* `docs/02`: "above the 35 GeV HLT threshold where trigger SFs are reliable" — see 3.7.
* `docs/05`: "no MC subtraction in the nominal (as requested)" — the request is not in the written prompt.
* `docs/08`: "Bins above 130 GeV are ~95 % fakes and fix the fake normalisation" — 92.5 % fakes, 7.5 %
  signal scaled by μ (3.1).
* `handoff.md` BLUE inputs: correctly warns that counting is useless; good. Add the DM- vs pT-binned
  SF swing to the "uncorrelated" list.

---

## 4. The fake background: why 81 %, and how to reduce it

The fakes are large because both legs use DeepTau Medium (jet misidentification ~0.5 %) and because the
selection is purely "two τh + trigger". Nothing kinematic is asked. Three handles were quantified on the
ntuples (`review/fake_studies.py`; signal = DY_NLO with all corrections, fakes = data − MC with genuine τ1;
"stat-only δμ" = 1/√Σ s²/(s+b) over the m_tt fit bins with fixed backgrounds, i.e. the same quantity as
`stat_only_sensitivity` in `results.json`, 1.83 % nominal).

### 4.1 DeepTau working point (both legs, pT > 40 GeV)

| WP | data | Z→ττ | fakes | fake fraction | S/B in 70–110 | signal kept | stat-only δμ |
|---|---:|---:|---:|---:|---:|---:|---:|
| Medium (used) | 47 586 | 7878 | 38 219 | 0.80 | 1.55 | 1.00 | 1.87 % |
| Tight | 21 160 | 5963 | 14 211 | 0.67 | 2.69 | 0.76 | 1.91 % |
| VTight | 9935 | 4075 | 5446 | 0.55 | 3.79 | 0.52 | 2.13 % |
| VVTight | 4642 | 2580 | 1902 | 0.41 | 4.60 | 0.33 | 2.46 % |

Reading: Tight removes 63 % of the fakes for 24 % of the signal *at no cost in statistical precision*. Since
every fake-related systematic (closure, OS/SS, FF statistics, W+jets composition) scales with B/S, and the
signal MC statistics (4.8 %) gets slightly worse, Tight is the better operating point for a systematics-
dominated measurement. Two caveats: the ID SFs and their uncertainties are WP-specific (the POG file has
all WPs; check that Tight's DM0/DM11 uncertainties are not larger than Medium's), and the fake factor must
be re-derived with T = Tight (the code takes it from one constant, `config.py:70`). Recommended: run the
full chain at Tight and compare the total uncertainty.

### 4.2 τ pT threshold (Medium)

| pT > | data | Z→ττ | fakes | fake fraction | S/B 70–110 | signal kept | stat-only δμ | ⟨trigger SF unc.⟩ |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 47 586 | 7878 | 38 219 | 0.80 | 1.55 | 1.00 | 1.87 % | 6.3 % |
| 45 | 26 937 | 4599 | 21 404 | 0.79 | 2.52 | 0.58 | 2.30 % | 6.5 % |
| 50 | 16 123 | 3142 | 12 248 | 0.76 | 3.02 | 0.40 | 2.71 % | 6.5 % |
| 60 | 6716 | 1645 | 4728 | 0.70 | 4.00 | 0.21 | 3.72 % | 6.3 % |

Reading: the fake fraction barely moves (QCD and Z→ττ have similar τ pT spectra once both legs are above
40 GeV), the signal halves by 50 GeV, and the trigger uncertainty does not improve. Do not raise the
threshold; 40 GeV is the right choice.

### 4.3 A k-fold BDT on mass-agnostic kinematics (proof of concept)

Set-up (all in `review/fake_studies.py`, 19 s):

* Signal: DY_NLO SR events with LHE ττ (10 437 raw, positive weights for training, signed weights for
  yields). Background: the nominal fake estimate itself, i.e. AR data events weighted by FF × C (187 973
  events), so that the classifier learns *fakes as the analysis models them*.
* 16 inputs, none of which can see the τ1 isolation that defines the FF regions: pT(τ1), pT(τ2), their
  ratio, |η| of both, ΔR, Δφ(τ1,τ2), MET, MET significance, pT of the visible pair, Δφ(MET, ττ),
  pT(ττ+MET), N_jets, leading-jet pT, DM1, DM2. No mass variable (m_vis, m_tt, m_col, m_T^tot), so the
  score can be used as a category and m_tt stays the fit variable inside each category.
* XGBoost, depth 4, 300 trees, 5 folds by `event mod 5`; every event (signal MC, AR, SS_L, SS_T, SR data,
  other MC) is scored with the model that did not see its fold.

Results:

| quantity | value |
|---|---|
| held-out AUC (signal vs fakes) | 0.888 |
| most important inputs | ΔR (0.28), pT(ττ+MET) (0.14), Δφ(τ1,τ2) (0.09), pT(vis ττ), pT(τ1), DM2, DM1 |
| top signal quartile (score > 0.948) | 1627 signal, 165 fakes: S/B = 9.8 (12.4 in 70–110 GeV) |
| second quartile | S/B = 2.1 |
| stat-only δμ, m_tt fit inclusive → in 4 score categories | 1.81 % → 1.65 % |

Same-sign closure of the FF *in the score* (SS_T observed / FF × SS_L), the test any ML selection on top
of a fake factor must pass:

| score | 0–0.1 | 0.1–0.2 | 0.2–0.3 | 0.3–0.4 | 0.4–0.5 | 0.5–0.6 | 0.6–0.7 | 0.7–0.8 | 0.8–0.9 | 0.9–1 |
|---|---|---|---|---|---|---|---|---|---|---|
| obs/pred | 0.97 | 0.98 | 1.04 | 1.03 | 1.04 | 1.03 | 0.99 | 1.02 | 0.97 | 1.12 ± 0.09 |

Signal region in the score, prefit, nominal FF (data / prediction): 0.93 at score < 0.1 (11 900 fakes, 180
signal), 0.98–1.01 in the middle, 1.05–1.15 above 0.7 where the signal dominates. See `review/review_bdt.png`.

Assessment:

* The separation is real and comes from topology (ΔR, Δφ, the pT balance of the ττ+MET system): Z→ττ
  at these pT is back-to-back with the neutrinos along the τ's, QCD dijets are back-to-back with the MET
  pointing elsewhere or nowhere. The classifier isolates a fifth of the signal at S/B = 10 and half of it
  at S/B ≥ 2, versus S/B = 0.2 inclusively.
* The closure of the FF in the score holds at the ±4 % level with a mild trend (the same pT(τ2)/ΔR
  correlation as in 3.5): usable, but a closure correction as a function of the score (or an FF binned in
  ΔR) would be needed before the categories are trusted at the percent level.
* The gain is where expected: not in statistics (1.81 → 1.65 %; the analysis is not stat-limited) but in
  the leverage of fake systematics, which in the high-score categories act on 165 events instead of
  38 000. The category structure also lets the fit measure the fake normalisation and shape where fakes
  dominate and the signal where signal dominates, which removes the sideband/signal degeneracy of 3.1.
* It does *not* touch the τh ID (12.5 %), signal modelling (9.8 %), trigger (5.1 %) or MC statistics (4.8 %,
  which gets a little worse). Anyone proposing the ML step should present it as the fix for the fake
  systematics and as a cleaner category structure, not as the path to a precision result.
* Practicalities for a real implementation: k-fold by event number exactly as here, and apply the fold
  models to AR, SS_L, SS_T, SR and all MC consistently (the AR events *are* the fake template, so a model
  trained on them must never score them); never use τ1 raw DeepTau, isolation sums or
  `leadTkPtOverTauPt`; train against AR × FF (fakes as modelled) and cross-check against SS_T data (fakes
  as observed); re-derive `FakeClosure` per category; with `mcsub` nominal, subtract MC from the AR
  in each category; keep m_tt as the fit variable per category (two or three categories are enough; the
  top quartile alone is a near-background-free Z peak, a good control of the τ ES too).

### 4.4 Decay-mode categories (no training needed)

The fake factors span 0.06 (DM11) to 0.44 (DM0) and the ID SF uncertainties span 5 % (DM1) to 16 %
(DM11). Splitting the SR by the leading (or both) τ decay modes therefore gives categories with very
different S/B and lets the fit constrain `TauID_DM*` from the relative yields once μ_Z is fixed by the
light-lepton channels. This costs nothing but histograms and is the natural companion of the combination
strategy in `handoff.md`. It should be done before, or together with, the BDT.

### 4.5 Recommended order

1. `mcsub` nominal (3.2); non-fiducial split (3.1); jet-binned + high-mass DY samples (3.6); rename or
   replace `SigModel` (3.4). These change the central value and the two largest uncertainties.
2. Decorrelate/smooth the closure NP, FF statistics bin-by-bin, |η| in the FF binning (3.5).
3. Tight working point on both legs (4.1) *or* the BDT categories (4.3) — the latter if the group wants
   the ML exercise, the former if it wants the smallest change. Both together are also fine.
4. Decay-mode categories (4.4) and the combined fit with `TauID_DM*` as measured parameters.

---

## 5. The η(τ1) non-closure: where it comes from

The observation (raised by the reviewer of this review, and visible in `step4_SS_T_t1_eta.png`,
`step4_OSAI_T_t1_eta.png`, `step4_SR_t1_eta.png`): data/prediction versus η(τ1) has the same symmetric
structure in every region, ≈ 0.88 at |η| < 0.4, rising to ≈ 1.15 at |η| ≈ 1.2–1.4, back to ≈ 0.9 beyond 1.8,
while η(τ2) closes to within a few percent. Since the pattern is identical in same-sign (where the FF is
measured) and in opposite-sign regions, it is not an OS/SS effect and not a signal effect: it is a property
of the fake factor of the leading τ. `review/eta_study.py` measures it directly.

**The fake factor depends on |η(τ1)| by ±15 %, and the binning does not know it.**
FF = N(SS, τ1 Medium) / N(SS, τ1 VVVLoose & !Medium), τ2 Medium, inclusive in era/DM/N_jets/pT:

| \|η(τ1)\| | 0–0.4 | 0.4–0.8 | 0.8–1.2 | 1.2–1.5 | 1.5–1.8 | 1.8–2.1 |
|---|---|---|---|---|---|---|
| FF | 0.177 | 0.192 | 0.203 | **0.219** | 0.186 | **0.169** |
| FF / mean | 0.93 | 1.00 | 1.07 | 1.15 | 0.97 | 0.89 |

Every η bin of a prediction uses the η-averaged FF of its (era, DM, N_jets, pT) cell, so the prediction is
15 % low where the true FF is 15 % high, which is exactly the ratio plot. The shape is the same as the
closure ratios (0.95 / 1.10 / 0.94 in three |η| bins, both eras).

**The η shape is universal**, i.e. it is not a hidden dependence on something already binned:

| split | FF(η)/mean in the six bins |
|---|---|
| Run2016G | 0.93 1.00 1.06 1.16 0.97 0.89 |
| Run2016H | 0.93 1.01 1.07 1.14 0.98 0.89 |
| 0 / 1 / ≥ 2 jets | same shape within statistics |
| pT 40–45 … > 80 GeV | same shape within statistics |
| DM0 | 0.88 0.95 **1.16 1.27** 1.02 0.82 (FF 0.31 → 0.44 → 0.29) |
| DM1 | 0.91 1.00 1.05 1.16 0.99 0.93 |
| DM10 | 0.99 1.03 1.07 1.01 0.91 0.93 |
| DM11 | 0.90 1.05 1.02 1.22 0.95 0.91 |

The two eras agree, so the different HLT τ isolation of G and H (the one era-dependent ingredient) is
not the cause. pT and N_jets binning do not touch it. The decay-mode composition of the fakes does change
with η (more DM0, fewer DM1 beyond |η| = 1.5), but the FF is binned in DM and the shape is present *inside*
each DM, strongest for 1-prong (DM0) and weakest for 3-prong (DM10).

**Two mechanisms, seen by tightening the loose definition** (FF shape vs |η| for L = X & !Medium):

| L | 0–0.4 | 0.4–0.8 | 0.8–1.2 | 1.2–1.5 | 1.5–1.8 | 1.8–2.1 |
|---|---|---|---|---|---|---|
| VVVLoose (used) | 0.93 | 1.00 | 1.07 | 1.15 | 0.97 | 0.89 |
| VVLoose | 0.95 | 1.03 | 1.06 | 1.11 | 0.95 | 0.87 |
| VLoose | 0.96 | 1.04 | 1.07 | 1.09 | 0.93 | 0.85 |
| Loose | 0.97 | 1.06 | 1.06 | 1.04 | 0.93 | 0.87 |

* The **bump at 1.2–1.5** shrinks from +15 % to +4 % as L approaches Medium: it is the *shape of the
  DeepTau VSjet score distribution just below the Medium threshold* that changes with η. In the barrel–
  endcap transition the loose jets are more τ-like (median raw score of the L τ's rises from 0.505 to 0.53
  beyond |η| = 1.2, the top quartile of the L population grows from 24 % to 27 %), so a larger fraction of
  them sits close to the threshold and the T/L ratio is higher. DeepTau uses η and the transition-region
  ECAL/tracker response as inputs; the misidentification probability is simply not η-flat, and a wide
  loose window (VVVLoose, raw score 0.26–0.88) integrates over the η-dependent part of the spectrum.
* The **drop beyond |η| ≈ 1.5** (−10 to −15 %) survives every loose definition, so it is a genuine
  η-dependence of the Medium efficiency for jets at the edge of the tracker/HLT acceptance (|η| < 2.1 is
  the HLT boundary, and both legs are trigger-matched), not a threshold effect. DM0 shows it most
  (single track, no π⁰).

**Does it matter for the result?** No. Applying a factorised closure correction f(|η(τ1)|) (three bins,
measured in SS) changes the SR fake yield by −0.03 % and every m_tt bin of the fake template by less than
0.4 %, and the SS closure in m_tt is unchanged, because m_tt is uncorrelated with η(τ1) inside the
application region. The correction restores the η(τ1) closure to within ±5 % (from ±15 %), see
`review/review_eta.png` and `review/eta_study.txt`.

**What to do.** Add the factorised correction f(|η(τ1)|) to `fakes.evaluate` (the universality means one
function for all eras/DMs/pT bins is enough; 6 bins have 2–4 % statistical precision), or add a 3-bin |η|
axis to the FF table (costs a factor 3 in bin statistics, 10 % → 17 %). Regenerate the η plots and the
same-sign closure plots. This is a prerequisite for any decay-mode or ML categorisation, where η enters the
classifier (|η| carried 4 % and 3 % of the BDT importance in section 4.3) and the current mis-modelling
would leak into the category yields. It is also a good reason to use a *narrower* loose window
(VLoose & !Medium halves the bump at no real cost: 53 700 denominator events instead of 140 000), which
makes the FF closer to a local property of the score distribution and less sensitive to all such
integrations, at the price of larger statistical errors per bin.

---

## 6. Reproduction log

* `python review/fake_studies.py` (after `source ../setup.sh`): loads all ntuples in 6 s, reproduces the
  nominal FF table and C_OS/SS to the committed values (fakes 38 758, C = 1.078), runs the three studies in
  19 s total on 8 cores.
* Committed numbers cross-checked against `output/results.json`: μ, impacts, pulls, contamination, C, cut
  flow all consistent with `docs/` and `handoff.md`. The only textual inconsistencies are listed in 3.9.
* `output/data/`, `fit/fitinputs/`, `fit/results/` are absent from this checkout (git-ignored); the fit was
  not re-run here. Nothing in the review depends on re-running it.

---

## 7. Addendum (15 September 2026): what v2 did about it

The analysis was redone after this review; the numbers below are the v2 result (`output/RESULTS.md`).

| finding | action in v2 | where |
|---|---|---|
| 3.1 non-fiducial signal | `DYtautau` = fiducial part only (μ_Z); `DYtautau_nonfid` = background with `XS_DYtautau_nonfid` 5 % and the theory variations acting on the non-fid/fid ratio | `analysis.mc_components`, `docs/06` |
| 3.2 no MC subtraction | `config.FF_SUBTRACT_MC = True`; `nosub` kept as a cross-check (μ_Z = 1.212 +0.144 -0.125 vs 1.212 +0.145 -0.125); W+jets subtracted with uniform weights after one event with weight −63 faked a 75 % non-closure in the top BDT bin | `fakes.py`, `analysis.subtraction_weights` |
| 3.3 τh ID prescriptions | documented (`docs/07`); the pT-binned cross-check and the DM categories remain open items | |
| 3.4 `SigModel` | renamed `SigModel_tautau`, computed for the fiducial C (C_LO/C_NLO = 0.867), reported but not fitted (`config.SIGMODEL_IN_FIT`); NLO scale/PS/PDF variations kept; visible-pT spectra checked against data in the signal-dominated category | `docs/07` |
| 3.5 fake uncertainty model | closure corrections f(\|η(τ1)\|), g(pT(τ2)); FF statistics per event into the `Fakes` variance; `FakeClosure_tautau_c<k>_lo\|hi` per category and mass region (sizes: c0/lo 15 %, c0/hi 1 %, c1/lo 5 %, c1/hi 7 %, c2/lo 12 %, c2/hi 55 %) | `fakes.closure_corrections`, `step4` |
| 3.6 MC statistics | DY 0J/1J/2J skimmed (with `LHE_NpNLO`) and stitched per jet bin with the inclusive sample's bin fractions; γ impact 4.7 % | `analysis.dy_norm` |
| 3.7 trigger | unchanged (external); documented as such | |
| 3.8 trigger SF on jet legs | SF = 1 for `genPartFlav = 0` legs | `corrections.trigger_sf` |
| 3.9 second TRExFitter build | the checkout's binary had been built against the StatAnalysis ROOT and crashed under LCG (`TExMap::Add key not unique`, `munmap_chunk`); rebuilt with `fitting/build_trexfitter.sh` | |
| 4.1 Tight WP | full chain rerun with `BND_TAUTAU_WP=Tight` (`variants/tight/`): μ_Z = 1.099 +0.116 −0.102, total +11.6/−10.2 % (Medium +14.5/−12.5 %), GoF p = 0.32 (Medium 0.01); recommended as the next nominal | `docs/08` |
| 4.3 BDT | `ztautau/bdt.py`, `scripts/step3b_bdt.py`: 5-fold XGBoost, mass-agnostic inputs, three categories fitted in m_tt; held-out AUC 0.970 (higher than the 0.89 of section 4.3 because the signal is now the fiducial part only); same-sign closure in the score verified | `docs/09` |
| 5 η non-closure | corrected; the residual is within ±5 % | `docs/05` |

Result: μ_Z = 1.212 +0.145 -0.125 (v1: 1.160 +0.192 −0.163), σ(60–120) = 2358 +282 −242 pb, grouped impacts in `output/RESULTS.md`.
Still open (in `CLAUDE.md`): decay-mode categories, in-situ trigger efficiency, HT-binned W+jets, the high-mass DY
sample, EWK Z.
