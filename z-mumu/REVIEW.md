# Review of the Z → μ⁺μ⁻ cross-section measurement (v2)

> **Historical document (15 Sep 2026).** The numbers below are the review's. The channel result was frozen on
> 17 Sep 2026 with two further changes that the review did not see: the muon reconstruction scale factor is
> measured (F8's 0.4 %/muon became 1.0001 ± 0.0013 per muon) and the acceptance uncertainties are evaluated for
> the 60–120 GeV denominator with the rows of `docs/16`. Frozen result: σ_fid = 790.2 ± 0.2 (stat) ± 6.1 (syst)
> ± 9.6 (lumi) pb, σ(60–120) = 1931 ± 30 pb (`handoff.md`).

Reviewed 15 Sep 2026 against the committed state of `z-mumu/` (result of 14 Sep: σ_fid = 791.8 ± 0.2 (stat)
± 6.2 (syst) ± 9.3 (lumi) pb, μ_Z = 0.990 ± 0.014). Everything below was re-derived from the repository
outputs and, where the outputs did not contain the information, from a dedicated pass over the v2 skims
(`review/review_pass.py` → `review/review_pass.pkl`) and from four additional TRExFitter fits
(`review/fitcheck/`). Figures: `review/figures/` (made by `review/deck_plots.py`); slides:
`review/deck/zmumu_review.pdf`.

The two questions asked explicitly — the large data/MC discrepancy in the e-μ control region and the
p_T^miss mismatch in the signal region — are answered in §2 and §3. Both are understood and neither
affects the cross section. The finding that *does* affect the result is in §4: the 30-bin mass fit is not
robust against the choice of shape systematics and binning (μ_Z moves by up to 1.5 %), because the
aMC@NLO template does not describe the 60–80 GeV tail and the fit uses over-constrained shape nuisance
parameters to fix that.

## 0. Status after the fixes (15 Sep 2026, same day)

Every finding below was acted on in the commit that follows fc9a998; the numbers in sections
1-8 describe the state *before*. What changed (details in `docs/14-fit-and-systematics.md`,
`docs/11-mc-weights.md`, `output/v2/RESULTS_v2.md`):

| finding | fix | effect |
|---|---|---|
| F1 e-μ region | all MC kept in the e-μ regions (+ W+jets, `XS_WJets` 30 %); same-sign e-μ region added (`SSemu`) | data/MC 0.992 in OS, flat in electron p_T, jets and mass; SS e-μ under-predicted by 24 % (W+jets MC statistics, no QCD) -- plots only |
| F2 p_T^miss | Puppi p_T^miss added to the plots; nothing else needed | none |
| F3/F5 fit not robust | 1 GeV input bins, **fit in 12 × 5 GeV bins**, no template smoothing, MINOS on every parameter (`UseMinos: all`; the Hesse covariance was ill-conditioned and gave the spurious 0.2-0.8 "constraints" on `Lumi`, `MuonReco`, ...), stability script `scripts/v2_5_fit_variants.py` | μ_Z = 0.9881 ± 0.0159, GoF p = 0.79; 0.9858-0.9901 over the binnings with p > 0.05; counting 0.9939; the no-`SigModel` and smoothed variants have p ≤ 0.01 and are rejected |
| F4 `SigModel` template | powheg / aMC@NLO ratio with both generators inside 50 < m_LHE < 120 GeV, normalised to the NLO fiducial prediction in the window, two-sided (mirrored) | no edge artefact (last-bin ratio 1.02 instead of 0.67); pull +0.66 σ, constraint 0.14; C(powheg)/C(aMC@NLO) = 1.0028 |
| 62-78 GeV deficit | generator-level comparison added (`sigmodel_lineshape.png`): powheg has 4-6 % less continuum than aMC@NLO at Born level; the data lie between | covered by the two-sided `SigModel`; listed as an open issue |
| F6 pileup | profile scale (1.035) and bunch-to-bunch smearing (9 %) fitted to the SR N_PV distribution (`scripts/v2_2_pileup.py`) | N_PV data/MC within 3 % over the bulk (χ²/ndf 104 414 → 4 729 / 50); `Pileup` pull −1.5 → +0.06 |
| F7 jets | `regions.clean_jet_count` (ΔR > 0.4 to tight/anti muons and selected electrons) | jet multiplicity peaks at 0 |
| F8 `MuonReco` | 0.4 %/muon fully correlated → 0.8 % per event (`fitting/CONVENTIONS.md` updated) | +0.77 % impact; total systematic 1.57 % instead of 1.40 % |
| F9 MC statistics | stated in `RESULTS_v2.md` and `handoff.md` | -- |
| §7 small items | `DYmumu_powheg`/`WJets` rows dropped from the yields table (note added); `err_decomp` removed from the JSON; `Fakes` template without Sumw2; luminosity quoted as the external 1.2 % | -- |

Result after the fixes: **σ_fid = 790.1 ± 0.2 (stat) ± 8.5 (syst) ± 9.6 (lumi) pb** (μ_Z = 0.988 ± 0.016),
against 791.8 ± 6.2 ± 9.3 pb before. The central value barely moved; the systematic grew because
`MuonReco` doubled and because the shape fit now runs in bins where the shape NPs are not
over-constrained. The "±0.7 % lineshape term" recommended below is no longer needed: the two-sided
`SigModel` carries it, and the fit is stable to ±0.2 % against the binning.

## 1. Summary of findings

| # | finding | severity | effect on σ_fid |
|---|---|---|---|
| F1 | e-μ region: the +9.2 % data excess is the non-prompt-electron contribution (W+jets jet→e, Z→μμ+γ conversions, Z→ττ τ_h→e, tt̄ b→e). The MC is filtered to prompt-prompt events and nothing replaces what is removed. Putting the removed MC back gives data/MC = 0.996, flat in every variable. | understood, cosmetic for the result (validation region only) | none |
| F2 | p_T^miss in the SR: data is broader than MC. Puppi p_T^miss in 0-jet events agrees to < 1 %; the PF p_T^miss residual grows with pileup (unclustered/pileup energy) and both grow with jet multiplicity (jets in MC are not JER-smeared). p_T^miss is not used in the selection or the fit. | understood, not a problem | none |
| F3 | The fitted μ_Z depends on the shape model: 0.990 (nominal 30 bins) / 0.996 (6 bins) / 0.9935 (1 bin = counting) / 1.005 (30 bins without the powheg `SigModel` template). Only the nominal configuration has an acceptable goodness of fit; it obtains that by pulling `SigModel` by 0.5σ and constraining it to 10 % of its prior. The data show a 3–4 % deficit relative to aMC@NLO at 62–78 GeV that no other nuisance parameter can describe. | **must fix before the combination** | up to ±0.8 % (not in the quoted 0.77 % non-lumi systematic) |
| F4 | The `SigModel` template is flawed: the powheg sample is generated with 50 < m_LHE < 120 GeV, so its last two reconstructed-mass bins are depleted by 10 % and 33 % (resolution migration across the generator cut), and its low-mass continuum differs from aMC@NLO by 6–11 %. The template mixes a phase-space artefact with the generator comparison. | must fix | part of F3 |
| F5 | Several nuisance parameters are constrained far below their prior by the 2 GeV shape (`MuonScale` 0.11, `SigModel` 0.10, `MuonRes` 0.21, `QCDScale` 0.51, `PS_FSR` 0.53, `Pileup` 0.69, `Lumi` 0.82, `MuonIso` 0.82) and some are pulled (`Pileup` −1.5σ, `PDF` +1.3σ, `MuonIso` −1.1σ). The luminosity is "measured" by the tt̄ sideband shape. These constraints are fit artefacts, not measurements. | should fix | the quoted systematic is ~0.1 % too small (lumi 1.16 % instead of 1.2 %) |
| F6 | Pileup profile: the data profile built from the luminosity CSV is ~1 interaction (4 %) too low; after reweighting, data/MC in N_PV rises from 0.8 to 1.2 across the distribution (era H is worse than era G). | should fix | ≤ 0.05 % (isolation SF measured in data; SF vs N_PV flat) |
| F7 | `njet` in the histograms is not lepton-cleaned: the two muons (or the muon and the electron) are counted as jets, so the "jets" plots peak at 2 and are meaningless. Lepton-cleaned multiplicities are used in the skim category C and in this review. | bug (plots only) | none |
| F8 | `MuonReco` is documented as 0.4 % per muon (docs/11, v1 used 0.8 % per event) but applied as 0.4 % per event in the fit. One of the two is wrong. | inconsistency | +0.35 % in quadrature if the per-muon reading is right |
| F9 | MC statistics: the per-bin MC uncertainty is 1.7× the data uncertainty in every mass bin (72 M generated DY events, 16 % negative weights, ~10 M selected); the 30 gamma parameters are constrained by the data and pulled up to 1.3σ, and their impact (0.35 %) is 10× the data statistics. | inherent, worth knowing | 0.35 % (already quoted) |
| F10 | Choices that are sound: tight ID + iso, dressed fiducial volume, normtag luminosity, trigger-object tag-and-probe with pass/fail fits, prompt subtraction in the fake factor, same-sign closure, Z-peak momentum calibration, NLO acceptance with 60 < m < 120 denominator recommended. See §6. | — | — |

Bottom line: the fiducial cross section is right at the 1 % level, but the quoted central value and
systematic come from a fit whose shape freedom is not under control. Until F3/F4 are fixed, the safer
number is the counting extraction, σ_fid = 794.4 pb, with an added ±0.7 % lineshape-model uncertainty.

## 2. The e-μ control region (F1)

### What is observed

`mumu_CRemu` (one tight trigger-matched muon, one medium electron, OS, 60 < m(eμ) < 120): 72 357 data
events versus 66 179 predicted, +9.3 %. The excess is not uniform:

| variable | where the excess sits | data / MC |
|---|---|---|
| electron p_T | 20–40 GeV (5 300 of the 6 150 excess events) | 1.27 at 20–25 GeV, 1.02 above 55 GeV |
| lepton-cleaned jets | 0 jets: +5 100 events; ≥ 2 jets (tt̄-dominated): 0.99–1.00 | 1.31 / 1.07 / 0.99 / 1.00 for 0/1/2/3 jets |
| p_T^miss | 10–50 GeV | 1.10–1.18 |
| m(eμ) | 75–90 GeV | 1.14–1.18 |
| N_PV | rising with N_PV | 0.98 → 1.2 (this part is F6, the pileup profile) |

Low electron p_T, no jets, moderate p_T^miss and a visible mass just above the Z→ττ→eμ region is the
signature of a **non-prompt electron** accompanying a real W/Z muon, not of a normalisation problem of
tt̄ or WW (which populate ≥ 2 jets and agree to 1 %).

### Why it happens

`histograms.fill_chunk` keeps MC events in every region only if both leptons are prompt
(`genPartFlav` ∈ {1, 15}). In the SR this is correct because the fake-factor estimate replaces the removed
non-prompt muons. In the e-μ region nothing replaces the removed events, and no data-driven
non-prompt-electron estimate exists (`handoff.md` open issue 2 already suspected this). The removal is
therefore an *asymmetric* treatment: the prediction is missing a whole background.

### Evidence

The review pass re-runs the e-μ selection on all skims, keeping the non-prompt MC events separately and
adding the same-sign e-μ region (`review/review_pass.py`, results in `review/review_numbers_emu.json`):

| OS e-μ | events |
|---|---:|
| data | 72 357 |
| prompt-prompt MC (as in v2) | 66 232 |
| **data − prompt MC** | **6 125 (+9.2 %)** |
| MC non-prompt, removed by the prompt requirement: W+jets (jet → e) 3 412, Z→μμ + γ→e (`genPartFlav` 22, FSR/conversion) 1 344, Z→μμ other 660, Z→ττ with τ_h → e 627, tt̄ 274, tW 54, dibosons 27 | **6 412** |
| data / (prompt + non-prompt MC) | **0.996** |

With the non-prompt MC added back, data/MC is 1.00 / 1.00 / 0.99 / 0.99 in the 0/1/2/3 lepton-cleaned
jet bins, 1.00 ± 0.01 in every electron-p_T bin below 60 GeV, and flat in m(eμ)
(`review/figures/emu_os_with_nonprompt.pdf`). The same-sign e-μ region, which was not built in v2, gives
an independent, data-driven cross-check: 6 718 data events versus 1 039 prompt MC, i.e. 5 679
charge-symmetric non-prompt events; adding the OS-only sources (Z→ττ τ_h→e, top b→e: 955 events from MC)
gives 6 630, consistent with the 6 125 excess (the MC underestimates the SS region by 26 %, as expected
with no QCD sample and a W+jets skim of 0.03 GB, so the SS number is the more reliable one).

### Consequences and what to do

- The e-μ region is `Type: VALIDATION` in the fit, so the excess has **no effect on σ_fid**. It does mean
  the region does not validate anything as plotted, and that promoting it to a control region with
  `--emu-control` (`mu_top`) would be wrong: `mu_top` would absorb 9 % of non-prompt events into tt̄.
- With the non-prompt component included, the flavour-symmetric backgrounds (tt̄, tW, WW, Z→ττ; 58 000
  events = 0.56 % of the SR) are validated to about 2 %, i.e. 0.01 % on σ_fid.
- Fix: either keep the non-prompt MC in the e-μ region (cheapest, and the numbers above show it is
  adequate for a validation region) or estimate it from the same-sign e-μ region
  (N_OS^fake = R × (N_SS^data − N_SS^prompt MC), R ≈ 1, plus the OS-only τ_h→e term from MC). Add the
  same-sign region to `regions.py` (the review's `emu_both` function does this) and to the plots.
- Electron scale factors are absent, but they would move the prediction *down* (2016 medium-ID SFs are
  0.97–0.99), i.e. in the wrong direction; they are not the explanation.

## 3. p_T^miss in the signal region (F2)

### What is observed

Data/MC of PF p_T^miss in the SR rises monotonically from 0.88 at 0–5 GeV to 1.4 at 70–90 GeV; the mean
is 26.5 GeV in data and 25.1 GeV in MC. In Z → μμ events there is no genuine p_T^miss, so this is a
*resolution* difference: the MC recoil is too well measured.

### Evidence (review pass, `review/figures/sr_met_by_jets.pdf`, `sr_met_by_pileup.pdf`)

| selection | PF p_T^miss data/MC (shape) 0–30 / 30–60 / 60–90 GeV | Puppi p_T^miss |
|---|---|---|
| 0 lepton-cleaned jets (81 % of events) | 0.95 / 1.10 / 1.33 | **1.00 / 1.00 / 0.99** |
| 1 jet (15 %) | 0.93 / 1.09 / 1.41 | 0.97 / 1.07 / 1.20 |
| ≥ 2 jets (5 %) | 0.89 / 1.06 / 1.57 | 0.91 / 1.08 / 1.49 |

| N_PV (0 jets, PF p_T^miss) | mean data / mean MC |
|---|---|
| < 15 | 1.03 |
| 15–22 | 1.04 |
| ≥ 22 | 1.08 |

Three effects, all expected:

1. **Jet energy resolution.** NanoAOD MC jets and the Type-1 p_T^miss are not JER-smeared; the analysis
   applies no smearing. The mismatch grows with the number of jets, for PF and Puppi alike.
2. **Pileup / unclustered energy.** In 0-jet events Puppi p_T^miss (which suppresses pileup) agrees
   perfectly while PF p_T^miss does not, and the PF discrepancy grows with N_PV. Part of this is the
   pileup profile (F6: MC has ~4 % less pileup than data), the rest is the usual unclustered-energy
   mismodelling.
3. **Z p_T.** aMC@NLO+Pythia8 CP5 mismodels the low Z p_T (data/MC 0.83 at 0–5 GeV, 1.13 at 10–15 GeV,
   0.95 at 30–60 GeV); the hadronic recoil that balances the Z inherits this.

### Consequences

None for the measurement: p_T^miss enters only through the MET filters (applied identically to data and
MC) and the skim category C (single-muon + jet fake region, which is documented as biased and unused).
If p_T^miss is ever cut on (e.g. to suppress tt̄, or in the ττ channel that shares these conventions),
apply JER smearing to MC jets and propagate it to p_T^miss, or use Puppi p_T^miss, and add the
`MET_Unclustered` nuisance parameter already foreseen in `fitting/CONVENTIONS.md`. The Z-p_T
mismodelling should be covered by a Z-p_T reweighting or by an explicit systematic; the powheg-vs-aMC@NLO
comparison gives a first estimate of its effect on C (0.17 %).

## 4. Signal extraction: the fit is not robust (F3, F4, F5)

### The mass shape is not described by the nominal template

Pre-fit, data/MC in the SR is 0.95–0.97 for 62 < m < 78 GeV, 0.98–0.99 up to 86 GeV, and 1.00 ± 0.005
elsewhere (`review/figures/sr_mass.pdf`, `fit_sigmodel_shape.pdf`). Backgrounds are only 7 % of the
60–64 GeV bin and are validated by the e-μ region, so this is the Z/γ* lineshape of the aMC@NLO sample
(continuum, FSR tail and momentum resolution together). Nothing in the systematic model describes a
3–4 % low-mass deficit except the one-sided `SigModel` template (powheg/aMC@NLO ratio: 0.89–0.94 at
60–66 GeV, 0.95–0.97 up to 80 GeV, +2.5 % at the peak, and −10 %/−33 % in the two highest bins).

### What the fit does with it

The nominal fit pulls `SigModel` by +0.51σ and constrains it to 0.10σ; `MuonScale` and `MuonRes` are
constrained to 0.11σ and 0.21σ; `QCDScale` and `PS_FSR` to 0.5σ. Post-fit the ratio is flat and the
goodness of fit is p = 0.33. Removing the powheg template makes the fit look for another shape handle:

| fit configuration (`review/fitcheck/`) | μ_Z | GoF p | notable pulls / constraints |
|---|---:|---:|---|
| nominal: 30 × 2 GeV bins | 0.9903 +0.0141 −0.0138 | 0.33 | SigModel +0.5 (0.10), MuonScale (0.11), MuonRes (0.21), Pileup −1.5, PDF +1.3, MuonIso −1.1 |
| 6 × 10 GeV bins | 0.9959 ± 0.0160 | 0.0009 | MuonIso −1.7, PS_ISR +1.3, PS_FSR +1.2 |
| 1 bin (= counting (N − B)/(C L) = 794.4 pb) | 0.9935 | — | Hesse errors degenerate; grouped total 1.49 % |
| 30 bins, no `SigModel` | 1.0050 +0.0141 −0.0137 | 0.0014 | **MuonIso −3.0σ**, QCDScale −0.7, MuonScale (0.06) |
| 30 bins, no `SigModel`, no smoothing | 1.0063 +0.0141 −0.0137 | 0.0003 | MuonIso −3.1σ, Pileup (0.48), Lumi (0.22) |

The spread of μ_Z is 0.990–1.006, i.e. ±0.8 %, larger than the entire non-luminosity systematic
(0.77 %). Without the powheg template the fit tilts the `MuonIso` variation (whose effect is +1.1 % at
60 GeV falling to +0.2 % at 120 GeV) by three standard deviations to mimic the tail deficit. In other
words, in the nominal fit the low-mass tail is "explained" by the generator difference, in the others by
a muon-isolation efficiency 1 % lower than the tag-and-probe measurement; the data cannot tell which,
and the choice moves σ_fid by 12 pb.

Why the constraints are artefacts: 10 M events in 2 GeV bins resolve shape differences of 10⁻³ per bin,
which is finer than the MC statistical precision of the variation templates (MC per-bin uncertainty is
1.7× the data one, F9) and finer than the physics content of the templates. A `Lumi` constraint of 0.82
means the fit has "measured" the luminosity via the 68 000 tt̄+diboson events' sideband shape: not
credible, and it lowers the quoted luminosity impact from 1.20 % to 1.16 %. The `Pileup` pull of −1.5σ
(towards *less* pileup) contradicts the N_PV distribution, which wants *more* (F6): the mass shape has no
pileup information, the NP is pulled by unrelated template noise.

### Why the `SigModel` template is wrong as built (F4)

The powheg sample `ZToMuMu_M-50To120` has a generator-level cut at m_LHE < 120 GeV; aMC@NLO M-50 has
none. After resolution smearing and FSR, reconstructed masses at 116–120 GeV in aMC@NLO come partly from
generated masses above 120 GeV, which the powheg sample lacks — hence −10 % and −33 % in the two highest
bins. That is a phase-space artefact, not generator modelling, and the fit treats it as a free shape.
In addition, a 6–11 % difference in the 60–70 GeV continuum between two NLO generators is larger than
usual and should be understood (PDF set of the powheg sample, photon-induced contribution, LHE
mass definition) before it is used as a one-sided systematic.

### Recommendation

1. Rebuild `SigModel` with the same LHE mass window applied to both generators (e.g. 50 < m_LHE < 120 on
   the aMC@NLO sample too), or restrict the comparison to a region without edge migration.
2. Investigate the 62–78 GeV deficit directly: compare the aMC@NLO and powheg *generator-level* dressed
   lineshapes with each other and with the data after unfolding the resolution (the momentum calibration
   already provides it). If the data prefer powheg, consider it as the nominal signal template or assign
   the full difference as a two-sided lineshape uncertainty.
3. Until then, quote the counting extraction (1-bin likelihood or (N − B)/(C L) = 794.4 pb) as the
   nominal, keep the 30-bin fit as a cross-check, and add half the spread of the table above (±0.7 %)
   as a lineshape-model systematic. The counting result is insensitive to the lineshape to first order
   (C is a ratio of the same MC inside the same mass window).
4. If the shape fit is kept: use coarser bins (5–10 GeV), remove the smoothing on the momentum
   templates (or smooth all), decorrelate `Lumi` from the backgrounds or quote the external 1.2 %, and
   report pre-fit (external) uncertainties for `Lumi`, `Pileup` and `MuonIso` rather than the profiled
   ones.

## 5. Pileup profile (F6)

`zmumu/pileup.py` builds the data true-pileup profile from the per-lumisection `avgpu` column of the
CMS Open Data luminosity CSV (σ_mb = 80 mb) scaled to 69.2 mb, weighted by recorded luminosity in the
certified lumisections of runs 278820–284044. This is a sensible replacement for the unreachable
`puWeights` file, but it is missing the bunch-to-bunch spread within a lumisection (the official
`pileupCalc.py` "true" mode adds it) and the resulting profile is both too narrow (rms 6.2) and too low:

| | mean N_PV | mean N_true |
|---|---:|---:|
| data SR | 18.22 (era G 17.72, era H 18.66) | 25.6 (from the MC relation N_PV = 2.16 + 0.63 N_true) |
| DY MC, no weight | 16.26 | 21.9 |
| DY MC, nominal weights | 17.72 | 24.6 |

Data/MC in N_PV rises from 0.80 at N_PV = 3 to 1.2 at N_PV = 28 in the SR and in the e-μ region
(`review/figures/sr_npv_pileup_check.pdf`; a 3 % higher profile flattens most of it). The fit's
`Pileup` NP (σ_mb ± 4.6 %) is pulled the *other* way, which shows that the NP is not doing its job.

Effect on σ_fid: small by construction — the isolation efficiency is measured in data by tag-and-probe
in (p_T, |η|) bins, and the scale factor changes by only 0.5 % over the full N_PV range
(`tnp_iso_vs_npv`), so a 4 % shift of the mean pileup changes the average SF by ~0.03 %. Fix by using
the official 2016 pileup JSON when it becomes reachable, or by scaling the profile to reproduce the
data N_PV distribution (an "N_PV-matched" profile) and assigning the difference as the systematic.

## 6. Analysis choices reviewed

| choice | verdict | comment |
|---|---|---|
| tight ID + `pfRelIso04_all` < 0.15, IP cuts | good | standard; the reviewer's reference uses the same; T&P SFs 0.97–0.99 (ID), 1.00–1.04 (iso) are in the expected range |
| fiducial volume: dressed muons (ΔR < 0.1), p_T > 26/20, \|η\| < 2.4, 60–120 GeV | good | dressed avoids the 1.8 % Born/dressed ambiguity; document that the reco FSR recovery uses ΔR < 0.5 while the gen dressing uses 0.1 (absorbed in C, not a problem) |
| own skims from the parent NanoAOD (TrigObj, prefiring, 1-muon events) | good | required for the trigger T&P and the e-μ region; laptop bundle provided |
| signal normalisation 6077.22/3 pb, L = 16393.381 pb⁻¹ | good | normtag value now used everywhere; the fit floats μ so the DY cross section only sets the reference |
| pileup from the luminosity CSV | acceptable stop-gap | see §5 |
| L1 prefiring from `L1PreFiringWeight_Nom` | good | 2.0 % correction, 0.5 % uncertainty; the T&P does not double count it (prefired events are lost entirely and never enter the T&P sample) |
| T&P: trigger-matched tag, both orientations, pass/fail simultaneous fits with MC templates, SS-subtracted counting cross-check | good | closure −0.1 %/−0.2 % (ID/iso); the ID fits have χ²/ndf ≈ 2 in data (template resolution), covered by the cmsshape/range/tag variants |
| trigger efficiency per muon, event weight 1 − (1 − ε₁)(1 − ε₂) | good | plateau 0.907 data / 0.923 MC; the L1-seed veto variant is in the data error |
| `MuonReco` 0.4 % | inconsistent | per muon or per event? see F8 |
| fake factor from the same-sign tag+probe region, prompt subtraction, anti-tight = tightId with 0.20 < iso < 1.0 | good | the single-muon + jet region is correctly rejected (online isolation bias); the two-template fit in the application region is a reasonable way to handle the Z-dominated anti-isolated sample; closure 2077 ± 20 vs 1922 ± 85 |
| fake-factor uncertainty floored at 30 % | fine | fakes are 0.035 % of the SR |
| prompt-prompt requirement on MC in SR/SS | good | 0.006 % of DY SR events fail the matching (checked) |
| prompt-prompt requirement on MC in the e-μ region without replacement | wrong | F1 |
| Z-peak momentum calibration (κ, extra smearing per \|η\|) | good | replaces the unavailable Rochester corrections; but the same Z peak then constrains `MuonScale`/`MuonRes` in the fit to 0.11σ/0.21σ — the calibration is used twice |
| theory variations renormalised to constant fiducial yield (vary C only) | good | correct for a fiducial measurement; the acceptance uncertainty (0.61 %) is added only for σ_tot |
| TRExFitter shape fit, 30 × 2 GeV bins, 24 NPs + 30 gammas | not robust | §4 |
| e-μ region as VALIDATION (not CONTROL) | good decision | and it must stay so until F1 is fixed |
| powheg `SigModel`, one-sided | wrong as built | F4 |
| acceptance: NLO, 60 < m_Born < 120 denominator recommended | good | A = 0.4092 reproduces the external reviewer's number; the LO/NLO 3.2 % difference is a convention issue for the combination, correctly not a systematic |
| jet counting without lepton cleaning | bug (plots only) | F7 |

## 7. Smaller items

- `RESULTS_v2.md` yields table lists `DYmumu_powheg` with raw generator weights (1.98 × 10⁹ "events");
  normalise or drop it from the table. `WJets` in the SR is −5 events (negative-weight fluctuation, not
  used in the fit) — say so or drop it.
- `err_decomp` in `zmumu_fit_result.json` (TOT_ERROR 0.0116, STAT_ERROR 0.0065) is inconsistent with the
  MINOS result (0.0140) and the stat-only fit (0.0003); the errDecomp step ran with a different setup.
  Remove it from the JSON or explain it.
- The `Fakes` template carries Sumw2 = (σ_b · F)², a fully correlated normalisation error, which the
  fit treats as per-bin gammas on top of `FakeStat_mumu` (double counting, negligible).
- The T&P fit closure and the FF closure outliers at 112 and 119 GeV
  (`ff_closure_and_template.png`: −45 and +52 events in single 0.5 GeV bins) are negative-weight
  fluctuations of the prompt subtraction; harmless but worth a caption.
- `handoff.md` quotes the 0.03 % statistical uncertainty; the effective statistical limit of the fit is
  the MC (0.35 %). Say so in the combination inputs so nobody plans on "more data".
- The luminosity NP should be quoted as the external 1.2 %, not the profiled 1.16 %.

## 8. Recommendations, in order of impact on the deliverable

1. **Signal extraction (F3, F4, F5)**: fix the `SigModel` template, understand the 62–78 GeV deficit,
   and either move to a counting-type extraction with an explicit lineshape systematic or to a coarse
   binning with the shape NPs decorrelated from the normalisation. Expected effect: central value
   between 792 and 805 pb, an additional ±0.7 % systematic until the tail is understood.
2. **e-μ region (F1)**: keep the non-prompt MC (or add the same-sign estimate) so the region validates
   the flavour-symmetric backgrounds; add the SS e-μ region to the plots.
3. **`MuonReco` (F8)**: decide per-muon or per-event and correct either the fit config or the documentation.
4. **Pileup (F6)**: N_PV-matched profile or the official JSON; re-check the `Pileup` pull afterwards.
5. **Lepton-cleaned jets (F7)** in `histograms.py` (ΔR > 0.4 to the selected leptons).
6. p_T^miss: nothing to do for μμ; JER smearing or Puppi p_T^miss and `MET_Unclustered` if p_T^miss is ever cut on.
7. Consider a Z-p_T reweighting (or an explicit NP) — the low-p_T mismodelling is visible and enters C
   through the p_T thresholds at the 0.1–0.2 % level.

## 9. How this review was done (reproducible)

```bash
source ../setup.sh
python review/review_pass.py 8            # e-mu OS/SS with MC flavour split, SR pileup/MET 2D histograms (~7 min, 8 workers)
python review/deck_plots.py               # all figures in review/figures/ (dark deck style, PDF + PNG)
# fit variants: review/fitcheck/*.config (rebin 5, rebin 30, no SigModel, no SigModel + no smoothing),
#   run with  trex-fitter h/w/f/i <config>  from a directory with the fitinputs path adjusted; the
#   resulting Fits/*.txt and GroupedImpact files are archived next to the configs.
cd review/deck && ./build.sh              # LaTeX slides (uses the cvmfs TeX Live 2025 binaries)
```

Numbers quoted above: `review/review_numbers_emu.json`, `review/fitcheck/*.txt`,
`output/v2/RESULTS_v2.md`, `fit/results/zmumu_fit_result.json`.
