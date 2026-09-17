# Systematic uncertainties compared with the published ATLAS and CMS measurements

The question is whether our three channels carry the same systematic uncertainties as the published
measurements of the same cross section, and which ones are missing.

**Our numbers** are the relative change of each channel's signal yield in the fitted bins (%), for a
±1σ variation as it enters the combined likelihood. They come from `checks/systematics.py` →
`checks/systematics.json`, and are the same quantity ATLAS tabulates as δC/C. Where the channel quotes
σ(60–120) = σ_fid / A, the acceptance uncertainty (δA/A) is listed separately. The μμ column is the frozen
result (17 Sep 2026): the reconstruction SF is measured and the acceptance block is computed for the 60–120 GeV
denominator (`z-mumu/zmumu/acceptance.py`, `z-mumu/docs/16`).

**The published numbers**, verified against the papers on 16 Sep 2026:

* **ATLAS ee / μμ:** arXiv:1603.09222 (PLB 759 (2016) 601), 81 pb⁻¹, 66–116 GeV, Table 1 (δC/C).
  δA/A is 1.8 % in total; its breakdown is not published.
* **CMS ℓℓ:** CMS-SMP-20-004, arXiv:2408.03744 (JHEP 04 (2025) 162), 206 pb⁻¹, 60–120 GeV, Table 7.
  These are uncertainties on σ_tot from one ee+μμ fit; the paper publishes no per-channel numbers.
* **CMS ττ:** arXiv:1801.03535 (EPJC 78 (2018) 708), 2.3 fb⁻¹, 60–120 GeV, Table 2. The first number
  is the change in yield or acceptance, the number in brackets the impact on σ.

> **17 Sep 2026 — the ττ column below is the frozen τhτh-only channel (v3).** The combination now uses the
> four-channel ττ (τhτh + μτh + eτh + eμ). What that changes in this comparison:
> * **Point 4 of "What is missing" is resolved.** Like CMS, z-tautau now fits several final states together and
>   measures the τh ID scale factors in situ (0.99 ± 0.04, 0.96 ± 0.03, 0.89 ± 0.03, 0.79 ± 0.05 for DM0/1/10/11 in
>   the combined fit; TauPOG: ± 0.05–0.15) and the τh energy scales to 0.6–2.2 %. Its own uncertainty went from
>   ~11 % to 3.6 %; CMS quotes 3.1 %.
> * **Point 5 no longer applies to ττ:** its theory templates are renormalised to a constant σ(60–120), so
>   `PDF`/`QCDScale`/`PS_*` move A × ε together inside the fit and there is no separate ττ `Acc_*` block.
> * **New, and the largest ττ term: `TauIDpT_tautau`, ± 6.3 %** — the p_T dependence of the τh ID scale factor that
>   z-tautau measured with a p_T-split fit and that its nominal model does not cover (`../README.md`). CMS covers
>   the same physics with p_T-binned tag-and-probe scale factors.
> * Input-level sizes on the ττ signal yield, summed over the 13 regions (`../checks/systematics.json`, in %):
>   `TauIDpT_tautau` 6.3, `TauES_DM1` 2.5, `MET_Unclustered` 2.1, `EmuTrigger` 1.8, `Pileup` 1.6, `QCDScale` 1.3,
>   `Lumi` 1.2, `TauES_DM0` 1.0, `TauES_DM10` 0.8, electron reco / ID 0.8 / 0.7, `PDF` 0.6, `PS_ISR` 0.6,
>   muon isolation / ID 0.4 / 0.3. The τh ID itself is a free factor per decay mode, not a prior.
> * New in ττ with no analogue in ee/μμ: `BTag`, `JES`, `TopPt` (the eμ tt̄ control region and b-jet vetoes), the
>   eμ and single-electron trigger efficiencies measured in situ, fake-factor parameters per final state.

## Side by side (%)

| source | ATLAS ee | ATLAS μμ | CMS ℓℓ | **ours ee** | **ours μμ** | CMS ττ | **ours τhτh** |
|---|---|---|---|---|---|---|---|
| luminosity | 2.1 | 2.1 | 2.3 | 1.2 | 1.2 | 2.3 (1.9) | 1.2 |
| lepton trigger | 0.1 | 0.2 | in efficiency | **none** | 0.19 | τh ID+trigger 6–12 (1.5) | τh trigger 0.8–3.1 per DM |
| lepton reconstruction / ID | 0.9 | 0.9 | 0.28 (stat) ⊕ 0.16 (syst) | reco 0.62, **ID 5.9** | reco 0.27 (measured, T&P), ID 0.21 | τh ID in the line above (constrained in situ to 2.2) | τh ID 1.9–8.0 per DM |
| isolation | 0.3 | 0.5 | in efficiency | in the cut-based ID | 0.36 | – | – |
| energy / momentum scale, resolution | 0.2 | 0.1 | shape NPs, "small" | shape only (0 on yield) | 0.02–0.03 | τh ES 2–17 (< 0.1), e/μ ES < 1 | τh ES 0.07–1.8 |
| charge misidentification | 0.1 | – | – | **none** | – | – | – |
| pileup | < 0.1 | < 0.1 | in efficiency | 0.79 | 0.26 | – | 1.0 |
| L1 prefiring | n/a | n/a | 0.34 | 0.47 | 0.53 | – | 0.29 |
| background cross sections | negligible | negligible | EW + tt̄ 0.04 | tt̄ 6, VV 10, W+jets 10, Z→ττ 5 (on 0.7 % of the yield) | tt̄, tW, VV, Z→ττ | tt̄ 7 (1.0), VV 15 (0.2), Z→ee/μμ (1.8) | as CMS, plus DY low mass, non-fiducial Z→ττ |
| multijet / non-prompt | – | – | 0.14 (syst) ⊕ 0.07 (stat) | **none** | fake factor (0.04) | jet→τh 6–16 (< 0.1) | fake factor, the dominant background |
| MC statistics | – | – | 0.16 | per-bin γ | per-bin γ | – | per-bin γ |
| PDF (+ α_s) | 0.1 (C) + in A | < 0.1 (C) + in A | 0.43 | +0.11 / +0.10 (A·C) | 0.06 (C) + 0.50 (A), α_s 0.01 + 0.03 | 1 (1.0) | 0.68 (C) + 0.49 (A), α_s 0.59 (A) |
| μ_R, μ_F scales | in A | in A | 0.66 | +0.31 / −0.10 (A·C) | 0.18 (C) + 0.29 (A) | < 6 (0.5) | 2.7 (C) + 3.4 (A) |
| parton shower / FSR / generator | in A | in A | resummation + FSR 0.12 | FSR +0.10 / −0.17 (A·C) | ISR 0.13, FSR 0.04, powheg vs aMC@NLO 0.28 (C); boson p_T 0.38, powheg vs aMC@NLO 0.04, FSR 0.03, QED FSR 0.10 (estimate) (A) | UE + PS 1 (1.0) | FSR 1.1–2.2, ISR 0.35 (C); ISR 0.93, FSR 0.25 (A) |
| acceptance MC statistics | – | – | – | in the template | 0.05 | – | 0.71 |
| total acceptance | 1.8 | 1.8 | (in the rows above) | in the template | 0.70 | ~2–6 (in the rows above) | 3.7 |
| total systematic on σ, without luminosity | 1.0 ⊕ 1.8 | 1.1 ⊕ 1.8 | 0.90 | 5.5 (standalone fit) | 1.0 (standalone fit) | 3.1 | 10.6 (standalone fit, incl. stat.) |

(A·C) marks templates normalised to the LHE 60–120 GeV cross section: they vary acceptance and
efficiency together, and are then fitted.

## What they have that we do not

1. **Electron trigger efficiency (z-ee).** ATLAS assigns 0.1 %; CMS includes the trigger in its
   efficiency measurement. z-ee applies no trigger scale factor and no uncertainty for
   `HLT_Ele27_WPTight_Gsf`: its event weight is `genWeight × prefiring × pileup × RecoSF × IDSF`. The
   data/MC efficiency ratio of a single-electron trigger in 2016 is typically a few per cent below one.
   A missing correction of that size shifts σ(ee) directly, and would be consistent with the ee channel
   measuring 4.7 % below μμ in its own fit.
2. **Charge misidentification (z-ee).** ATLAS assigns 0.1 % for ee and nothing for μμ. z-ee requires
   opposite charges and assigns nothing.
3. **Multijet / non-prompt background (z-ee).** CMS: 0.14 % ⊕ 0.07 %. z-mumu and z-tautau estimate it
   from data; z-ee has neither an estimate nor an uncertainty. At the ee peak it is small.
4. **An in-situ τh identification constraint (z-tautau).** CMS fits five ττ final states together
   and constrains the τh ID scale factor to 2.2 %, so the τh ID contributes 1.5 % to its σ. Our τhτh-only
   channel carries the TauPOG uncertainties (2–8 % per decay mode). This is most of why our ττ channel
   is at ~11 % and CMS's at 3.1 %.
5. **A common treatment of acceptance and efficiency theory uncertainties.** CMS fits them jointly;
   ATLAS puts them all in A. Our μμ and ττ channels vary C inside the fit (`PDF`, `QCDScale`, …) and A
   through separate `Acc_*` parameters, so a scale variation is not forced to move A and C together. At
   0.3 % (μμ) that is immaterial; at 3.4 % (ττ) it is not, but ττ carries only a few per cent of the
   combination's weight (its uncertainty is ~7× larger than the combined one).

## What we carry that they do not, or with a different size

* **Electron ID ±5.9 % (z-ee): an artefact, not a physics uncertainty.** The official UL2016postVFP
  scale-factor map (`datasets/corrections/.../POG/EGM/2016postVFP_UL/electron.json.gz`,
  `UL-Electron-ID-SF`, `Medium`) gives 0.3–0.5 % per electron at Z-peak p_T and ~1.5 % at 20–30 GeV,
  i.e. about 1.2 % per Z → ee event. In the barrel–endcap gap (1.444 < |η_SC| < 1.566) the same map
  returns a placeholder, sf = 1.0 with sfup = 2.0. CMS analyses veto that region; z-ee does not. A
  toy with Z-like electron kinematics gives 8.3 % with the gap and 1.2 % without. The 5.9 % in the z-ee
  templates is therefore the ~5 % of events with a gap electron varying by ±100 %. What this does to the
  fit, and why the combination splits the ee shape and normalisation parameters, is in `README.md`.
* **L1 prefiring** (0.3–0.5 %). ATLAS has no such effect; CMS assigns 0.34 %. Consistent.
* **Pileup** (0.3–1.0 %) is larger than ATLAS's < 0.1 %. The 2016 G+H pileup (~25 interactions per crossing) is much higher
  than in the 2015 ATLAS data or the CMS low-pileup run.
* **Generator comparison in μμ** (`SigModel_mumu`, powheg vs aMC@NLO, 0.28 % on the yield in the fit; 0.04 %
  on the acceptance, `Acc_Generator_mumu`). Together with the boson-p_T row (`Acc_PTZ_mumu`, 0.38 %: the
  generator p_T(Z) reweighted to the measured p_T(μμ)) and the QED FSR estimate (`Acc_QEDFSR_mumu`, 0.10 %)
  it plays the role of CMS's "resummation + FSR" (0.12 %).
* **MET unclustered energy** (ττ only). There is no analogue in a dilepton counting measurement.
* **Luminosity 1.2 %** is smaller than any published number: this is the final 2016 legacy
  calibration, against 2.1–2.3 % for the early 2015 and low-pileup 2017 datasets.

## Summary

The μμ channel has the same list of systematics as ATLAS and CMS, with comparable sizes. The ττ channel
is complete but not constrained in situ. The **ee channel is missing three systematics** (trigger, charge
misidentification, multijet), each small in the published analyses. It **also has one systematic ~5×
too large** because the ECAL gap is not vetoed. That last one, not a missing uncertainty, dominates how
much ee can contribute to the combination.
