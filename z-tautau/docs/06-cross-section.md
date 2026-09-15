# 06 — From a signal strength to a cross section

## What is measured

The fit (`08-fit-and-results.md`) returns μ_Z, the signal strength of the **fiducial** Z/γ*→ττ template
relative to the NNLO prediction σ(Z/γ*→ℓℓ, m > 50 GeV) = 6077.22 pb / 3 per flavour. The template is
aMC@NLO FxFx (inclusive + jet-binned samples, `01-data-and-samples.md`) with all TauPOG corrections. The
POI name `mu_Z` and the normalisation are those of every channel (`fitting/CONVENTIONS.md`), so the
combination fits one μ_Z.

## Fiducial volume

Generator level, in LHE Z/γ*→ττ events (`ztautau/gen.py`):

* 60 < m_LHE(ττ) < 120 GeV,
* **both τ decay hadronically**,
* both visible τ (`GenVisTau`, the τ four-momentum minus neutrinos) with pT > 40 GeV and |η| < 2.1.

It follows the reconstruction-level cuts, so the extrapolation from the measured phase space is small.
m_LHE (pre-shower, Born-like) is used for the mass window because the visible di-τ mass is not a physical
Z mass; this is the Born-level 60–120 GeV denominator the combination agreed on (the ee/μμ channels must
use a Born-level window for σ(60–120) as well; their *fiducial* volumes use dressed leptons, which is a
different quantity).

| quantity | value |
|---|---:|
| σ(Z/γ*→ττ, m > 50 GeV), prediction | 2019.7 pb |
| σ(Z/γ*→ττ, 60 < m < 120 GeV), prediction | 1944.9 pb |
| σ_fid, prediction | 4.501 pb |
| A = σ_fid / σ(60–120) | 0.002314 |
| C = N_SR(fiducial signal, all categories) / (σ_fid · L) | 0.0654 |

A = 0.23 % is B(τhτh) ≈ 42 % times the chance that both visible τ carry pT > 40 GeV (≈ 0.55 %). A τ of
45 GeV gives its visible products on average ~half its momentum, so only strongly boosted or asymmetric
decays pass. This is why a τhτh channel with 16 fb⁻¹ has "only" ~5 000 fiducial signal events.

Uncertainties on A (from the generator-weight sums over all generated events): QCD scales 3.4 %, PDF 0.5 %,
α_s 0.6 %, ISR 0.9 %, FSR 0.2 %, MC statistics 0.7 % → 3.7 %. The scale dependence is large because the
double 40 GeV cut on the visible τ selects the Z pT tail.

## The non-fiducial part of the selected Z/γ*→ττ (new in v2)

The selected Z/γ*→ττ events are only 62 % fiducial:

| component of the selected Z/γ*→ττ | fraction |
|---|---:|
| inside the fiducial volume → sample `DYtautau`, scaled by μ_Z | 62 % |
| m_LHE > 120 GeV | 30 % |
| 60–120 GeV, both genuine τh, a visible τ below 40 GeV or beyond \|η\| = 2.1 at generator level (migration) | 6 % |
| 60–120 GeV, a leg that is not a hadronic τ (τ→e/μ misidentified as τh) | 1.5 % |
| m_LHE < 60 GeV | 0.2 % |

The two 40 GeV visible-pT cuts have an acceptance of 0.23 % at the Z peak but tens of percent at
m > 200 GeV, so the γ* continuum is enriched by two orders of magnitude relative to the peak, and 88 % of
the Z/γ*→ττ events with m_tt > 130 GeV — the bins that fix the fake normalisation — have m_LHE > 120 GeV.
In v1 all of this was one template scaled by μ_Z, i.e. the measured "σ(60–120)" was a measurement of the
whole DY spectrum weighted by the SM shape, and the fake sideband contained μ-dependent signal
(`REVIEW.md` 3.1).

In v2 the non-fiducial part is the sample **`DYtautau_nonfid`**, a background with

* a 5 % normalisation uncertainty (`XS_DYtautau_nonfid`: the NNLO/NLO and electroweak knowledge of the
  high-mass Drell-Yan continuum),
* the theory variations `QCDScale`, `PDF`, `PS_ISR`, `PS_FSR` with the *same* per-event multipliers as
  the signal (renormalised to a constant fiducial yield), so that they move the non-fiducial / fiducial
  ratio: that is precisely the theory uncertainty of extrapolating from the fiducial to the non-fiducial
  phase space,
* all experimental nuisance parameters, correlated with the signal.

The BDT (`09-bdt.md`) is trained on the fiducial part only; 84 % of the non-fiducial events land in the
fake-dominated category.

## Formulae

```
σ_fid                     = μ̂ · σ_fid^pred                          (theory NPs in the fit move only C)
σ(Z/γ*→ττ, 60<m<120 GeV)  = μ̂ · σ_pred(60–120)  = σ_fid / A        (+ the A uncertainty above)
```

The theory variations (QCD scales, PDF, PS) are applied to the signal template **renormalised to a
constant fiducial yield**, so inside the fit they change only the reconstruction efficiency and the
migrations (C). Their effect on the extrapolation (A) is quoted separately. The generator comparison
(`SigModel_tautau`: madgraph LO vs aMC@NLO, both normalised to the NLO fiducial cross section) gives
C_LO / C_NLO = 0.867: the LO sample's softer visible-τ spectrum inside the fiducial volume meets the
trigger turn-on. It is reported and not used as an uncertainty (`07-corrections-and-systematics.md`).

## Counting cross-check

With prefit numbers, σ = (N_obs − N_bkg) / (A·ε · L), with A·ε = N_SR^sig / (σ_pred(60–120) · L): see
`output/RESULTS.md` ("Inputs for the combination"). The counting number is useless here (B/S ≈ 8 over the
whole signal region); only the shape fit in the categories constrains the fake normalisation.
