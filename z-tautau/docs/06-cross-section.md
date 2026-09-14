# 06 — From a signal strength to a cross section

## What is measured

The fit (`08-fit-and-results.md`) returns μ_Z, the signal strength of the Z/γ*→ττ template relative to the
NNLO prediction σ(Z/γ*→ℓℓ, m > 50 GeV) = 6077.22 pb / 3 per flavour. The template is aMC@NLO FxFx with all
TauPOG corrections. μ_Z multiplies the *whole* ττ template, including events outside any fiducial volume.
The POI name `mu_Z` and the normalisation are those of every channel (`fitting/CONVENTIONS.md`), so the
combination fits one μ_Z.

## Fiducial volume

Generator level, in LHE Z/γ*→ττ events (`ztautau/gen.py`):

* 60 < m_LHE(ττ) < 120 GeV,
* **both τ decay hadronically**,
* both visible τ (`GenVisTau`, the τ four-momentum minus neutrinos) with pT > 40 GeV and |η| < 2.1.

It follows the reconstruction-level cuts, so the extrapolation from the measured phase space is small.
m_LHE (pre-shower, Born-like) is used for the mass window because the visible di-τ mass is not a physical
Z mass.

| quantity | value |
|---|---:|
| σ(Z/γ*→ττ, m > 50 GeV), prediction | 2019.7 pb |
| σ(Z/γ*→ττ, 60 < m < 120 GeV), prediction | 1944.9 pb |
| σ_fid, prediction | 4.501 pb |
| A = σ_fid / σ(60–120) | 0.002314 |
| C = N_SR(signal) / (σ_fid · L) | 0.107 |
| fraction of SR signal outside the fiducial volume | 38 % |

A = 0.23 % is B(τhτh) ≈ 42 % times the chance that both visible τ carry pT > 40 GeV (≈ 0.55 %). A τ of
45 GeV gives its visible products on average ~half its momentum, so only strongly boosted or asymmetric
decays pass. This is why a τhτh channel with 16 fb⁻¹ has "only" ~8 000 signal events.

Uncertainties on A (from the generator-weight sums over all generated events): QCD scales 3.4 %, PDF 0.5 %,
α_s 0.6 %, ISR 0.9 %, FSR 0.2 %, MC statistics 0.7 % → 3.7 %. The scale dependence is large because the
double 40 GeV cut on the visible τ selects the Z pT tail.

## Formulae

```
σ_fid                     = μ̂ · σ_fid^pred                          (theory NPs in the fit move only C)
σ(Z/γ*→ττ, 60<m<120 GeV)  = μ̂ · σ_pred(60–120)  = σ_fid / A        (+ the A uncertainty above)
```

The theory variations (QCD scales, PDF, PS) are applied to the signal template **renormalised to a
constant fiducial yield**, so inside the fit they change only the reconstruction efficiency and the
migrations (C). Their effect on the extrapolation (A) is quoted separately. The generator comparison
(`SigModel`, madgraph LO vs aMC@NLO, both normalised to the NLO fiducial cross section) is a
normalisation-only nuisance: C differs by 7.3 %, mostly through the out-of-fiducial migration, which depends
on the slope of the visible-τ pT spectrum at 40 GeV.

## Counting cross-check

With prefit numbers, σ = (N_obs − N_bkg) / (A·ε · L), with A·ε = N_SR^sig / (σ_pred(60–120) · L) =
2.47·10⁻⁴: (47 586 − 40 235) / (2.47·10⁻⁴ · 16 393.4) = 1815 pb. The fit gives more (μ̂ ≈ 1.16) because the
high-m_tt sideband pulls the fake normalisation down by ~3 %. That moves ~1 100 events from the fake
template (38 700 events) to the signal (7 900 events). The ratio of the two templates makes the counting
method useless here; the shape fit is essential.
