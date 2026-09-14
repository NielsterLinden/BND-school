# 04 — The di-τ mass: correcting the visible mass with MET

Code: `ztautau/mass.py`. In Z→τhτh at least two neutrinos escape, so the visible mass m_vis of the two τh
peaks at ~0.8 m_Z with a long tail. The missing transverse momentum carries the neutrinos' information.

## Three estimators

**Visible mass** m_vis — the invariant mass of the two τh four-vectors.

**Collinear approximation** m_col. The τ are boosted, so each neutrino system is collinear with its τh.
With x_i the visible momentum fraction of τ i, MET = (1/x₁ − 1) pT^vis₁ + (1/x₂ − 1) pT^vis₂ is two equations
for two unknowns and m_col = m_vis / √(x₁x₂). It is exact in the collinear limit, but the system is singular
for back-to-back τ — the typical Z topology — and any MET mismeasurement pushes x out of (0, 1]:
**m_col is defined in only 38 % of signal events**.

**MET-likelihood mass** m_tt (the fit variable; same idea as CMS FastMTT / "SVfit-lite", ATLAS MMC).
Instead of solving, scan (x₁, x₂) on a 60 × 60 grid over the physical range x_i ∈ [m_vis,i²/m_τ², 1] and weigh
every point with

* the **MET transfer function** exp(−½ rᵀV⁻¹r), r = MET_measured − MET_predicted(x₁, x₂), with V the
  per-event MET covariance matrix from NanoAOD (`MET_covXX/XY/YY`);
* the **two-body phase space** of τ → τh ν, which is flat in x;
* a **1/m² prior** on m = m_vis/√(x₁x₂), mimicking the falling Drell–Yan spectrum. Without it the weakly
  constrained low-x corner (large neutrino momenta hidden in the MET resolution) pulls the estimate up.

The estimate is the posterior **median** of m. It is always defined and uses the MET resolution event by
event (vectorised numpy). When the τ energy scale or the MET is varied for a
systematic, m_tt is recomputed (`analysis.kinematics`).

## Performance on Z→ττ simulation

Genuine τhτh in the signal region, 70 < m_LHE < 110 GeV (m_LHE = generator ττ mass):

| estimator | median m / m_LHE | core resolution (IQR/2) | 68 % half-width | defined |
|---|---:|---:|---:|---:|
| m_vis | 0.80 | 13.1 % | 18.3 % | 100 % |
| m_col (defined events only) | 1.05 | 16 % | — | 38 % |
| m_tt, posterior mean, no prior | 1.04 | 11.9 % | 19.1 % | 100 % |
| **m_tt, posterior median, 1/m² prior (used)** | **0.99** | **10.9 %** | **16.0 %** | 100 % |
| m_tt, maximum a posteriori, 1/m² prior | 0.91 | 9.1 % | 14.1 % | 100 % |

Relative to m_vis the MET-corrected mass has the correct scale and ~15 % better resolution for the signal.
For **separating the signal from the dominant fakes** the two are, honestly, equivalent here. The best
windows give Z→ττ / fakes = 2980 / 1165 (S/B 2.6, 38 % of the signal) for 70 < m_tt < 100 GeV and
2822 / 1212 (S/B 2.3, 36 %) for 50 < m_vis < 80 GeV. The stat-only expected precision on μ with fixed
backgrounds is 1.83 % (m_tt) and 1.82 % (m_vis). With pT > 40 GeV on both τh, the selected Z→ττ events
already have hard visible decay products and little MET, so the MET correction adds less than in a
low-threshold ℓτh selection. m_tt is used because it is the physical ττ mass the prompt asked for
(peak at m_Z). Its fake template is also well separated: fakes peak at m_tt ≈ 130–170 GeV, which gives a
clean normalisation region for the fit (`08-fit-and-results.md`). Computing it costs ~1 800 events/s per
core.

MAP has the best resolution but a 9 % scale bias and grid artefacts; the median was preferred for a
peak that sits at m_Z. The choice does not change the result: the templates come from the same estimator in
data and simulation.
