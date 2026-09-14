# 07 — Corrections to the simulation and systematic uncertainties

## Corrections (`ztautau/corrections.py`, `analysis.weights`)

The official TauPOG inputs for UL2016 postVFP, DeepTau2017v2p1, are taken from the public GitHub
repositories `cms-tau-pog/TauIDSFs` and `cms-tau-pog/TauTriggerSFs` (branch `run2_SFs`). CMS cvmfs is not
mounted here and the CERN GitLab `jsonpog-integration` needs a login. `scripts/step0_external.py`
converts them to the committed `external/tau_pog_UL2016postVFP.json`.

| correction | applied to | value (UL2016 postVFP) |
|---|---|---|
| τh ID SF, VSjet Medium, **per decay mode** (TauPOG recommendation for τhτh, pT > 40) | genuine τh (`genPartFlav` 5) | DM0 0.923 ± 0.142, DM1 0.880 ± 0.054, DM10 0.868 ± 0.085, DM11 0.898 ± 0.163 |
| τh energy scale | genuine τh; propagated to MET and m_tt | DM0 0.993 ± 0.009, DM1 0.991 ± 0.007, DM10 1.001 ± 0.007, DM11 0.997 ± 0.016 |
| e→τh VSe VVLoose SF | prompt e / τ→e (flav 1, 3) | 1.06 ± 0.05 (barrel), 0.95 ± 0.06 (endcap) |
| μ→τh VSmu VLoose SF | prompt μ / τ→μ (flav 2, 4) | 0.88–1.04 per \|η\| bin |
| di-τ trigger leg SF(pT, DM) | every leg | ratio of the fitted data and MC efficiencies, e.g. 0.94–0.97 at 50 GeV |
| pileup | all MC | data profile from the Open Data lumi-by-LS table (record 1059), σ_mb = 69.2 mb; same construction as z-mumu v2 |
| L1 prefiring | all MC | `L1PreFiringWeight_Nom` (mean 0.990 for the signal) |

The pT-binned ID SF for pT > 40 GeV (0.823 +0.086 −0.036) was checked as an alternative. It is not
recommended for τhτh and its uncertainty is larger.

## Systematic uncertainties (nuisance parameters in the fit)

Names follow `fitting/CONVENTIONS.md`: identical strings are correlated between channels; `_tautau`
marks channel-specific parameters. "Prefit effect" = change of the Z→ττ (or fake) yield in the SR.

| NP | type | prefit effect | notes |
|---|---|---|---|
| `Lumi` | OVERALL, all MC | 1.2 % | correlated with ee, μμ |
| `Pileup` | HISTO | 1.3 % | σ_mb ± 4.6 %; correlated |
| `L1Prefiring` | HISTO | 0.2 % | correlated |
| `TauID_DM0/1/10/11` | HISTO | 8.9 / 5.7 / 3.6 / 2.5 % | TauPOG uncertainty, uncorrelated between DMs, correlated between the two legs |
| `TauTrigger_DM*` | HISTO | 2.5 / 3.2 / 2.1 / 0.9 % | fitted SF uncertainty vs pT |
| `TauES_DM*` | HISTO, shape (smoothed) | 1–2 % + shape | τ pT, MET and m_tt recomputed |
| `TauFakeEle`, `TauFakeMu` | HISTO | < 1.3 % on Z→ee | |
| `MET_Unclustered` | HISTO, shape | < 0.1 % norm. | NanoAOD unclustered-energy shift, m_tt recomputed |
| `QCDScale` | HISTO on Z→ττ, fiducial-normalised | ± 4.1 % | 6-point envelope per bin (no (0.5,2),(2,0.5)) |
| `PDF` | HISTO, fiducial-normalised | ± 1.0 % | 100 Hessian eigenvectors in quadrature per bin, + α_s |
| `PS_ISR`, `PS_FSR` | HISTO, fiducial-normalised | < 1 % | PSWeight 2 / 0.5 |
| `SigModel` | normalisation only, symmetrised | 7.3 % | madgraph LO vs aMC@NLO C factor |
| `XS_DYll` 5 %, `XS_DYlowmass` 10 %, `XS_TTbar` 6 %, `XS_SingleTop`/`XS_WW`/`XS_WZ`/`XS_ZZ`/`XS_WJets` 10 % | OVERALL | | |
| `MCStatNorm_WJets_tautau` | OVERALL | 50 % | W+jets SR yield = 514 events from very few simulated events with weights up to ~200; shape taken from τ2 ≥ VVVLoose |
| `FakeStat_tautau_DM*` | HISTO on fakes | 0.7–1.7 % | FF statistical uncertainty, per DM |
| `FakeOSSS_tautau` | HISTO on fakes | 3.0 % | C_OS/SS stat ⊕ 3 % extrapolation |
| `FakeClosure_tautau` | HISTO, shape | shape | same-sign non-closure in m_tt |
| MC statistics | γ per bin (`MCstatThreshold 0`) | | aMC@NLO signal (negative weights), AR data statistics |

Two pitfalls, both hit and fixed during this analysis (`CLAUDE.md` lists them):
envelope- and Hessian-type theory variations must be combined **on histograms**, never event by event.
The event-wise versions gave ±15 % (scales) and ±4 % (PDF) instead of ±4 % and ±1 %. And a sample with a
handful of large-weight events (W+jets) must not enter bin-by-bin: its "MC stat" of ±100 events in empty
bins inflated the γ parameters and the μ uncertainty.

## Not included (known limitations)

* e→τh and μ→τh energy scales (Z→ee/μμ is 2.6 % of the SR prediction and far from the Z→ττ peak in m_tt).
* Jet energy scale/resolution: jets only enter through N_jets (fake-factor binning, data-driven) and the
  unclustered MET. A JES effect on the MET of the simulation is not propagated.
* Jet→τh fake-rate scale factors for the τ2 leg of W+jets/tt̄ in simulation (covered by the cross-section
  and W+jets normalisation parameters, 1.2 % of the SR).
* Z pT reweighting of the aMC@NLO sample (not needed at NLO; covered by `QCDScale`/`SigModel`).
* τh ID SF dependence on the VSe/VSmu working points.
