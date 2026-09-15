# 07 — Corrections to the simulation and systematic uncertainties

## Corrections (`ztautau/corrections.py`, `analysis.weights`)

The official TauPOG inputs for UL2016 postVFP, DeepTau2017v2p1, are taken from the public GitHub
repositories `cms-tau-pog/TauIDSFs` and `cms-tau-pog/TauTriggerSFs` (branch `run2_SFs`). CMS cvmfs is not
mounted here and the CERN GitLab `jsonpog-integration` needs a login. `scripts/step0_external.py`
converts them to the committed `external/tau_pog_UL2016postVFP.json`.

| correction | applied to | value (UL2016 postVFP) |
|---|---|---|
| τh ID SF, VSjet Medium, **per decay mode** (TauPOG recommendation for τhτh, pT > 40) | genuine τh (`genPartFlav` 5) | DM0 0.923 ± 0.142, DM1 0.880 ± 0.054, DM10 0.868 ± 0.085, DM11 0.898 ± 0.163 |
| τh energy scale | genuine τh; propagated to MET, m_tt and the BDT inputs | DM0 0.993 ± 0.009, DM1 0.991 ± 0.007, DM10 1.001 ± 0.007, DM11 0.997 ± 0.016 |
| e→τh VSe VVLoose SF | prompt e / τ→e (flav 1, 3) | 1.06 ± 0.05 (barrel), 0.95 ± 0.06 (endcap) |
| μ→τh VSmu VLoose SF | prompt μ / τ→μ (flav 2, 4) | 0.88–1.04 per \|η\| bin |
| di-τ trigger leg SF(pT, DM) | every **genuine** leg (SF = 1 for a jet leg, `genPartFlav` 0, since v2) | ratio of the fitted data and MC efficiencies, e.g. 0.94–0.97 at 50 GeV; data efficiency 0.3–0.6 at 40 GeV (turn-on) |
| pileup | all MC | data profile from the Open Data lumi-by-LS table (record 1059), σ_mb = 69.2 mb; same construction as z-mumu v2 |
| L1 prefiring | all MC | `L1PreFiringWeight_Nom` (mean 0.990 for the signal) |

**The two TauPOG ID prescriptions disagree.** The pT-binned SF for pT > 40 GeV is 0.823 +0.086 −0.036 for
every decay mode, 7 % per leg below the DM-binned values used here, i.e. ~14 % on the τhτh yield (μ_Z would
be ~14 % higher with it). The DM-binned set is the POG's recommendation for selections with a
decay-mode-dependent trigger like the di-τ trigger, and it is what is used; the pT-binned result is a
cross-check to be quoted, not a nuisance parameter (`REVIEW.md` 3.3). In the combination with ee/μμ the
`TauID_DM*` parameters are measured in situ, which settles the question.

## Systematic uncertainties (nuisance parameters in the fit)

Names follow `fitting/CONVENTIONS.md`: identical strings are correlated between channels; `_tautau`
marks channel-specific parameters. "Prefit effect" = change of the fiducial Z→ττ (or fake) yield summed
over the three categories. The complete list with per-category sizes is in
`output/data/yields.json["prefit_norm_effects"]`.

| NP | type | prefit effect | notes |
|---|---|---|---|
| `Lumi` | OVERALL, all MC | 1.2 % | correlated with ee, μμ |
| `Pileup` | HISTO | ~1 % | σ_mb ± 4.6 %; correlated |
| `L1Prefiring` | HISTO | 0.2 % | correlated |
| `TauID_DM0/1/10/11` | HISTO | 9 / 6 / 4 / 3 % | TauPOG uncertainty, uncorrelated between DMs, correlated between the two legs |
| `TauTrigger_DM*` | HISTO | 1–3 % | fitted SF uncertainty vs pT |
| `TauES_DM*` | HISTO, shape (smoothed) | 1–2 % + shape | τ pT, MET, m_tt and the BDT score recomputed (events migrate between categories) |
| `TauFakeEle`, `TauFakeMu` | HISTO | < 1.3 % on Z→ee | |
| `MET_Unclustered` | HISTO, shape | < 0.1 % norm. | NanoAOD unclustered-energy shift, m_tt and the BDT score recomputed |
| `QCDScale` | HISTO on `DYtautau` and `DYtautau_nonfid`, fiducial-normalised | few % | 6-point envelope per bin; on the non-fiducial part it varies the non-fid / fid ratio |
| `PDF` | HISTO, fiducial-normalised | ~1 % | 100 Hessian eigenvectors in quadrature per bin, + α_s |
| `PS_ISR`, `PS_FSR` | HISTO, fiducial-normalised | < 1 % | PSWeight 2 / 0.5 |
| `XS_DYtautau_nonfid` 5 %, `XS_DYll` 5 %, `XS_DYlowmass` 10 %, `XS_TTbar` 6 %, `XS_SingleTop`/`XS_WW`/`XS_WZ`/`XS_ZZ`/`XS_WJets` 10 % | OVERALL | | |
| `MCStatNorm_WJets_tautau` | OVERALL | 50 % | W+jets SR yield from very few simulated events with weights up to ~200; shape taken from τ2 ≥ VVVLoose per category |
| `FakeOSSS_tautau` | HISTO on fakes | 3.0 % | C_OS/SS stat ⊕ 3 % extrapolation |
| `FakeClosure_tautau_c<k>_lo`, `_hi` | HISTO on fakes, per category k and m_tt region | c0/lo 15 %, c0/hi 1 %, c1/lo 5 %, c1/hi 7 %, c2/lo 12 %, c2/hi 55 % | residual same-sign non-closure ⊕ its statistical error (`05-fake-factors.md`) |
| MC statistics | γ per bin (`MCstatThreshold 0`) | | signal and background templates; the `Fakes` template carries the AR data statistics **and** the FF statistics per event |

**Not a nuisance parameter (reported only): `SigModel_tautau`.** The fiducial C factor with madgraph LO
instead of aMC@NLO, both normalised to the NLO fiducial cross section, is 0.867 (v1 used the LO/NLO
difference of the whole template, 7.3 %, as a one-sided NP). Reasoning (`REVIEW.md` 3.4): the LO sample is
the worse model of the Z pT and hence of the visible-τ pT spectrum, whose slope at the 40 GeV threshold
and on the trigger turn-on is what the number measures; an uncertainty built as "nominal minus a worse
model" is a placeholder, it double counts `QCDScale`/`PS_ISR`, and its former name `SigModel` would have
been correlated in the combination with z-mumu's powheg-vs-aMC@NLO shape NP (0.2 %), a different quantity.
The NLO prescription (scales, PS, PDF) is used, and the visible-pT spectra are checked against data in the
signal-dominated category (`step4_tautau_SR2_t1_pt.png`, `step4_tautau_SR2_t2_pt.png`). Set
`config.SIGMODEL_IN_FIT = True` to include it as a one-sided normalisation NP.

Pitfalls, all hit and fixed during this analysis (`CLAUDE.md` lists them): theory envelopes / Hessians must
be combined **on histograms**, never event by event; a sample with a handful of large-weight events
(W+jets) must not enter bin-by-bin nor be subtracted event by event from the fake-factor regions.

## Not included (known limitations)

* e→τh and μ→τh energy scales (Z→ee/μμ is 0.5 % of the SR prediction and far from the Z→ττ peak in m_tt).
* Jet energy scale/resolution: jets only enter through N_jets (fake-factor binning, data-driven), the
  unclustered MET and two BDT inputs (N_jets, leading-jet pT). A JES effect on the MET of the simulation
  is not propagated.
* Jet→τh fake-rate scale factors for the τ2 leg of W+jets/tt̄ in simulation (covered by the cross-section
  and W+jets normalisation parameters, ~1 % of the SR).
* Z pT reweighting of the aMC@NLO sample (not needed at NLO; covered by `QCDScale`).
* τh ID SF dependence on the VSe/VSmu working points; the DM- vs pT-binned prescription difference (above).
* Trigger efficiency in situ: the POG SFs sit on the turn-on for both legs (data efficiency 0.3–0.6 at
  40 GeV) and their uncertainty does not shrink with a higher threshold (`REVIEW.md` 3.7).
