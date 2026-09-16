# MC event weights

`zmumu/weights.py`, `zmumu/pileup.py`.

```
w = genWeight x sigma L / sum genWeight x w_PU(nTrueInt) x L1PreFiringWeight_Nom
    x SF_ID(mu1) SF_ID(mu2) x SF_iso(mu1) SF_iso(mu2) x SF_trigger(event)
```

- **Normalisation** `sigma L / sum genWeight`: `L = 16393.381 pb^-1`; the sum runs over the
  `GenSums` of the processed skim files (equal to `Runs.genEventSumw`). aMC@NLO weights are
  +-constant (16% negative).
- **Pileup.** No official `puWeights` file is reachable from this cluster (the jsonpog
  GitLab needs CERN credentials). The data profile is built from the by-lumisection table of
  the CMS Open Data luminosity record (recid 1059, `pp_2016lumibyls.csv`, PHYSICS normtag):
  the certified lumisections of runs 278820-284044 weighted by recorded luminosity (their
  sum is exactly 16393.381 pb^-1). brilcalc's `avgpu` uses sigma_minbias = 80 mb; the CMS
  UL recommendation is 69.2 mb, so the profile is `avgpu x 69.2/80` (mean 24.7 for G+H,
  28.5 with 80 mb). That raw profile is ~4% too low and too narrow (the CSV has one average
  per lumisection and no bunch-to-bunch spread, which the official `pileupCalc.py` adds):
  with it, data/MC in N_PV rises from 0.8 to 1.2 across the distribution (REVIEW.md F6).
  `scripts/v2_2_pileup.py` therefore fits two parameters -- a scale of the per-lumisection
  mean and a relative Gaussian smearing per lumisection -- to the N_PV distribution of the
  Z -> mu mu signal region (`pileup.match_npv`); the **N_PV-matched profile** is the nominal
  (parameters and chi2 in `pileup_weights.json`, plots `pu_profile_data_mc.png`,
  `pu_npv_matching.png`). The MC profile is the `Pileup_nTrueInt` histogram of all generated
  events (`GenSums`, mean 21.9). Up/Down: sigma_minbias +- 4.6% on top of the matched scale.
  Weights are renormalised to unit mean over the MC profile and clipped at 10.
- **L1 pre-firing.** 2016 ECAL and muon-system prefiring; MC is weighted by
  `L1PreFiringWeight_Nom`, the systematic uses `_Up/_Dn`. Mean weight of DY -> mu mu events in
  the signal region: **0.9803** (muon system 0.9818, ECAL 0.9985) -- exactly the reviewer's
  0.980, i.e. the 2% effect that v1 could not correct. The Up/Down variation is the largest
  non-luminosity uncertainty of the result (0.51%).
- **Scale factors** from our own tag-and-probe (docs/12): per-muon ID and isolation, per-event
  trigger `[1-(1-e1^data)(1-e2^data)] / [1-(1-e1^MC)(1-e2^MC)]` with the measured per-muon
  efficiencies including the turn-on bins. Reconstruction: SF = 1 +- 0.4%/muon (assigned; the
  simulation truth is 0.9986. It *is* measurable on the unskimmed NanoAOD: 1.0001 +- 0.0013 per muon,
  docs/16, applied only in the tagged fit `zmumu_recosf`), taken as fully correlated between the two muons,
  i.e. the `MuonReco` nuisance parameter is **0.8% per event** (the first v2 fit applied 0.4% per
  event, REVIEW.md F8).
- **Theory variations on the signal** (PDF Hessian members, alpha_s, 7-point muR/muF, PS ISR/FSR)
  are renormalised so that the *fiducial* generator yield stays constant
  (`Weighter.theory_renorm`): a fiducial cross section only carries their effect on the
  C factor, not on the total cross section.
