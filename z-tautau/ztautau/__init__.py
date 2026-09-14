"""Z -> tau_h tau_h cross-section measurement with CMS 2016 Open Data (Run2016G+H, Tau dataset).

Modules (every physics choice is in `config`):
  config       selection, fake-factor binning, fit binning, paths
  samples      data / simulation registry, cross sections, file lists
  batch        subprocess-per-file runner (resumable, stall-safe)
  io           golden JSON, trigger OR, MET filters
  gen          LHE flavour split and the generator-level fiducial volume
  skim         NanoAOD -> skim (preselection, generator sums)
  objects      tau / lepton / jet selection, trigger matching, pair choice
  mass         visible, collinear and MET-likelihood di-tau masses
  corrections  TauPOG scale factors, tau energy scale, pileup
  fakes        fake-factor measurement and application
  analysis     ntuple loading, regions, event weights and systematic variations
  plotting     CMS-style stacked data / prediction plots
"""
