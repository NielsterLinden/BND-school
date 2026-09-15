"""Central configuration of the Z -> tau_h tau_h measurement: every physics choice lives here.

Nothing in this file is computed at runtime, so it is also the written record of the analysis
definition. The physics motivation for each block is in docs/ (the section is quoted).
If you change a number here, rerun from the step that uses it (see README.md).
"""

from __future__ import annotations

import os
from pathlib import Path

# ------------------------------------------------------------------------------ paths
CHANNEL_DIR = Path(__file__).resolve().parents[1]          # z-tautau/
REPO_DIR = CHANNEL_DIR.parent                               # BND-school/
GRL_PATH = REPO_DIR / "datasets" / "GRL" / "GRL.txt"
FILELIST_DIR = CHANNEL_DIR / "filelists"
EXTERNAL_DIR = CHANNEL_DIR / "external"                     # small, committed POG inputs (JSON)
# Working point of both tau legs. v3 nominal: DeepTau VSjet Tight (v1/v2 used Medium: 80% fakes; Tight
# removes 63% of them for 24% of the signal, and gave the more precise result with the same chain,
# docs/08). BND_TAUTAU_WP=<other WP> runs the whole chain from step 3 with that working point and the
# matching TauPOG ID and trigger scale factors into variants/<wp>/{output,fit} and a separate BDT.
NOMINAL_WP = "Tight"
VERSION = "v3"                        # analysis version stamped on every plot (with the working point)
TAU_WP = os.environ.get("BND_TAUTAU_WP", NOMINAL_WP)
PLOT_TAG = f"{VERSION}: DeepTau {TAU_WP} " + r"$\tau_h$"   # drawn by plotting.label on every figure
_VARIANT = CHANNEL_DIR if TAU_WP == NOMINAL_WP else CHANNEL_DIR / "variants" / TAU_WP.lower()
OUTPUT_DIR = _VARIANT / "output"
PLOT_DIR = OUTPUT_DIR / "plots"
DATA_DIR = OUTPUT_DIR / "data"                              # git-ignored intermediates
FIT_DIR = _VARIANT / "fit"

# Bulk storage (not in git). Override with BND_TAUTAU_CACHE=/somewhere on a laptop.
CACHE_DIR = Path(os.environ.get("BND_TAUTAU_CACHE", "/data/atlas/users/sjankovy/BND-school-cache/ztautau"))
SKIM_DIR = CACHE_DIR / "skims_v1"          # NanoAOD-format skims, one file per parent file
NTUPLE_DIR = CACHE_DIR / "ntuples_v1"      # flat analysis ntuples, one file per sample (laptop bundle)
DOWNLOAD_DIR = CACHE_DIR / "downloads"     # raw POG ROOT files before conversion
BDT_DIR = CACHE_DIR / ("bdt" if TAU_WP == NOMINAL_WP else f"bdt_{TAU_WP.lower()}")   # k-fold BDT models (not committed)

# Luminosity-by-lumisection table of CMS Open Data record 1059 (PHYSICS normtag), used for the
# pileup profile. Same file the z-mumu v2 analysis uses.
LUMIBYLS_CSV = Path("/data/atlas/users/nterlind/BND-school-cache/external/pp_2016lumibyls.csv")
LUMIBYLS_URL = "root://eospublic.cern.ch//eos/opendata/cms/luminosity/2016/pp_2016lumibyls.csv"

# TRExFitter v1.8.0: the build of the git submodule. A worktree without its own build falls back to
# the main checkout, so every channel uses literally the same binary (fitting/CONVENTIONS.md).
TREX_FALLBACK_HOME = Path("/project/atlas/Users/nterlind/BND-school/TRExFitter-v1.8.0")

# ------------------------------------------------------------------------------ data set
# Tau primary dataset, Run2016G + Run2016H (UL2016 NanoAODv9), docs/01-data-and-samples.md.
RUN_MIN, RUN_MAX = 278820, 284044
ERA_H_FIRST_RUN = 281613
# Normtag luminosity of the certified lumisections (CMS Open Data record 1059), identical in all
# three channels (fitting/CONVENTIONS.md section 2).
LUMI_PB = 16393.381
LUMI_REL_UNC = 0.012

# ------------------------------------------------------------------------------ trigger
# The 2016 di-tau paths (docs/02-selection.md section "Trigger"). Run2016G has only the first,
# Run2016H only the second ("CombinedIso"); both exist in the UL16 simulation, where the OR is
# used, as in the TauPOG trigger scale-factor derivation.
DITAU_TRIGGERS = ["HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
                  "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg"]
# Trigger-object matching of both offline taus: TrigObj id 15, filterBits bit 1 (value 2) =
# "Medium(Comb)Iso" tau leg, above the path threshold, within dR < 0.5.
TRIGOBJ_TAU_ID = 15
TRIGOBJ_TAU_BITS = 2
TRIGOBJ_PT_MIN = 35.0
TRIG_MATCH_DR = 0.5

# ------------------------------------------------------------------------------ taus
# DeepTau2017v2p1 bitmasks in NanoAODv9: VSjet/VSe 1=VVVLoose 2=VVLoose 4=VLoose 8=Loose 16=Medium
# 32=Tight 64=VTight 128=VVTight; VSmu 1=VLoose 2=Loose 4=Medium 8=Tight.
TAU_DMS = (0, 1, 10, 11)                 # the decay modes with TauPOG scale factors
TAU_PT_MIN = 40.0                        # GeV, on the trigger plateau side of the 35 GeV threshold
TAU_ETA_MAX = 2.1                        # HLT eta2p1
TAU_DZ_MAX = 0.2                         # cm
TAU_VSE_BIT = 2                          # VVLoose  (TauPOG recommendation for tau_h tau_h)
TAU_VSMU_BIT = 1                         # VLoose
VSJET_BITS = {"VVVLoose": 1, "VVLoose": 2, "VLoose": 4, "Loose": 8, "Medium": 16, "Tight": 32, "VTight": 64, "VVTight": 128}
TAU_VSJET_TIGHT_BIT = VSJET_BITS[TAU_WP]  # the signal-region ("tight") working point: Tight nominally (v3)
TAU_VSJET_LOOSE_BIT = 1                  # VVVLoose: "loose" for the fake-factor regions
PAIR_DR_MIN = 0.5
# Candidates entering the pair choice (and the ntuples) are a little looser than the final cut so the
# tau energy-scale variations can move taus across TAU_PT_MIN without changing the chosen pair.
TAU_PT_NTUPLE = 38.0

# Skim preselection (docs/03-skims.md): at least two such candidates.
SKIM_TAU_PT = 35.0
SKIM_TAU_ETA = 2.3
SKIM_TAU_VSJET_BIT = 1
SKIM_TAU_VSE_BIT = 1
SKIM_TAU_VSMU_BIT = 1

# ------------------------------------------------------------------------------ other objects
# Extra-lepton vetoes keep the channel orthogonal to e tau_h, mu tau_h, ee and mumu.
VETO_MU = dict(pt=10.0, eta=2.4, dxy=0.045, dz=0.2, iso=0.3)        # + Muon_mediumId
VETO_EL = dict(pt=10.0, eta=2.5, dxy=0.045, dz=0.2, iso=0.3)        # + mvaFall17V2noIso_WP90, convVeto, lostHits <= 1
JET_PT, JET_ETA, JET_ID_BIT, JET_TAU_DR = 30.0, 4.7, 2, 0.5          # jetId bit 2 = tight
BJET_PT, BJET_ETA = 20.0, 2.4
BTAG_DEEPJET_MEDIUM = 0.2489                                        # UL2016 postVFP DeepJet medium WP

MET_FILTERS = ["Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
               "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
               "Flag_BadPFMuonFilter", "Flag_BadPFMuonDzFilter", "Flag_eeBadScFilter"]

# ------------------------------------------------------------------------------ di-tau mass
# docs/04-ditau-mass.md. The likelihood mass scans the visible energy fractions x1, x2 of the two
# taus on a grid, weighs each point with the MET transfer function (MET covariance from NanoAOD) and
# the two-body hadronic phase space (flat in x on [m_vis^2/m_tau^2, 1]) times a 1/m^2 prior, and returns the
# posterior median. On DY simulation (70 < m_LHE < 110 GeV, genuine tau_h tau_h in the SR): scale m/m_LHE
# 0.99 and core resolution 10.9% (IQR/2), vs 0.80 and 13.1% for m_vis; posterior mean without prior:
# 1.04 / 11.9%; collinear mass defined in only 38% of the events.
M_TAU = 1.77686
MASS_GRID_N = 60                 # grid points per x
MASS_CHUNK = 2000                # events per vectorised block
MASS_PRIOR_POW = 2.0             # prior m^-2
COV_MIN = 1.0                    # GeV^2 floor on the MET covariance diagonal (regularisation)

# ------------------------------------------------------------------------------ fake factors
# docs/05-fake-factors.md. FF(pT, DM) of the leading tau, measured in same-sign events whose
# subleading tau passes the tight working point; applied to opposite-sign events with the leading tau
# loose-not-tight. Events whose leading tau is a jet in simulation are dropped (the FF covers them).
# Binned in (DM, number of jets 0 / 1 / >= 2, tau1 pT): the first iteration without the jet binning
# over-predicted same-sign events with jets by 11% (1 jet) to 25% (>= 3 jets) and under-predicted 0-jet
# events by 9%, which moved the m_tt shape (docs/05-fake-factors.md, "Closure").
FF_PT_BINS = [40.0, 45.0, 50.0, 60.0, 80.0, 1000.0]
FF_NJET_BINS = [0, 1, 2]                 # lower edges; the last bin is inclusive
# ... and per era: Run2016G and Run2016H use different di-tau trigger paths (Iso / CombinedIso), and an
# era-inclusive FF over-predicted same-sign G by 4.6% and under-predicted H by 5.1%. Simulation (no era)
# uses the luminosity-weighted average. Recorded luminosity per era (normtag, record 1059):
FF_BY_ERA = True
ERA_LUMI_PB = (7653.261, 8740.119)       # (Run2016G, Run2016H), sum = LUMI_PB
FF_DMS = TAU_DMS
# Simulated events with a genuine leading tau are subtracted from the determination and application
# regions (the classic fake factor; without it the application region double counts ~6% of the signal
# and C_OS/SS is biased by +2.5%, REVIEW.md 3.2). The unsubtracted variant is only computed on request
# (step 4/5 --ff-variant nosub), it is not part of the nominal chain (v3).
FF_SUBTRACT_MC = True
# Factorised closure corrections of the FF, measured in same-sign data after the era x DM x N_jets x pT
# table: the FF varies by +-15% with |eta(tau1)| (universal shape: barrel-endcap transition, tracker
# edge, REVIEW.md section 5) and by -7% with pT(tau2) at 60-100 GeV (isolation correlation of the two
# jets). Applied multiplicatively in fakes.evaluate; each is obs/pred in same-sign events.
FF_CLOSURE_ETA_BINS = [0.0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1]
FF_CLOSURE_PT2_BINS = [40.0, 45.0, 50.0, 60.0, 80.0, 1000.0]
# The residual non-closure in the fit variable becomes one normalisation-type nuisance parameter per
# (BDT category, mass region below / above this split), docs/05-fake-factors.md.
FF_CLOSURE_MASS_SPLIT = 110.0
# Extra (non-statistical) relative uncertainty on the OS/SS extrapolation factors C_OS/SS (one per jet
# category), added in quadrature to their statistical uncertainty. C rises from 1.05 to 1.08-1.09 as
# the tau2 sideband is tightened towards Medium (after subtracting genuine taus), so the extrapolation
# to tau2 = Medium is uncertain at the few-percent level (docs/05-fake-factors.md, "OS/SS extrapolation").
FF_OSSS_SYST = 0.03

# ------------------------------------------------------------------------------ fiducial volume
# docs/06-cross-section.md. Generator level, in LHE Z/gamma* -> tau tau events:
# 60 < m(tau tau, LHE) < 120 GeV, both taus decay hadronically, both visible taus (GenVisTau)
# with pT > 40 GeV and |eta| < 2.1.
MASS_LO, MASS_HI = 60.0, 120.0
FID_VIS_PT, FID_VIS_ETA = 40.0, 2.1

# ------------------------------------------------------------------------------ BDT (docs/09-bdt.md)
# k-fold gradient-boosted classifier, Z -> tautau (fiducial, simulation) against the fake estimate
# (application-region data x FF). Mass-agnostic inputs only: the score defines categories and m_tt stays
# the fit variable in each. Nothing that enters the definition of the fake-factor regions (tau1 isolation)
# may be an input. Fold = event number mod BDT_K; every event is scored by the model that never saw it.
BDT_K = 5
BDT_FEATURES = ["t1_pt", "t2_pt", "pt_ratio", "t1_abseta", "t2_abseta", "dr_tt", "dphi_tt", "met", "met_sig",
                "pt_vis", "dphi_met_tt", "pt_tt", "njets", "jet1_pt", "t1_dm", "t2_dm"]
BDT_PARAMS = dict(n_estimators=300, max_depth=4, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                  min_child_weight=5, tree_method="hist", n_jobs=1, random_state=20260915)
BDT_TRAIN_ON_FIDUCIAL = True
# category edges in the score: [0, 0.55) fake dominated (fixes the fake normalisation and shape),
# [0.55, 0.90) mixed, [0.90, 1] signal dominated (S/B ~ 5)
BDT_CATEGORY_EDGES = [0.0, 0.55, 0.90, 1.0]

# Generator comparison for the fiducial C factor: madgraph LO (MLM) vs aMC@NLO FxFx, normalised to the
# same fiducial cross section. The LO sample's softer visible-tau pT spectrum inside the fiducial volume
# lowers C by ~15% (trigger turn-on), which is not a credible uncertainty of the NLO prediction: the LO
# sample is the worse model, the NLO scale / PS / PDF variations already move the spectrum, and the pT
# spectra are checked against data in the signal-dominated category (docs/07). The number is computed and
# reported (SigModel_tautau); set SIGMODEL_IN_FIT = True to include it as a one-sided normalisation NP.
SIGMODEL_IN_FIT = False

# ------------------------------------------------------------------------------ fit
FIT_VARIABLE = "m_tt"
FIT_BINS = [0.0, 40.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 150.0, 175.0, 200.0, 250.0, 350.0]
REGION_PREFIX = "tautau_SR"
REGIONS = [f"{REGION_PREFIX}{i}" for i in range(len(BDT_CATEGORY_EDGES) - 1)]     # one per BDT category
REGION_LABELS = ["BDT < 0.55 (fake dominated)", "0.55 < BDT < 0.90", "BDT > 0.90 (signal dominated)"]
# The fake-dominated category is the fake sideband: its bins below this m_tt are dropped from the fit
# (S/B = 0.09 there with a 15% closure prior; left in, its OS/SS parameter was 49% correlated with mu_Z).
SIDEBAND_REGION_MTT_MIN = 110.0
JOB = "ztautau"
DY_XSEC_PB = 6077.22             # sigma(Z/gamma* -> ll, m > 50) summed over flavours, NNLO (all channels)

# ------------------------------------------------------------------------------ processing
N_WORKERS = int(os.environ.get("BND_TAUTAU_WORKERS", "10"))
CHUNK_SIZE = "300 MB"
