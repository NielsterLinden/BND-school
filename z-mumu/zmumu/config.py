"""Central configuration: every physics choice in the measurement lives here.

Nothing in this file is derived at runtime, so it doubles as a written record
of the analysis definition. If you change a number here, rerun the whole chain
(`python run_all.py`) -- the steps cache nothing across runs.
"""

from pathlib import Path

# --------------------------------------------------------------------------
# Input data
# --------------------------------------------------------------------------
# Skimmed CMS 2016 SingleMuon Open Data (NanoAODv9, UL2016). The skim keeps
# every event with >= 2 entries in the Muon collection and ~300 branches.
# Provenance: /dcache/atlas/sjankovy/BND/collision_data/SingleMuon
SKIM_DIR = Path("/dcache/atlas/kdevries/BND2026/DoubleMuonSkimmed")
ERAS = ["Run2016G__30530", "Run2016H__30563"]

# Golden JSON (certified runs/lumisections). Events outside it are not usable
# for a cross-section measurement because they are not in the luminosity sum.
GRL_PATH = Path(__file__).resolve().parents[2] / "datasets" / "GRL" / "GRL.txt"

REPO_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_DIR / "output"
PLOT_DIR = OUTPUT_DIR / "plots"
DATA_DIR = OUTPUT_DIR / "data"

# --------------------------------------------------------------------------
# Luminosity
# --------------------------------------------------------------------------
# Run2016G + Run2016H, run range 278820-284044, certified lumisections of
# datasets/GRL/GRL.txt. The value is the normtag_PHYSICS recorded luminosity from the
# CMS Open Data luminosity record for 2016 (recid 1059, pp_2016lumibyls.csv):
# Run2016G 7653.261 + Run2016H 8740.119 pb^-1 (zmumu/pileup.py reproduces it exactly).
# The first iteration used 16290.713420 pb^-1 (brilcalc without --normtag), 0.63% low;
# the committed v1 outputs (output/RESULTS.md) were produced with that value.
# Uncertainty: 1.2% (CMS 2016 legacy luminosity calibration, CMS-LUM-17-003).
LUMI_PB = 16393.381
LUMI_REL_UNC = 0.012
LUMI_PB_NORMTAG = 16393.381        # kept for the review scripts
LUMI_PB_V1 = 16290.713420          # what steps 1-6 were run with

# --------------------------------------------------------------------------
# Trigger
# --------------------------------------------------------------------------
# The analysis trigger. Both paths exist in all 152 skim files (the rest of the
# 2016 muon menu does not -- see docs/01-selection.md).
TRIGGERS = ["HLT_IsoMu24", "HLT_IsoTkMu24"]

# Reference paths used to measure the trigger efficiency (docs/04-efficiency.md).
# These are single-muon paths that are *not* IsoMu24/IsoTkMu24 and are present
# in every file of both eras.
REFERENCE_TRIGGERS = ["HLT_Mu50", "HLT_IsoMu27", "HLT_Mu27", "HLT_Mu45_eta2p1"]

# --------------------------------------------------------------------------
# Muon selection (defines the fiducial volume)
# --------------------------------------------------------------------------
MU_ID_BRANCH = "Muon_mediumId"   # Muon POG medium working point
MU_PT_LEAD = 26.0                # GeV -- on the IsoMu24 efficiency plateau
MU_PT_SUBLEAD = 20.0             # GeV
MU_ETA_MAX = 2.4                 # muon system acceptance
MU_ISO_MAX = 0.15                # pfRelIso04_all, "tight" working point
MU_DXY_MAX = 0.2                 # cm, prompt-muon requirement
MU_DZ_MAX = 0.5                  # cm

# Exactly this many selected muons are required (suppresses WZ/ZZ and makes the
# fiducial volume unambiguous). See docs/01-selection.md.
N_MUONS_REQUIRED = 2

# --------------------------------------------------------------------------
# Event cleaning
# --------------------------------------------------------------------------
MET_FILTERS = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_eeBadScFilter",
]

# --------------------------------------------------------------------------
# FSR recovery (docs/02-fsr-photons.md)
# --------------------------------------------------------------------------
FSR_ENABLED = True
FSR_REL_ISO_MAX = 1.8      # FsrPhoton_relIso03
FSR_DR_OVER_ET2_MAX = 0.012  # FsrPhoton_dROverEt2, the standard CMS H->4l cut
FSR_PT_MIN = 2.0           # GeV

# --------------------------------------------------------------------------
# Signal region and binning
# --------------------------------------------------------------------------
MASS_LO, MASS_HI = 60.0, 120.0     # GeV, the fiducial mass window
MASS_NBINS = 120                   # 0.5 GeV bins

# Binning used for the tag-and-probe efficiency maps.
EFF_PT_BINS = [20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 120.0, 200.0]
EFF_ETA_BINS = [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, 0.0, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]

# --------------------------------------------------------------------------
# External inputs that cannot be measured from NanoAOD
# --------------------------------------------------------------------------
# Muon reconstruction efficiency (muon track -> reco Muon object). NanoAOD does
# not store generalTracks, so this cannot be measured here. Value and
# uncertainty from the CMS Muon POG for 2016 UL. See docs/04-efficiency.md.
MU_RECO_EFF = 0.9960
MU_RECO_EFF_UNC = 0.0040

# --------------------------------------------------------------------------
# Tag-and-probe definitions (docs/04-efficiency.md)
# --------------------------------------------------------------------------
TAG_PT_MIN = 26.0
TAG_ISO_MAX = 0.15
PROBE_PT_MIN = 20.0
TP_MASS_LO, TP_MASS_HI = 70.0, 110.0   # tighter window: purer tag-and-probe pairs

# --------------------------------------------------------------------------
# Background estimation (docs/03-fake-leptons.md, docs/05-backgrounds.md)
# --------------------------------------------------------------------------
# Transfer factor from the same-sign control region to the opposite-sign signal
# region for non-prompt ("fake") muons. 1.0 with a 50% uncertainty is the
# standard conservative choice when it cannot be measured in situ.
R_OS_SS = 1.0
R_OS_SS_UNC = 0.50

# Electron selection for the flavour-symmetric (e-mu) background estimate.
EL_PT_MIN = 20.0
EL_ETA_MAX = 2.5
EL_ID_MIN = 3          # Electron_cutBased >= 3 (medium)
EL_ISO_MAX = 0.15      # Electron_pfRelIso03_all

# Flavour-symmetric transfer factor: N(mumu) = FS_TRANSFER * k * N(emu).
# 0.5 is the combinatorial factor (ttbar/WW/Ztautau give e-mu twice as often as
# mu-mu). `k = eff_mu / eff_e` cannot be measured in the SingleMuon dataset --
# there is no Z->ee here -- so it is taken as 1.0 with a 50% uncertainty. The
# background is ~0.1% of the signal, so this dominates nothing.
FS_TRANSFER = 0.5
FS_EFF_RATIO = 1.0
FS_REL_UNC = 0.50

# --------------------------------------------------------------------------
# Systematic uncertainties that are not derived from a measurement
# --------------------------------------------------------------------------
# Bias of the reference-trigger method relative to true trigger-object
# matching (docs/04-efficiency.md).
TRIG_METHOD_REL_UNC = 0.005
# Muon momentum scale, propagated to the mass-window acceptance.
MUON_SCALE_REL_UNC = 0.002

# --------------------------------------------------------------------------
# Cross-checks on the parent NanoAOD (scripts/xcheck_efficiency.py)
# --------------------------------------------------------------------------
# The unskimmed SingleMuon NanoAOD. Same era sub-directories as the skim, and it
# still has the TrigObj_* branches the skim dropped.
PARENT_DIR = Path("/dcache/atlas/sjankovy/BND/collision_data/SingleMuon")

# HLT muon objects of the analysis paths. TrigObj_filterBits for muons:
# 2 = Iso (hltL3crIso*, IsoMu*), 8 = IsoTkMu. Objects above the path threshold
# are required; the matching is only trusted in events that fired the OR.
TRIG_OBJ_BITS = 2 | 8
TRIG_OBJ_PT_MIN = 24.0
TRIG_MATCH_DR = 0.1
# Systematic variation: drop objects whose L1 seed is in [20, 22) GeV. Those come
# from L1_SingleMu20-seeded paths sharing the Iso filter, not from IsoMu24.
TRIG_L1_VETO_BAND = (20.0, 22.0)

# Finer pT binning than EFF_PT_BINS around the 24 GeV turn-on.
TRIG_PT_BINS = [20.0, 22.0, 24.0, 25.0, 26.0, 27.0, 28.0, 30.0, 35.0, 40.0,
                50.0, 60.0, 80.0, 120.0, 200.0]

# --------------------------------------------------------------------------
# Simulation: acceptance of the fiducial volume (scripts/mc_acceptance.py)
# --------------------------------------------------------------------------
# Streamed from CERN Open Data EOS over xrootd; the file lists come from
# `cernopendata-client get-file-locations --recid <id> --protocol xrootd --verbose`.
FILELIST_DIR = REPO_DIR / "filelists"
DY_SAMPLES = {
    "nlo": {"recid": 35669, "filelist": "DYJetsToLL_M-50_NLO_amcatnloFXFX_35669.sizes.tsv",
            "label": "DYJetsToLL_M-50 amcatnloFXFX (NLO)"},
    "lo":  {"recid": 35671, "filelist": "DYJetsToLL_M-50_LO_madgraphMLM_35671.sizes.tsv",
            "label": "DYJetsToLL_M-50 madgraphMLM (LO)"},
}
# sigma(pp -> Z/gamma* -> ll, m_ll > 50 GeV) summed over e, mu, tau, NNLO. The
# same number the z-ee subgroup uses; the combination expects 6077.22/3 per flavour.
DY_XSEC_PB = 6077.22
# Generator-level fiducial muons: GenDressedLepton (muon + photons within
# dR < 0.1), |pdgId| == 13, no tau ancestor. The kinematic cuts and the mass
# window are the reco-level ones above.
GEN_MATCH_DR = 0.1

# --------------------------------------------------------------------------
# Processing
# --------------------------------------------------------------------------
N_WORKERS = 12
CHUNK_SIZE = "400 MB"
