"""Sample registry: every data set and simulation sample used by the Z -> tau_h tau_h analysis.

Each entry: CERN Open Data record id, whether a dCache mirror exists, cross section (pb) with its
reference, how the sample enters the fit (`fit_sample`; `None` = split by LHE flavour into
DYtautau / DYee / DYmumu), and an optional `max_files` (large samples whose contribution is tiny are
only partly processed; the normalisation uses the generator-weight sum of the processed files, so a
subset is unbiased, only statistically poorer).

QCD multijet simulation is deliberately absent: all jet -> tau_h fakes come from the fake-factor
method in data (docs/05-fake-factors.md).

File lists: the dCache manifest `/dcache/atlas/sjankovy/BND/_manifest/manifest.json` already has
(xrootd URI, size, adler32) for every mirrored record; the others come from the Open Data API.
Both are cached as filelists/<key>_<recid>.tsv.
"""

from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path

from . import config

API = "https://opendata.cern.ch/api/records/{recid}?format=json"
MANIFEST = Path("/dcache/atlas/sjankovy/BND/_manifest/manifest.json")
DCACHE_MC = Path("/dcache/atlas/sjankovy/BND/mc")
DCACHE_DATA = Path("/dcache/atlas/sjankovy/BND/collision_data")
DCACHE_TAU = DCACHE_DATA / "Tau"

# Cross-section references (pb, 13 TeV):
#   [SMP]  CMS StandardModelCrossSectionsat13TeV TWiki (NNLO where available)
#   [TOP]  CMS TtbarNNLO TWiki: sigma(ttbar) = 831.76 pb; BR(W -> l nu, l = e/mu/tau) = 0.3259
#   [XSDB] CMS GenXSecAnalyzer database / values used by CMS H->tautau 2016 for inclusive dibosons
SAMPLES = {
    # ---- data ---------------------------------------------------------------------------------
    "data_2016G": dict(recid=30532, is_mc=False, group="data", fit_sample="Data",
                       dataset="/Tau/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                       dcache=DCACHE_TAU / "Run2016G__30532"),
    "data_2016H": dict(recid=30565, is_mc=False, group="data", fit_sample="Data",
                       dataset="/Tau/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                       dcache=DCACHE_TAU / "Run2016H__30565"),
    # ---- v4 lepton channels: SingleMuon (mu tau_h), SingleElectron (e tau_h, EOS only), MuonEG (e mu) ----
    "data_mu_2016G": dict(recid=30530, is_mc=False, group="data_mu", fit_sample="Data", stream="SingleMuon", v4=True,
                          dataset="/SingleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                          dcache=DCACHE_DATA / "SingleMuon/Run2016G__30530"),
    "data_mu_2016H": dict(recid=30563, is_mc=False, group="data_mu", fit_sample="Data", stream="SingleMuon", v4=True,
                          dataset="/SingleMuon/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                          dcache=DCACHE_DATA / "SingleMuon/Run2016H__30563"),
    "data_el_2016G": dict(recid=30529, is_mc=False, group="data_el", fit_sample="Data", stream="SingleElectron", v4=True,
                          dataset="/SingleElectron/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD", dcache=None),
    "data_el_2016H": dict(recid=30562, is_mc=False, group="data_el", fit_sample="Data", stream="SingleElectron", v4=True,
                          dataset="/SingleElectron/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD", dcache=None),
    "data_emu_2016G": dict(recid=30528, is_mc=False, group="data_emu", fit_sample="Data", stream="MuonEG", v4=True,
                           dataset="/MuonEG/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                           dcache=DCACHE_DATA / "MuonEG/Run2016G__30528"),
    "data_emu_2016H": dict(recid=30561, is_mc=False, group="data_emu", fit_sample="Data", stream="MuonEG", v4=True,
                           dataset="/MuonEG/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                           dcache=DCACHE_DATA / "MuonEG/Run2016H__30561"),
    # ---- Drell-Yan: signal (LHE tau tau) and Z -> ee / mumu with lepton -> tau_h fakes ----------
    "DY_NLO": dict(recid=35669, is_mc=True, group="DY", split_lhe=True, keep_pdf=True,
                   dataset="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                   dcache=DCACHE_MC / "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__35669",
                   xsec_pb=6077.22, xsec_ref="[SMP] FEWZ NNLO, m_ll > 50 GeV, sum of e/mu/tau (all channels)",
                   fit_sample=None),
    # jet-binned aMC@NLO FxFx samples (0, 1, 2 partons in the NLO matrix element, LHE_NpNLO exclusive):
    # extra signal statistics, stitched with the inclusive sample per LHE_NpNLO bin (analysis.dy_norm);
    # the jet-bin fractions come from the inclusive sample, so no separate cross sections are needed.
    "DY_0J": dict(recid=35577, is_mc=True, group="DY", split_lhe=True, keep_pdf=True, optional=True,
                  dataset="DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                  dcache=DCACHE_MC / "DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__35577",
                  xsec_pb=6077.22, xsec_ref="stitched to the inclusive sample per LHE_NpNLO bin", fit_sample=None),
    "DY_1J": dict(recid=35595, is_mc=True, group="DY", split_lhe=True, keep_pdf=True, optional=True,
                  dataset="DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                  dcache=DCACHE_MC / "DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8__35595",
                  xsec_pb=6077.22, xsec_ref="stitched to the inclusive sample per LHE_NpNLO bin", fit_sample=None),
    "DY_2J": dict(recid=35613, is_mc=True, group="DY", split_lhe=True, keep_pdf=True, optional=True,
                  dataset="DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                  dcache=DCACHE_MC / "DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8__35613",
                  xsec_pb=6077.22, xsec_ref="stitched to the inclusive sample per LHE_NpNLO bin", fit_sample=None),
    "DY_LO": dict(recid=35671, is_mc=True, group="DY_alt", split_lhe=True, keep_pdf=False, max_files=24,
                  dataset="DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
                  dcache=DCACHE_MC / "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__35671",
                  xsec_pb=6077.22, xsec_ref="normalised to the same NNLO; alternative generator for C only",
                  fit_sample=None),
    "DY_lowmass": dict(recid=35631, is_mc=True, group="DY_lowmass", fit_sample="DYlowmass",
                       dataset="DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                       dcache=DCACHE_MC / "DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8__35631",
                       xsec_pb=18610.0, xsec_ref="[XSDB] aMC@NLO, 10 < m_ll < 50 GeV, sum of flavours"),
    # ---- W + jets: only events whose leading tau is genuine (the rest is in the fake factor) ------
    "WJets": dict(recid=69745, is_mc=True, group="wjets", fit_sample="WJets",
                  dataset="WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                  dcache=DCACHE_MC / "WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8__69745",
                  xsec_pb=61526.7, xsec_ref="[SMP] NNLO 3 x 20508.9"),
    # W+jets for the lepton channels (v4): the madgraph MLM sample has unit weights and 81M events, which the
    # fake-factor determination regions (m_T > 70 GeV) and the AR fractions need; the tau_h tau_h channel keeps
    # the aMC@NLO sample above (v1-v3 chain unchanged).
    "WJets_LO": dict(recid=69747, is_mc=True, group="wjets_lo", fit_sample="WJets", v4_only=True,
                     dataset="WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
                     dcache=DCACHE_MC / "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8__69747",
                     xsec_pb=61526.7, xsec_ref="[SMP] NNLO 3 x 20508.9"),
    # ---- top (EOS only) ----------------------------------------------------------------------------
    "TTTo2L2Nu": dict(recid=67801, is_mc=True, group="top", fit_sample="TTbar", max_files=20,
                      dataset="TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8", dcache=None,
                      xsec_pb=88.29, xsec_ref="[TOP] 831.76 x 0.3259^2 x ... = 88.29"),
    "TTToSemiLeptonic": dict(recid=67993, is_mc=True, group="top", fit_sample="TTbar", max_files=14,
                             dataset="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8", dcache=None,
                             xsec_pb=365.34, xsec_ref="[TOP] 831.76 x 2 x 0.3259 x 0.6741"),
    "ST_tW_top": dict(recid=64895, is_mc=True, group="top", fit_sample="SingleTop",
                      dataset="ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8", dcache=None,
                      xsec_pb=19.47, xsec_ref="[XSDB] 35.85 x (1 - 0.676^2), as z-mumu"),
    "ST_tW_antitop": dict(recid=64839, is_mc=True, group="top", fit_sample="SingleTop",
                          dataset="ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8", dcache=None,
                          xsec_pb=19.47, xsec_ref="[XSDB] 35.85 x (1 - 0.676^2), as z-mumu"),
    # ---- diboson: inclusive pythia8 samples (all decays, so tau_h from any boson is included) --------
    "WW": dict(recid=72696, is_mc=True, group="diboson", fit_sample="WW",
               dataset="WW_TuneCP5_13TeV-pythia8",
               dcache=DCACHE_MC / "WW_TuneCP5_13TeV-pythia8__72696",
               xsec_pb=118.7, xsec_ref="[SMP] NNLO"),
    "WZ": dict(recid=72754, is_mc=True, group="diboson", fit_sample="WZ",
               dataset="WZ_TuneCP5_13TeV-pythia8",
               dcache=DCACHE_MC / "WZ_TuneCP5_13TeV-pythia8__72754",
               xsec_pb=47.13, xsec_ref="[XSDB] MCFM NLO"),
    "ZZ": dict(recid=75593, is_mc=True, group="diboson", fit_sample="ZZ",
               dataset="ZZ_TuneCP5_13TeV-pythia8",
               dcache=DCACHE_MC / "ZZ_TuneCP5_13TeV-pythia8__75593",
               xsec_pb=16.523, xsec_ref="[XSDB] MCFM NLO"),
}
for _k, _s in SAMPLES.items():
    _s.setdefault("stream", "Tau" if not _s["is_mc"] else None)
    _s.setdefault("v4", True)            # part of the v4 (four-channel) skim / ntuple production
    _s.setdefault("v4_only", False)      # not part of the v1-v3 tau_h tau_h chain
    _s.setdefault("split_lhe", False)
    _s.setdefault("keep_pdf", False)
    _s.setdefault("max_files", None)
    _s.setdefault("xsec_pb", None)
    _s.setdefault("optional", False)
    _s["key"] = _k

# the tau_h tau_h chain (v1-v3): Tau stream data and the simulation registered before v4
DATA_KEYS = [k for k, s in SAMPLES.items() if not s["is_mc"] and s["stream"] == "Tau"]
MC_KEYS = [k for k, s in SAMPLES.items() if s["is_mc"] and not s["v4_only"]]
NOMINAL_MC_KEYS = [k for k in MC_KEYS if SAMPLES[k]["group"] != "DY_alt"]
# v4: the three lepton-channel data streams and every simulation sample (skims_v4)
DATA_KEYS_V4 = {stream: [k for k, s in SAMPLES.items() if not s["is_mc"] and s["stream"] == stream]
                for stream in ("SingleMuon", "SingleElectron", "MuonEG")}
V4_SKIM_KEYS = [k for k, s in SAMPLES.items() if s["v4"] and (s["is_mc"] or s["stream"] != "Tau")]
MC_KEYS_V4 = [k for k, s in SAMPLES.items() if s["is_mc"] and s["v4"]]
# W+jets: NLO sample in tau_h tau_h, madgraph in the lepton channels
MC_KEYS_LEPTON = [k for k in MC_KEYS_V4 if k != "WJets" and SAMPLES[k]["group"] != "DY_alt"]
# the Drell-Yan signal group: inclusive sample (normalisation, acceptance) + jet-binned samples (statistics)
DY_INCLUSIVE = "DY_NLO"
DY_STITCHED = ["DY_NLO", "DY_0J", "DY_1J", "DY_2J"]

# Relative normalisation uncertainties of the backgrounds (OVERALL nuisance parameters; the names
# follow fitting/CONVENTIONS.md, XS_DYll and XS_WJets are specific to this channel).
XSEC_UNC = {"DYlowmass": 0.10, "TTbar": 0.06, "SingleTop": 0.10, "WW": 0.10, "WZ": 0.10, "ZZ": 0.10, "WJets": 0.10,
            "DYee": 0.05, "DYmumu": 0.05, "DYtautau_nonfid": 0.05}
# Z/gamma* -> tautau outside the fiducial volume (mostly m_LHE > 120 GeV, docs/06-cross-section.md) is a
# background normalised to theory, not part of the signal strength.
SIGNAL = "DYtautau"
SIGNAL_NONFID = "DYtautau_nonfid"
FIT_BACKGROUNDS = [SIGNAL_NONFID, "DYee", "DYmumu", "DYlowmass", "WJets", "TTbar", "SingleTop", "WW", "WZ", "ZZ"]
LHE_SPLIT = {15: "DYtautau", 11: "DYee", 13: "DYmumu"}


# ------------------------------------------------------------------------------ file lists
def filelist_path(key: str) -> Path:
    return config.FILELIST_DIR / f"{key}_{SAMPLES[key]['recid']}.tsv"


def _from_manifest(recid: int):
    if not MANIFEST.exists():
        return None
    rec = json.load(open(MANIFEST))["recids"].get(str(recid))
    if not rec:
        return None
    return sorted((f["uri"], int(f["size"]), f["checksum"]) for f in rec["files"])


def _from_api(recid: int, timeout: float = 60.0):
    with urllib.request.urlopen(API.format(recid=recid), timeout=timeout) as resp:
        meta = json.load(resp)["metadata"]
    rows = [(f["uri"], int(f["size"]), f["checksum"])
            for idx in meta.get("_file_indices", []) for f in idx.get("files", [])]
    if not rows:
        raise RuntimeError(f"record {recid}: no files in the Open Data API response")
    return sorted(rows)


def write_filelist(key: str) -> Path:
    recid = SAMPLES[key]["recid"]
    rows = _from_manifest(recid) or _from_api(recid)
    path = filelist_path(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{u}\t{n}\t{c}\n" for u, n, c in rows))
    return path


def read_filelist(key: str):
    path = filelist_path(key)
    if not path.exists():
        write_filelist(key)
    rows = []
    for line in path.read_text().splitlines():
        if line.strip():
            u, n, c = line.split("\t")
            rows.append((u, int(n), c))
    return rows


def sources(key: str, prefer: str = "dcache", max_files: int | None = None):
    """Per parent file: dict(name, stem, eos, local, size, primary, fallback).

    The subset for `max_files` is the first N of the sorted list, so it is reproducible.
    """
    s = SAMPLES[key]
    out = []
    # one directory listing instead of a stat per file: stats on dCache NFS can hang for minutes
    mirrored = set(os.listdir(s["dcache"])) if s["dcache"] and s["dcache"].is_dir() else set()
    for url, size, checksum in read_filelist(key):
        name = url.rsplit("/", 1)[1]
        local = s["dcache"] / name if name in mirrored else None
        if prefer == "dcache" and local is not None:
            primary, fallback = str(local), url
        else:
            primary, fallback = url, (str(local) if local else None)
        out.append(dict(name=name, stem=Path(name).stem, eos=url, local=str(local) if local else None,
                        size=size, checksum=checksum, primary=primary, fallback=fallback))
    limit = max_files or s["max_files"]
    return out[:limit] if limit else out
