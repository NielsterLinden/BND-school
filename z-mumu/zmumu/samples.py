"""Sample registry for the v2 (MC-based) measurement: datasets, cross sections, file lists.

Every sample used anywhere in v2 is defined here once. `xsec_pb` normalises the MC to
`LUMI_PB` (`config.LUMI_PB_NORMTAG`); the measured sigma_fid does not depend on the DY
value (the signal strength floats in the fit and the C factor is a ratio), and the
backgrounds are < 0.5% of the selected sample, so the 5-10% uncertainties assigned to
the background cross sections in the fit are more than enough.

File lists come from the CERN Open Data API (`https://opendata.cern.ch/api/records/<recid>`),
which embeds every file's xrootd URI, size and adler32. They are cached as
`filelists/<key>_<recid>.sizes.tsv` (url, size, checksum; readable by `batch.read_catalogue`).
Samples that were mirrored to dCache are read locally with the EOS URL as fallback.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

from . import config

API = "https://opendata.cern.ch/api/records/{recid}?format=json"
DCACHE_MC = Path("/dcache/atlas/sjankovy/BND/mc")
DCACHE_DATA = Path("/dcache/atlas/sjankovy/BND/collision_data/SingleMuon")

# Cross-section references (all in pb, 13 TeV):
#   [SMP]  CMS StandardModelCrossSectionsat13TeV TWiki (NNLO where available)
#   [TOP]  CMS TopSystematics/TtbarNNLO TWiki: sigma(ttbar) = 831.76 pb (NNLO+NNLL, m_t = 172.5)
#   [XSDB] CMS cross-section database (GenXsecAnalyzer, NLO) for the exclusive diboson samples
SAMPLES = {
    # ---- data -------------------------------------------------------------------
    "data_2016G": dict(recid=30530, is_mc=False, era="G", group="data",
                       dataset="/SingleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                       dcache=DCACHE_DATA / "Run2016G__30530", xsec_pb=None, fit_sample="Data"),
    "data_2016H": dict(recid=30563, is_mc=False, era="H", group="data",
                       dataset="/SingleMuon/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD",
                       dcache=DCACHE_DATA / "Run2016H__30563", xsec_pb=None, fit_sample="Data"),
    # ---- Drell-Yan (signal + Z->tautau, Z->ee via the LHE flavour split) ------
    "DY_NLO": dict(recid=35669, is_mc=True, group="DY", split_lhe=True, keep_pdf=True,
                   dataset="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                   dcache=DCACHE_MC / "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__35669",
                   xsec_pb=6077.22, xsec_ref="[SMP] NNLO FEWZ, m_ll > 50 GeV, sum over e/mu/tau; same as z-ee",
                   fit_sample=None),   # -> DYmumu / DYee / DYtautau
    "DY_powheg": dict(recid=75482, is_mc=True, group="DY_alt", split_lhe=False, keep_pdf=True,
                      dataset="ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8",
                      dcache=DCACHE_MC / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8__75482",
                      xsec_pb=None, xsec_ref="alternative generator for the C factor only (normalised to the NLO fiducial prediction)",
                      fit_sample="DYmumu_powheg"),
    # ---- top ----------------------------------------------------------------------
    "TTbar": dict(recid=67801, is_mc=True, group="top", dcache=None,
                  dataset="TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
                  xsec_pb=88.29, xsec_ref="[TOP] 831.76 pb x BR(WW -> l l nu nu, incl. tau) = 0.1062",
                  fit_sample="TTbar"),
    "ST_tW_top": dict(recid=64895, is_mc=True, group="top", dcache=None,
                      dataset="ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
                      xsec_pb=19.47, xsec_ref="[XSDB] 35.85 pb (NLO+NNLL) x (1 - 0.6741^2)",
                      fit_sample="SingleTop"),
    "ST_tW_antitop": dict(recid=64839, is_mc=True, group="top", dcache=None,
                          dataset="ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
                          xsec_pb=19.47, xsec_ref="[XSDB] 35.85 pb (NLO+NNLL) x (1 - 0.6741^2)",
                          fit_sample="SingleTop"),
    # ---- diboson ------------------------------------------------------------------
    "WW": dict(recid=72676, is_mc=True, group="diboson",
               dataset="WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
               dcache=DCACHE_MC / "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8__72676",
               xsec_pb=12.178, xsec_ref="[SMP] 118.7 pb (NNLO) x BR(2l2nu)", fit_sample="WW"),
    "WZ_3LNu": dict(recid=72752, is_mc=True, group="diboson",
                    dataset="WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                    dcache=DCACHE_MC / "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8__72752",
                    xsec_pb=4.43, xsec_ref="[XSDB] NLO", fit_sample="WZ"),
    "WZ_2Q2L": dict(recid=72742, is_mc=True, group="diboson",
                    dataset="WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                    dcache=DCACHE_MC / "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8__72742",
                    xsec_pb=5.60, xsec_ref="[XSDB] NLO", fit_sample="WZ"),
    "ZZ_4L": dict(recid=75589, is_mc=True, group="diboson", max_files=25,
                  dataset="ZZTo4L_TuneCP5_13TeV_powheg_pythia8",
                  dcache=DCACHE_MC / "ZZTo4L_TuneCP5_13TeV_powheg_pythia8__75589",
                  xsec_pb=1.256, xsec_ref="[XSDB] powheg NLO", fit_sample="ZZ"),
    "ZZ_2L2Nu": dict(recid=75567, is_mc=True, group="diboson",
                     dataset="ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8",
                     dcache=DCACHE_MC / "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8__75567",
                     xsec_pb=0.564, xsec_ref="[XSDB] powheg NLO", fit_sample="ZZ"),
    "ZZ_2Q2L": dict(recid=75573, is_mc=True, group="diboson",
                    dataset="ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                    dcache=DCACHE_MC / "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8__75573",
                    xsec_pb=3.22, xsec_ref="[XSDB] NLO", fit_sample="ZZ"),
    # ---- W+jets: prompt subtraction in the fake-factor regions only, never in the SR ------
    "WJets": dict(recid=69745, is_mc=True, group="wjets",
                  dataset="WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                  dcache=DCACHE_MC / "WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8__69745",
                  xsec_pb=61526.7, xsec_ref="[SMP] NNLO, 3 x 20508.9 pb", fit_sample="WJets"),
}

# Defaults for keys not given above.
for _k, _s in SAMPLES.items():
    _s.setdefault("split_lhe", False)
    _s.setdefault("keep_pdf", False)
    _s.setdefault("max_files", None)
    _s.setdefault("dcache", None)
    _s.setdefault("era", None)

# Relative normalisation uncertainties used as OVERALL nuisance parameters in the fit.
XSEC_UNC = {"TTbar": 0.06, "SingleTop": 0.10, "WW": 0.10, "WZ": 0.10, "ZZ": 0.10, "DYtautau": 0.05,
            "WJets": 0.30}     # WJets: non-prompt-electron rate from MC, e-mu validation region only

DATA_KEYS = [k for k, s in SAMPLES.items() if not s["is_mc"]]
MC_KEYS = [k for k, s in SAMPLES.items() if s["is_mc"]]


# --------------------------------------------------------------------------
# File lists
# --------------------------------------------------------------------------
def filelist_path(key: str) -> Path:
    s = SAMPLES[key]
    return config.FILELIST_DIR / f"{key}_{s['recid']}.sizes.tsv"


def fetch_filelist(recid: int, timeout: float = 60.0) -> list[tuple[str, int, str]]:
    """[(xrootd URI, size, checksum)] for a record, from the Open Data API."""
    with urllib.request.urlopen(API.format(recid=recid), timeout=timeout) as resp:
        meta = json.load(resp)["metadata"]
    files = []
    for index in meta.get("_file_indices", []):
        for f in index.get("files", []):
            files.append((f["uri"], int(f["size"]), f["checksum"]))
    if not files:
        raise RuntimeError(f"record {recid}: no files in the API response")
    expected = meta.get("distribution", {}).get("number_files")
    if expected and expected != len(files):
        raise RuntimeError(f"record {recid}: API lists {len(files)} files, catalogue says {expected}")
    return sorted(files)


def write_filelist(key: str) -> Path:
    path = filelist_path(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = fetch_filelist(SAMPLES[key]["recid"])
    path.write_text("".join(f"{u}\t{n}\t{c}\n" for u, n, c in rows))
    return path


def read_filelist(key: str) -> list[tuple[str, int, str]]:
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
    """Per input file: dict(name, eos, local, size, checksum, primary, fallback).

    `primary` is the path to read first (dCache copy when it exists and `prefer` is
    'dcache'), `fallback` the other one (or None).
    """
    s = SAMPLES[key]
    out = []
    for url, size, checksum in read_filelist(key):
        name = url.rsplit("/", 1)[1]
        local = (s["dcache"] / name) if s["dcache"] else None
        if local is not None and not local.exists():
            local = None
        if prefer == "dcache" and local is not None:
            primary, fallback = str(local), url
        else:
            primary, fallback = url, (str(local) if local is not None else None)
        out.append(dict(name=name, eos=url, local=str(local) if local else None, size=size,
                        checksum=checksum, primary=primary, fallback=fallback))
    limit = max_files or s["max_files"]
    return out[:limit] if limit else out
