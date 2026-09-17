"""Official CMS POG corrections for the lepton channels, read with correctionlib (docs/10-v4-plan.md section 7).

The files are the public CMS Open Data copy of jsonpog-integration for UL2016 postVFP
(`datasets/corrections/POG/`, unpacked from opendata.cern.ch/eos/opendata/cms/corrections/
jsonpog-integration-2016post.tar by the z-ee group). Everything here is a thin, vectorised wrapper:
numpy arrays in, numpy arrays out; `syst` selects the variation. The tau_h VSjet ID scale factor is *not*
here on purpose: in v4 it is a free parameter of the fit (docs/10, section 8). The di-tau trigger and the tau
energy scale keep the v3 inputs (`corrections.py`, from the TauPOG GitHub repositories) so the tau_h tau_h
channel is unchanged.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from . import config

POG_DIR = config.REPO_DIR / "datasets" / "corrections" / "POG"
FILES = {
    "muon": POG_DIR / "MUO/2016postVFP_UL/muon_Z.json.gz",
    "electron": POG_DIR / "EGM/2016postVFP_UL/electron.json.gz",
    "tau": POG_DIR / "TAU/2016postVFP_UL/tau.json.gz",
    "btag": POG_DIR / "BTV/2016postVFP_UL/btagging.json.gz",
    "jerc": POG_DIR / "JME/2016postVFP_UL/jet_jerc.json.gz",
}
MUON_KEYS = {"id": "NUM_TightID_DEN_TrackerMuons", "iso": "NUM_TightRelIso_DEN_TightIDandIPCut",
             "trig": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"}
MUON_PT_RANGE = {"id": (15.0, 119.9), "iso": (15.0, 119.9), "trig": (26.0, 199.9)}
JES_TOTAL = "Summer19UL16_V7_MC_Total_AK4PFchs"


@lru_cache(maxsize=None)
def cset(name: str):
    import correctionlib
    return correctionlib.CorrectionSet.from_file(str(FILES[name]))


def _f(x):
    return np.asarray(x, dtype=np.float64)


# ------------------------------------------------------------------------------ muons
def muon_sf(kind: str, pt, abseta, syst: str = "nominal") -> np.ndarray:
    """Tight-ID / tight-iso / IsoMu24-or-IsoTkMu24 scale factor per muon (kind = 'id', 'iso', 'trig');
    syst 'nominal', 'systup', 'systdown'. pT clipped into the measured range."""
    lo, hi = MUON_PT_RANGE[kind]
    pt = np.clip(_f(pt), lo, hi)
    eta = np.clip(np.abs(_f(abseta)), 0.0, 2.399)
    if len(pt) == 0:
        return np.ones(0)
    return np.asarray(cset("muon")[MUON_KEYS[kind]].evaluate(eta, pt, syst), dtype=np.float64)


# ------------------------------------------------------------------------------ electrons
def electron_sf(kind: str, pt, sceta, syst: str = "sf") -> np.ndarray:
    """Reconstruction ('reco': RecoAbove20 / RecoBelow20) or ID ('id': wp90noiso) scale factor per electron,
    syst 'sf', 'sfup', 'sfdown'. `sceta` is the supercluster eta (eta + deltaEtaSC)."""
    pt = _f(pt)
    eta = np.clip(_f(sceta), -2.499, 2.499)
    out = np.ones(len(pt))
    if len(pt) == 0:
        return out
    c = cset("electron")["UL-Electron-ID-SF"]
    if kind == "reco":
        lo = pt < 20.0
        if lo.any():
            out[lo] = c.evaluate("2016postVFP", syst, "RecoBelow20", eta[lo], np.clip(pt[lo], 10.0, 19.99))
        if (~lo).any():
            out[~lo] = c.evaluate("2016postVFP", syst, "RecoAbove20", eta[~lo], np.clip(pt[~lo], 20.0, 499.0))
        return out
    return np.asarray(c.evaluate("2016postVFP", syst, "wp90noiso", eta, np.clip(pt, 10.0, 499.0)), dtype=np.float64)


# ------------------------------------------------------------------------------ tau_h: e / mu -> tau_h
def tau_vse_sf(wp: str, eta, genmatch, syst: str = "nom") -> np.ndarray:
    """DeepTau VSe scale factor at working point `wp` for e -> tau_h (genmatch 1, 3); 1 elsewhere."""
    eta = np.clip(np.abs(_f(eta)), 0.0, 2.299)
    gm = np.asarray(genmatch, dtype=np.int64)
    out = np.ones(len(eta))
    sel = (gm == 1) | (gm == 3)
    if sel.any():
        out[sel] = cset("tau")["DeepTau2017v2p1VSe"].evaluate(eta[sel], gm[sel], wp, syst)
    return out


def tau_vsmu_sf(wp: str, eta, genmatch, syst: str = "nom") -> np.ndarray:
    """DeepTau VSmu scale factor at working point `wp` for mu -> tau_h (genmatch 2, 4); 1 elsewhere."""
    eta = np.clip(np.abs(_f(eta)), 0.0, 2.299)
    gm = np.asarray(genmatch, dtype=np.int64)
    out = np.ones(len(eta))
    sel = (gm == 2) | (gm == 4)
    if sel.any():
        out[sel] = cset("tau")["DeepTau2017v2p1VSmu"].evaluate(eta[sel], gm[sel], wp, syst)
    return out


# ------------------------------------------------------------------------------ b tagging
def btag_sf(flavour, abseta, pt, wp: str = "M", syst: str = "central") -> np.ndarray:
    """DeepJet fixed-WP scale factor per jet: 'deepJet_comb' for b/c (hadronFlavour 5/4), 'deepJet_incl' for
    light jets. syst 'central', 'up', 'down' (fully correlated between flavours here: one NP `BTag`)."""
    fl = np.asarray(flavour, dtype=np.int64)
    pt = np.clip(_f(pt), 20.0, 999.0)
    eta = np.clip(np.abs(_f(abseta)), 0.0, 2.399)
    out = np.ones(len(pt))
    if len(pt) == 0:
        return out
    heavy = fl > 0
    if heavy.any():
        out[heavy] = cset("btag")["deepJet_comb"].evaluate(syst, wp, fl[heavy], eta[heavy], pt[heavy])
    if (~heavy).any():
        out[~heavy] = cset("btag")["deepJet_incl"].evaluate(syst, wp, fl[~heavy], eta[~heavy], pt[~heavy])
    return out


# ------------------------------------------------------------------------------ jets
def jes_total_unc(eta, pt) -> np.ndarray:
    """Relative total JES uncertainty per jet (Summer19UL16_V7_MC)."""
    pt = np.clip(_f(pt), 9.0, 4999.0)
    eta = np.clip(_f(eta), -5.19, 5.19)
    if len(pt) == 0:
        return np.zeros(0)
    return np.asarray(cset("jerc")[JES_TOTAL].evaluate(eta, pt), dtype=np.float64)


# ------------------------------------------------------------------------------ top pT
def top_pt_weight(pt_top, pt_antitop, power: float = 1.0) -> np.ndarray:
    """CMS Run-2 top pT reweighting, w = sqrt(f(pt_t) f(pt_tbar)), f = exp(0.0615 - 0.0005 pT) with the pT
    capped at 500 GeV; `power` 0 (no reweighting) / 1 (nominal) / 2 (twice) for the shape uncertainty."""
    a = np.clip(_f(pt_top), 0.0, 500.0)
    b = np.clip(_f(pt_antitop), 0.0, 500.0)
    w = np.sqrt(np.exp(0.0615 - 0.0005 * a) * np.exp(0.0615 - 0.0005 * b))
    w = np.where((a > 0) & (b > 0), w, 1.0)
    return w ** power
