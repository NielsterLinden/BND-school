"""Generator-level information: LHE flavour split and the tau_h tau_h fiducial volume.

* LHE flavour: exactly two status-1 LHE leptons of one flavour with opposite charge -> 11 / 13 / 15
  (0 otherwise), and their invariant mass m_LHE (pre-shower, i.e. Born level).
* Fiducial volume (docs/06-cross-section.md): LHE flavour 15, 60 < m_LHE < 120 GeV, and at least two
  visible hadronic taus (`GenVisTau`, which only exists for hadronic decays) with pT > 40 GeV and
  |eta| < 2.1.
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config

LHE_BRANCHES = ["LHEPart_pt", "LHEPart_eta", "LHEPart_phi", "LHEPart_mass", "LHEPart_pdgId", "LHEPart_status"]
VISTAU_BRANCHES = ["GenVisTau_pt", "GenVisTau_eta", "GenVisTau_phi", "GenVisTau_mass", "GenVisTau_charge",
                   "GenVisTau_status", "GenVisTau_genPartIdxMother"]


def _pair_mass(pt, eta, phi, mass, mask):
    rec = ak.zip({"pt": pt, "eta": eta, "phi": phi, "mass": mass})[mask]
    n = ak.to_numpy(ak.num(rec, axis=1))
    rec = ak.pad_none(rec, 2, axis=1)

    def f(x):
        return ak.to_numpy(ak.fill_none(x, 0.0)).astype(np.float64)

    a, b = rec[:, 0], rec[:, 1]
    px = f(a.pt) * np.cos(f(a.phi)) + f(b.pt) * np.cos(f(b.phi))
    py = f(a.pt) * np.sin(f(a.phi)) + f(b.pt) * np.sin(f(b.phi))
    pz = f(a.pt) * np.sinh(f(a.eta)) + f(b.pt) * np.sinh(f(b.eta))
    ea = np.sqrt((f(a.pt) * np.cosh(f(a.eta))) ** 2 + f(a.mass) ** 2)
    eb = np.sqrt((f(b.pt) * np.cosh(f(b.eta))) ** 2 + f(b.mass) ** 2)
    m = np.sqrt(np.maximum((ea + eb) ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))
    return n, m


def lhe_flavour(ev):
    """(flavour per event 11/13/15/0, LHE dilepton mass)."""
    n_ev = len(ev)
    if "LHEPart_pdgId" not in ak.fields(ev):
        return np.zeros(n_ev, dtype=np.uint8), np.zeros(n_ev)
    flav = np.zeros(n_ev, dtype=np.uint8)
    mll = np.zeros(n_ev)
    st1 = ev.LHEPart_status == 1
    for code in (11, 13, 15):
        lep = st1 & (abs(ev.LHEPart_pdgId) == code)
        n, m = _pair_mass(ev.LHEPart_pt, ev.LHEPart_eta, ev.LHEPart_phi, ev.LHEPart_mass, lep)
        charge_sum = ak.to_numpy(ak.sum(np.sign(ev.LHEPart_pdgId[lep]), axis=1))
        sel = (n == 2) & (charge_sum == 0)
        flav[sel] = code
        mll[sel] = m[sel]
    return flav, mll


def fiducial(ev, flav, mll):
    """dict(fid: bool, n_vis: number of visible taus passing the cuts, vis1_pt, vis2_pt)."""
    n_ev = len(ev)
    out = {"fid": np.zeros(n_ev, dtype=bool), "n_vis": np.zeros(n_ev, dtype=np.int32),
           "vis1_pt": np.zeros(n_ev, dtype=np.float32), "vis2_pt": np.zeros(n_ev, dtype=np.float32)}
    if "GenVisTau_pt" not in ak.fields(ev):
        return out
    good = (ev.GenVisTau_pt > config.FID_VIS_PT) & (abs(ev.GenVisTau_eta) < config.FID_VIS_ETA)
    n_vis = ak.to_numpy(ak.sum(good, axis=1)).astype(np.int32)
    pts = ak.sort(ev.GenVisTau_pt, axis=1, ascending=False)
    pts = ak.fill_none(ak.pad_none(pts, 2, axis=1), 0.0)
    out["n_vis"] = n_vis
    out["vis1_pt"] = ak.to_numpy(pts[:, 0]).astype(np.float32)
    out["vis2_pt"] = ak.to_numpy(pts[:, 1]).astype(np.float32)
    out["fid"] = (flav == 15) & (mll > config.MASS_LO) & (mll < config.MASS_HI) & (n_vis >= 2)
    return out
