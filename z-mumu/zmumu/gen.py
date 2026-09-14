"""Generator-level quantities for the MC skims: LHE flavour split and the fiducial volume.

The definitions are those of the review (`scripts/mc_acceptance.py`) so that v2 numbers
are directly comparable:

  * LHE flavour: exactly two status-1 LHE leptons of the same flavour -> 11 / 13 / 15, else 0.
  * Fiducial (dressed): exactly two `GenDressedLepton` muons (dR < 0.1 dressing, NanoAOD),
    no tau ancestor, pT > 20, |eta| < 2.4, leading pT > 26, opposite charge,
    60 < m < 120 GeV, in an LHE Z -> mu mu event.
  * Fiducial (Born): the same kinematic cuts on the two LHE muons.

`GenDressedLepton` is only stored above a pT threshold, so the fiducial *denominators*
use the LHE (Born) dimuon mass, never a "dressed 60-120" requirement.
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config, objects

LHE_BRANCHES = ["LHEPart_pt", "LHEPart_eta", "LHEPart_phi", "LHEPart_mass", "LHEPart_pdgId",
                "LHEPart_status"]
DRESSED_BRANCHES = ["GenDressedLepton_pt", "GenDressedLepton_eta", "GenDressedLepton_phi",
                    "GenDressedLepton_mass", "GenDressedLepton_pdgId", "GenDressedLepton_hasTauAnc"]


def leading_pair(pt, eta, phi, mass, pdg, mask):
    """(n passing, leading pT, subleading pT, pair mass, pdg product) per event, 0 where absent."""
    rec = ak.zip({"pt": pt, "eta": eta, "phi": phi, "mass": mass, "pdg": pdg})[mask]
    rec = rec[ak.argsort(rec.pt, axis=1, ascending=False)]
    n = ak.to_numpy(ak.num(rec, axis=1))
    rec = ak.pad_none(rec, 2, axis=1)
    a, b = rec[:, 0], rec[:, 1]

    def f(x):
        return ak.to_numpy(ak.fill_none(x, 0.0)).astype(float)

    px1, py1, pz1, e1 = objects.p4(f(a.pt), f(a.eta), f(a.phi), f(a.mass))
    px2, py2, pz2, e2 = objects.p4(f(b.pt), f(b.eta), f(b.phi), f(b.mass))
    m = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
    return n, f(a.pt), f(b.pt), m, f(a.pdg) * f(b.pdg)


def lhe_flavour(ev):
    """(flavour per event: 11/13/15/0, LHE dilepton mass, n LHE leptons)."""
    fields = set(ak.fields(ev))
    n_ev = len(ev)
    if "LHEPart_pdgId" not in fields:
        return np.zeros(n_ev, dtype=np.uint8), np.zeros(n_ev), np.zeros(n_ev, dtype=np.int32)
    status1 = ev.LHEPart_status == 1
    pdg = abs(ev.LHEPart_pdgId)
    flav = np.zeros(n_ev, dtype=np.uint8)
    mll = np.zeros(n_ev)
    nlep = ak.to_numpy(ak.sum(status1 & ((pdg == 11) | (pdg == 13) | (pdg == 15)), axis=1)).astype(np.int32)
    for code in (11, 13, 15):
        mask = status1 & (pdg == code)
        n, _, _, m, prod = leading_pair(ev.LHEPart_pt, ev.LHEPart_eta, ev.LHEPart_phi,
                                        ev.LHEPart_mass, ev.LHEPart_pdgId, mask)
        sel = (n == 2) & (prod < 0)
        flav[sel] = code
        mll[sel] = m[sel]
    return flav, mll, nlep


def fiducial_flags(ev, flav, mll):
    """Dressed and Born fiducial flags and the dressed dimuon kinematics (per event)."""
    n_ev = len(ev)
    zeros = np.zeros(n_ev)
    out = {"fid_dressed": np.zeros(n_ev, dtype=bool), "fid_born": np.zeros(n_ev, dtype=bool),
           "m_dressed": zeros.copy(), "pt1_dressed": zeros.copy(), "pt2_dressed": zeros.copy(),
           "n_dressed_mu": np.zeros(n_ev, dtype=np.int32)}
    is_mumu = flav == 13
    fields = set(ak.fields(ev))
    if "GenDressedLepton_pt" in fields:
        dressed = ((abs(ev.GenDressedLepton_pdgId) == 13) & ~ev.GenDressedLepton_hasTauAnc
                   & (ev.GenDressedLepton_pt > config.MU_PT_SUBLEAD)
                   & (abs(ev.GenDressedLepton_eta) < config.MU_ETA_MAX))
        n_d, pt1, pt2, m_d, prod = leading_pair(ev.GenDressedLepton_pt, ev.GenDressedLepton_eta,
                                                ev.GenDressedLepton_phi, ev.GenDressedLepton_mass,
                                                ev.GenDressedLepton_pdgId, dressed)
        fid = (is_mumu & (n_d == config.N_MUONS_REQUIRED) & (pt1 > config.MU_PT_LEAD) & (prod < 0)
               & (m_d > config.MASS_LO) & (m_d < config.MASS_HI))
        out.update(fid_dressed=fid, m_dressed=m_d, pt1_dressed=pt1, pt2_dressed=pt2, n_dressed_mu=n_d)
    if "LHEPart_pt" in fields:
        lhe_mu = ((ev.LHEPart_status == 1) & (abs(ev.LHEPart_pdgId) == 13)
                  & (ev.LHEPart_pt > config.MU_PT_SUBLEAD) & (abs(ev.LHEPart_eta) < config.MU_ETA_MAX))
        n_b, pt1b, _, m_b, prodb = leading_pair(ev.LHEPart_pt, ev.LHEPart_eta, ev.LHEPart_phi,
                                                ev.LHEPart_mass, ev.LHEPart_pdgId, lhe_mu)
        out["fid_born"] = (is_mumu & (n_b == 2) & (pt1b > config.MU_PT_LEAD) & (prodb < 0)
                           & (m_b > config.MASS_LO) & (m_b < config.MASS_HI))
    return out
