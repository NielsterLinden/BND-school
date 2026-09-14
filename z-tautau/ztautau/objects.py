"""Object selection and the tau_h tau_h pair choice (docs/02-selection.md).

Pure functions on one awkward chunk of the skim. The chain used by step 2 is

    taus  = tau_candidates(ev)             # ID/kinematic preselection, loose isolation (VVVLoose)
    pair  = best_pair(taus)                # most isolated pair with dR > 0.5, then ordered in pT
    match = trigger_matched(ev, pair.t1)   # both legs matched to a di-tau trigger object
    veto  = extra_lepton_veto(ev)          # orthogonal to e tau_h / mu tau_h / ee / mumu

Pitfall worth knowing (it cost the z-mumu group real time): indexing a jagged array with a flat
integer array selects *events*, not elements. The pair choice therefore uses `ak.argmax(...,
keepdims=True)` and `ak.firsts`, which pick one element per event.
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config


# ------------------------------------------------------------------------------ kinematics
def delta_phi(a, b):
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    return np.sqrt((eta1 - eta2) ** 2 + delta_phi(phi1, phi2) ** 2)


def p4(pt, eta, phi, mass):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e = np.sqrt(px ** 2 + py ** 2 + pz ** 2 + mass ** 2)
    return px, py, pz, e


def inv_mass(px, py, pz, e):
    return np.sqrt(np.maximum(e ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))


# ------------------------------------------------------------------------------ taus
TAU_FIELDS = ["pt", "eta", "phi", "mass", "charge", "decayMode", "dz", "idDeepTau2017v2p1VSjet",
              "idDeepTau2017v2p1VSe", "idDeepTau2017v2p1VSmu", "rawDeepTau2017v2p1VSjet"]


def tau_collection(ev, is_mc: bool, tes: dict | None = None):
    """Tau record array with the nominal energy scale applied to genuine simulated tau_h.

    `tes` = {dm: scale}; genuine = genPartFlav 5. Adds `pt_raw` (uncorrected) and `genflav`
    (-1 in data).
    """
    cols = {f: ev[f"Tau_{f}"] for f in TAU_FIELDS}
    cols["idx"] = ak.local_index(ev.Tau_pt, axis=1)
    cols["pt_raw"] = ev.Tau_pt
    if is_mc:
        flav = ev.Tau_genPartFlav
        cols["genflav"] = ak.values_astype(flav, np.int16)
        if tes:
            scale = ak.ones_like(ev.Tau_pt)
            for dm, s in tes.items():
                scale = ak.where((flav == 5) & (ev.Tau_decayMode == dm), s, scale)
            cols["pt"] = ev.Tau_pt * scale
            cols["mass"] = ev.Tau_mass * scale
    else:
        cols["genflav"] = ak.values_astype(ak.full_like(ev.Tau_pt, -1), np.int16)
    return ak.zip(cols)


def tau_candidates(taus):
    dm = taus.decayMode
    good_dm = (dm == 0) | (dm == 1) | (dm == 10) | (dm == 11)
    return taus[(taus.pt > config.TAU_PT_NTUPLE) & (abs(taus.eta) < config.TAU_ETA_MAX)
                & (abs(taus.dz) < config.TAU_DZ_MAX) & good_dm
                & ((taus.idDeepTau2017v2p1VSe & config.TAU_VSE_BIT) > 0)
                & ((taus.idDeepTau2017v2p1VSmu & config.TAU_VSMU_BIT) > 0)
                & ((taus.idDeepTau2017v2p1VSjet & config.TAU_VSJET_LOOSE_BIT) > 0)]


def best_pair(cands):
    """(t1, t2, has_pair): the pair with the most isolated tau (then the more isolated second tau, then
    the larger scalar pT sum) among pairs with dR > PAIR_DR_MIN; t1 is the leading-pT tau of the pair.
    t1/t2 are option-type records (None where there is no pair)."""
    pairs = ak.combinations(cands, 2, fields=["a", "b"])
    pairs = pairs[delta_r(pairs.a.eta, pairs.a.phi, pairs.b.eta, pairs.b.phi) > config.PAIR_DR_MIN]
    ia, ib = pairs.a.rawDeepTau2017v2p1VSjet, pairs.b.rawDeepTau2017v2p1VSjet
    hi = np.maximum(ia, ib)
    lo = np.minimum(ia, ib)
    score = ak.values_astype(hi, np.float64) + 1e-4 * lo + 1e-9 * (pairs.a.pt + pairs.b.pt)
    best = ak.firsts(pairs[ak.argmax(score, axis=1, keepdims=True)])
    has_pair = ~ak.is_none(best)
    a_lead = best.a.pt >= best.b.pt
    t1 = ak.where(a_lead, best.a, best.b)
    t2 = ak.where(a_lead, best.b, best.a)
    return t1, t2, ak.to_numpy(ak.fill_none(has_pair, False))


def trigger_matched(ev, tau):
    """True where `tau` (one option record per event) has a di-tau trigger object within dR."""
    to = ev.TrigObj
    ok = (to.id == config.TRIGOBJ_TAU_ID) & ((to.filterBits & config.TRIGOBJ_TAU_BITS) > 0) & (to.pt > config.TRIGOBJ_PT_MIN)
    to = to[ok]
    eta = ak.fill_none(tau.eta, 99.0)
    phi = ak.fill_none(tau.phi, 0.0)
    m = ak.any(delta_r(to.eta, to.phi, eta, phi) < config.TRIG_MATCH_DR, axis=1)
    return ak.to_numpy(ak.fill_none(m, False))


# ------------------------------------------------------------------------------ vetoes and jets
def extra_lepton_veto(ev):
    """True where the event has no veto muon and no veto electron."""
    vm, ve = config.VETO_MU, config.VETO_EL
    mu = (ev.Muon_mediumId & (ev.Muon_pt > vm["pt"]) & (abs(ev.Muon_eta) < vm["eta"])
          & (abs(ev.Muon_dxy) < vm["dxy"]) & (abs(ev.Muon_dz) < vm["dz"]) & (ev.Muon_pfRelIso04_all < vm["iso"]))
    el = (ev.Electron_mvaFall17V2noIso_WP90 & (ev.Electron_pt > ve["pt"]) & (abs(ev.Electron_eta) < ve["eta"])
          & (abs(ev.Electron_dxy) < ve["dxy"]) & (abs(ev.Electron_dz) < ve["dz"])
          & (ev.Electron_pfRelIso03_all < ve["iso"]) & ev.Electron_convVeto & (ev.Electron_lostHits <= 1))
    return ak.to_numpy((ak.sum(mu, axis=1) == 0) & (ak.sum(el, axis=1) == 0))


def jets(ev, t1, t2):
    """(n jets pT > 30, n b-tagged jets pT > 20, leading jet pT), jets cleaned against the pair."""
    e1, p1 = ak.fill_none(t1.eta, 99.0), ak.fill_none(t1.phi, 0.0)
    e2, p2 = ak.fill_none(t2.eta, 99.0), ak.fill_none(t2.phi, 0.0)
    clean = ((delta_r(ev.Jet_eta, ev.Jet_phi, e1, p1) > config.JET_TAU_DR)
             & (delta_r(ev.Jet_eta, ev.Jet_phi, e2, p2) > config.JET_TAU_DR)
             & ((ev.Jet_jetId & config.JET_ID_BIT) > 0))
    j = clean & (ev.Jet_pt > config.JET_PT) & (abs(ev.Jet_eta) < config.JET_ETA)
    b = clean & (ev.Jet_pt > config.BJET_PT) & (abs(ev.Jet_eta) < config.BJET_ETA) & (ev.Jet_btagDeepFlavB > config.BTAG_DEEPJET_MEDIUM)
    lead = ak.fill_none(ak.max(ev.Jet_pt[j], axis=1), 0.0)
    return ak.to_numpy(ak.sum(j, axis=1)), ak.to_numpy(ak.sum(b, axis=1)), ak.to_numpy(lead)
