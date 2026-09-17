"""Objects and pair selection of the lepton channels mu tau_h, e tau_h, e mu (docs/10-v4-plan.md sections 2-4).

Pure functions on one awkward chunk of a v4 skim, in the style of objects.py. The ntuple step stores the
selected objects with their isolation values and identification bits, so the signal region, the anti-isolated
lepton sidebands and the tau_h application region are all cut later (analysis_v4.regions).

    mu  = muons(ev, is_mc)            el = electrons(ev, is_mc)
    tau = ltau_taus(ev, is_mc, tes, "mutau")
    sel = select_ltau(ev, mu_or_el, tau, is_mu=True)          # lepton, tau_h, veto counts, trigger matching
    sel = select_emu(ev, mu, el, role="emu" | "emu_mu" | "emu_el")
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config
from .objects import delta_r

MU_FIELDS = ["pt", "eta", "phi", "mass", "charge", "pfRelIso04_all", "tightId", "mediumId", "looseId", "dxy", "dz"]
EL_FIELDS = ["pt", "eta", "phi", "mass", "charge", "pfRelIso03_all", "mvaFall17V2noIso_WP90", "mvaFall17V2noIso_WP80",
             "cutBased", "deltaEtaSC", "dxy", "dz", "convVeto", "lostHits", "dEscaleUp", "dEscaleDown", "dEsigmaUp", "dEsigmaDown"]
TAU_FIELDS = ["pt", "eta", "phi", "mass", "charge", "decayMode", "dz", "idDeepTau2017v2p1VSjet",
              "idDeepTau2017v2p1VSe", "idDeepTau2017v2p1VSmu", "rawDeepTau2017v2p1VSjet"]


def _match(objs, ev, trig_id, bits, pt_min, dr=config.LEP_TRIG_MATCH_DR):
    """(events, objs) boolean: object within dR of a trigger object of `trig_id` with any of `bits` set."""
    to = ev.TrigObj
    ok = (to.id == trig_id) & ((to.filterBits & bits) > 0) & (to.pt > pt_min)
    to = to[ok]
    o, t = ak.unzip(ak.cartesian([objs, to], axis=1, nested=True))
    return ak.any(delta_r(o.eta, o.phi, t.eta, t.phi) < dr, axis=2)


def muons(ev, is_mc: bool):
    cols = {f: ev[f"Muon_{f}"] for f in MU_FIELDS if f"Muon_{f}" in ak.fields(ev)}
    cols["idx"] = ak.local_index(ev.Muon_pt, axis=1)
    cols["genflav"] = ak.values_astype(ev.Muon_genPartFlav, np.int16) if is_mc else ak.values_astype(ak.full_like(ev.Muon_pt, -1), np.int16)
    mu = ak.zip(cols)
    ip = (abs(mu.dxy) < config.MU_DXY) & (abs(mu.dz) < config.MU_DZ)
    mu["good"] = mu.tightId & ip & (abs(mu.eta) < config.MU_ETA_MAX) & (mu.pfRelIso04_all < config.LEP_ANTIISO[1])
    vm = config.VETO_MU_V4
    mu["veto"] = (mu.looseId & (mu.pt > vm["pt"]) & (abs(mu.eta) < vm["eta"]) & (abs(mu.dxy) < vm["dxy"])
                  & (abs(mu.dz) < vm["dz"]) & (mu.pfRelIso04_all < vm["iso"]))
    mu["match_iso24"] = _match(mu, ev, config.TRIGOBJ_MU_ID, config.TRIGOBJ_MU_ISO_BITS, 24.0)
    mu["match_trkiso8"] = _match(mu, ev, config.TRIGOBJ_MU_ID, config.TRIGOBJ_MU_TRKISOVVL_BIT, 8.0)
    mu["match_trkiso23"] = _match(mu, ev, config.TRIGOBJ_MU_ID, config.TRIGOBJ_MU_TRKISOVVL_BIT, 23.0)
    return mu


def electrons(ev, is_mc: bool):
    cols = {f: ev[f"Electron_{f}"] for f in EL_FIELDS if f"Electron_{f}" in ak.fields(ev)}
    cols["idx"] = ak.local_index(ev.Electron_pt, axis=1)
    cols["genflav"] = ak.values_astype(ev.Electron_genPartFlav, np.int16) if is_mc else ak.values_astype(ak.full_like(ev.Electron_pt, -1), np.int16)
    el = ak.zip(cols)
    sceta = abs(el.eta + el.deltaEtaSC)
    el["sceta"] = el.eta + el.deltaEtaSC
    gap = (sceta > config.EL_GAP[0]) & (sceta < config.EL_GAP[1])
    ip = (abs(el.dxy) < config.EL_DXY) & (abs(el.dz) < config.EL_DZ)
    base = el.mvaFall17V2noIso_WP90 & ip & el.convVeto & (el.lostHits <= 1)
    el["good"] = base & ~gap & (abs(el.eta) < 2.5) & (el.pfRelIso03_all < config.LEP_ANTIISO[1])
    ve = config.VETO_EL_V4
    el["veto"] = (base & (el.pt > ve["pt"]) & (abs(el.eta) < ve["eta"]) & (el.pfRelIso03_all < ve["iso"]))
    el["match_ele27"] = _match(el, ev, config.TRIGOBJ_EL_ID, config.TRIGOBJ_EL_WPTIGHT_BIT, 27.0)
    el["match_emu12"] = _match(el, ev, config.TRIGOBJ_EL_ID, config.TRIGOBJ_EL_EMU_BIT, 12.0)
    el["match_emu23"] = _match(el, ev, config.TRIGOBJ_EL_ID, config.TRIGOBJ_EL_EMU_BIT, 23.0)
    return el


def ltau_taus(ev, is_mc: bool, tes: dict | None, channel: str):
    """tau_h candidates of the lepton channels (nominal energy scale applied to genuine tau_h), loose
    isolation (VVVLoose) and the channel's VSe / VSmu working points."""
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
    t = ak.zip(cols)
    vse, vsmu = (config.MUTAU_VSE_BIT, config.MUTAU_VSMU_BIT) if channel == "mutau" else (config.ETAU_VSE_BIT, config.ETAU_VSMU_BIT)
    dm = t.decayMode
    good = ((t.pt > config.LTAU_TAU_PT_NTUPLE) & (abs(t.eta) < config.LTAU_TAU_ETA_MAX) & (abs(t.dz) < config.TAU_DZ_MAX)
            & ((dm == 0) | (dm == 1) | (dm == 10) | (dm == 11))
            & ((t.idDeepTau2017v2p1VSjet & config.TAU_VSJET_LOOSE_BIT) > 0)
            & ((t.idDeepTau2017v2p1VSe & vse) > 0) & ((t.idDeepTau2017v2p1VSmu & vsmu) > 0))
    return t[good]


def _first(objs, mask):
    """The first object of `objs[mask]` per event (option type) and whether there is one."""
    sel = objs[mask]
    first = ak.firsts(sel)
    return first, ak.to_numpy(ak.fill_none(~ak.is_none(first), False))


def select_ltau(ev, leptons, taus, is_mu: bool):
    """Lepton + tau_h pair: the highest-pT good lepton above the plateau threshold (any isolation < 0.5,
    trigger-matched) and, among tau_h candidates with dR > 0.5 to it, the one with the highest raw VSjet score.
    Returns (lepton record, tau record, has_pair, n_veto_mu_other, n_veto_el_other)."""
    pt_min = (config.MU_PT_MIN_MUTAU if is_mu else config.EL_PT_MIN_ETAU) - 1.0      # 1 GeV margin for scale variations
    eta_max = config.MU_ETA_MAX if is_mu else config.EL_ETA_MAX
    matched = leptons.match_iso24 if is_mu else leptons.match_ele27
    cand = leptons.good & (leptons.pt > pt_min) & (abs(leptons.eta) < eta_max) & matched
    # highest pT first: NanoAOD collections are pT-ordered
    lep, has_lep = _first(leptons, cand)
    l_eta = ak.fill_none(lep.eta, 99.0)
    l_phi = ak.fill_none(lep.phi, 0.0)
    far = delta_r(taus.eta, taus.phi, l_eta, l_phi) > config.LTAU_DR_MIN
    tc = taus[far]
    best = ak.firsts(tc[ak.argmax(tc.rawDeepTau2017v2p1VSjet, axis=1, keepdims=True)])
    has_tau = ak.to_numpy(ak.fill_none(~ak.is_none(best), False))
    l_idx = ak.fill_none(lep.idx, -1)
    n_mu = ak.to_numpy(ak.sum(leptons.veto & (leptons.idx != l_idx), axis=1)) if is_mu else None
    n_el = ak.to_numpy(ak.sum(leptons.veto & (leptons.idx != l_idx), axis=1)) if not is_mu else None
    return lep, best, has_lep & has_tau, n_mu, n_el


def select_emu(ev, mu, el, role: str):
    """e mu pair for the signal channel (`role` 'emu': cross triggers) or the trigger-efficiency side samples
    ('emu_mu': IsoMu24-triggered, muon on the plateau and matched; 'emu_el': Ele27-triggered, electron on the
    plateau and matched). The highest-pT good muon and electron (isolation < 0.5 stored, cut later).
    Returns (muon, electron, has_pair, n_veto_mu_other, n_veto_el_other)."""
    mcand = mu.good & (mu.pt > config.EMU_MU_PT_MIN - 1.0) & (abs(mu.eta) < config.MU_ETA_MAX)
    ecand = el.good & (el.pt > config.EMU_EL_PT_MIN - 1.0) & (abs(el.eta) < 2.5)
    if role == "emu_mu":
        mcand = mcand & (mu.pt > config.MU_PT_MIN_MUTAU) & mu.match_iso24
    if role == "emu_el":
        ecand = ecand & (el.pt > config.EL_PT_MIN_ETAU) & el.match_ele27
    m, has_m = _first(mu, mcand)
    e, has_e = _first(el, ecand)
    dr = delta_r(ak.fill_none(m.eta, 99.0), ak.fill_none(m.phi, 0.0), ak.fill_none(e.eta, -99.0), ak.fill_none(e.phi, 0.0))
    has = has_m & has_e & ak.to_numpy(ak.fill_none(dr > 0.3, False))
    n_mu = ak.to_numpy(ak.sum(mu.veto & (mu.idx != ak.fill_none(m.idx, -1)), axis=1))
    n_el = ak.to_numpy(ak.sum(el.veto & (el.idx != ak.fill_none(e.idx, -1)), axis=1))
    return m, e, has, n_mu, n_el


def count_veto(mu, el):
    """(n veto muons, n veto electrons) of the whole event (used with the tau_h tau_h-like accounting)."""
    return ak.to_numpy(ak.sum(mu.veto, axis=1)), ak.to_numpy(ak.sum(el.veto, axis=1))


def jets_v4(ev, clean_eta, clean_phi, n_lep):
    """Jets cleaned against the selected objects: (njets, nbjets, jet1_pt, JES up/down variants of the counts
    and of the MET shift, and the b-tag inputs of up to 4 jets with pT > 20, |eta| < 2.4).

    clean_eta/phi: lists of per-event (option) arrays of the objects to clean against.
    Returns a dict of flat numpy arrays.
    """
    from . import pog
    j = ak.zip({"pt": ev.Jet_pt, "eta": ev.Jet_eta, "phi": ev.Jet_phi, "id": ev.Jet_jetId, "btag": ev.Jet_btagDeepFlavB,
                "flav": ev.Jet_hadronFlavour if "Jet_hadronFlavour" in ak.fields(ev) else ak.zeros_like(ev.Jet_pt, dtype=np.int32)})
    clean = (j.id & config.JET_ID_BIT) > 0
    for eta, phi in zip(clean_eta, clean_phi):
        clean = clean & (delta_r(j.eta, j.phi, ak.fill_none(eta, 99.0), ak.fill_none(phi, 0.0)) > config.JET_TAU_DR)
    j = j[clean]
    # relative JES uncertainty per jet (flattened evaluation)
    counts = ak.num(j.pt, axis=1)
    flat_unc = pog.jes_total_unc(ak.to_numpy(ak.flatten(j.eta)), ak.to_numpy(ak.flatten(j.pt)))
    j["unc"] = ak.unflatten(flat_unc, counts)
    out = {}
    for tag, scale in (("", 1.0), ("_jesUp", None), ("_jesDown", None)):
        pt = j.pt if scale else (j.pt * (1 + j.unc) if tag == "_jesUp" else j.pt * (1 - j.unc))
        sel = (pt > config.JET_PT) & (abs(j.eta) < config.JET_ETA)
        b = (pt > config.BJET_PT) & (abs(j.eta) < config.BJET_ETA) & (j.btag > config.BTAG_DEEPJET_MEDIUM)
        out[f"njets{tag}"] = ak.to_numpy(ak.sum(sel, axis=1)).astype(np.int8)
        out[f"nbjets{tag}"] = ak.to_numpy(ak.sum(b, axis=1)).astype(np.int8)
        out[f"jet1_pt{tag}"] = ak.to_numpy(ak.fill_none(ak.max(pt[sel], axis=1), 0.0)).astype(np.float32)
        if tag:
            dpt = pt - j.pt
            out[f"met_dx{tag}"] = -ak.to_numpy(ak.sum(dpt * np.cos(j.phi), axis=1)).astype(np.float32)
            out[f"met_dy{tag}"] = -ak.to_numpy(ak.sum(dpt * np.sin(j.phi), axis=1)).astype(np.float32)
    # b-tag weight inputs: jets in the b-tag acceptance, up to 4 per event
    bj = j[(j.pt > config.BJET_PT) & (abs(j.eta) < config.BJET_ETA)]
    bj = ak.pad_none(bj, 4, axis=1, clip=True)
    for i in range(4):
        out[f"bj{i + 1}_pt"] = ak.to_numpy(ak.fill_none(bj.pt[:, i], 0.0)).astype(np.float32)
        out[f"bj{i + 1}_eta"] = ak.to_numpy(ak.fill_none(bj.eta[:, i], 0.0)).astype(np.float32)
        out[f"bj{i + 1}_flav"] = ak.to_numpy(ak.fill_none(bj.flav[:, i], -1)).astype(np.int8)
        out[f"bj{i + 1}_tag"] = ak.to_numpy(ak.fill_none(bj.btag[:, i] > config.BTAG_DEEPJET_MEDIUM, False))
    out["nbjets_acc"] = ak.to_numpy(ak.num(j[(j.pt > config.BJET_PT) & (abs(j.eta) < config.BJET_ETA)].pt, axis=1)).astype(np.int8)
    return out
