"""Object and region definitions of the v2 analysis, evaluated on a skim chunk.

Single source of truth for everything that is selected: the tag-and-probe, the fake
factor, the histogramming and the fit all import from here. See docs/09-v2-overview.md.

Muon classes (all: pT > 20, |eta| < 2.4, |dxy| < 0.2 cm, |dz| < 0.5 cm)
    loose   (isGlobal | isTracker)                      T&P probes, skim category
    tight   tightId & pfRelIso04_all < 0.15              signal muons
    anti    tightId & 0.20 < pfRelIso04_all < 1.0        fake-factor "anti-tight" muons
Trigger matching: dR < 0.1 to a TrigObj (id 13, filterBits & (2|8), pT > 24), only in events
that fired HLT_IsoMu24 || HLT_IsoTkMu24.

Regions (all after MET filters and PV_npvsGood >= 1; masses are FSR-recovered):
    SR      exactly two tight muons, leading pT > 26, >= 1 trigger-matched, OS, 60 < m < 120
    SS      same, same-sign (fake-factor closure, plots)
    CRemu   exactly one tight muon (pT > 26, matched), no other tight/anti muon, exactly one
            electron (pT > 20, |eta| < 2.5, cutBased >= 3, not in the barrel-endcap gap), OS,
            60 < m(e mu) < 120
    SSemu   same, same-sign (non-prompt-electron check: W+jets, conversions)
    FFapp   exactly one tight (pT > 26, matched) and exactly one anti muon, no other tight, OS,
            60 < m < 120    (fake-factor application region; SS variant for the closure)
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config, objects

MU_PT_MIN, MU_ETA = 20.0, 2.4
ISO_TIGHT = 0.15
ANTI_ISO_LO, ANTI_ISO_HI = 0.20, 1.0
TRIG_BITS, TRIG_OBJ_PT, TRIG_DR = 2 | 8, 24.0, 0.1
EL_PT, EL_ETA, EL_ID = 20.0, 2.5, 3
EL_GAP = (1.444, 1.566)


# ----------------------------------------------------------------------------- objects
def with_muon_pt(ev, pt):
    """Copy of the chunk with a replaced Muon_pt (momentum-scale variations)."""
    return ak.with_field(ev, pt, "Muon_pt")


def muon_masks(ev):
    ip = (abs(ev.Muon_dxy) < config.MU_DXY_MAX) & (abs(ev.Muon_dz) < config.MU_DZ_MAX)
    acc = (ev.Muon_pt > MU_PT_MIN) & (abs(ev.Muon_eta) < MU_ETA) & ip
    loose = (ev.Muon_isGlobal | ev.Muon_isTracker) & (ev.Muon_pt > MU_PT_MIN) & (abs(ev.Muon_eta) < MU_ETA)
    tight = acc & ev.Muon_tightId & (ev.Muon_pfRelIso04_all < ISO_TIGHT)
    anti = (acc & ev.Muon_tightId & (ev.Muon_pfRelIso04_all > ANTI_ISO_LO)
            & (ev.Muon_pfRelIso04_all < ANTI_ISO_HI))
    return {"loose": loose, "tight": tight, "anti": anti}


def trigger_match(ev, veto_l1_band=False):
    """(events, muons) boolean: muon matched to an IsoMu24/IsoTkMu24 trigger object."""
    obj = (ev.TrigObj_id == 13) & ((ev.TrigObj_filterBits & TRIG_BITS) > 0) & (ev.TrigObj_pt > TRIG_OBJ_PT)
    if veto_l1_band:
        lo, hi = config.TRIG_L1_VETO_BAND
        obj = obj & ~((ev.TrigObj_l1pt >= lo) & (ev.TrigObj_l1pt < hi))
    mu = ak.zip({"eta": ev.Muon_eta, "phi": ev.Muon_phi})
    to = ak.zip({"eta": ev.TrigObj_eta[obj], "phi": ev.TrigObj_phi[obj]})
    m, t = ak.unzip(ak.cartesian([mu, to], axis=1, nested=True))
    dr = objects.delta_r(m.eta, m.phi, t.eta, t.phi)
    fired = ak.to_numpy(ev.HLT_IsoMu24 | ev.HLT_IsoTkMu24)
    return ak.any(dr < TRIG_DR, axis=2) & fired


def electron_mask(ev):
    sc_eta = abs(ev.Electron_eta + ev.Electron_deltaEtaSC)
    return ((ev.Electron_pt > EL_PT) & (abs(ev.Electron_eta) < EL_ETA) & (ev.Electron_cutBased >= EL_ID)
            & ~((sc_eta > EL_GAP[0]) & (sc_eta < EL_GAP[1])))


def event_clean(ev):
    flags = [f for f in config.MET_FILTERS + ["Flag_BadPFMuonDzFilter"] if f in ak.fields(ev)]
    ok = ak.to_numpy(ev.PV_npvsGood >= 1)
    for f in flags:
        ok = ok & ak.to_numpy(ev[f])
    return ok


def fired(ev):
    return ak.to_numpy(ev.HLT_IsoMu24 | ev.HLT_IsoTkMu24)


# ----------------------------------------------------------------------------- dimuon regions
def _pair(ev, mask, matched):
    """Events with exactly two muons passing `mask`: indices, kinematics, charge product."""
    n = ak.to_numpy(ak.sum(mask, axis=1))
    idx = np.nonzero(n == 2)[0]
    if idx.size == 0:
        return None
    sub, m = ev[idx], mask[idx]
    i1, i2 = objects.leading_two(sub, m)
    out = {
        "idx": idx,
        "pt1": ak.to_numpy(objects.take(sub.Muon_pt, i1)), "pt2": ak.to_numpy(objects.take(sub.Muon_pt, i2)),
        "eta1": ak.to_numpy(objects.take(sub.Muon_eta, i1)), "eta2": ak.to_numpy(objects.take(sub.Muon_eta, i2)),
        "phi1": ak.to_numpy(objects.take(sub.Muon_phi, i1)), "phi2": ak.to_numpy(objects.take(sub.Muon_phi, i2)),
        "q1": ak.to_numpy(objects.take(sub.Muon_charge, i1)), "q2": ak.to_numpy(objects.take(sub.Muon_charge, i2)),
        "iso1": ak.to_numpy(objects.take(sub.Muon_pfRelIso04_all, i1)),
        "iso2": ak.to_numpy(objects.take(sub.Muon_pfRelIso04_all, i2)),
        "match1": ak.to_numpy(objects.take(matched[idx], i1)), "match2": ak.to_numpy(objects.take(matched[idx], i2)),
        "mass": ak.to_numpy(objects.dimuon_mass(sub, i1, i2, with_fsr=config.FSR_ENABLED)),
        "mass_bare": ak.to_numpy(objects.dimuon_mass(sub, i1, i2, with_fsr=False)),
    }
    zpt, zy = objects.dimuon_pt_y(sub, i1, i2)
    out["zpt"], out["zy"] = ak.to_numpy(zpt), ak.to_numpy(zy)
    if "Muon_genPartFlav" in ak.fields(ev):
        out["flav1"] = ak.to_numpy(objects.take(sub.Muon_genPartFlav, i1))
        out["flav2"] = ak.to_numpy(objects.take(sub.Muon_genPartFlav, i2))
    return out


def dimuon_regions(ev, masks, matched, clean):
    """SR / SS selections on tight muons. Returns dict region -> pair dict (indices into ev)."""
    tight = masks["tight"]
    pair = _pair(ev, tight, matched)
    out = {"SR": None, "SS": None}
    if pair is None:
        return out
    base = (clean[pair["idx"]] & (pair["pt1"] > config.MU_PT_LEAD) & (pair["match1"] | pair["match2"])
            & (pair["mass"] > config.MASS_LO) & (pair["mass"] < config.MASS_HI))
    os_ = pair["q1"] * pair["q2"] < 0
    for name, sel in (("SR", base & os_), ("SS", base & ~os_)):
        out[name] = {k: (v[sel] if isinstance(v, np.ndarray) else v) for k, v in pair.items()}
    return out


def ff_application(ev, masks, matched, clean):
    """One tight (leading, matched) + one anti-tight muon, no other tight; OS and SS."""
    tight, anti = masks["tight"], masks["anti"]
    n_t = ak.to_numpy(ak.sum(tight, axis=1))
    n_a = ak.to_numpy(ak.sum(anti, axis=1))
    both = tight | anti
    pair = _pair(ev, both, matched)
    out = {"FFapp_OS": None, "FFapp_SS": None}
    if pair is None:
        return out
    idx = pair["idx"]
    ok = (n_t[idx] == 1) & (n_a[idx] == 1) & clean[idx]
    # the tight muon must be the trigger-matched one with pT > 26; identify which of the two it is
    t1 = ak.to_numpy(objects.take(tight[idx], objects.leading_two(ev[idx], both[idx])[0]))
    pt_t = np.where(t1, pair["pt1"], pair["pt2"])
    match_t = np.where(t1, pair["match1"], pair["match2"])
    pt_a, eta_a = np.where(t1, pair["pt2"], pair["pt1"]), np.where(t1, pair["eta2"], pair["eta1"])
    ok &= (pt_t > config.MU_PT_LEAD) & match_t & (pair["mass"] > config.MASS_LO) & (pair["mass"] < config.MASS_HI)
    os_ = pair["q1"] * pair["q2"] < 0
    for name, sel in (("FFapp_OS", ok & os_), ("FFapp_SS", ok & ~os_)):
        d = {k: (v[sel] if isinstance(v, np.ndarray) else v) for k, v in pair.items()}
        d["pt_anti"], d["eta_anti"] = pt_a[sel], eta_a[sel]
        if "flav1" in pair:
            d["flav_anti"] = np.where(t1, pair["flav2"], pair["flav1"])[sel]
            d["flav_tight"] = np.where(t1, pair["flav1"], pair["flav2"])[sel]
        out[name] = d
    return out


def emu_region(ev, masks, matched, clean, same_sign=False):
    """One tight trigger-matched muon (pT > 26), no other tight/anti muon, one electron, OS
    (or SS with `same_sign`)."""
    tight, anti = masks["tight"], masks["anti"]
    el = electron_mask(ev)
    n_t = ak.to_numpy(ak.sum(tight, axis=1))
    n_a = ak.to_numpy(ak.sum(anti, axis=1))
    n_e = ak.to_numpy(ak.sum(el, axis=1))
    idx = np.nonzero((n_t == 1) & (n_a == 0) & (n_e == 1) & clean)[0]
    if idx.size == 0:
        return None
    sub = ev[idx]
    t = tight[idx]
    e = el[idx]
    mu_pt = ak.to_numpy(ak.flatten(sub.Muon_pt[t]))
    mu_eta = ak.to_numpy(ak.flatten(sub.Muon_eta[t]))
    mu_phi = ak.to_numpy(ak.flatten(sub.Muon_phi[t]))
    mu_m = ak.to_numpy(ak.flatten(sub.Muon_mass[t]))
    mu_q = ak.to_numpy(ak.flatten(sub.Muon_charge[t]))
    mu_match = ak.to_numpy(ak.flatten(matched[idx][t]))
    el_pt = ak.to_numpy(ak.flatten(sub.Electron_pt[e]))
    el_eta = ak.to_numpy(ak.flatten(sub.Electron_eta[e]))
    el_phi = ak.to_numpy(ak.flatten(sub.Electron_phi[e]))
    el_m = ak.to_numpy(ak.flatten(sub.Electron_mass[e]))
    el_q = ak.to_numpy(ak.flatten(sub.Electron_charge[e]))
    px1, py1, pz1, e1 = objects.p4(mu_pt, mu_eta, mu_phi, mu_m)
    px2, py2, pz2, e2 = objects.p4(el_pt, el_eta, el_phi, el_m)
    mass = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
    charge = (mu_q * el_q > 0) if same_sign else (mu_q * el_q < 0)
    sel = ((mu_pt > config.MU_PT_LEAD) & mu_match & charge
           & (mass > config.MASS_LO) & (mass < config.MASS_HI))
    out = {"idx": idx[sel], "pt1": mu_pt[sel], "eta1": mu_eta[sel], "pt_el": el_pt[sel], "eta_el": el_eta[sel],
           "mass": mass[sel], "match1": mu_match[sel]}
    if "Muon_genPartFlav" in ak.fields(ev):
        out["flav1"] = ak.to_numpy(ak.flatten(sub.Muon_genPartFlav[t]))[sel]
        out["flav_el"] = ak.to_numpy(ak.flatten(sub.Electron_genPartFlav[e]))[sel]
    return out


def clean_jet_count(ev, masks, pt_min=30.0, eta_max=2.4, dr=0.4):
    """Per-event number of jets (pT > 30, |eta| < 2.4, tight jet ID) with dR > 0.4 to every
    tight/anti-tight muon and every selected electron (the leptons of all regions)."""
    jet = (ev.Jet_pt > pt_min) & (abs(ev.Jet_eta) < eta_max) & (ev.Jet_jetId >= 2)
    lep = masks["tight"] | masks["anti"]
    el = electron_mask(ev)
    lep_eta = ak.concatenate([ev.Muon_eta[lep], ev.Electron_eta[el]], axis=1)
    lep_phi = ak.concatenate([ev.Muon_phi[lep], ev.Electron_phi[el]], axis=1)
    j, l = ak.unzip(ak.cartesian([ak.zip({"eta": ev.Jet_eta, "phi": ev.Jet_phi}),
                                  ak.zip({"eta": lep_eta, "phi": lep_phi})], axis=1, nested=True))
    near = ak.any(objects.delta_r(j.eta, j.phi, l.eta, l.phi) < dr, axis=2)
    return np.asarray(ak.sum(jet & ~near, axis=1), dtype=float)


def is_prompt(flav):
    """NanoAOD genPartFlav: 1 = prompt, 15 = from tau decay -> treated as prompt."""
    return (flav == 1) | (flav == 15)
