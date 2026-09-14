"""Tag-and-probe efficiencies with simultaneous pass/fail fits (docs/12-tag-and-probe-fits.md).

Sequential Muon-POG-style definitions, measured identically in data and simulation:

    ID      probe: loose muon (isGlobal|isTracker, pT > 20, |eta| < 2.4)   pass: tightId + IP cuts
    ISO     probe: tightId + IP                                            pass: pfRelIso04_all < 0.15
    ANTIISO probe: tightId + IP                                            pass: 0.20 < pfRelIso04_all < 1.0
    TRIG    probe: tight (ID + iso)                                        pass: matched to an
            IsoMu24/IsoTkMu24 trigger object (events that fired the OR only)

Tag: tight, pT > 26, trigger-matched (alternative for systematics: pT > 30, iso < 0.10).
Both orientations of every opposite-sign pair with 60 < m < 120 GeV are used; same-sign
pairs are kept for the counting cross-checks.

ID and ISO efficiencies come from a simultaneous extended binned likelihood fit of the pass
and fail mass spectra in every (pT, |eta|) cell: signal = MC template (gen-matched probes,
pass and fail separately) convolved with a Gaussian (shift, width), background = exponential.
Alternatives (background shape, fit range, tag definition, MC-truth counting) give the
systematic uncertainty of the scale factor. The trigger efficiency is a same-sign-subtracted
count in fine pT bins across the turn-on.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import awkward as ak
import numpy as np

from . import config, objects, regions

PT_EDGES = np.array([20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 80.0, 120.0, 200.0])
ETA_EDGES = np.array([0.0, 0.9, 1.2, 2.1, 2.4])
TRIG_PT_EDGES = np.array(config.TRIG_PT_BINS)
NPV_EDGES = np.array([0.0, 15.0, 22.0, 30.0, 100.0])
MASS_LO, MASS_HI, MASS_NBINS = 60.0, 120.0, 120
MASS_EDGES = np.linspace(MASS_LO, MASS_HI, MASS_NBINS + 1)
MASS_CENTRES = 0.5 * (MASS_EDGES[1:] + MASS_EDGES[:-1])
COUNT_WINDOW = (MASS_CENTRES > config.TP_MASS_LO) & (MASS_CENTRES < config.TP_MASS_HI)

EFFS = ("id", "iso", "antiiso")
TAGS = ("nom", "alt")
TAG_ALT_PT, TAG_ALT_ISO = 30.0, 0.10


# ----------------------------------------------------------------------------- filling
def blank():
    npt, neta = len(PT_EDGES) - 1, len(ETA_EDGES) - 1
    out = {"n_files": 0, "n_events": 0, "n_pairs": 0}
    for eff in EFFS:
        for tag in TAGS:
            for charge in ("os", "ss"):
                for pf in ("pass", "fail"):
                    for suffix in ("", "_w2", "_gen"):
                        out[f"{eff}_{tag}_{charge}_{pf}{suffix}"] = np.zeros((npt, neta, MASS_NBINS))
    for variant in ("nominal", "l1veto"):
        for charge in ("os", "ss"):
            for pf in ("pass", "tot"):
                out[f"trig_{variant}_{charge}_{pf}"] = np.zeros((len(TRIG_PT_EDGES) - 1, neta))
                out[f"trig_{variant}_{charge}_{pf}_w2"] = np.zeros((len(TRIG_PT_EDGES) - 1, neta))
    for charge in ("os", "ss"):
        for pf in ("pass", "tot"):
            out[f"isonpv_{charge}_{pf}"] = np.zeros(len(NPV_EDGES) - 1)
    return out


def _pairs(ev):
    idx = ak.local_index(ev.Muon_pt, axis=1)
    p = ak.combinations(idx, 2, fields=["a", "b"])
    tag = ak.concatenate([p.a, p.b], axis=1)
    probe = ak.concatenate([p.b, p.a], axis=1)
    pt, eta, phi, mass = ev.Muon_pt, ev.Muon_eta, ev.Muon_phi, ev.Muon_mass
    px1, py1, pz1, e1 = objects.p4(pt[tag], eta[tag], phi[tag], mass[tag])
    px2, py2, pz2, e2 = objects.p4(pt[probe], eta[probe], phi[probe], mass[probe])
    m = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
    opposite = (ev.Muon_charge[tag] * ev.Muon_charge[probe]) < 0
    return tag, probe, m, opposite


def _flat(x):
    return ak.to_numpy(ak.flatten(x))


def _cell_hist(pt, eta, mass, w, pt_edges):
    """3D histogram (pT, |eta|, mass) with the pT clipped into the last bin."""
    pt = np.clip(pt, pt_edges[0], pt_edges[-1] - 1e-3)
    h, _ = np.histogramdd(np.stack([pt, np.abs(eta), mass], axis=1),
                          bins=[pt_edges, ETA_EDGES, MASS_EDGES], weights=w)
    return h


def fill(ev, out, weight=None, is_mc=False):
    """Fill the pass/fail spectra of one chunk of skim events (events that fired the trigger)."""
    keep = regions.fired(ev) & regions.event_clean(ev)
    ev = ev[keep]
    w_ev = np.ones(len(ev)) if weight is None else np.asarray(weight, dtype=float)[keep]
    if len(ev) == 0:
        return
    masks = regions.muon_masks(ev)
    matched = regions.trigger_match(ev)
    matched_l1 = regions.trigger_match(ev, veto_l1_band=True)
    tight, loose = masks["tight"], masks["loose"]
    ip = (abs(ev.Muon_dxy) < config.MU_DXY_MAX) & (abs(ev.Muon_dz) < config.MU_DZ_MAX)
    id_pass = ev.Muon_tightId & ip
    iso = ev.Muon_pfRelIso04_all
    gen = (ev.Muon_genPartFlav == 1) if is_mc else ak.ones_like(ev.Muon_pt, dtype=bool)

    tag, probe, mass, opposite = _pairs(ev)
    in_range = (mass > MASS_LO) & (mass < MASS_HI)
    out["n_events"] += len(ev)
    tags = {
        "nom": tight & (ev.Muon_pt > config.TAG_PT_MIN) & matched,
        "alt": tight & (ev.Muon_pt > TAG_ALT_PT) & (iso < TAG_ALT_ISO) & matched,
    }
    probes = {"id": loose, "iso": id_pass & (ev.Muon_pt > regions.MU_PT_MIN) & (abs(ev.Muon_eta) < regions.MU_ETA),
              "antiiso": id_pass & (ev.Muon_pt > regions.MU_PT_MIN) & (abs(ev.Muon_eta) < regions.MU_ETA)}
    passes = {"id": id_pass, "iso": iso < regions.ISO_TIGHT,
              "antiiso": (iso > regions.ANTI_ISO_LO) & (iso < regions.ANTI_ISO_HI)}
    n_pairs = ak.to_numpy(ak.num(tag, axis=1))
    w_pair = np.repeat(w_ev, n_pairs)
    for eff in EFFS:
        for tname, is_tag in tags.items():
            base = is_tag[tag] & probes[eff][probe] & in_range
            base_flat = _flat(base)
            for cname, sign in (("os", opposite), ("ss", ~opposite)):
                valid = _flat(base & sign)
                if not valid.any():
                    continue
                ppt = _flat(ev.Muon_pt[probe])[valid]
                peta = _flat(ev.Muon_eta[probe])[valid]
                pm = _flat(mass)[valid]
                pw = w_pair[valid]
                ppass = _flat(passes[eff][probe])[valid]
                pgen = _flat((gen[tag] & gen[probe]))[valid]
                for pf, sel in (("pass", ppass), ("fail", ~ppass)):
                    key = f"{eff}_{tname}_{cname}_{pf}"
                    out[key] += _cell_hist(ppt[sel], peta[sel], pm[sel], pw[sel], PT_EDGES)
                    out[key + "_w2"] += _cell_hist(ppt[sel], peta[sel], pm[sel], pw[sel] ** 2, PT_EDGES)
                    if is_mc:
                        out[key + "_gen"] += _cell_hist(ppt[sel & pgen], peta[sel & pgen], pm[sel & pgen],
                                                        pw[sel & pgen], PT_EDGES)
                if eff == "id" and tname == "nom" and cname == "os":
                    out["n_pairs"] += int(valid.sum())
    # trigger: probe = tight muon, pass = matched; window 70-110 to keep the pairs pure
    window = in_range & (mass > config.TP_MASS_LO) & (mass < config.TP_MASS_HI)
    for variant, match in (("nominal", matched), ("l1veto", matched_l1)):
        is_tag = tight & (ev.Muon_pt > config.TAG_PT_MIN) & match
        base = is_tag[tag] & tight[probe] & window
        for cname, sign in (("os", opposite), ("ss", ~opposite)):
            valid = _flat(base & sign)
            if not valid.any():
                continue
            ppt = np.clip(_flat(ev.Muon_pt[probe])[valid], TRIG_PT_EDGES[0], TRIG_PT_EDGES[-1] - 1e-3)
            peta = np.abs(_flat(ev.Muon_eta[probe])[valid])
            pw = w_pair[valid]
            ppass = _flat(match[probe])[valid]
            for pf, sel in (("tot", np.ones_like(ppass)), ("pass", ppass)):
                h, _, _ = np.histogram2d(ppt[sel], peta[sel], bins=[TRIG_PT_EDGES, ETA_EDGES], weights=pw[sel])
                out[f"trig_{variant}_{cname}_{pf}"] += h
                h2, _, _ = np.histogram2d(ppt[sel], peta[sel], bins=[TRIG_PT_EDGES, ETA_EDGES], weights=pw[sel] ** 2)
                out[f"trig_{variant}_{cname}_{pf}_w2"] += h2
    # isolation efficiency vs pileup (counting, probes with pT > 26, window 70-110)
    is_tag = tags["nom"]
    base = is_tag[tag] & probes["iso"][probe] & window & (ev.Muon_pt[probe] > config.TAG_PT_MIN)
    npv = ak.to_numpy(ev.PV_npvsGood).astype(float)
    npv_pair = np.repeat(npv, n_pairs)
    for cname, sign in (("os", opposite), ("ss", ~opposite)):
        valid = _flat(base & sign)
        if not valid.any():
            continue
        ppass = _flat(passes["iso"][probe])[valid]
        pw = w_pair[valid]
        out[f"isonpv_{cname}_tot"] += np.histogram(npv_pair[valid], bins=NPV_EDGES, weights=pw)[0]
        out[f"isonpv_{cname}_pass"] += np.histogram(npv_pair[valid][ppass], bins=NPV_EDGES, weights=pw[ppass])[0]


# ----------------------------------------------------------------------------- fitting
def _shift_smear(template, shift, sigma):
    """Template (0.5 GeV bins) shifted by `shift` GeV and smeared by a Gaussian of `sigma` GeV."""
    from scipy.ndimage import gaussian_filter1d
    step = MASS_EDGES[1] - MASS_EDGES[0]
    t = gaussian_filter1d(template, max(sigma / step, 1e-3), mode="nearest")
    x = MASS_CENTRES - shift
    t = np.interp(x, MASS_CENTRES, t, left=t[0], right=t[-1])
    s = t.sum()
    return t / s if s > 0 else np.full_like(t, 1.0 / len(t))


def _bkg(kind, params, x, lo):
    if kind == "expo":
        (lam,) = params
        b = np.exp(-lam * (x - lo))
    elif kind == "cmsshape":
        alpha, beta, gamma = params
        from scipy.special import erfc
        b = erfc((alpha - x) * beta) * np.exp(-gamma * (x - lo))
    else:
        raise ValueError(kind)
    s = b.sum()
    return b / s if s > 0 else np.full_like(b, 1.0 / len(b))


def fit_cell(pass_h, fail_h, tmpl_pass, tmpl_fail, pass_var=None, fail_var=None, bkg="expo",
             mass_range=(MASS_LO, MASS_HI)):
    """Simultaneous pass/fail fit; returns dict(eps, err, N, chi2, ok, ...).

    Weighted histograms (MC) are handled with a global scale k = sum(var)/sum(n) per spectrum
    (Bohm-Zech), so the Poisson likelihood applies to n/k.
    """
    from iminuit import Minuit

    sel = (MASS_CENTRES > mass_range[0]) & (MASS_CENTRES < mass_range[1])
    x = MASS_CENTRES[sel]
    lo = mass_range[0]
    tp = np.clip(tmpl_pass, 0, None)[sel]
    tf = np.clip(tmpl_fail, 0, None)[sel]
    tp = tp / tp.sum() if tp.sum() > 0 else np.ones_like(tp) / len(tp)
    tf = tf / tf.sum() if tf.sum() > 0 else tp
    n_p, n_f = pass_h[sel], fail_h[sel]
    k_p = (pass_var[sel].sum() / n_p.sum()) if (pass_var is not None and n_p.sum() > 0) else 1.0
    k_f = (fail_var[sel].sum() / n_f.sum()) if (fail_var is not None and n_f.sum() > 0) else 1.0
    k_p, k_f = max(k_p, 1e-9), max(k_f, 1e-9)
    n_bkg = 1 if bkg == "expo" else 3

    def model(N, eps, dm_p, sg_p, dm_f, sg_f, Bp, Bf, *bp):
        bp_p, bp_f = bp[:n_bkg], bp[n_bkg:]
        sp = _shift_smear(tp, dm_p, sg_p)
        sf = _shift_smear(tf, dm_f, sg_f)
        mu_p = N * eps * sp + Bp * _bkg(bkg, bp_p, x, lo)
        mu_f = N * (1 - eps) * sf + Bf * _bkg(bkg, bp_f, x, lo)
        return mu_p, mu_f

    def nll(*args):
        mu_p, mu_f = model(*args)
        mu_p = np.clip(mu_p, 1e-9, None) / k_p
        mu_f = np.clip(mu_f, 1e-9, None) / k_f
        a = n_p / k_p
        b = n_f / k_f
        return float(np.sum(mu_p - a * np.log(mu_p)) + np.sum(mu_f - b * np.log(mu_f)))

    tot_p, tot_f = float(n_p.sum()), float(n_f.sum())
    eps0 = tot_p / max(tot_p + tot_f, 1.0)
    names = ["N", "eps", "dm_p", "sg_p", "dm_f", "sg_f", "Bp", "Bf"]
    start = [max(tot_p + tot_f, 1.0) * 0.9, min(max(eps0, 0.02), 0.995), 0.0, 0.5, 0.0, 1.0,
             0.05 * max(tot_p, 1.0), 0.3 * max(tot_f, 1.0)]
    if bkg == "expo":
        names += ["lam_p", "lam_f"]
        start += [0.02, 0.02]
        limits = {"lam_p": (-0.1, 0.3), "lam_f": (-0.1, 0.3)}
    else:
        names += ["a_p", "b_p", "g_p", "a_f", "b_f", "g_f"]
        start += [60.0, 0.1, 0.02, 60.0, 0.1, 0.02]
        limits = {k: v for n in ("p", "f") for k, v in
                  ((f"a_{n}", (40.0, 90.0)), (f"b_{n}", (0.01, 1.0)), (f"g_{n}", (-0.1, 0.3)))}
    m = Minuit(nll, *start, name=names)
    m.errordef = Minuit.LIKELIHOOD
    m.limits["N"] = (0, None)
    m.limits["eps"] = (0.0, 1.0)
    m.limits["dm_p"] = (-3.0, 3.0)
    m.limits["dm_f"] = (-5.0, 5.0)
    m.limits["sg_p"] = (0.05, 4.0)
    m.limits["sg_f"] = (0.05, 8.0)
    m.limits["Bp"] = (0, None)
    m.limits["Bf"] = (0, None)
    for k, v in limits.items():
        m.limits[k] = v
    if tot_f < 50:                        # too few failing probes to constrain the fail shape
        m.fixed["dm_f"] = True
        m.fixed["sg_f"] = True
    m.strategy = 1
    m.migrad(ncall=20000)
    if not m.valid:
        m.simplex()
        m.migrad(ncall=20000)
    m.hesse()
    eps, err = m.values["eps"], m.errors["eps"]
    mu_p, mu_f = model(*m.values)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = float(np.nansum((n_p - mu_p) ** 2 / np.where(mu_p > 0, mu_p * k_p, 1))
                     + np.nansum((n_f - mu_f) ** 2 / np.where(mu_f > 0, mu_f * k_f, 1)))
    return {"eps": float(eps), "err": float(err), "N": float(m.values["N"]), "valid": bool(m.valid),
            "chi2": chi2, "ndf": int(2 * sel.sum() - m.nfit), "params": dict(zip(names, m.values)),
            "model_pass": mu_p, "model_fail": mu_f, "x": x,
            "bkg_frac_fail": float((mu_f - m.values["N"] * (1 - eps) * _shift_smear(tf, m.values["dm_f"], m.values["sg_f"])).sum() / max(mu_f.sum(), 1e-9))}


def count_eff(pass_h, fail_h, pass_ss=None, fail_ss=None, window=COUNT_WINDOW):
    """Counting efficiency in the 70-110 window, optionally same-sign subtracted."""
    p, f = pass_h[window].sum(), fail_h[window].sum()
    if pass_ss is not None:
        p -= pass_ss[window].sum()
        f -= fail_ss[window].sum()
    tot = p + f
    eps = p / tot if tot > 0 else 0.0
    err = np.sqrt(max(eps * (1 - eps) / tot, 0.0)) if tot > 0 else 0.0
    return float(eps), float(err)


# ----------------------------------------------------------------------------- scale factors
class ScaleFactors:
    """Per-muon ID / iso scale factors and per-muon trigger efficiencies from the T&P JSON."""

    def __init__(self, path):
        with open(path) as fh:
            d = json.load(fh)
        self.pt_edges = np.array(d["pt_edges"])
        self.eta_edges = np.array(d["eta_edges"])
        self.trig_pt_edges = np.array(d["trig_pt_edges"])
        self.sf = {k: np.array(v) for k, v in d["sf"].items()}         # id, iso: (npt, neta)
        self.sf_err = {k: np.array(v) for k, v in d["sf_err"].items()}
        self.eff = {k: np.array(v) for k, v in d["eff"].items()}       # trig_data, trig_mc, ...
        self.eff_err = {k: np.array(v) for k, v in d["eff_err"].items()}
        self.meta = d.get("meta", {})

    def _cell(self, pt, eta, pt_edges):
        ipt = np.clip(np.digitize(np.clip(pt, pt_edges[0], pt_edges[-1] - 1e-3), pt_edges) - 1, 0, len(pt_edges) - 2)
        ieta = np.clip(np.digitize(np.abs(eta), self.eta_edges) - 1, 0, len(self.eta_edges) - 2)
        return ipt, ieta

    def muon_sf(self, kind, pt, eta, shift=0.0):
        ipt, ieta = self._cell(pt, eta, self.pt_edges)
        return self.sf[kind][ipt, ieta] + shift * self.sf_err[kind][ipt, ieta]

    def trig_eff(self, sample, pt, eta, shift=0.0):
        ipt, ieta = self._cell(pt, eta, self.trig_pt_edges)
        e = self.eff[f"trig_{sample}"][ipt, ieta] + shift * self.eff_err[f"trig_{sample}"][ipt, ieta]
        return np.clip(e, 0.0, 1.0)

    def event_weights(self, pt1, eta1, pt2, eta2):
        """Dict variation -> per-event SF weight for two selected muons."""
        out = {}
        id1, id2 = self.muon_sf("id", pt1, eta1), self.muon_sf("id", pt2, eta2)
        iso1, iso2 = self.muon_sf("iso", pt1, eta1), self.muon_sf("iso", pt2, eta2)

        def trig(shift):
            d1, d2 = self.trig_eff("data", pt1, eta1, shift), self.trig_eff("data", pt2, eta2, shift)
            m1, m2 = self.trig_eff("mc", pt1, eta1), self.trig_eff("mc", pt2, eta2)
            num = 1.0 - (1.0 - d1) * (1.0 - d2)
            den = 1.0 - (1.0 - m1) * (1.0 - m2)
            return np.where(den > 0, num / np.where(den > 0, den, 1.0), 1.0)

        t0 = trig(0.0)
        out["nominal"] = id1 * id2 * iso1 * iso2 * t0
        for name, s in (("Up", +1.0), ("Down", -1.0)):
            out[f"MuonID{name}"] = self.muon_sf("id", pt1, eta1, s) * self.muon_sf("id", pt2, eta2, s) * iso1 * iso2 * t0
            out[f"MuonIso{name}"] = id1 * id2 * self.muon_sf("iso", pt1, eta1, s) * self.muon_sf("iso", pt2, eta2, s) * t0
            out[f"MuonTrigger{name}"] = id1 * id2 * iso1 * iso2 * trig(s)
        return out

    def single_muon_weights(self, pt, eta):
        """For the e-mu region: one muon, trigger SF = eff_data/eff_mc."""
        out = {}
        idw, isow = self.muon_sf("id", pt, eta), self.muon_sf("iso", pt, eta)
        m = self.trig_eff("mc", pt, eta)
        tr = lambda s: np.where(m > 0, self.trig_eff("data", pt, eta, s) / np.where(m > 0, m, 1), 1.0)
        out["nominal"] = idw * isow * tr(0.0)
        for name, s in (("Up", 1.0), ("Down", -1.0)):
            out[f"MuonID{name}"] = self.muon_sf("id", pt, eta, s) * isow * tr(0.0)
            out[f"MuonIso{name}"] = idw * self.muon_sf("iso", pt, eta, s) * tr(0.0)
            out[f"MuonTrigger{name}"] = idw * isow * tr(s)
        return out
