"""Muon reconstruction efficiency by tag-and-probe on the unskimmed NanoAOD (docs/16).

Until 16 Sep 2026 the analysis assigned SF_reco = 1 +- 0.4 % per muon ("not measurable in
NanoAOD", docs/11). It is measurable: NanoAODv9 stores stand-alone-only muons
(`Muon_isStandalone` without `isGlobal`/`isTracker`, pT > 15 GeV) and isolated tracks
(`IsoTrack`, cleaned of loose muons). The probability that a real muon becomes a loose muon
(`isGlobal | isTracker`, the denominator of our ID tag-and-probe, docs/12) factorises as

    eps_reco = eps_trk (inner track | muon-system track) x eps_mu (global/tracker muon | inner track)

and both factors have an unbiased probe in NanoAOD, as in the CMS Muon POG and in
CMS-SMP-20-004 (arXiv:2408.03744, section 7):

    trk    probe: stand-alone muon (isStandalone, pT > 20, |eta| < 2.4)
           pass:  the probe is global or tracker, or a global/tracker muon without its own
                  stand-alone track lies within dR < 0.3 (track found but not linked).
                  A stand-alone-only probe with a *stand-alone-carrying* global/tracker muon
                  within dR < 0.3 is a duplicate of that probe and is dropped.
    mutrk  probe: isolated track, pT > 20, |eta| < 2.4, |dxy| < 0.2, |dz| < 0.5,
                  pfRelIso03_chg < 0.10, from either collection:
           pass:  global/tracker muon with those cuts (its own track)
           fail:  IsoTrack with those cuts and no global/tracker muon within dR < 0.05
    mutrkT the same with a purer probe (the nominal for the scale factor): pfRelIso03_all < 0.05,
           |dxy| < 0.05, |dz| < 0.2, in events with PuppiMET < 40 GeV. The loose probe leaves a
           W+jets / ttbar background under the Z in the failing probes that the fit models only to
           ~0.3 %; `mutrk` is kept as a cross-check.

Tag: tight muon (tightId, pfRelIso04_all < 0.15, IP cuts), pT > 26, matched to an
IsoMu24/IsoTkMu24 trigger object; alternative tag pT > 30 and iso < 0.10. dR(tag, probe) > 0.3.
The spectra use the same 0.5 GeV mass bins as `zmumu.tnp`, so `tnp.fit_cell` fits them
unchanged (pass and fail templates from generator-matched simulation).
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config, io, objects, regions, tnp

PT_EDGES = np.array([20.0, 30.0, 40.0, 50.0, 60.0, 200.0])
ETA_EDGES = tnp.ETA_EDGES
EFFS = ("trk", "mutrk", "mutrkT")
TRUTH_OF = {"trk": "trk", "mutrk": "mutrk", "mutrkT": "mutrk"}
TAGS = ("nom", "alt")
TAG_ALT_PT, TAG_ALT_ISO = 30.0, 0.10
TAG_PROBE_DR = 0.3
STA_MATCH_DR = 0.3
TRK_MATCH_DR = 0.05
TRK_ISO_MAX = 0.10
TIGHT_ISO_ALL_MAX, TIGHT_DXY_MAX, TIGHT_DZ_MAX, TIGHT_MET_MAX = 0.05, 0.05, 0.2, 40.0
GEN_DR = {"trk": 0.3, "mutrk": 0.1, "mutrkT": 0.1}

BRANCHES_DATA = (["run", "luminosityBlock", "PV_npvsGood", "HLT_IsoMu24", "HLT_IsoTkMu24", "PuppiMET_pt"] + config.MET_FILTERS
                 + ["Flag_BadPFMuonDzFilter"]
                 + [f"Muon_{f}" for f in ("pt", "eta", "phi", "mass", "charge", "tightId", "pfRelIso04_all",
                                          "pfRelIso03_chg", "pfRelIso03_all", "dxy", "dz", "isGlobal", "isTracker", "isStandalone")]
                 + [f"IsoTrack_{f}" for f in ("pt", "eta", "phi", "charge", "dxy", "dz", "pfRelIso03_chg",
                                              "pfRelIso03_all")]
                 + [f"TrigObj_{f}" for f in ("pt", "eta", "phi", "id", "filterBits", "l1pt")])
BRANCHES_MC = (BRANCHES_DATA + ["genWeight", "Pileup_nTrueInt", "Muon_genPartFlav"]
               + [f"GenPart_{f}" for f in ("pt", "eta", "phi", "pdgId", "status", "statusFlags")])


def blank():
    npt, neta = len(PT_EDGES) - 1, len(ETA_EDGES) - 1
    out = {"n_files": 0, "n_events": 0, "n_probes": {e: 0 for e in EFFS}}
    for eff in EFFS:
        for tag in TAGS:
            for charge in ("os", "ss"):
                for pf in ("pass", "fail"):
                    for suffix in ("", "_w2", "_gen"):
                        out[f"{eff}_{tag}_{charge}_{pf}{suffix}"] = np.zeros((npt, neta, tnp.MASS_NBINS))
    # truth for the closure: generator-level prompt muons from Z (pT > 20, |eta| < 2.4) matched
    # to a stand-alone / an isolated track, and how many of those are global or tracker muons
    for eff in EFFS:
        for pf in ("pass", "tot"):
            out[f"truth_{eff}_{pf}"] = np.zeros((npt, neta))
    return out


def _cross_dr(eta1, phi1, eta2, phi2):
    a = ak.zip({"eta": eta1, "phi": phi1})
    b = ak.zip({"eta": eta2, "phi": phi2})
    x, y = ak.unzip(ak.cartesian([a, b], axis=1, nested=True))
    return objects.delta_r(x.eta, x.phi, y.eta, y.phi)


def _prompt_gen_muons(ev):
    g = ((abs(ev.GenPart_pdgId) == 13) & (ev.GenPart_status == 1) & ((ev.GenPart_statusFlags & 1) == 1)
         & (ev.GenPart_pt > 15) & (abs(ev.GenPart_eta) < 2.5))
    return ev.GenPart_eta[g], ev.GenPart_phi[g], ev.GenPart_pt[g]


def probes(ev):
    """Dict eff -> record array (events x probes) with pt, eta, phi, charge, passed, index of the
    muon it is (-1 for IsoTrack probes)."""
    glob_trk = ev.Muon_isGlobal | ev.Muon_isTracker
    muon_idx = ak.local_index(ev.Muon_pt, axis=1)
    kin = (ev.Muon_pt > 20) & (abs(ev.Muon_eta) < 2.4)
    out = {}

    # ---- trk: stand-alone probes ------------------------------------------------------------
    sta = kin & ev.Muon_isStandalone
    dr_mm = _cross_dr(ev.Muon_eta, ev.Muon_phi, ev.Muon_eta, ev.Muon_phi)
    other = muon_idx[:, :, None] != muon_idx[:, None, :]
    # the partner index j runs over all muons: broadcast its flags along axis 2
    near = (dr_mm < STA_MATCH_DR) & other
    near_gt_sta = ak.any(near & (glob_trk & ev.Muon_isStandalone)[:, None, :], axis=2)
    near_gt_nosta = ak.any(near & (glob_trk & ~ev.Muon_isStandalone)[:, None, :], axis=2)
    sta_only = ~glob_trk
    keep = sta & ~(sta_only & near_gt_sta)
    passed = glob_trk | (sta_only & near_gt_nosta)
    out["trk"] = ak.zip({"pt": ev.Muon_pt[keep], "eta": ev.Muon_eta[keep], "phi": ev.Muon_phi[keep],
                         "charge": ev.Muon_charge[keep], "passed": passed[keep], "muidx": muon_idx[keep]})

    # ---- mutrk: isolated-track probes -------------------------------------------------------
    dr_tm = _cross_dr(ev.IsoTrack_eta, ev.IsoTrack_phi, ev.Muon_eta, ev.Muon_phi)
    dup = ak.any((dr_tm < TRK_MATCH_DR) & glob_trk[:, None, :], axis=2)
    loose_mu = (kin & glob_trk & (abs(ev.Muon_dxy) < config.MU_DXY_MAX) & (abs(ev.Muon_dz) < config.MU_DZ_MAX)
                & (ev.Muon_pfRelIso03_chg < TRK_ISO_MAX))
    loose_it = ((ev.IsoTrack_pt > 20) & (abs(ev.IsoTrack_eta) < 2.4) & (abs(ev.IsoTrack_dxy) < config.MU_DXY_MAX)
                & (abs(ev.IsoTrack_dz) < config.MU_DZ_MAX) & (ev.IsoTrack_pfRelIso03_chg < TRK_ISO_MAX) & ~dup)
    tight_mu = (loose_mu & (abs(ev.Muon_dxy) < TIGHT_DXY_MAX) & (abs(ev.Muon_dz) < TIGHT_DZ_MAX)
                & (ev.Muon_pfRelIso03_all < TIGHT_ISO_ALL_MAX))
    tight_it = (loose_it & (abs(ev.IsoTrack_dxy) < TIGHT_DXY_MAX) & (abs(ev.IsoTrack_dz) < TIGHT_DZ_MAX)
                & (ev.IsoTrack_pfRelIso03_all < TIGHT_ISO_ALL_MAX))
    for name, mu_sel, it in (("mutrk", loose_mu, loose_it), ("mutrkT", tight_mu, tight_it)):
        mu_part = ak.zip({"pt": ev.Muon_pt[mu_sel], "eta": ev.Muon_eta[mu_sel],
                          "phi": ev.Muon_phi[mu_sel], "charge": ev.Muon_charge[mu_sel],
                          "passed": ak.ones_like(ev.Muon_pt[mu_sel], dtype=bool),
                          "muidx": muon_idx[mu_sel]})
        trk_part = ak.zip({"pt": ev.IsoTrack_pt[it], "eta": ev.IsoTrack_eta[it], "phi": ev.IsoTrack_phi[it],
                           "charge": ev.IsoTrack_charge[it],
                           "passed": ak.zeros_like(ev.IsoTrack_pt[it], dtype=bool),
                           "muidx": ak.values_astype(ak.zeros_like(ev.IsoTrack_pt[it]) - 1, np.int64)})
        out[name] = ak.concatenate([mu_part, trk_part], axis=1)
    return out


def fill(ev, out, weight=None, is_mc=False, grl=None):
    """Fill the pass/fail spectra of one chunk of unskimmed NanoAOD."""
    keep = regions.fired(ev) & regions.event_clean(ev)
    if grl is not None:
        keep &= io.lumi_mask(ev.run, ev.luminosityBlock, grl)
    ev = ev[keep]
    if len(ev) == 0:
        return
    w_ev = np.ones(len(ev)) if weight is None else np.asarray(weight, dtype=float)[keep]
    out["n_events"] += len(ev)

    masks = regions.muon_masks(ev)
    matched = regions.trigger_match(ev)
    tight = masks["tight"]
    tags = {"nom": tight & (ev.Muon_pt > config.TAG_PT_MIN) & matched,
            "alt": tight & (ev.Muon_pt > TAG_ALT_PT) & (ev.Muon_pfRelIso04_all < TAG_ALT_ISO) & matched}
    tag_gen = (ev.Muon_genPartFlav == 1) if is_mc else None
    if is_mc:
        g_eta, g_phi, _ = _prompt_gen_muons(ev)

    muon_idx = ak.local_index(ev.Muon_pt, axis=1)
    met_ok = np.asarray(ak.to_numpy(ev.PuppiMET_pt)) < TIGHT_MET_MAX
    for eff, pr in probes(ev).items():
        if is_mc:
            dr_g = _cross_dr(pr.eta, pr.phi, g_eta, g_phi)
            probe_gen = ak.any(dr_g < GEN_DR[eff], axis=2)
        for tname, tmask in tags.items():
            t_idx = muon_idx[tmask]
            tag = ak.zip({"pt": ev.Muon_pt[tmask], "eta": ev.Muon_eta[tmask], "phi": ev.Muon_phi[tmask],
                          "mass": ev.Muon_mass[tmask], "charge": ev.Muon_charge[tmask], "idx": t_idx,
                          "gen": tag_gen[tmask] if is_mc else ak.ones_like(ev.Muon_pt[tmask], dtype=bool)})
            pairs = ak.cartesian({"t": tag, "p": ak.zip({"pt": pr.pt, "eta": pr.eta, "phi": pr.phi,
                                                         "charge": pr.charge, "passed": pr.passed,
                                                         "muidx": pr.muidx,
                                                         "gen": probe_gen if is_mc else pr.passed})},
                                 axis=1, nested=False)
            t, p = pairs.t, pairs.p
            dr = objects.delta_r(t.eta, t.phi, p.eta, p.phi)
            px1, py1, pz1, e1 = objects.p4(t.pt, t.eta, t.phi, t.mass)
            px2, py2, pz2, e2 = objects.p4(p.pt, p.eta, p.phi, 0.1056583745)
            m = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
            ok = (p.muidx != t.idx) & (dr > TAG_PROBE_DR) & (m > tnp.MASS_LO) & (m < tnp.MASS_HI)
            if eff == "mutrkT":
                ok = ok & ak.Array(met_ok)
            n_pairs = ak.to_numpy(ak.num(ok, axis=1))
            w_pair = np.repeat(w_ev, n_pairs)
            flat = lambda x: ak.to_numpy(ak.flatten(x))
            okf = flat(ok)
            q = flat(t.charge * p.charge)[okf]
            ppt, peta, pm, pw = flat(p.pt)[okf], flat(p.eta)[okf], flat(m)[okf], w_pair[okf]
            ppass = flat(p.passed)[okf]
            pgen = (flat(t.gen)[okf] & flat(p.gen)[okf]) if is_mc else np.zeros_like(ppass)
            for cname, sign in (("os", q < 0), ("ss", q > 0)):
                for pf, sel in (("pass", ppass), ("fail", ~ppass)):
                    s = sign & sel
                    key = f"{eff}_{tname}_{cname}_{pf}"
                    out[key] += tnp._cell_hist(ppt[s], peta[s], pm[s], pw[s], PT_EDGES)
                    out[key + "_w2"] += tnp._cell_hist(ppt[s], peta[s], pm[s], pw[s] ** 2, PT_EDGES)
                    if is_mc:
                        out[key + "_gen"] += tnp._cell_hist(ppt[s & pgen], peta[s & pgen], pm[s & pgen],
                                                            pw[s & pgen], PT_EDGES)
            if tname == "nom":
                out["n_probes"][eff] += int((q < 0).sum())

    if is_mc:
        _truth(ev, w_ev, out)


def _truth(ev, w_ev, out):
    """Generator-level closure: prompt muons (pT > 20, |eta| < 2.4) in events with a nominal tag."""
    g = ((abs(ev.GenPart_pdgId) == 13) & (ev.GenPart_status == 1) & ((ev.GenPart_statusFlags & 1) == 1)
         & (ev.GenPart_pt > 20) & (abs(ev.GenPart_eta) < 2.4))
    geta, gphi, gpt = ev.GenPart_eta[g], ev.GenPart_phi[g], ev.GenPart_pt[g]
    glob_trk = ev.Muon_isGlobal | ev.Muon_isTracker
    # trk: gen muons with a stand-alone track within 0.3; passed if a global/tracker muon within 0.1
    dr_gm = _cross_dr(geta, gphi, ev.Muon_eta, ev.Muon_phi)
    has_sta = ak.any((dr_gm < STA_MATCH_DR) & ev.Muon_isStandalone[:, None, :], axis=2)
    has_gt = ak.any((dr_gm < 0.1) & glob_trk[:, None, :], axis=2)
    # mutrk: gen muons with an isolated track (IsoTrack or the track of a global/tracker muon) within 0.1
    dr_gt = _cross_dr(geta, gphi, ev.IsoTrack_eta, ev.IsoTrack_phi)
    has_track = ak.any(dr_gt < 0.1, axis=2) | has_gt
    w = ak.broadcast_arrays(ak.Array(w_ev), gpt)[0]
    flat = lambda x: ak.to_numpy(ak.flatten(x))
    pt, eta, ww = np.clip(flat(gpt), PT_EDGES[0], PT_EDGES[-1] - 1e-3), np.abs(flat(geta)), flat(w)
    for eff in EFFS:
        den, num = {"trk": (flat(has_sta), flat(has_sta & has_gt)),
                    "mutrk": (flat(has_track), flat(has_track & has_gt))}[TRUTH_OF[eff]]
        for pf, sel in (("tot", den), ("pass", num)):
            h, _, _ = np.histogram2d(pt[sel], eta[sel], bins=[PT_EDGES, ETA_EDGES], weights=ww[sel])
            out[f"truth_{eff}_{pf}"] += h


# ----------------------------------------------------------------------------- fitting
def ss_template(ss_h, fallback=None, min_entries=30.0, smooth_bins=2.0):
    """Normalised background shape from the same-sign spectrum of the same cell.

    The failing probes' background (a tag plus a random isolated track or stand-alone muon) is sculpted by
    the kinematic cuts into a broad peak near the Z mass, which an exponential cannot describe; the
    same-sign pairs have the same kinematics and no Z. Smoothed with a Gaussian of `smooth_bins` bins; a
    cell with fewer than `min_entries` same-sign entries takes `fallback` (the eta- or cell-inclusive
    shape)."""
    from scipy.ndimage import gaussian_filter1d
    h = np.clip(np.asarray(ss_h, dtype=float), 0, None)
    if h.sum() < min_entries and fallback is not None:
        h = np.clip(np.asarray(fallback, dtype=float), 0, None)
    h = gaussian_filter1d(h, smooth_bins, mode="nearest") + 1e-4 * max(h.mean(), 1e-9)
    return h / h.sum()


def _eps_error(m, n_eff):
    """Uncertainty of a near-unit efficiency. HESSE is meaningless at the eps = 1 limit (it returns
    anything up to 1): use MINOS, its lower error when the upper one reaches the limit, and a binomial
    floor 1/n_eff when MINOS fails as well."""
    err = float(m.errors["eps"])
    eps = float(m.values["eps"])
    try:
        m.minos("eps")
        me = m.merrors["eps"]
        lo, hi = abs(me.lower), abs(me.upper)
        if me.is_valid or lo > 0:
            err = lo if (me.at_upper_limit or eps + hi >= 1.0 - 1e-6) else 0.5 * (lo + hi)
    except Exception:
        pass
    floor = 1.0 / max(n_eff, 1.0)
    if not np.isfinite(err) or err > max(10 * np.sqrt(max(eps * (1 - eps), floor) / max(n_eff, 1.0)), 0.01):
        err = np.sqrt(max(eps * (1 - eps), floor) / max(n_eff, 1.0))
    return max(err, floor)


def fit_cell_template(pass_h, fail_h, tmpl_pass, tmpl_fail, bkg_pass, bkg_fail, pass_var=None, fail_var=None,
                      mass_range=(tnp.MASS_LO, tnp.MASS_HI), fail_sigma_max=3.0, fail_shift_max=2.0):
    """Simultaneous pass/fail fit like `tnp.fit_cell`, with the background *shapes* fixed to same-sign
    templates (`ss_template`) and only their normalisations free.

    `fail_sigma_max`/`fail_shift_max` bound the extra smearing and shift of the failing-probe signal: a
    failing isolated track is measured by the tracker like a passing one (3 GeV / 2 GeV), a failing
    stand-alone muon is not (use ~15 / 8 GeV). Without the bound the failing "signal" can widen into the
    background and take part of it."""
    from iminuit import Minuit

    sel = (tnp.MASS_CENTRES > mass_range[0]) & (tnp.MASS_CENTRES < mass_range[1])
    x = tnp.MASS_CENTRES[sel]
    tp = np.clip(tmpl_pass, 0, None)[sel]
    tf = np.clip(tmpl_fail, 0, None)[sel]
    tp = tp / tp.sum() if tp.sum() > 0 else np.ones_like(tp) / len(tp)
    tf = tf / tf.sum() if tf.sum() > 0 else tp
    bp = np.asarray(bkg_pass)[sel]
    bf = np.asarray(bkg_fail)[sel]
    bp, bf = bp / bp.sum(), bf / bf.sum()
    n_p, n_f = pass_h[sel], fail_h[sel]
    k_p = (pass_var[sel].sum() / n_p.sum()) if (pass_var is not None and n_p.sum() > 0) else 1.0
    k_f = (fail_var[sel].sum() / n_f.sum()) if (fail_var is not None and n_f.sum() > 0) else 1.0
    k_p, k_f = max(k_p, 1e-9), max(k_f, 1e-9)

    def model(N, eps, dm_p, sg_p, dm_f, sg_f, Bp, Bf):
        mu_p = N * eps * tnp._shift_smear(tp, dm_p, sg_p, x) + Bp * bp
        mu_f = N * (1 - eps) * tnp._shift_smear(tf, dm_f, sg_f, x) + Bf * bf
        return mu_p, mu_f

    def nll(*args):
        mu_p, mu_f = model(*args)
        mu_p = np.clip(mu_p, 1e-9, None) / k_p
        mu_f = np.clip(mu_f, 1e-9, None) / k_f
        return float(np.sum(mu_p - (n_p / k_p) * np.log(mu_p)) + np.sum(mu_f - (n_f / k_f) * np.log(mu_f)))

    tot_p, tot_f = float(n_p.sum()), float(n_f.sum())
    eps0 = tot_p / max(tot_p + tot_f, 1.0)
    names = ["N", "eps", "dm_p", "sg_p", "dm_f", "sg_f", "Bp", "Bf"]
    start = [max(tot_p + tot_f, 1.0) * 0.99, min(max(eps0, 0.02), 0.9995), 0.0, 0.5, 0.0, min(1.5, fail_sigma_max),
             0.002 * max(tot_p, 1.0), 0.5 * max(tot_f, 1.0)]
    m = Minuit(nll, *start, name=names)
    m.errordef = Minuit.LIKELIHOOD
    m.limits["N"] = (0, None)
    m.limits["eps"] = (0.0, 1.0)
    m.limits["dm_p"] = (-3.0, 3.0)
    m.limits["dm_f"] = (-fail_shift_max, fail_shift_max)
    m.limits["sg_p"] = (0.05, 4.0)
    m.limits["sg_f"] = (0.05, fail_sigma_max)
    m.limits["Bp"] = (0, None)
    m.limits["Bf"] = (0, None)
    m.strategy = 1
    m.migrad(ncall=20000)
    if not m.valid:
        m.simplex()
        m.migrad(ncall=20000)
    m.hesse()
    mu_p, mu_f = model(*m.values)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = float(np.nansum((n_p - mu_p) ** 2 / np.where(mu_p > 0, mu_p * k_p, 1))
                     + np.nansum((n_f - mu_f) ** 2 / np.where(mu_f > 0, mu_f * k_f, 1)))
    eps = m.values["eps"]
    err = _eps_error(m, n_p.sum() / k_p + n_f.sum() / k_f)
    sig_f = m.values["N"] * (1 - eps) * tnp._shift_smear(tf, m.values["dm_f"], m.values["sg_f"], x)
    return {"eps": float(eps), "err": float(err), "N": float(m.values["N"]), "valid": bool(m.valid),
            "chi2": chi2, "ndf": int(2 * sel.sum() - m.nfit), "params": dict(zip(names, m.values)),
            "model_pass": mu_p, "model_fail": mu_f, "x": x,
            "bkg_frac_fail": float(1 - sig_f.sum() / max(mu_f.sum(), 1e-9))}
