"""Fake-factor method for jet -> tau_h fakes (docs/05-fake-factors.md).

The classic ("old-time normal") fake factor, applied to the leading tau:

    FF(era, DM, N_jets, pT) = N(SS, tau1 T, tau2 T) / N(SS, tau1 L, tau2 T)       measured in data
    N_fake(SR)              = sum over data events in AR of C_OS/SS(era, N_jets) x FF(era, DM_1, N_jets, pT_1)

where T = Medium, L = VVVLoose-and-not-Medium, SS = same-sign pair (QCD multijet enriched, no Z).
Requiring tau2 T in both the measurement and the application makes the FF conditional on an isolated
second tau, which absorbs the correlation of the two taus' isolation in dijet events. The jet
multiplicity binning (0 / 1 / >= 2 jets) is needed for closure: quark- and gluon-initiated jets fake
tau_h at different rates and the jet multiplicity changes that mixture.

Everything is binned per data-taking era as well (Run2016G / Run2016H trigger paths differ).

C_OS/SS corrects the extrapolation from same-sign to opposite-sign pairs. It is measured in the tau2
anti-isolated sideband (tau2 L), where both charges are fake-dominated: FF measured in SS-AI and
applied to OS-AI_L predicts OS-AI_T; C = observed / predicted.

Genuine-tau contamination: the SS and AR regions contain some genuine taus (mostly Z -> tau tau with
a leading tau failing Medium in AR, 1.3% of it but ~6% of the signal). The nominal subtracts the
simulation with a genuine leading tau (weighted like the SR prediction) everywhere (`subtract`); the
unsubtracted variant is kept as a cross-check.

Closure corrections (`closure_corrections`): after the table, the FF still depends on |eta(tau1)| (+-15%,
the same shape in every era / DM / pT / N_jets bin: DeepTau's jet rejection is not eta-flat) and on
pT(tau2) (-7% at 60-100 GeV). Both are corrected multiplicatively by their same-sign obs/pred ratios.

Statistics: the statistical uncertainty of every FF bin is carried per event (`fake_weights(with_err)`)
into the sum of squares of the fake template, i.e. into the per-bin gamma parameters of the fit.
"""

from __future__ import annotations

import numpy as np

from . import config

PT_BINS = np.asarray(config.FF_PT_BINS)
NJ_BINS = list(config.FF_NJET_BINS)
DMS = list(config.FF_DMS)
N_ERA = 2 if config.FF_BY_ERA else 1
ERA_FRAC = (np.asarray(config.ERA_LUMI_PB) / np.sum(config.ERA_LUMI_PB)) if config.FF_BY_ERA else np.ones(1)
SHAPE = (N_ERA, len(DMS), len(NJ_BINS), len(PT_BINS) - 1)


def era_index(arr, mask):
    """Era index per selected event (0 = Run2016G, 1 = Run2016H); None for simulation (no `era` field)."""
    if "era" not in arr:
        return None
    return np.asarray(arr["era"][mask]).astype(int) if N_ERA > 1 else np.zeros(int(np.sum(mask)), dtype=int)


def _bin_index(pt, dm, nj):
    ipt = np.clip(np.searchsorted(PT_BINS, pt, side="right") - 1, 0, len(PT_BINS) - 2)
    inj = np.clip(np.searchsorted(np.asarray(NJ_BINS), nj, side="right") - 1, 0, len(NJ_BINS) - 1)
    idm = np.full(len(dm), -1)
    for i, d in enumerate(DMS):
        idm[dm == d] = i
    return idm, inj, ipt


def _fill(arr, mask, w):
    """Weighted counts per (era, DM, N_jet, pT). Simulation is shared out over the eras by luminosity."""
    h = np.zeros(SHAPE)
    idm, inj, ipt = _bin_index(arr["t1_pt"][mask], arr["t1_dm"][mask], arr["njets"][mask])
    ok = idm >= 0
    ww = np.asarray(w)[mask][ok]
    era = era_index(arr, mask)
    if era is None:
        for e in range(N_ERA):
            np.add.at(h, (np.full(int(ok.sum()), e), idm[ok], inj[ok], ipt[ok]), ww * ERA_FRAC[e])
    else:
        np.add.at(h, (era[ok], idm[ok], inj[ok], ipt[ok]), ww)
    return h


def measure(data, num_mask, den_mask, subtract=()):
    """FF table {'ff', 'err', 'num', 'den', 'num_mc', 'den_mc'} of shape (era, DM, N_jet, pT).

    subtract: iterable of (arrays, num_mask, den_mask, weights) of simulation removed from both.
    """
    one = np.ones(len(data["t1_pt"]))
    num = _fill(data, num_mask, one)
    den = _fill(data, den_mask, one)
    num_mc = np.zeros(SHAPE)
    den_mc = np.zeros(SHAPE)
    for arr, nm, dmask, w in subtract:
        num_mc += _fill(arr, nm, w)
        den_mc += _fill(arr, dmask, w)
    n = np.maximum(num - num_mc, 0.0)
    dd = np.maximum(den - den_mc, 1e-9)
    ff = n / dd
    err = ff * np.sqrt(1.0 / np.maximum(num, 1) + 1.0 / np.maximum(den, 1))
    err = np.where(num > 0, err, 1.0 / dd)          # empty numerator: one-event scale
    return {"ff": ff, "err": err, "num": num, "den": den, "num_mc": num_mc, "den_mc": den_mc,
            "pt_bins": PT_BINS.tolist(), "njet_bins": NJ_BINS, "dms": DMS, "eras": ["G", "H"][:N_ERA]}


def evaluate(table, arr, mask, shift_dm: int | None = None, direction: int = 0, closure: bool = True,
             with_err: bool = False):
    """FF of the leading tau for the events in `mask` (array of length mask.sum()): per era for data, the
    luminosity-weighted era average for simulation. Optional +-1 sigma (stat.) shift of one decay mode.
    `closure`: multiply by the factorised closure corrections stored in the table (`closure_corrections`).
    `with_err`: also return the relative statistical uncertainty of the FF bin of every event."""
    idm, inj, ipt = _bin_index(arr["t1_pt"][mask], arr["t1_dm"][mask], arr["njets"][mask])
    ok = idm >= 0
    vals = np.asarray(table["ff"], dtype=float).copy()
    errs = np.asarray(table["err"], dtype=float)
    if shift_dm is not None and direction:
        i = DMS.index(shift_dm)
        vals[:, i] = np.maximum(vals[:, i] + direction * errs[:, i], 0.0)
    ff = np.zeros(int(np.sum(mask)))
    rel = np.zeros(int(np.sum(mask)))
    era = era_index(arr, mask)
    with np.errstate(divide="ignore", invalid="ignore"):
        relerr = np.where(vals > 0, errs / vals, 0.0)
    if era is None:
        avg = np.tensordot(ERA_FRAC, vals, axes=1)
        ff[ok] = avg[idm[ok], inj[ok], ipt[ok]]
        rel[ok] = np.tensordot(ERA_FRAC, relerr, axes=1)[idm[ok], inj[ok], ipt[ok]]
    else:
        ff[ok] = vals[era[ok], idm[ok], inj[ok], ipt[ok]]
        rel[ok] = relerr[era[ok], idm[ok], inj[ok], ipt[ok]]
    if closure and table.get("closure"):
        for var, corr in table["closure"].items():
            x = np.abs(arr["t1_eta"][mask]) if var == "eta" else arr["t2_pt"][mask]
            edges = np.asarray(corr["edges"], dtype=float)
            ib = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
            ff = ff * np.asarray(corr["values"], dtype=float)[ib]
    return (ff, rel) if with_err else ff


def closure_corrections(data, regions_data, table, subtract=()):
    """Factorised closure corrections f(|eta(tau1)|) then g(pT(tau2)), each obs/pred in same-sign events
    after the previous one (config.FF_CLOSURE_*_BINS). Stored in table["closure"] and applied by `evaluate`.
    subtract: iterable of (arrays, regions, weights) of simulation with a genuine leading tau."""
    table["closure"] = {}
    for var, edges in (("eta", config.FF_CLOSURE_ETA_BINS), ("pt2", config.FF_CLOSURE_PT2_BINS)):
        edges = np.asarray(edges, dtype=float)

        def xval(a, m):
            return np.abs(a["t1_eta"][m]) if var == "eta" else np.clip(a["t2_pt"][m], edges[0], edges[-1] - 1e-6)

        obs = np.histogram(xval(data, regions_data["SS_T"]), bins=edges)[0].astype(float)
        wp = evaluate(table, data, regions_data["SS_L"])
        pred = np.histogram(xval(data, regions_data["SS_L"]), bins=edges, weights=wp)[0]
        pred2 = np.histogram(xval(data, regions_data["SS_L"]), bins=edges, weights=wp ** 2)[0]
        for a, r, w in subtract:
            obs -= np.histogram(xval(a, r["SS_T"]), bins=edges, weights=w[r["SS_T"]])[0]
            pred -= np.histogram(xval(a, r["SS_L"]), bins=edges, weights=w[r["SS_L"]] * evaluate(table, a, r["SS_L"]))[0]
        ratio = np.where(pred > 0, obs / np.maximum(pred, 1e-9), 1.0)
        err = ratio * np.sqrt(1.0 / np.maximum(obs, 1.0) + pred2 / np.maximum(pred, 1e-9) ** 2)
        table["closure"][var] = {"edges": edges.tolist(), "values": ratio.tolist(), "err": err.tolist(),
                                 "obs": obs.tolist(), "pred": pred.tolist()}
    return table


def njet_category(arr, mask=None):
    nj = arr["njets"] if mask is None else arr["njets"][mask]
    return np.clip(np.searchsorted(np.asarray(NJ_BINS), nj, side="right") - 1, 0, len(NJ_BINS) - 1)


def per_event(values, arr, mask=None):
    """Look up a number, a per-N_jet array, an (era, N_jet) array or an (era, N_jet, BDT category) array for
    every event (in `mask`); the category comes from arr["_cat"] (set by the caller)."""
    v = np.asarray(values, dtype=float)
    n = len(arr["t1_pt"]) if mask is None else int(np.sum(mask))
    if v.ndim == 0:
        return np.full(n, float(v))
    m = np.ones(len(arr["t1_pt"]), dtype=bool) if mask is None else mask
    nj = njet_category(arr, m)
    if v.ndim == 1:
        return v[nj]
    era = era_index(arr, m)
    if v.ndim == 2:
        return np.tensordot(ERA_FRAC, v, axes=1)[nj] if era is None else v[era, nj]
    cat = np.asarray(arr["_cat"])[m]
    return np.tensordot(ERA_FRAC, v, axes=1)[nj, cat] if era is None else v[era, nj, cat]


def fake_weights(arr, mask, table, c_osss, shift_dm=None, direction=0, with_err=False):
    """Per-event fake weight FF x C for events in `mask` (zeros elsewhere); `c_osss` may be a number,
    one value per N_jet category or an (era, N_jet) array. `with_err`: also the weight to use for the
    sum of squares, w x sqrt(1 + relerr_FF^2), so the FF statistics enter the template variance bin by bin."""
    w = np.zeros(len(arr["t1_pt"]))
    w2 = np.zeros(len(arr["t1_pt"]))
    ff, rel = evaluate(table, arr, mask, shift_dm, direction, with_err=True)
    w[mask] = ff * per_event(c_osss, arr, mask)
    w2[mask] = w[mask] * np.sqrt(1.0 + rel ** 2)
    return (w, w2) if with_err else w


def _osss(data, regions_data, table_ai, subtract, sel=None):
    sel_d = np.ones(len(data["t1_pt"]), dtype=bool) if sel is None else sel(data, True)
    obs = float((regions_data["OSAI_T"] & sel_d).sum())
    ff = evaluate(table_ai, data, regions_data["OSAI_L"] & sel_d)
    pred, pred_var = float(ff.sum()), float((ff ** 2).sum())
    obs_mc = pred_mc = 0.0
    for a, r, w in subtract:
        s_a, frac = (np.ones(len(a["t1_pt"]), dtype=bool), 1.0) if sel is None else sel(a, False)
        obs_mc += frac * float(w[r["OSAI_T"] & s_a].sum())
        m = r["OSAI_L"] & s_a
        pred_mc += frac * float((w[m] * evaluate(table_ai, a, m)).sum())
    c = (obs - obs_mc) / (pred - pred_mc)
    stat = c * np.sqrt(1.0 / max(obs, 1) + pred_var / pred ** 2)
    return {"C": c, "stat": float(stat), "obs": obs, "pred": pred, "obs_mc": obs_mc, "pred_mc": pred_mc}


def osss_correction(data, regions_data, subtract=(), n_cat: int | None = None, max_rel_stat: float = 0.10):
    """C_OS/SS per (era, N_jet) [and BDT category] and inclusive, from the tau2 anti-isolated sideband: FF
    measured in SS-AI and applied to OS-AI_L predicts OS-AI_T; C = observed / predicted.
    subtract: iterable of (arrays, regions, weights) of simulation (shared out over eras by luminosity).
    n_cat: if given, C is also measured per BDT category, read from arr["_cat"] of the data and of every
    subtracted sample (the OS/SS charge correlation of the two jets depends on the topology the BDT selects:
    with the inclusive C the signal-like categories were under-predicted by 3-6%, docs/05).
    """
    sub_ai = [(a, r["SSAI_T"], r["SSAI_L"], w) for a, r, w in subtract]
    table_ai = measure(data, regions_data["SSAI_T"], regions_data["SSAI_L"], sub_ai)
    out = {"inclusive": _osss(data, regions_data, table_ai, subtract)}
    cats = [None] if n_cat is None else list(range(n_cat))
    shape = (N_ERA, len(NJ_BINS)) + (() if n_cat is None else (n_cat,))
    C = np.zeros(shape)
    stat = np.zeros_like(C)
    for e in range(N_ERA):
        for j in range(len(NJ_BINS)):
            for k in cats:
                def sel(a, is_data, e=e, j=j, k=k):
                    m = njet_category(a) == j
                    if k is not None:
                        m = m & (np.asarray(a["_cat"]) == k)
                    if is_data:
                        return m & ((a["era"] == e) if N_ERA > 1 else True)
                    return m, float(ERA_FRAC[e])
                res = _osss(data, regions_data, table_ai, subtract, sel)
                idx = (e, j) if k is None else (e, j, k)
                C[idx], stat[idx] = res["C"], res["stat"]
    if n_cat is not None:
        # a (era, N_jets, category) bin with a poor measurement (relative statistical error above
        # `max_rel_stat`, e.g. 0-jet events in the signal-like category) takes the N_jets-inclusive value of
        # its era and category instead
        C_incl = np.zeros((N_ERA, n_cat)); stat_incl = np.zeros_like(C_incl)
        for e in range(N_ERA):
            for k in range(n_cat):
                def sel(a, is_data, e=e, k=k):
                    m = np.asarray(a["_cat"]) == k
                    if is_data:
                        return m & ((a["era"] == e) if N_ERA > 1 else True)
                    return m, float(ERA_FRAC[e])
                res = _osss(data, regions_data, table_ai, subtract, sel)
                C_incl[e, k], stat_incl[e, k] = res["C"], res["stat"]
        poor = stat / np.maximum(C, 1e-9) > max_rel_stat
        out["replaced_by_njet_inclusive"] = [[int(e), int(j), int(k)] for e, j, k in zip(*np.where(poor))]
        for e, j, k in zip(*np.where(poor)):
            C[e, j, k], stat[e, j, k] = C_incl[e, k], stat_incl[e, k]
        out["C_njet_inclusive"] = C_incl.tolist()
        out["per_category"] = True
    out["C"] = C.tolist()
    out["stat"] = stat.tolist()
    out["table_ai"] = table_ai
    return out


def to_json(table):
    return {k: (np.asarray(v).tolist() if isinstance(v, np.ndarray) else v) for k, v in table.items()}


def from_json(js):
    return {k: (np.asarray(v) if k in ("ff", "err", "num", "den", "num_mc", "den_mc") else v) for k, v in js.items()}
