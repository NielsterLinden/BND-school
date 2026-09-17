"""Fake-factor method of the mu tau_h and e tau_h channels (docs/10-v4-plan.md section 5; CMS arXiv:1801.03535 6.1).

    FF_p(DM, N_jets, pT)   p = multijet (same-sign DR), W+jets (m_T > 70 DR), ttbar (simulation)
    w(AR event) = sum_p R_p(N_jets, m_T) x FF_p           R_p: fractions of the jet-fake processes in the AR
    N_fake(SR) = sum_AR data x w  -  sum_AR simulation with a genuine / lepton-faked tau_h x w

Corrections: the multijet FF is measured in same-sign events and extrapolated to opposite sign with C_OS/SS
from the anti-isolated-lepton sidebands; the W+jets FF is measured at high m_T and extrapolated to m_T < 40
with r_W(DM) from the W+jets simulation. Statistical uncertainties of every FF bin enter the template
variance per event (as in fakes.py); the corrections and the fractions are nuisance parameters.
"""

from __future__ import annotations

import numpy as np

from . import config

PT_BINS = np.asarray(config.LTAU_FF_PT_BINS, dtype=float)
NJ_BINS = list(config.LTAU_FF_NJET_BINS)
MT_BINS = np.asarray(config.LTAU_FF_MT_BINS, dtype=float)
DMS = list(config.TAU_DMS)
DM_GROUPS = {"1prong": (0, 1), "3prong": (10, 11)}
TT_PT_BINS = np.asarray([30.0, 40.0, 50.0, 70.0, 1000.0])
SHAPE = (len(DMS), len(NJ_BINS), len(PT_BINS) - 1)
PROCESSES = ("qcd", "w", "tt")


def _idx(pt, dm, nj, pt_bins=PT_BINS):
    ipt = np.clip(np.searchsorted(pt_bins, pt, side="right") - 1, 0, len(pt_bins) - 2)
    inj = np.clip(np.searchsorted(np.asarray(NJ_BINS), nj, side="right") - 1, 0, len(NJ_BINS) - 1)
    idm = np.full(len(dm), -1)
    for i, d in enumerate(DMS):
        idm[dm == d] = i
    return idm, inj, ipt


def _fill(arr, mask, w, kin=None):
    kin = kin or arr
    h = np.zeros(SHAPE)
    idm, inj, ipt = _idx(kin["t_pt"][mask], arr["t_dm"][mask], kin["njets"][mask])
    ok = idm >= 0
    np.add.at(h, (idm[ok], inj[ok], ipt[ok]), np.asarray(w)[mask][ok])
    return h


def measure(data, num_mask, den_mask, subtract=()):
    """FF table (DM, N_jets, pT) = (N_T - MC_T) / (N_L - MC_L); subtract: [(arr, num_mask, den_mask, weights)]."""
    one = np.ones(len(data["t_pt"]))
    num, den = _fill(data, num_mask, one), _fill(data, den_mask, one)
    num_mc, den_mc = np.zeros(SHAPE), np.zeros(SHAPE)
    for arr, nm, dm, w in subtract:
        num_mc += _fill(arr, nm, w)
        den_mc += _fill(arr, dm, w)
    n = np.maximum(num - num_mc, 0.0)
    dd = np.maximum(den - den_mc, 1e-9)
    ff = n / dd
    err = ff * np.sqrt(1.0 / np.maximum(num, 1) + 1.0 / np.maximum(den, 1))
    err = np.where(num > 0, err, 1.0 / np.maximum(dd, 1.0))
    # a bin whose denominator is (almost) entirely subtracted simulation, or empty, carries no information
    bad = (den - den_mc < 1.0) | (den < 2)
    ff = np.where(bad, 0.0, ff)
    err = np.where(bad, 0.5, err)
    return {"ff": ff, "err": err, "num": num, "den": den, "num_mc": num_mc, "den_mc": den_mc,
            "pt_bins": PT_BINS.tolist(), "njet_bins": NJ_BINS, "dms": DMS}


def measure_mc(arrs_masks_weights, num=True):
    """FF table from weighted simulation only: [(arr, T mask, L mask, weights)] -> (DM, N_jets, pT)."""
    num_h, den_h, num2, den2 = np.zeros(SHAPE), np.zeros(SHAPE), np.zeros(SHAPE), np.zeros(SHAPE)
    for arr, tm, lm, w in arrs_masks_weights:
        num_h += _fill(arr, tm, w)
        den_h += _fill(arr, lm, w)
        num2 += _fill(arr, tm, np.asarray(w) ** 2)
        den2 += _fill(arr, lm, np.asarray(w) ** 2)
    ff = np.where(den_h > 0, num_h / np.maximum(den_h, 1e-9), 0.0)
    rel = np.sqrt(np.where(num_h > 0, num2 / np.maximum(num_h, 1e-9) ** 2, 1.0) + np.where(den_h > 0, den2 / np.maximum(den_h, 1e-9) ** 2, 1.0))
    return {"ff": ff, "err": ff * rel, "num": num_h, "den": den_h, "pt_bins": PT_BINS.tolist(), "njet_bins": NJ_BINS, "dms": DMS}


def coarsen_tt(table):
    """ttbar FF: merge the decay modes into 1-prong / 3-prong and the jet bins (few simulated events)."""
    ff, num, den = np.asarray(table["ff"]), np.asarray(table["num"]), np.asarray(table["den"])
    out_ff, out_err = np.zeros(SHAPE), np.zeros(SHAPE)
    for dms in DM_GROUPS.values():
        sel = [DMS.index(d) for d in dms]
        n = num[sel].sum(axis=(0, 1))
        dd = den[sel].sum(axis=(0, 1))
        f = np.where(dd > 0, n / np.maximum(dd, 1e-9), 0.0)
        e = f * np.sqrt(1.0 / np.maximum(n, 1e-9) + 1.0 / np.maximum(dd, 1e-9)) * 0 + f * 0.3   # 30% (paper) on the ttbar FF
        for i in sel:
            out_ff[i, :, :] = f[None, :]
            out_err[i, :, :] = e[None, :]
    return {**table, "ff": out_ff, "err": out_err}


def evaluate(table, arr, mask, kin=None, shift: int = 0):
    """(FF, relative error) per event in `mask`; shift = +-1 moves every bin by its error (stat)."""
    kin = kin or arr
    idm, inj, ipt = _idx(kin["t_pt"][mask], arr["t_dm"][mask], kin["njets"][mask])
    ok = idm >= 0
    vals = np.asarray(table["ff"], dtype=float)
    errs = np.asarray(table["err"], dtype=float)
    if shift:
        vals = np.maximum(vals + shift * errs, 0.0)
    ff = np.zeros(int(mask.sum()))
    rel = np.zeros(int(mask.sum()))
    with np.errstate(divide="ignore", invalid="ignore"):
        relerr = np.where(vals > 0, errs / vals, 0.0)
    ff[ok] = vals[idm[ok], inj[ok], ipt[ok]]
    rel[ok] = relerr[idm[ok], inj[ok], ipt[ok]]
    return ff, rel


# ------------------------------------------------------------------------------ fractions
def _frac_idx(arr, mask, kin=None):
    kin = kin or arr
    inj = np.clip(np.searchsorted(np.asarray(NJ_BINS), kin["njets"][mask], side="right") - 1, 0, len(NJ_BINS) - 1)
    imt = np.clip(np.searchsorted(MT_BINS, kin["mt_1"][mask], side="right") - 1, 0, len(MT_BINS) - 2)
    return inj, imt


def fractions(data, ar_mask, mc_parts):
    """R_p(N_jets, m_T) in the AR. mc_parts: [(arr, mask_all_origins, mask_jet_fake, weights, process)] with
    process 'w' (W+jets and Drell-Yan jet fakes), 'tt' (ttbar, single top) or 'other' (genuine / lepton-faked
    tau_h, any sample; subtracted from the data). Multijet = data - all simulation."""
    shape = (len(NJ_BINS), len(MT_BINS) - 1)
    n_data = np.zeros(shape)
    inj, imt = _frac_idx(data, ar_mask)
    np.add.at(n_data, (inj, imt), 1.0)
    n = {"w": np.zeros(shape), "tt": np.zeros(shape), "other": np.zeros(shape)}
    for arr, m_all, m_jet, w, proc in mc_parts:
        for m, key in ((m_all & m_jet, proc), (m_all & ~m_jet, "other")):
            if key == "other" and proc == "other":
                m = m_all
            inj, imt = _frac_idx(arr, m)
            np.add.at(n[key], (inj, imt), np.asarray(w)[m])
    qcd = np.maximum(n_data - n["w"] - n["tt"] - n["other"], 0.0)
    tot = np.maximum(qcd + n["w"] + n["tt"], 1e-9)
    f = {"qcd": qcd / tot, "w": n["w"] / tot, "tt": n["tt"] / tot}
    return {"fractions": {k: v.tolist() for k, v in f.items()}, "counts": {"data": n_data.tolist(), "w": n["w"].tolist(),
            "tt": n["tt"].tolist(), "other": n["other"].tolist(), "qcd": qcd.tolist()},
            "njet_bins": NJ_BINS, "mt_bins": MT_BINS.tolist()}


def fraction_per_event(fr, arr, mask, kin=None, shift_w: float = 0.0):
    """(f_qcd, f_w, f_tt) per event; shift_w moves f_w by the relative amount (compensated by f_qcd)."""
    inj, imt = _frac_idx(arr, mask, kin)
    f = {p: np.asarray(fr["fractions"][p])[inj, imt] for p in PROCESSES}
    if shift_w:
        d = np.minimum(f["w"] * shift_w, f["qcd"]) if shift_w > 0 else np.maximum(f["w"] * shift_w, -f["w"])
        f["w"] = f["w"] + d
        f["qcd"] = f["qcd"] - d
    return f["qcd"], f["w"], f["tt"]


# ------------------------------------------------------------------------------ corrections
def osss_correction(data, reg, subtract=()):
    """C_OS/SS of the multijet FF from the anti-isolated-lepton sidebands: FF(AI, OS) / FF(AI, SS), inclusive,
    both with the simulation of genuine / lepton-faked tau_h subtracted; subtract: [(arr, regions, weights)]."""
    out = {}
    for tag in ("OS", "SS"):
        t, lo = float(reg[f"AI_{tag}_T"].sum()), float(reg[f"AI_{tag}_L"].sum())
        t_mc = sum(float(w[r[f"AI_{tag}_T"]].sum()) for _, r, w in subtract)
        l_mc = sum(float(w[r[f"AI_{tag}_L"]].sum()) for _, r, w in subtract)
        ff = (t - t_mc) / max(lo - l_mc, 1e-9)
        out[tag] = {"ff": ff, "T": t, "L": lo, "T_mc": t_mc, "L_mc": l_mc, "stat": ff * np.sqrt(1 / max(t, 1) + 1 / max(lo, 1))}
    if out["OS"]["ff"] <= 0 or out["SS"]["ff"] <= 0 or out["OS"]["L"] - out["OS"]["L_mc"] < 10 or out["SS"]["L"] - out["SS"]["L_mc"] < 10:
        return {"C": 1.0, "stat": 0.5, "syst": config.LTAU_FF_OSSS_SYST, "rel_unc": 0.5, "sidebands": out, "note": "no measurement (statistics)"}
    c = out["OS"]["ff"] / out["SS"]["ff"]
    stat = c * np.hypot(out["OS"]["stat"] / out["OS"]["ff"], out["SS"]["stat"] / out["SS"]["ff"])
    return {"C": c, "stat": float(stat), "syst": config.LTAU_FF_OSSS_SYST, "rel_unc": float(np.hypot(stat / c, config.LTAU_FF_OSSS_SYST)), "sidebands": out}


def w_mt_correction(w_parts):
    """r_W per decay-mode group from the W+jets simulation (jet-faked tau_h): FF(OS, m_T < 40) / FF(OS, m_T > 70).
    w_parts: [(arr, regions with is_mc=False, jet-fake mask, weights)]."""
    out = {}
    for group, dms in DM_GROUPS.items():
        n = {k: 0.0 for k in ("lo_T", "lo_L", "hi_T", "hi_L")}
        n2 = {k: 0.0 for k in n}
        for arr, r, jet, w in w_parts:
            dm_sel = np.isin(arr["t_dm"], dms) & jet
            for k, m in (("lo_T", r["SR"]), ("lo_L", r["AR"]), ("hi_T", r["W_T"]), ("hi_L", r["W_L"])):
                ww = np.asarray(w)[m & dm_sel]
                n[k] += float(ww.sum())
                n2[k] += float((ww ** 2).sum())
        ff_lo = n["lo_T"] / max(n["lo_L"], 1e-9)
        ff_hi = n["hi_T"] / max(n["hi_L"], 1e-9)
        r = ff_lo / max(ff_hi, 1e-9) if ff_hi > 0 else 1.0
        rel = np.sqrt(sum(n2[k] / max(n[k], 1e-9) ** 2 for k in n))
        out[group] = {"r": float(r), "rel_stat": float(rel), "ff_lo": ff_lo, "ff_hi": ff_hi, "counts": n}
    return out


def r_w_per_event(rw, arr, mask, shift: int = 0):
    dm = arr["t_dm"][mask]
    r = np.ones(len(dm))
    for group, dms in DM_GROUPS.items():
        v = rw[group]["r"] * (1 + shift * np.hypot(rw[group]["rel_stat"], 0.5 * abs(rw[group]["r"] - 1)))
        r[np.isin(dm, dms)] = v
    return r


# ------------------------------------------------------------------------------ application
def fake_weights(arr, mask, tables, fr, osss, rw, kin=None, variation: str | None = None, direction: int = 0, with_err=False):
    """Per-event fake weight w = sum_p R_p FF_p (corrected) for the events in `mask` (zeros elsewhere).

    variation: None | 'osss' (C +- rel_unc) | 'wmt' (r_W shifted) | 'frac' (f_W +- 20% relative) | 'tt' (tt FF +- 30%)
               | 'qcdstat' / 'wstat' (every bin of that table by its statistical error).
    with_err: also the weight for the sum of squares, w sqrt(1 + (delta w_stat / w)^2)."""
    n = len(arr["t_pt"])
    w = np.zeros(n)
    w2 = np.zeros(n)
    if not mask.any():
        return (w, w2) if with_err else w
    f_qcd, f_w, f_tt = fraction_per_event(fr, arr, mask, kin, shift_w=0.2 * direction if variation == "frac" else 0.0)
    ff_qcd, rel_qcd = evaluate(tables["qcd"], arr, mask, kin, shift=direction if variation == "qcdstat" else 0)
    ff_w, rel_w = evaluate(tables["w"], arr, mask, kin, shift=direction if variation == "wstat" else 0)
    ff_tt, rel_tt = evaluate(tables["tt"], arr, mask, kin)
    c = osss["C"] * (1 + direction * osss["rel_unc"] if variation == "osss" else 1.0)
    r = r_w_per_event(rw, arr, mask, shift=direction if variation == "wmt" else 0)
    tt_scale = 1 + 0.3 * direction if variation == "tt" else 1.0
    parts = (f_qcd * c * ff_qcd, f_w * r * ff_w, f_tt * ff_tt * tt_scale)
    tot = parts[0] + parts[1] + parts[2]
    var = (parts[0] * rel_qcd) ** 2 + (parts[1] * rel_w) ** 2
    w[mask] = tot
    w2[mask] = np.sqrt(tot ** 2 + var)
    return (w, w2) if with_err else w


def to_json(obj):
    if isinstance(obj, dict):
        return {k: to_json(v) for k, v in obj.items()}
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    return obj


def from_json(tables):
    return {k: {kk: (np.asarray(vv) if kk in ("ff", "err", "num", "den", "num_mc", "den_mc") else vv) for kk, vv in t.items()}
            for k, t in tables.items()}
