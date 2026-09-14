"""Fake-factor estimate of the non-prompt muon background (docs/13-fake-factor.md).

Muon classes: tight (signal) and anti-tight (tightId, 0.20 < iso < 1.0); alternative window
(0.15 < iso < 0.6) for the systematic. The fake factor FF = N_tight / N_anti-tight is measured
in a fake-enriched region with the prompt contribution subtracted using simulation (events in
which both muons are prompt: `genPartFlav` in {1, 15}):

  * nominal: same-sign tag + probe pairs (tag tight, pT > 26, trigger-matched; probe tight or
    anti-tight, pT > 20; any mass above 12 GeV) -- the probe did not fire the trigger, so its
    isolation is unbiased, and the fake composition (W+jets, QCD) is the one of the signal region;
  * alternative: single muon + jet events (skim category C: exactly one muon, MET < 30, mT < 30).

Application: opposite-sign events with exactly one tight (pT > 26, matched) and one anti-tight
muon, 60 < m < 120 GeV. N_fake(m) = sum FF(pT, |eta| of the anti-tight muon) over data minus the
prompt-prompt simulation (calibrated with the anti-isolation scale factor from tag-and-probe).
Closure: the same in the same-sign application region against the same-sign tight-tight yield.
"""

from __future__ import annotations

import json

import awkward as ak
import numpy as np

from . import config, objects, regions

PT_EDGES = np.array([20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 200.0])
ETA_EDGES = np.array([0.0, 0.9, 1.2, 2.1, 2.4])
MASS_EDGES = np.linspace(config.MASS_LO, config.MASS_HI, config.MASS_NBINS + 1)
MIN_PAIR_MASS = 12.0
ANTI_ALT = (0.15, 0.60)
REGIONS = ("ss", "c")            # FF measurement regions
CLASSES = ("tight", "anti", "anti_alt")


def blank():
    shape = (len(PT_EDGES) - 1, len(ETA_EDGES) - 1)
    out = {"n_files": 0}
    for reg in REGIONS:
        for cls in CLASSES:
            for suf in ("", "_w2"):
                out[f"ff_{reg}_{cls}{suf}"] = np.zeros(shape)
    for charge in ("os", "ss"):
        for cls in ("anti", "anti_alt"):
            for suf in ("", "_w2"):
                out[f"app_{charge}_{cls}{suf}"] = np.zeros(shape + (config.MASS_NBINS,))
        out[f"tt_{charge}"] = np.zeros(config.MASS_NBINS)      # tight-tight (SR / SS) mass, for closure
        out[f"tt_{charge}_w2"] = np.zeros(config.MASS_NBINS)
    return out


def _cell_hist2d(pt, eta, w):
    pt = np.clip(pt, PT_EDGES[0], PT_EDGES[-1] - 1e-3)
    h, _, _ = np.histogram2d(pt, np.abs(eta), bins=[PT_EDGES, ETA_EDGES], weights=w)
    return h


def _cell_hist3d(pt, eta, mass, w):
    pt = np.clip(pt, PT_EDGES[0], PT_EDGES[-1] - 1e-3)
    h, _ = np.histogramdd(np.stack([pt, np.abs(eta), mass], axis=1), bins=[PT_EDGES, ETA_EDGES, MASS_EDGES], weights=w)
    return h


def muon_classes(ev):
    masks = regions.muon_masks(ev)
    ip = (abs(ev.Muon_dxy) < config.MU_DXY_MAX) & (abs(ev.Muon_dz) < config.MU_DZ_MAX)
    acc = (ev.Muon_pt > regions.MU_PT_MIN) & (abs(ev.Muon_eta) < regions.MU_ETA) & ip & ev.Muon_tightId
    iso = ev.Muon_pfRelIso04_all
    masks["anti_alt"] = acc & (iso > ANTI_ALT[0]) & (iso < ANTI_ALT[1])
    return masks


def fill(ev, out, weight, is_mc, prompt_only=True):
    """Fill FF measurement and application histograms for one chunk (weight per event)."""
    keep = regions.fired(ev) & regions.event_clean(ev)
    ev, weight = ev[keep], np.asarray(weight, dtype=float)[keep]
    if len(ev) == 0:
        return
    masks = muon_classes(ev)
    matched = regions.trigger_match(ev)
    tight = masks["tight"]
    prompt = regions.is_prompt(ev.Muon_genPartFlav) if is_mc else ak.ones_like(ev.Muon_pt, dtype=bool)

    # --- same-sign tag + probe measurement region ------------------------------------
    is_tag = tight & (ev.Muon_pt > config.MU_PT_LEAD) & matched
    idx = ak.local_index(ev.Muon_pt, axis=1)
    p = ak.combinations(idx, 2, fields=["a", "b"])
    tag = ak.concatenate([p.a, p.b], axis=1)
    probe = ak.concatenate([p.b, p.a], axis=1)
    px1, py1, pz1, e1 = objects.p4(ev.Muon_pt[tag], ev.Muon_eta[tag], ev.Muon_phi[tag], ev.Muon_mass[tag])
    px2, py2, pz2, e2 = objects.p4(ev.Muon_pt[probe], ev.Muon_eta[probe], ev.Muon_phi[probe], ev.Muon_mass[probe])
    m = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
    same = (ev.Muon_charge[tag] * ev.Muon_charge[probe]) > 0
    n_pairs = ak.to_numpy(ak.num(tag, axis=1))
    w_pair = np.repeat(weight, n_pairs)
    base = is_tag[tag] & same & (m > MIN_PAIR_MASS)
    if is_mc and prompt_only:
        base = base & prompt[tag] & prompt[probe]
    for cls in CLASSES:
        valid = ak.to_numpy(ak.flatten(base & masks[cls][probe]))
        if not valid.any():
            continue
        ppt = ak.to_numpy(ak.flatten(ev.Muon_pt[probe]))[valid]
        peta = ak.to_numpy(ak.flatten(ev.Muon_eta[probe]))[valid]
        out[f"ff_ss_{cls}"] += _cell_hist2d(ppt, peta, w_pair[valid])
        out[f"ff_ss_{cls}_w2"] += _cell_hist2d(ppt, peta, w_pair[valid] ** 2)

    # --- single-muon + jet region (skim category C) ---------------------------------
    in_c = ak.to_numpy((ev.skim_cat & 4) > 0)
    n_loose = ak.to_numpy(ak.sum(masks["loose"], axis=1))
    sel_c = in_c & (n_loose == 1)
    if sel_c.any():
        evc, wc = ev[sel_c], weight[sel_c]
        mc_masks = {cls: masks[cls][sel_c] for cls in CLASSES}
        lead = masks["loose"][sel_c]
        matched_c = matched[sel_c]
        pr_c = prompt[sel_c]
        for cls in CLASSES:
            sel = lead & mc_masks[cls] & matched_c & (evc.Muon_pt > config.MU_PT_LEAD)
            if is_mc and prompt_only:
                sel = sel & pr_c
            flat = ak.to_numpy(ak.flatten(sel))
            if not flat.any():
                continue
            npm = ak.to_numpy(ak.num(evc.Muon_pt, axis=1))
            wm = np.repeat(wc, npm)[flat]
            ppt = ak.to_numpy(ak.flatten(evc.Muon_pt))[flat]
            peta = ak.to_numpy(ak.flatten(evc.Muon_eta))[flat]
            out[f"ff_c_{cls}"] += _cell_hist2d(ppt, peta, wm)
            out[f"ff_c_{cls}_w2"] += _cell_hist2d(ppt, peta, wm ** 2)

    # --- application region (OS and SS), both anti-tight definitions ------------------
    clean = np.ones(len(ev), dtype=bool)
    for cls in ("anti", "anti_alt"):
        mk = dict(masks); mk["anti"] = masks[cls]
        app = regions.ff_application(ev, mk, matched, clean)
        for charge in ("os", "ss"):
            d = app[f"FFapp_{charge.upper()}"]
            if d is None or len(d["idx"]) == 0:
                continue
            w = weight[d["idx"]]
            if is_mc and prompt_only:
                keep_p = regions.is_prompt(d["flav_anti"]) & regions.is_prompt(d["flav_tight"])
                w = np.where(keep_p, w, 0.0)
            out[f"app_{charge}_{cls}"] += _cell_hist3d(d["pt_anti"], d["eta_anti"], d["mass"], w)
            out[f"app_{charge}_{cls}_w2"] += _cell_hist3d(d["pt_anti"], d["eta_anti"], d["mass"], w ** 2)
    # tight-tight (for the SS closure and the SR comparison)
    dimu = regions.dimuon_regions(ev, masks, matched, clean)
    for charge, name in (("os", "SR"), ("ss", "SS")):
        d = dimu[name]
        if d is None or len(d["idx"]) == 0:
            continue
        w = weight[d["idx"]]
        if is_mc and prompt_only:
            w = np.where(regions.is_prompt(d["flav1"]) & regions.is_prompt(d["flav2"]), w, 0.0)
        out[f"tt_{charge}"] += np.histogram(d["mass"], bins=MASS_EDGES, weights=w)[0]
        out[f"tt_{charge}_w2"] += np.histogram(d["mass"], bins=MASS_EDGES, weights=w ** 2)[0]


# ----------------------------------------------------------------------------- estimate
def fake_factor(data, mc, region="ss", cls="anti", mc_scale=1.0):
    """FF map and its statistical error from data minus prompt MC, numerator tight / denominator anti."""
    num = data[f"ff_{region}_tight"] - mc_scale * mc[f"ff_{region}_tight"]
    den = data[f"ff_{region}_{cls}"] - mc_scale * mc[f"ff_{region}_{cls}"]
    num_var = data[f"ff_{region}_tight_w2"] + mc_scale ** 2 * mc[f"ff_{region}_tight_w2"]
    den_var = data[f"ff_{region}_{cls}_w2"] + mc_scale ** 2 * mc[f"ff_{region}_{cls}_w2"]
    with np.errstate(divide="ignore", invalid="ignore"):
        ff = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
        err = np.where((den > 0) & (num > 0), ff * np.sqrt(num_var / np.maximum(num, 1e-9) ** 2 + den_var / np.maximum(den, 1e-9) ** 2), ff)
    return np.clip(ff, 0, None), err


def ff_weighted(hists, ff, charge="os", cls="anti"):
    """FF-weighted mass spectrum and its variance from a (pT, |eta|, mass) application histogram."""
    h = hists[f"app_{charge}_{cls}"]
    v = hists[f"app_{charge}_{cls}_w2"]
    return np.einsum("ij,ijk->k", ff, h), np.einsum("ij,ijk->k", ff ** 2, v)


def fit_application(data, mc, ff, cls="anti", fake_shape="ss"):
    """Non-prompt yield in the opposite-sign application region from a two-template fit.

    D(m) = a P(m) + b F(m): P = FF-weighted prompt-prompt MC (Z peak), F = non-prompt shape
    (FF-weighted same-sign data minus prompt MC, or an exponential). The prompt
    normalisation `a` is fitted (replacing an assumed anti-isolation scale factor), so the
    estimate does not rely on the absolute MC rate of Z events with a non-isolated muon.
    Returns dict(n_fake, n_fake_err, a, b, template, template_var, chi2).
    """
    from scipy.optimize import nnls
    d, d_var = ff_weighted(data, ff, "os", cls)
    p, _ = ff_weighted(mc, ff, "os", cls)
    if fake_shape == "ss":
        f_ss, f_ss_var = ff_weighted(data, ff, "ss", cls)
        pm, _ = ff_weighted(mc, ff, "ss", cls)
        f = np.clip(f_ss - pm, 0, None)
    else:
        c = 0.5 * (MASS_EDGES[1:] + MASS_EDGES[:-1])
        f = np.exp(-0.02 * (c - 60.0))
    # 2 GeV bins for the fit (less noise in F), weighted least squares with data variances
    k = 4
    D, P, F = d.reshape(-1, k).sum(1), p.reshape(-1, k).sum(1), f.reshape(-1, k).sum(1)
    sig = np.sqrt(np.maximum(d_var.reshape(-1, k).sum(1), 1.0))
    A = np.stack([P / sig, F / sig], axis=1)
    (a, b), _ = nnls(A, D / sig)
    resid = D - a * P - b * F
    chi2 = float(np.sum((resid / sig) ** 2))
    # error on b from the normal equations
    cov = np.linalg.pinv(A.T @ A)
    b_err = float(np.sqrt(max(cov[1, 1], 0.0)))
    template = b * f
    return {"n_fake": float(template.sum()), "n_fake_err": float(b_err * f.sum()), "a": float(a), "b": float(b),
            "template": template, "template_var": (b_err * f) ** 2, "chi2": chi2, "ndf": int(len(D) - 2),
            "n_prompt_fitted": float(a * p.sum()), "n_prompt_mc": float(p.sum())}


def summarise(data, mc, antiiso_sf=None):
    """All fake-factor maps, the SR template with its variations, and closure numbers."""
    res = {"pt_edges": PT_EDGES.tolist(), "eta_edges": ETA_EDGES.tolist(), "maps": {}, "templates": {}, "fits": {}}
    ff, ff_err = fake_factor(data, mc, "ss", "anti")
    res["maps"]["nominal"] = ff.tolist(); res["maps"]["nominal_err"] = ff_err.tolist()
    variants = {
        "ff_stat_up": (ff + ff_err, "anti", "ss"),
        "ff_stat_down": (np.clip(ff - ff_err, 0, None), "anti", "ss"),
        "region_c": (fake_factor(data, mc, "c", "anti")[0], "anti", "ss"),
        "anti_alt": (fake_factor(data, mc, "ss", "anti_alt")[0], "anti_alt", "ss"),
        "mcsub_up": (fake_factor(data, mc, "ss", "anti", 1.3)[0], "anti", "ss"),
        "mcsub_down": (fake_factor(data, mc, "ss", "anti", 0.7)[0], "anti", "ss"),
        "shape_expo": (ff, "anti", "expo"),
    }
    res["maps"]["region_c"] = variants["region_c"][0].tolist()
    res["maps"]["anti_alt"] = variants["anti_alt"][0].tolist()
    nom = fit_application(data, mc, ff, "anti", "ss")
    res["templates"]["nominal"] = nom["template"].tolist()
    res["templates"]["nominal_var"] = nom["template_var"].tolist()
    res["fits"]["nominal"] = {k: v for k, v in nom.items() if not isinstance(v, np.ndarray)}
    for name, (f, cls, shape) in variants.items():
        r = fit_application(data, mc, f, cls, shape)
        res["templates"][name] = r["template"].tolist()
        res["fits"][name] = {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}
    # plain subtraction with the anti-iso SF, for the record
    d, _ = ff_weighted(data, ff, "os", "anti"); p, _ = ff_weighted(mc, ff, "os", "anti")
    sf = antiiso_sf if antiiso_sf is not None else np.ones_like(ff)
    p_sf, _ = ff_weighted({k: (v * sf[:, :, None] if k.startswith("app_os") else v) for k, v in mc.items()}, ff, "os", "anti")
    res["fits"]["plain_subtraction"] = {"n_fake": float((d - p_sf).sum()), "n_prompt_mc_sf": float(p_sf.sum()), "n_data": float(d.sum())}
    # closure in the same-sign application region (plain subtraction: prompt is small there)
    pred, pred_var = ff_weighted(data, ff, "ss", "anti")
    pm, _ = ff_weighted(mc, ff, "ss", "anti")
    pred = pred - pm
    obs = data["tt_ss"] - mc["tt_ss"]
    res["closure"] = {"predicted_ss": float(pred.sum()), "predicted_ss_err": float(np.sqrt(pred_var.sum())),
                      "observed_ss": float(obs.sum()), "observed_ss_err": float(np.sqrt(data["tt_ss_w2"].sum() + mc["tt_ss_w2"].sum())),
                      "predicted_hist": pred.tolist(), "observed_hist": obs.tolist()}
    n_nom = nom["n_fake"]
    res["yields"] = {"sr_fakes": float(n_nom), "sr_fakes_stat": float(nom["n_fake_err"]),
                     "sr_data_tight_tight": float(data["tt_os"].sum()), "ss_data_tight_tight": float(data["tt_ss"].sum()),
                     "variants": {k: float(np.sum(v)) for k, v in res["templates"].items() if not k.endswith("_var")},
                     "prompt_norm_fitted": nom["a"], "os_over_ss_fakes": float(n_nom / max(res["closure"]["predicted_ss"], 1e-9))}
    # region_c (single-muon + jet) is kept as a documented cross-check only: the IsoMu24 trigger
    # depletes its anti-isolated denominator, so its FF is biased high by an order of magnitude.
    devs = {k: abs(res["yields"]["variants"][k] - n_nom) for k in ("anti_alt", "mcsub_up", "mcsub_down", "shape_expo")}
    closure_dev = abs(res["closure"]["observed_ss"] - res["closure"]["predicted_ss"]) / max(res["closure"]["predicted_ss"], 1e-9) * abs(n_nom)
    devs["closure"] = float(closure_dev)
    res["yields"]["method_deviations"] = devs
    worst = max(devs, key=devs.get)
    res["yields"]["method_rel_unc"] = float(devs[worst] / max(abs(n_nom), 1e-9))
    res["yields"]["method_worst"] = worst
    return res
