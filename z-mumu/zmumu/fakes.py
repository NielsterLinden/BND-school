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


def apply_ff(data, mc, ff, charge="os", cls="anti", mc_scale=1.0, antiiso_sf=None):
    """Fake mass template: sum over cells of FF x (data - prompt MC) in the application region."""
    d = data[f"app_{charge}_{cls}"]
    m = mc[f"app_{charge}_{cls}"] * mc_scale
    if antiiso_sf is not None:
        m = m * antiiso_sf[:, :, None]
    diff = d - m
    template = np.einsum("ij,ijk->k", ff, diff)
    var = np.einsum("ij,ijk->k", ff ** 2, data[f"app_{charge}_{cls}_w2"] + mc_scale ** 2 * mc[f"app_{charge}_{cls}_w2"])
    return template, var, diff


def summarise(data, mc, antiiso_sf=None):
    """All fake-factor maps, the SR template with its variations, and closure numbers."""
    res = {"pt_edges": PT_EDGES.tolist(), "eta_edges": ETA_EDGES.tolist(), "maps": {}, "templates": {}}
    ff, ff_err = fake_factor(data, mc, "ss", "anti")
    res["maps"]["nominal"] = ff.tolist(); res["maps"]["nominal_err"] = ff_err.tolist()
    variants = {
        "ff_stat_up": (ff + ff_err, "anti", 1.0),
        "ff_stat_down": (np.clip(ff - ff_err, 0, None), "anti", 1.0),
        "region_c": (fake_factor(data, mc, "c", "anti")[0], "anti", 1.0),
        "anti_alt": (fake_factor(data, mc, "ss", "anti_alt")[0], "anti_alt", 1.0),
        "mcsub_up": (fake_factor(data, mc, "ss", "anti", 1.3)[0], "anti", 1.3),
        "mcsub_down": (fake_factor(data, mc, "ss", "anti", 0.7)[0], "anti", 0.7),
    }
    res["maps"]["region_c"] = variants["region_c"][0].tolist()
    res["maps"]["anti_alt"] = variants["anti_alt"][0].tolist()
    nom, nom_var, _ = apply_ff(data, mc, ff, "os", "anti", 1.0, antiiso_sf)
    res["templates"]["nominal"] = nom.tolist(); res["templates"]["nominal_var"] = nom_var.tolist()
    for name, (f, cls, scale) in variants.items():
        t, _, _ = apply_ff(data, mc, f, "os", cls, scale, antiiso_sf)
        res["templates"][name] = t.tolist()
    # closure in the same-sign application region
    pred, pred_var, _ = apply_ff(data, mc, ff, "ss", "anti", 1.0, antiiso_sf)
    obs = data["tt_ss"] - mc["tt_ss"]
    win = slice(0, None)
    res["closure"] = {"predicted_ss": float(pred.sum()), "predicted_ss_err": float(np.sqrt(pred_var.sum())),
                      "observed_ss": float(obs.sum()), "observed_ss_err": float(np.sqrt(data["tt_ss_w2"].sum() + mc["tt_ss_w2"].sum())),
                      "predicted_hist": pred.tolist(), "observed_hist": obs.tolist()}
    n_nom = nom.sum()
    res["yields"] = {"sr_fakes": float(n_nom), "sr_fakes_stat": float(np.sqrt(nom_var.sum())),
                     "sr_data_tight_tight": float(data["tt_os"].sum()), "ss_data_tight_tight": float(data["tt_ss"].sum()),
                     "variants": {k: float(np.sum(v)) for k, v in res["templates"].items() if not k.endswith("_var")}}
    # FakeMethod = largest deviation among the method variants (region, window, MC subtraction, closure)
    devs = {k: abs(res["yields"]["variants"][k] - n_nom) for k in ("region_c", "anti_alt", "mcsub_up", "mcsub_down")}
    closure_dev = abs(res["closure"]["observed_ss"] - res["closure"]["predicted_ss"]) / max(res["closure"]["predicted_ss"], 1e-9) * n_nom
    devs["closure"] = float(closure_dev)
    res["yields"]["method_deviations"] = devs
    worst = max(devs, key=devs.get)
    res["yields"]["method_rel_unc"] = float(devs[worst] / max(n_nom, 1e-9))
    res["yields"]["method_worst"] = worst
    return res
