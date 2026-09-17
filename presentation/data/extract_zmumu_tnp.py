#!/usr/bin/env python
"""Freeze the Z -> mumu tag-and-probe inputs of the section-4 T&P clips (mumu_t1 ... mumu_t5).

    source setup.sh && python presentation/data/extract_zmumu_tnp.py [--json PATH] [--check-only]

Reads (read-only): z-mumu/output/v2/tnp/tnp_fits.pkl (per-cell pass/fail m(mumu) spectra of data and simulation and the
nominal fit of every cell), tnp_result.json, and the first data skim file of $BND_SKIM_DIR/data_2016G (two display events).
Writes presentation/data/zmumu_tnp.json:

    cell     the tight-ID cell with the most data probes: pass / fail spectra (0.5 GeV, 60-120 GeV) of data and simulation,
             the nominal simultaneous fit (model = signal template (x) Gauss + exponential background, both components),
             eps +- err, N (docs/12-tag-and-probe-fits.md)
    id       eps_data, eps_mc (+ errors) and SF = eps_data / eps_mc maps over (pT, |eta|); map mean (docs/12: 0.980)
    iso      SF map and map mean (docs/12: 1.005); trigger per-muon plateaus (docs/12: 0.907 / 0.923)
    events   two real opposite-sign tag + probe pairs from data (tnp.fill definitions): the probe passes tight ID
             (+ IP cuts) / the probe fails it. Chosen for DISPLAY by documented rules -- not a physics result.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import V2, Checker, add_zmumu_path, finalize, load_json, provenance, standard_args  # noqa: E402

add_zmumu_path()
import awkward as ak  # noqa: E402
from zmumu import config, objects, regions, skim, tnp  # noqa: E402
from extract_zmumu_events import BRANCHES, _event_head, _jets, _muon, r5, r6  # noqa: E402

PKL = V2 / "tnp" / "tnp_fits.pkl"
RES = V2 / "tnp" / "tnp_result.json"
RECO = V2 / "tnp" / "reco_result.json"
APPLY_RULES = {
    "barrel": "both muons |eta| < 0.9 and 40 < pT < 50 GeV; 88 < m < 94 GeV",
    "forward": "one muon 2.1 < |eta| < 2.4, the other |eta| < 0.9; 86 < m < 96 GeV",
    "overlap_endcap": "one muon 0.9 < |eta| < 1.2, the other 1.2 < |eta| < 2.1; 86 < m < 96 GeV",
    "common": "DY aMC@NLO (first skim file, first 150 MB chunk), signal-region selection of zmumu/regions.py with the momentum "
              "calibration, both muons generator-matched prompt, LHE flavour mu; first match in file order",
}
EVENT_RULES = {
    "pass": "tag tight + pT > 26 + trigger-matched; probe loose, passes tightId + |dxy| < 0.2 + |dz| < 0.5; opposite sign; "
            "88 < m < 94 GeV; both |eta| < 0.9; |dphi| > 2.8; probe 30 < pT < 60; first match in file order",
    "fail": "same tag and pair; probe fails tightId + IP; 86 < m < 96 GeV; both |eta| < 0.9; |dphi| > 2.5; tracker-only probe "
            "(isGlobal false) preferred, else any failing probe; first match in file order",
}


def _cell_fit(det, tmpl_p, tmpl_f, x):
    """Model components of a stored nominal fit (tnp.fit_cell, expo background) at the fit's mass points."""
    p = det["params"]
    sig_p = p["N"] * p["eps"] * tnp._shift_smear(tmpl_p, p["dm_p"], p["sg_p"], x)
    sig_f = p["N"] * (1 - p["eps"]) * tnp._shift_smear(tmpl_f, p["dm_f"], p["sg_f"], x)
    bkg_p = p["Bp"] * tnp._bkg("expo", [p["lam_p"]], x, tnp.MASS_LO)
    bkg_f = p["Bf"] * tnp._bkg("expo", [p["lam_f"]], x, tnp.MASS_LO)
    return sig_p, bkg_p, sig_f, bkg_f


def _templates(mc, i, j):
    tp = np.clip(mc["id_nom_os_pass_gen"][i, j], 0, None)
    tf = np.clip(mc["id_nom_os_fail_gen"][i, j], 0, None)
    tp = tp / tp.sum()
    tf = tf / tf.sum() if tf.sum() > 0 else tp
    return tp, tf


def pick_events():
    import uproot
    f0 = skim.skim_files("data_2016G")[0]
    with uproot.open(f0) as f:
        ev = f["Events"].arrays(filter_name=BRANCHES + ["Muon_nStations", "Muon_nTrackerLayers"], library="ak")
    ev = ev[regions.fired(ev) & regions.event_clean(ev)]
    masks = regions.muon_masks(ev)
    matched = regions.trigger_match(ev)
    ip = (abs(ev.Muon_dxy) < config.MU_DXY_MAX) & (abs(ev.Muon_dz) < config.MU_DZ_MAX)
    id_pass = ev.Muon_tightId & ip
    is_tag = masks["tight"] & (ev.Muon_pt > config.TAG_PT_MIN) & matched
    tag, probe, mass, opposite = tnp._pairs(ev)
    flat = lambda a: ak.to_numpy(ak.flatten(a))
    n_pairs = ak.to_numpy(ak.num(tag, axis=1))
    evi = np.repeat(np.arange(len(ev)), n_pairs)
    t, p, m = flat(tag), flat(probe), flat(mass)
    base = flat(is_tag[tag] & masks["loose"][probe] & opposite)
    teta, peta = np.abs(flat(ev.Muon_eta[tag])), np.abs(flat(ev.Muon_eta[probe]))
    dphi = np.abs(objects.delta_phi(flat(ev.Muon_phi[tag]), flat(ev.Muon_phi[probe])))
    ppt, ppass, pglob = flat(ev.Muon_pt[probe]), flat(id_pass[probe]), flat(ev.Muon_isGlobal[probe])
    barrel = (teta < 0.9) & (peta < 0.9)
    c_pass = np.flatnonzero(base & ppass & barrel & (m > 88) & (m < 94) & (dphi > 2.8) & (ppt > 30) & (ppt < 60))
    c_fail_all = base & ~ppass & barrel & (m > 86) & (m < 96) & (dphi > 2.5)
    c_fail = np.flatnonzero(c_fail_all & ~pglob)
    fail_tracker_only = len(c_fail) > 0
    if not fail_tracker_only:
        c_fail = np.flatnonzero(c_fail_all)
    if not len(c_pass) or not len(c_fail):
        sys.exit(f"no display pair found in {f0} (pass {len(c_pass)}, fail {len(c_fail)})")

    def record(jp, rule):
        k = int(evi[jp])
        sub = ev[k:k + 1]
        it, iq = np.array([int(t[jp])]), np.array([int(p[jp])])
        rt, rp = _muon(sub, it, matched[k:k + 1], 0), _muon(sub, iq, matched[k:k + 1], 0)
        for r, i in ((rt, it), (rp, iq)):
            r["nStations"] = int(objects.take(sub["Muon_nStations"], i)[0])
            r["nTrackerLayers"] = int(objects.take(sub["Muon_nTrackerLayers"], i)[0])
        rp["passes_id"] = bool(objects.take(id_pass[k:k + 1], iq)[0])
        rec = _event_head(sub, 0)
        rec.update({"rule": rule, "tag": rt, "probe": rp, "muons": [rt, rp], "fsr_photons": [], "mass_bare": r5(m[jp]),
                    "dphi": r6(objects.delta_phi(rt["phi"], rp["phi"])), "jets": _jets(sub, 0, {"tag": rt, "probe": rp}),
                    "skim_file": Path(f0).name})
        return rec

    return {"rules": EVENT_RULES, "fail_tracker_only": fail_tracker_only,
            "n_candidates": {"pass": int(len(c_pass)), "fail": int(len(c_fail))},
            "pass": record(int(c_pass[0]), "pass"), "fail": record(int(c_fail[0]), "fail")}


def _rebin1(a):
    """120 x 0.5 GeV fit bins -> 60 x 1 GeV (adjacent pairs summed), the binning the clips draw."""
    a = np.asarray(a, dtype=float)
    assert a.shape == (120,), a.shape
    return a.reshape(60, 2).sum(axis=1)


def barrel_cells(det, data, mc, j=0):
    """Every pT cell of one |eta| column: spectra, stored fit model and its background part, in 1 GeV bins, data and simulation."""
    cells = []
    for i in range(len(tnp.PT_EDGES) - 1):
        tp, tf = _templates(mc, i, j)
        c = {"ipt": i, "ieta": j, "pt_range": [float(tnp.PT_EDGES[i]), float(tnp.PT_EDGES[i + 1])],
             "eta_range": [float(tnp.ETA_EDGES[j]), float(tnp.ETA_EDGES[j + 1])]}
        for sample, h in (("data", data), ("mc", mc)):
            d = det[(sample, i, j)]
            _, bkg_p, _, bkg_f = _cell_fit(d, tp, tf, np.asarray(d["x"]))
            c[sample] = {"pass": _rebin1(h["id_nom_os_pass"][i, j]), "fail": _rebin1(h["id_nom_os_fail"][i, j]),
                         "model_pass": _rebin1(d["model_pass"]), "model_fail": _rebin1(d["model_fail"]),
                         "bkg_pass": _rebin1(bkg_p), "bkg_fail": _rebin1(bkg_f),
                         "eps": d["eps"], "err": d["err"], "N": d["N"], "n_pass_signal": float(d["N"] * d["eps"]),
                         "n_fail_signal": float(d["N"] * (1 - d["eps"])), "bkg_frac_fail": d["bkg_frac_fail"]}
        cells.append(c)
    return cells


def pick_apply_examples():
    """Three real simulated Z -> mumu events with the scale factors the analysis gives them (tnp.ScaleFactors, the fit input weights)."""
    import uproot
    from v2_4_histograms import BRANCHES as HIST_BRANCHES          # z-mumu/scripts, read-only import
    from zmumu import histograms as H, momentum
    f0 = skim.skim_files("DY_NLO")[0]
    sf = tnp.ScaleFactors(RES)
    calib = momentum.MomentumCalibration(V2 / "momentum.json")
    with uproot.open(f0) as f:
        ev = next(f["Events"].iterate(filter_name=HIST_BRANCHES, step_size="150 MB", library="ak"))
    ev = ev[regions.fired(ev) & regions.event_clean(ev)]
    e = regions.with_muon_pt(ev, calib.scaled_pt(ev, "nominal"))
    d = regions.dimuon_regions(e, regions.muon_masks(e), regions.trigger_match(e), np.ones(len(e), dtype=bool))["SR"]
    idx = d["idx"]
    flav = np.asarray(ev.gen_lhe_flavour)[idx]
    ok = regions.is_prompt(d["flav1"]) & regions.is_prompt(d["flav2"]) & np.array([H.FLAVOUR_SAMPLE.get(int(v)) == "DYmumu" for v in flav])
    pt1, pt2, m = np.asarray(d["pt1"]), np.asarray(d["pt2"]), np.asarray(d["mass"])
    a1, a2 = np.abs(np.asarray(d["eta1"])), np.abs(np.asarray(d["eta2"]))
    zwin = (m > 86) & (m < 96)
    sel = {
        "barrel": ok & (a1 < 0.9) & (a2 < 0.9) & (pt1 > 40) & (pt1 < 50) & (pt2 > 40) & (pt2 < 50) & (m > 88) & (m < 94),
        "forward": ok & zwin & (((a1 > 2.1) & (a2 < 0.9)) | ((a2 > 2.1) & (a1 < 0.9))),
        "overlap_endcap": ok & zwin & (((a1 > 0.9) & (a1 < 1.2) & (a2 > 1.2) & (a2 < 2.1)) | ((a2 > 0.9) & (a2 < 1.2) & (a1 > 1.2) & (a1 < 2.1))),
    }
    npt, neta = len(tnp.PT_EDGES) - 1, len(tnp.ETA_EDGES) - 1

    def muon(pt, eta):
        p, a = np.array([float(pt)]), np.array([float(eta)])
        return {"pt": r5(pt), "eta": r6(eta),
                "ipt": int(np.clip(np.searchsorted(tnp.PT_EDGES, float(pt), side="right") - 1, 0, npt - 1)),
                "ieta": int(np.clip(np.searchsorted(tnp.ETA_EDGES, abs(float(eta)), side="right") - 1, 0, neta - 1)),
                "sf_id": float(sf.muon_sf("id", p, a)[0]), "sf_iso": float(sf.muon_sf("iso", p, a)[0])}

    out = {"rules": APPLY_RULES, "skim_file": Path(f0).name, "n_candidates": {k: int(np.count_nonzero(v)) for k, v in sel.items()}, "events": []}
    for rule, mask in sel.items():
        hits = np.flatnonzero(mask)
        if not len(hits):
            sys.exit(f"no simulated example event for rule {rule!r} in the first chunk of {f0}")
        k = int(hits[0])
        mu1, mu2 = muon(pt1[k], d["eta1"][k]), muon(pt2[k], d["eta2"][k])
        w = sf.event_weights(pt1[k:k + 1], np.asarray(d["eta1"])[k:k + 1], pt2[k:k + 1], np.asarray(d["eta2"])[k:k + 1])["nominal"]
        out["events"].append({"rule": rule, "mass": r5(m[k]), "muons": [mu1, mu2], "w_id": mu1["sf_id"] * mu2["sf_id"],
                              "w_iso": mu1["sf_iso"] * mu2["sf_iso"], "w_event": float(w[0])})
    return out


def build():
    o = pickle.load(open(PKL, "rb"))
    res, data, mc = o["result"], o["data"], o["mc"]
    det = res["details"]["id"]
    npt, neta = len(tnp.PT_EDGES) - 1, len(tnp.ETA_EDGES) - 1
    n_data = np.array([[det[("data", i, j)]["N"] for j in range(neta)] for i in range(npt)])
    i, j = (int(v) for v in np.unravel_index(np.argmax(n_data), n_data.shape))
    tp, tf = _templates(mc, i, j)
    cell = {"ipt": i, "ieta": j, "pt_range": [float(tnp.PT_EDGES[i]), float(tnp.PT_EDGES[i + 1])],
            "eta_range": [float(tnp.ETA_EDGES[j]), float(tnp.ETA_EDGES[j + 1])],
            "mass_edges": tnp.MASS_EDGES.tolist(), "definition": "probe loose muon; pass = tightId + |dxy| < 0.2 + |dz| < 0.5 cm; "
            "opposite-sign pairs, 60 < m < 120 GeV, both orientations; simultaneous extended binned Poisson fit (docs/12)"}
    for sample, h in (("data", data), ("mc", mc)):
        d = det[(sample, i, j)]
        x = np.asarray(d["x"])
        sig_p, bkg_p, sig_f, bkg_f = _cell_fit(d, tp, tf, x)
        cell[sample] = {"pass": np.asarray(h["id_nom_os_pass"][i, j]), "fail": np.asarray(h["id_nom_os_fail"][i, j]),
                        "eps": d["eps"], "err": d["err"], "N": d["N"], "chi2": d["chi2"], "ndf": d["ndf"],
                        "fit_x": x, "model_pass": np.asarray(d["model_pass"]), "model_fail": np.asarray(d["model_fail"]),
                        "signal_pass": sig_p, "bkg_pass": bkg_p, "signal_fail": sig_f, "bkg_fail": bkg_f,
                        "n_pass_signal": float(d["N"] * d["eps"]), "n_fail_signal": float(d["N"] * (1 - d["eps"])),
                        "bkg_frac_fail": d["bkg_frac_fail"]}
    sf_id, sf_iso = np.asarray(res["sf"]["id"]), np.asarray(res["sf"]["iso"])
    out = {
        "provenance": provenance("extract_zmumu_tnp.py", [PKL, RES, RECO, "$BND_SKIM_DIR/data_2016G (first file)", "$BND_SKIM_DIR/DY_NLO (first file)"],
                                 binning="`cell` keeps the 120 x 0.5 GeV fit bins; `cells_barrel` is summed to 60 x 1 GeV here (adjacent pairs), "
                                         "so a scene never rebins it",
                                 anchors=["z-mumu/docs/12-tag-and-probe-fits.md", "z-mumu/output/v2/RESULTS_v2.md"]),
        "pt_edges": tnp.PT_EDGES.tolist(), "eta_edges": tnp.ETA_EDGES.tolist(),
        "n_pairs_data": int(res["meta"]["n_pairs_data"]), "n_pairs_mc": float(res["meta"]["n_pairs_mc"]),
        "cell": cell,
        "cells_barrel": barrel_cells(det, data, mc, 0),
        "reco": {k: load_json(RECO)["meta"][k] for k in ("sf_reco_per_muon", "sf_reco_per_muon_err", "sf_reco_per_event", "sf_reco_per_event_err",
                                                          "sf_trk_per_muon", "sf_trk_per_muon_err", "sf_mutrk_per_muon", "sf_mutrk_per_muon_err")},
        "apply_examples": pick_apply_examples(),
        "id": {"eff_data": res["eff"]["id_data_nominal"], "eff_mc": res["eff"]["id_mc_nominal"],
               "eff_err_data": res["eff_err"]["id_data_nominal"], "eff_err_mc": res["eff_err"]["id_mc_nominal"],
               "sf": sf_id, "sf_err": res["sf_err"]["id"], "map_mean": float(np.mean(sf_id)), "n_data": n_data},
        "iso": {"sf": sf_iso, "sf_err": res["sf_err"]["iso"], "map_mean": float(np.mean(sf_iso))},
        "trigger": {"pt_edges": res["trig_pt_edges"], "eff_data": res["eff"]["trig_data"], "eff_mc": res["eff"]["trig_mc"],
                    "plateau_data": res["meta"]["trig_plateau_data"], "plateau_mc": res["meta"]["trig_plateau_mc"]},
        "events": pick_events(),
    }
    return out


def verify(d, ck: Checker):
    r = load_json(RES)
    c = d["cell"]
    i, j = c["ipt"], c["ieta"]
    D = "z-mumu/docs/12-tag-and-probe-fits.md"
    ck.check("pairs in data 20 050 129", d["n_pairs_data"], 20050129, 0, D)
    for s in ("data", "mc"):
        cs = c[s]
        ck.check(f"cell eps_{s} == tnp_result.json", cs["eps"], r["eff"][f"id_{s}_nominal"][i][j], 1e-12, "tnp_result.json")
        ck.check(f"cell {s}: signal + bkg == stored fit model (pass)", np.add(cs["signal_pass"], cs["bkg_pass"]), cs["model_pass"], 1e-6, rel=True)
        ck.check(f"cell {s}: signal + bkg == stored fit model (fail)", np.add(cs["signal_fail"], cs["bkg_fail"]), cs["model_fail"], 1e-6, rel=True)
        ck.check_true(f"cell {s}: 120 bins pass / fail", len(cs["pass"]) == 120 and len(cs["fail"]) == 120)
    ck.check_true("cell: the most populated data cell", c["data"]["N"] == max(max(row) for row in d["id"]["n_data"]))
    ck.check_true("cell: data below simulation", c["data"]["eps"] < c["mc"]["eps"])
    ck.check("ID SF map mean 0.980", round(d["id"]["map_mean"], 3), 0.980, 1e-9, D)
    ck.check("iso SF map mean 1.005", round(d["iso"]["map_mean"], 3), 1.005, 1.5e-3, D)
    ck.check("trigger plateau data 0.907", round(d["trigger"]["plateau_data"], 3), 0.907, 1e-9, D)
    ck.check("trigger plateau simulation 0.923", round(d["trigger"]["plateau_mc"], 3), 0.923, 1e-9, D)
    ck.check("SF map == eps_data / eps_mc", d["id"]["sf"], np.divide(d["id"]["eff_data"], d["id"]["eff_mc"]), 1e-12, "tnp_result.json")
    ck.check_true("data below simulation in every ID cell", bool(np.all(np.less(d["id"]["eff_data"], d["id"]["eff_mc"]))),
                  D, "docs/12: 'The data ID efficiency is 1-2% below the simulation everywhere'")
    cb = d["cells_barrel"]
    ck.check_true("cells_barrel: the ten barrel pT cells in order", [x["ipt"] for x in cb] == list(range(10)) and all(x["ieta"] == 0 for x in cb))
    for x in cb:
        for smp in ("data", "mc"):
            ck.check(f"cells_barrel[{x['ipt']}] eps_{smp} == tnp_result.json", x[smp]["eps"], r["eff"][f"id_{smp}_nominal"][x["ipt"]][0], 1e-12, "tnp_result.json")
            ck.check_true(f"cells_barrel[{x['ipt']}] {smp}: 60 bins, background below the model",
                          all(len(x[smp][k]) == 60 for k in ("pass", "fail", "model_pass", "model_fail", "bkg_pass", "bkg_fail"))
                          and bool(np.all(np.asarray(x[smp]["bkg_fail"]) <= np.asarray(x[smp]["model_fail"]) * (1 + 1e-9))))
    for smp in ("data", "mc"):
        ck.check(f"cells_barrel[4] {smp} pass == `cell` pass summed in pairs", cb[i][smp]["pass"], np.asarray(c[smp]["pass"]).reshape(60, 2).sum(1), 1e-9, rel=True)
        ck.check(f"cells_barrel[4] {smp} model_fail == `cell` model_fail summed in pairs", cb[i][smp]["model_fail"],
                 np.asarray(c[smp]["model_fail"]).reshape(60, 2).sum(1), 1e-9, rel=True)
    ck.check("background share of the failing probes, barrel 20-25 GeV: 61 %", cb[0]["data"]["bkg_frac_fail"], 0.610, 5e-4, "z-mumu/docs/17 section 5")
    ck.check("background share of the failing probes, barrel 40-45 GeV: 7.4 %", cb[4]["data"]["bkg_frac_fail"], 0.074, 5e-4, "z-mumu/docs/17 section 5")
    ck.check("reconstruction SF per muon 1.0001 +- 0.0013", [d["reco"]["sf_reco_per_muon"], d["reco"]["sf_reco_per_muon_err"]], [1.0001, 0.0013], 5e-5, "FREEZE.md")
    ck.check("reconstruction SF per event 1.0002", d["reco"]["sf_reco_per_event"], 1.0002, 5e-5, "FREEZE.md")
    ax = d["apply_examples"]["events"]
    ck.check_true("apply examples: one event per rule", [x["rule"] for x in ax] == ["barrel", "forward", "overlap_endcap"])
    for x in ax:
        for mu in x["muons"]:
            ck.check(f"apply example {x['rule']}: SF_ID of the muon == its map cell [{mu['ipt']}][{mu['ieta']}]", mu["sf_id"], d["id"]["sf"][mu["ipt"]][mu["ieta"]], 1e-12, "tnp_result.json")
            ck.check(f"apply example {x['rule']}: SF_iso of the muon == its map cell", mu["sf_iso"], d["iso"]["sf"][mu["ipt"]][mu["ieta"]], 1e-12, "tnp_result.json")
        ck.check_true(f"apply example {x['rule']}: event weight (ID pair x iso pair x trigger) within 6 % of 1, mass in the window", abs(x["w_event"] - 1) < 0.06 and 60 < x["mass"] < 120,
                      detail=f"w_id {x['w_id']:.4f}, w_event {x['w_event']:.4f}, m {x['mass']:.2f}")
    e = d["events"]
    ck.check_true("event pass: probe passes ID", e["pass"]["probe"]["passes_id"] and e["pass"]["tag"]["matched"])
    ck.check_true("event fail: probe fails ID", (not e["fail"]["probe"]["passes_id"]) and e["fail"]["tag"]["matched"])
    for k in ("pass", "fail"):
        ev = e[k]
        ck.check_true(f"event {k}: opposite sign", ev["tag"]["charge"] * ev["probe"]["charge"] < 0)
        ck.check_true(f"event {k}: tag pT > 26", ev["tag"]["pt"] > 26.0)
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "zmumu_tnp.json")
    finalize(ap.parse_args(), build, verify, "zmumu_tnp")


if __name__ == "__main__":
    main()
