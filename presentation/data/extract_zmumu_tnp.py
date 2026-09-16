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
        "provenance": provenance("extract_zmumu_tnp.py", [PKL, RES, "$BND_SKIM_DIR/data_2016G (first file)"],
                                 anchors=["z-mumu/docs/12-tag-and-probe-fits.md", "z-mumu/output/v2/RESULTS_v2.md"]),
        "pt_edges": tnp.PT_EDGES.tolist(), "eta_edges": tnp.ETA_EDGES.tolist(),
        "n_pairs_data": int(res["meta"]["n_pairs_data"]), "n_pairs_mc": float(res["meta"]["n_pairs_mc"]),
        "cell": cell,
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
