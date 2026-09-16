#!/usr/bin/env python
"""Freeze the Z -> tautau jet -> tau_h fake-factor estimate: the era x DM x N_jets x pT table, the |eta(tau1)| and pT(tau2)
closure corrections (before / after), C_OS/SS, the region counts and the fake yields per BDT category.

    source setup.sh && python presentation/data/extract_ztautau_fakes.py [--json PATH] [--check-only] [--no-recompute]

Reads (read-only): z-tautau/output/data/fakefactors.json (mcsub: fakes.measure + closure_corrections + osss_correction),
z-tautau/fit/fitinputs/ztautau.root.meta.json (C per era x N_jets x BDT category), z-tautau/output/data/yields.json, and --
for the closure "after" the corrections -- the data and MC ntuples through ztautau.analysis (same-sign obs / pred with the
final, closure-corrected table, MC subtracted; ~1-2 min). Writes presentation/data/ztautau_fakes.json.
Anchors: output/RESULTS.md (13592 fakes), docs/05-fake-factors.md (closure table, C table, region counts).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATASET_TAUTAU, VERSION_TAUTAU, ZTAUTAU, Checker, add_ztautau_path, finalize, load_json,  # noqa: E402
                             provenance, standard_args)

FF_JSON = ZTAUTAU / "output" / "data" / "fakefactors.json"
META = ZTAUTAU / "fit" / "fitinputs" / "ztautau.root.meta.json"
YIELDS = ZTAUTAU / "output" / "data" / "yields.json"
REGIONS = ["tautau_SR0", "tautau_SR1", "tautau_SR2"]
DOCS05_ETA_BEFORE = [0.92, 0.99, 1.08, 1.16, 0.98, 0.93]
DOCS05_PT2_BEFORE = [1.02, 0.99, 1.00, 0.94, 0.95]


def closure_after(table):
    """Same-sign obs/pred of the fake factor in |eta(tau1)| and pT(tau2) bins, MC subtracted, with the raw table
    (closure=False), with the final table (closure=True) -- the re-implementation of fakes.closure_corrections
    with the stored table instead of measuring it. Returns {var: {"edges", "before_raw", "after", "obs", "pred_after", "pred_raw"}}."""
    add_ztautau_path()
    from ztautau import analysis, config, fakes
    data = analysis.load_data()
    reg = analysis.regions(data)
    mc = []
    for key in analysis.available_mc():
        d, _ = analysis.load(key)
        if not len(d):
            continue
        r = analysis.regions(d, is_mc=True)
        w = analysis.subtraction_weights(key, analysis.weights(d, key))
        mc.append((d, r, w))
    out = {}
    for var, edges in (("eta", config.FF_CLOSURE_ETA_BINS), ("pt2", config.FF_CLOSURE_PT2_BINS)):
        edges = np.asarray(edges, dtype=float)

        def xval(a, m):
            return np.abs(a["t1_eta"][m]) if var == "eta" else np.clip(a["t2_pt"][m], edges[0], edges[-1] - 1e-6)

        obs = np.histogram(xval(data, reg["SS_T"]), bins=edges)[0].astype(float)
        res = {"edges": edges.tolist(), "obs_data": obs.tolist()}
        for tag, closure in (("raw", False), ("after", True)):
            wp = fakes.evaluate(table, data, reg["SS_L"], closure=closure)
            pred = np.histogram(xval(data, reg["SS_L"]), bins=edges, weights=wp)[0]
            pred2 = np.histogram(xval(data, reg["SS_L"]), bins=edges, weights=wp ** 2)[0]
            o = obs.copy()
            for a, r, w in mc:
                o -= np.histogram(xval(a, r["SS_T"]), bins=edges, weights=w[r["SS_T"]])[0]
                pred -= np.histogram(xval(a, r["SS_L"]), bins=edges, weights=w[r["SS_L"]] * fakes.evaluate(table, a, r["SS_L"], closure=closure))[0]
            ratio = np.where(pred > 0, o / np.maximum(pred, 1e-9), 1.0)
            err = ratio * np.sqrt(1.0 / np.maximum(o, 1.0) + pred2 / np.maximum(pred, 1e-9) ** 2)
            res[f"{'before_raw' if tag == 'raw' else 'after'}"] = ratio.tolist()
            res[f"{'before_raw' if tag == 'raw' else 'after'}_err"] = err.tolist()
            res[f"obs_mcsub_{tag}"] = o.tolist()
            res[f"pred_{tag}"] = pred.tolist()
        out[var] = res
    out["n_mc_samples_subtracted"] = len(mc)
    return out


def fold_signed_eta(block):
    """Fold the signed eta(tau1) closure histogram (14 bins of 0.3 from -2.1 to 2.1) into 7 |eta| bins of 0.3."""
    edges = np.asarray(block["edges"], dtype=float)
    obs = np.asarray(block["obs"], dtype=float) - np.asarray(block["obs_mc"], dtype=float)
    pred, pred_raw = np.asarray(block["pred"], dtype=float), np.asarray(block["pred_nocorr"], dtype=float)
    n = len(edges) - 1
    assert n % 2 == 0 and np.isclose(edges[n // 2], 0.0)
    fold = lambda a: a[n // 2:] + a[: n // 2][::-1]
    o, p, pr = fold(obs), fold(pred), fold(pred_raw)
    return {"edges": edges[n // 2:].tolist(), "before": (o / pr).tolist(), "after": (o / p).tolist(), "obs_mcsub": o.tolist(),
            "pred_after": p.tolist(), "pred_before": pr.tolist(),
            "note": "fakefactors.json mcsub.closure.t1_eta (signed, 0.3 bins, MC subtracted) folded into |eta|; the 0.3 binning does not "
                    "align with the 6 correction bins, so this is the cross-check, `closure.eta` the correction bins"}


def build(recompute=True):
    ff = load_json(FF_JSON)
    meta = load_json(META)
    y = load_json(YIELDS)
    m = ff["mcsub"]
    t = m["ff"]
    eras, dms, njb, pte = t["eras"], t["dms"], t["njet_bins"], t["pt_bins"]
    vals, errs = np.asarray(t["ff"]), np.asarray(t["err"])
    assert vals.shape == (2, 4, 3, 5), vals.shape
    clo = t["closure"]
    after = None
    if recompute:
        add_ztautau_path()
        from ztautau import fakes
        after = closure_after(fakes.from_json(t))
    c = m["osss"]
    C_cat = np.asarray(meta["C_osss"])
    out = {
        "provenance": provenance("extract_ztautau_fakes.py", [FF_JSON, META, YIELDS], dataset=DATASET_TAUTAU, version=VERSION_TAUTAU,
                                 method=("FF(era, DM, N_jets, pT) = [N(SS, tau1 Tight, tau2 Tight) - MC] / [N(SS, tau1 VVVLoose&!Tight, tau2 Tight) - MC] "
                                         "(fakes.measure), times closure corrections f(|eta(tau1)|) g(pT(tau2)) (same-sign obs/pred, fakes.closure_corrections), "
                                         "applied to OS events with tau1 failing Tight (AR, MC genuine-tau1 subtracted), times C_OS/SS(era, N_jets, BDT category) "
                                         "from the tau2 anti-isolated sideband (fakes.osss_correction)"),
                                 closure_after=("same-sign obs/pred recomputed from the ntuples with the final table through ztautau.analysis / fakes.evaluate "
                                                "(MC subtracted with analysis.subtraction_weights)" if recompute else "not recomputed (--no-recompute)"),
                                 anchors=["z-tautau/output/RESULTS.md", "z-tautau/docs/05-fake-factors.md", "z-tautau/handoff.md"]),
        "ff": {"eras": eras, "dms": dms, "njet_bins": njb, "pt_edges": pte, "values": vals.tolist(), "err": errs.tolist(),
               "num": t["num"], "den": t["den"], "num_mc": t["num_mc"], "den_mc": t["den_mc"],
               "shape": "[era G/H][DM 0,1,10,11][N_jets 0,1,>=2][pT 40-45,45-50,50-60,60-80,>80]"},
        "ff_map": {"era": "G", "njet": 0, "row_labels_dm": dms, "pt_edges": pte, "values": vals[0, :, 0, :].tolist(), "err": errs[0, :, 0, :].tolist(),
                   "num": np.asarray(t["num"])[0, :, 0, :].tolist(), "den": np.asarray(t["den"])[0, :, 0, :].tolist()},
        "closure": {
            "eta": {"edges": clo["eta"]["edges"], "before": clo["eta"]["values"], "err": clo["eta"]["err"],
                    "obs": clo["eta"]["obs"], "pred": clo["eta"]["pred"],
                    "after": after["eta"]["after"] if after else None, "after_err": after["eta"]["after_err"] if after else None,
                    "before_recomputed": after["eta"]["before_raw"] if after else None,
                    "docs05_before": DOCS05_ETA_BEFORE,
                    "note": "before = the correction f itself (same-sign obs/pred with the era x DM x N_jets x pT table, MC subtracted); "
                            "after = obs/pred with f and g applied (recomputed); docs/05 quotes the v2 (Medium) pass"},
            "pt2": {"edges": clo["pt2"]["edges"], "before": clo["pt2"]["values"], "err": clo["pt2"]["err"],
                    "obs": clo["pt2"]["obs"], "pred": clo["pt2"]["pred"],
                    "after": after["pt2"]["after"] if after else None, "after_err": after["pt2"]["after_err"] if after else None,
                    "before_raw": after["pt2"]["before_raw"] if after else None,
                    "docs05_before": DOCS05_PT2_BEFORE,
                    "note": "before = g measured after f (the stored correction); before_raw = without f and g; after = with both"},
            "eta_folded_0p3": fold_signed_eta(m["closure"]["t1_eta"]),
        },
        "osss": {"inclusive": c["inclusive"]["C"], "inclusive_stat": c["inclusive"]["stat"], "inclusive_detail": dict(c["inclusive"]),
                 "per_era_njet": c["C"], "per_era_njet_stat": c["stat"],
                 "per_category": C_cat.tolist(), "per_category_stat": meta["C_osss_stat"], "rel_unc": meta["C_osss_rel_unc"],
                 "per_category_njet_inclusive": meta["C_osss_njet_inclusive"], "replaced_by_njet_inclusive": meta["C_osss_replaced"],
                 "axes": "[era G/H][N_jets 0,1,>=2][BDT category 0,1,2]", "syst_extrapolation": 0.03,
                 "range_per_category": [float(C_cat.min()), float(C_cat.max())],
                 "osss_nps": meta["osss_nps"]},
        "yields": {"n_fake": {**{r: y["regions"][r]["Fakes"]["value"] for r in REGIONS}, "total": y["yields"]["Fakes"]["value"]},
                   "n_fake_stat": {**{r: y["regions"][r]["Fakes"]["stat"] for r in REGIONS}, "total": y["yields"]["Fakes"]["stat"]},
                   "fraction_of_sr": y["yields"]["Fakes"]["value"] / y["yields"]["Data"]["value"],
                   "regions": {"SS_T": ff["contamination"]["SS_T"]["data"], "SS_L": ff["contamination"]["SS_L"]["data"],
                               "AR": ff["contamination"]["AR"]["data"], "SR": ff["contamination"]["SR"]["data"],
                               "OSAI_T": ff["contamination"]["OSAI_T"]["data"], "OSAI_L": ff["contamination"]["OSAI_L"]["data"]},
                   "contamination": ff["contamination"],
                   "closure_nps": meta["closure_nps"]},
        "nonclosure_m_tt": m["nonclosure_m_tt"],
    }
    return out


def verify(d, ck: Checker):
    D5 = "z-tautau/docs/05-fake-factors.md"
    R = "z-tautau/output/RESULTS.md"
    ck.check("n_fake total 13592 +- 1", d["yields"]["n_fake"]["total"], 13592, 1, R)
    ck.check("n_fake per region 12351 / 1044 / 197", [d["yields"]["n_fake"][r] for r in REGIONS], [12351, 1044, 197], 0.5, R)
    ck.check("n_fake total == sum of regions", d["yields"]["n_fake"]["total"], sum(d["yields"]["n_fake"][r] for r in REGIONS), 1e-6, "consistency")
    ck.check("fraction of SR 0.642", d["yields"]["fraction_of_sr"], 0.642, 0.0005, "docs/02, docs/05 (64 %)")
    ck.check("regions SS_T 9554 / SS_L 86929 / AR 120805", [d["yields"]["regions"][k] for k in ("SS_T", "SS_L", "AR")], [9554, 86929, 120805], 0, D5 + " contamination table")
    ck.check("AR genuine-tau1 fraction 2.5 %", 100 * d["yields"]["contamination"]["AR"]["fraction"], 2.5, 0.05, D5)
    ck.check("FF table shape 2 x 4 x 3 x 5", np.asarray(d["ff"]["values"]).shape, (2, 4, 3, 5), 0, D5 + " (120 bins)")
    ck.check("ff_map = era G, 0 jets slice", d["ff_map"]["values"], np.asarray(d["ff"]["values"])[0, :, 0, :], 0, "consistency")
    ck.check_true("FF values in (0, 1)", bool(np.all((np.asarray(d["ff"]["values"]) > 0) & (np.asarray(d["ff"]["values"]) < 1))))
    ck.check("closure eta edges", d["closure"]["eta"]["edges"], [0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1], 0, "config.FF_CLOSURE_ETA_BINS")
    ck.check("closure pt2 edges", d["closure"]["pt2"]["edges"], [40, 45, 50, 60, 80, 1000], 0, "config.FF_CLOSURE_PT2_BINS")
    ck.check("closure eta before within 0.03 of docs/05 (0.92, 0.99, 1.08, 1.16, 0.98, 0.93)", d["closure"]["eta"]["before"], DOCS05_ETA_BEFORE, 0.03,
             D5 + " closure table (quoted from the v2 Medium pass; v3 Tight stored values differ in bins 3 and 6)", soft=True)
    ck.check("closure pt2 before within 0.03 of docs/05 (1.02, 0.99, 1.00, 0.94, 0.95)", d["closure"]["pt2"]["before"], DOCS05_PT2_BEFORE, 0.03,
             D5 + " (v2 Medium pass)", soft=True)
    ck.check("closure eta before == obs/pred stored", d["closure"]["eta"]["before"], np.asarray(d["closure"]["eta"]["obs"]) / np.asarray(d["closure"]["eta"]["pred"]), 1e-9, "consistency", rel=True)
    if d["closure"]["eta"]["after"] is not None:
        ck.check("closure eta before recomputed == stored", d["closure"]["eta"]["before_recomputed"], d["closure"]["eta"]["before"], 1e-6, "ntuples vs fakefactors.json", rel=True)
        ck.check_true("closure eta after within 3 % of 1", bool(np.all(np.abs(np.asarray(d["closure"]["eta"]["after"]) - 1) < 0.03)), detail=str(np.round(d["closure"]["eta"]["after"], 4).tolist()))
        ck.check_true("closure pt2 after within 3 % of 1", bool(np.all(np.abs(np.asarray(d["closure"]["pt2"]["after"]) - 1) < 0.03)), detail=str(np.round(d["closure"]["pt2"]["after"], 4).tolist()))
    f = d["closure"]["eta_folded_0p3"]
    ck.check_true("folded eta after closer to 1 than before", float(np.abs(np.asarray(f["after"]) - 1).max()) < float(np.abs(np.asarray(f["before"]) - 1).max()),
                  detail=f"before max|r-1| {np.abs(np.asarray(f['before']) - 1).max():.3f}, after {np.abs(np.asarray(f['after']) - 1).max():.3f}")
    ck.check("C inclusive within 0.01 of 1.052", d["osss"]["inclusive"], 1.052, 0.01, D5 + " (stored 1.0557; docs quote 1.052)")
    ck.check("C per era x N_jets (docs/05 table)", d["osss"]["per_era_njet"], [[1.057, 1.086, 1.077], [1.042, 1.037, 1.070]], 5e-4, D5)
    ck.check("C per category (docs/05 table)", np.transpose(np.asarray(d["osss"]["per_category"]), (2, 0, 1)),
             [[[1.056, 1.079, 1.059], [1.039, 1.029, 1.058]], [[1.109, 1.172, 1.097], [1.163, 1.149, 1.113]], [[1.187, 1.045, 1.305], [1.072, 0.960, 1.151]]],
             5e-4, D5 + " (category, era, N_jets)")
    # docs/05 per-category table: 0.960 (cat 2, H, 1 jet) .. 1.305 (cat 2, G, >= 2 jets); the brief's "1.04-1.31" was a
    # rounded reading and is corrected from this value
    ck.check("C range per category 0.960 .. 1.305 (docs/05 table)", d["osss"]["range_per_category"], [0.960, 1.305], 0.005,
             "z-tautau/docs/05-fake-factors.md, C per (era, N_jets, BDT category)")
    ck.check("C rel_unc shape and >= 3 %", np.asarray(d["osss"]["rel_unc"]).min() >= 0.03, True, 0, "config.FF_OSSS_SYST")
    ck.check("closure NPs c0_hi 1 %, c1_lo 7 %, c1_hi 21 %, c2_lo 15 %, c2_hi 44 %",
             [100 * d["yields"]["closure_nps"][k]["delta"] for k in ("FakeClosure_tautau_c0_hi", "FakeClosure_tautau_c1_lo", "FakeClosure_tautau_c1_hi", "FakeClosure_tautau_c2_lo", "FakeClosure_tautau_c2_hi")],
             [1, 7, 21, 15, 44], 0.6, "docs/07 table")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_fakes.json")
    ap.add_argument("--no-recompute", action="store_true", help="skip the ntuple pass (closure 'after' = null)")
    args = ap.parse_args()
    finalize(args, lambda: build(recompute=not args.no_recompute), verify, "ztautau_fakes")


if __name__ == "__main__":
    main()
