#!/usr/bin/env python
"""Freeze the normalised 1D distributions of the 16 BDT inputs (config.BDT_FEATURES) for the two classes the Z -> tautau
BDT was trained on: the fiducial Z -> tautau simulation in the signal region and the fake estimate (application-region data
weighted by the nominal fake factor).

    source setup.sh && python presentation/data/extract_ztautau_bdt_inputs.py [--json PATH] [--check-only]

Reads (read-only, uproot through ztautau.analysis): the inclusive aMC@NLO ntuple DY_NLO of $BND_TAUTAU_CACHE/ntuples_v1
(SR, LHE tautau, gen_fid, both legs genPartFlav 5; weights sigma x L / sumw x genWeight, sign kept, no scale factors), the
Tau Run2016G+H data ntuples (AR = OS, tau1 VVVLoose-and-not-Tight, tau2 Tight, both pT > 40; weights FF x C with the
MC-subtracted table of z-tautau/output/data/fakefactors.json and the inclusive C_OS/SS, the genuine-tau1 simulation NOT
subtracted), the feature matrix of ztautau.bdt.features (the training definitions), and the feature order of
presentation/data/ztautau_bdt.json (importance_sorted). Writes presentation/data/ztautau_bdt_inputs.json.
Anchors: docs/09-bdt.md (signal at small dR), config.BDT_FEATURES, ztautau_bdt.json (importance).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATA_DIR, DATASET_TAUTAU, TAUTAU_NTUPLES, VERSION_TAUTAU, ZTAUTAU, Checker, add_ztautau_path,  # noqa: E402
                             finalize, load_json, provenance, standard_args)

FF_JSON = ZTAUTAU / "output" / "data" / "fakefactors.json"
BDT_FROZEN = DATA_DIR / "ztautau_bdt.json"
INCLUSIVE = "DY_NLO"
N_BINS = 20
# MathTex symbols (no words) and the binning of every input: (lo, hi) -> N_BINS equal bins, overflow into the last bin;
# integers: one bin per value; categories: one bin per decay mode
SYMBOL = {"t1_pt": r"p_T(\tau_1)", "t2_pt": r"p_T(\tau_2)", "pt_ratio": r"p_T(\tau_2)/p_T(\tau_1)", "t1_abseta": r"|\eta(\tau_1)|",
          "t2_abseta": r"|\eta(\tau_2)|", "dr_tt": r"\Delta R(\tau_1,\tau_2)", "dphi_tt": r"\Delta\phi(\tau_1,\tau_2)",
          "met": r"p_T^{\mathrm{miss}}", "met_sig": r"S(p_T^{\mathrm{miss}})", "pt_vis": r"p_T(\tau\tau)",
          "dphi_met_tt": r"\Delta\phi(p_T^{\mathrm{miss}},\tau\tau)", "pt_tt": r"p_T(\tau\tau+p_T^{\mathrm{miss}})",
          "njets": r"N_{\mathrm{jets}}", "jet1_pt": r"p_T(\mathrm{jet}_1)", "t1_dm": r"\mathrm{DM}(\tau_1)", "t2_dm": r"\mathrm{DM}(\tau_2)"}
UNIT = {"t1_pt": "GeV", "t2_pt": "GeV", "met": "GeV", "pt_vis": "GeV", "pt_tt": "GeV", "jet1_pt": "GeV"}
RANGE = {"t1_pt": (40.0, 200.0), "t2_pt": (40.0, 200.0), "pt_ratio": (0.0, 1.0), "t1_abseta": (0.0, 2.1), "t2_abseta": (0.0, 2.1),
         "dr_tt": (0.5, 5.0), "dphi_tt": (0.0, np.pi), "met": (0.0, 150.0), "met_sig": (0.0, 10.0), "pt_vis": (0.0, 300.0),
         "dphi_met_tt": (0.0, np.pi), "pt_tt": (0.0, 300.0), "jet1_pt": (0.0, 300.0)}
INTEGER = {"njets": (0, 6)}                        # values 0..6, >= 6 in the last bin
CATEGORIES = {"t1_dm": [0, 1, 10, 11], "t2_dm": [0, 1, 10, 11]}
DM_LABELS = {"0": r"1\text{-}\mathrm{prong}", "1": r"1\text{-}\mathrm{prong}+\pi^0", "10": r"3\text{-}\mathrm{prong}", "11": r"3\text{-}\mathrm{prong}+\pi^0"}


def binning(name):
    """('edges' | 'categories', values) of a feature."""
    if name in CATEGORIES:
        return "categories", list(CATEGORIES[name])
    if name in INTEGER:
        lo, hi = INTEGER[name]
        return "edges", [float(v) for v in range(lo, hi + 2)]
    lo, hi = RANGE[name]
    return "edges", np.linspace(lo, hi, N_BINS + 1).tolist()


def fill(name, x, w):
    """Weighted counts (and unweighted entries) of one feature: overflow into the last bin, underflow into the first."""
    kind, b = binning(name)
    x = np.asarray(x, dtype=np.float64)
    if kind == "categories":
        cats = np.asarray(b)
        h = np.array([w[np.isclose(x, c)].sum() for c in cats])
        n = np.array([int(np.isclose(x, c).sum()) for c in cats])
        return h, n
    e = np.asarray(b)
    xc = np.clip(x, e[0], e[-1] - 1e-9 * max(1.0, abs(e[-1])))
    h, _ = np.histogram(xc, bins=e, weights=w)
    n, _ = np.histogram(xc, bins=e)
    return h, n


def wmedian(x, w):
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    cw /= cw[-1]
    return float(x[o][min(int(np.argmax(cw >= 0.5)), len(x) - 1)])


def build():
    add_ztautau_path()
    from ztautau import analysis, bdt, config, fakes
    feats = list(config.BDT_FEATURES)
    frozen = load_json(BDT_FROZEN)
    order = [k for k, _ in frozen["importance_sorted"]]
    importance = dict((k, float(v)) for k, v in frozen["importance_sorted"])
    assert set(order) == set(feats) and bdt.FEATURES == feats
    # signal: fiducial Z -> tautau of the inclusive sample in the SR, both legs genuine tau_h
    d, _ = analysis.load(INCLUSIVE)
    r = analysis.regions(d, is_mc=True)
    sel = r["SR"] & (d["gen_lhe_flavour"] == 15) & (d["gen_fid"] > 0) & (d["t1_genflav"] == 5) & (d["t2_genflav"] == 5)
    Xs = bdt.features(d)[sel].astype(np.float64)
    ws = analysis.normalisation(INCLUSIVE) * d["genWeight"][sel].astype(np.float64)
    # fakes: AR data x FF(era, DM, N_jets, pT) x closure corrections x inclusive C_OS/SS
    ff = load_json(FF_JSON)
    table = fakes.from_json(ff["mcsub"]["ff"])
    C = float(ff["mcsub"]["osss"]["inclusive"]["C"])
    data = analysis.load_data()
    ar = analysis.regions(data)["AR"]
    Xb = bdt.features(data)[ar].astype(np.float64)
    wb = fakes.fake_weights(data, ar, table, C)[ar]
    print(f"[ztautau_bdt_inputs] signal {len(Xs)} events (sum w {ws.sum():.1f}, {np.mean(ws < 0):.1%} negative), fakes {len(Xb)} AR events (sum w {wb.sum():.1f}, C = {C:.4f})")
    features = []
    for name in order:
        i = feats.index(name)
        kind, b = binning(name)
        hs, ns = fill(name, Xs[:, i], ws)
        hb, nb = fill(name, Xb[:, i], wb)
        row = {"name": name, "symbol": SYMBOL[name], kind: b,
               "signal": (hs / hs.sum()).tolist(), "fakes": (hb / hb.sum()).tolist(),
               "signal_raw": hs.tolist(), "fakes_raw": hb.tolist(), "signal_n": ns.tolist(), "fakes_n": nb.tolist(),
               "signal_unit": (ns / ns.sum()).tolist(),
               "importance": importance[name],
               "median_signal": wmedian(Xs[:, i], ws), "median_fakes": wmedian(Xb[:, i], wb),
               "overflow_frac_signal": float(ws[Xs[:, i] >= b[-1]].sum() / ws.sum()) if kind == "edges" else 0.0,
               "overflow_frac_fakes": float(wb[Xb[:, i] >= b[-1]].sum() / wb.sum()) if kind == "edges" else 0.0}
        if name in UNIT:
            row["unit"] = UNIT[name]
        if name in INTEGER:
            row["integer"] = True
        features.append(row)
    out = {
        "provenance": provenance("extract_ztautau_bdt_inputs.py", [TAUTAU_NTUPLES / f"{INCLUSIVE}.root", TAUTAU_NTUPLES / "data_2016G.root", TAUTAU_NTUPLES / "data_2016H.root",
                                                                   FF_JSON, BDT_FROZEN, ZTAUTAU / "ztautau" / "config.py", ZTAUTAU / "docs" / "09-bdt.md"],
                                 dataset=DATASET_TAUTAU, version=VERSION_TAUTAU,
                                 signal={"sample": INCLUSIVE, "selection": "SR (OS, both DeepTau Tight, pT > 40, |eta| < 2.1), LHE tautau, gen_fid (60 < m_LHE < 120, "
                                                                          "both GenVisTau pT > 40, |eta| < 2.1), t1_genflav == 5 and t2_genflav == 5",
                                         "weights": "sigma x L / sumw(DY_NLO) x genWeight, sign kept (analysis.normalisation); no scale factors, no pileup, no prefiring",
                                         "note": "the BDT was trained on the stitched DY_NLO + DY_0J/1J/2J fiducial signal with the full event weights (positive part), "
                                                 "step3b_bdt.py; the inclusive sample alone with the uniform normalisation is drawn for the shapes"},
                                 fakes={"selection": "AR: OS, tau1 VVVLoose and not Tight, tau2 Tight, both pT > 40 (Tau Run2016G+H data)",
                                        "weights": f"FF(era, DM, N_jets, pT) x f(|eta(tau1)|) x g(pT(tau2)) (fakefactors.json mcsub.ff, fakes.evaluate) x inclusive C_OS/SS = {C:.4f}; "
                                                   "the genuine-tau1 simulation in the AR is NOT subtracted (2.5 % of the AR)"},
                                 features="ztautau.bdt.features (the training matrix): pt_ratio = pT2/pT1, |eta|, dphi from the tau phis, met from met_x/y, "
                                          "met_sig = met / sqrt(covxx + covyy), pt_vis = |pT1 + pT2|, pt_tt = |pT1 + pT2 + MET|, dphi_met_tt between MET and the visible pair",
                                 order="importance_sorted of ztautau_bdt.json (xgboost gain importance, mean over the 5 folds)",
                                 binning=f"{N_BINS} equal bins per continuous input over the stated range, overflow in the last bin (underflow, if any, in the first); "
                                         "njets one bin per value 0..6 (>= 6 in the last); decay modes one bin per category [0, 1, 10, 11]",
                                 normalisation="`signal` / `fakes` = weighted counts / sum (each sums to 1); `*_raw` weighted counts; `*_n` unweighted entries; "
                                               "`signal_unit` = signal_n / sum (unit weights, never negative: the signed genWeight leaves a sparsely populated "
                                               "tail bin of `signal` slightly negative where the negative-weight events outnumber the positive ones)",
                                 anchors=["z-tautau/docs/09-bdt.md (signal at small dR)", "z-tautau/ztautau/config.py BDT_FEATURES", "presentation/data/ztautau_bdt.json importance_sorted"]),
        "features": features,
        "dm_labels": DM_LABELS,
        "n_signal": int(len(Xs)), "sumw_signal": float(ws.sum()),
        "n_fakes": int(len(Xb)), "sumw_fakes": float(wb.sum()),
        "C_osss_inclusive": C,
        "feature_order": order,
    }
    return out


def verify(d, ck: Checker):
    D9 = "z-tautau/docs/09-bdt.md"
    add_ztautau_path()
    from ztautau import config
    f = {row["name"]: row for row in d["features"]}
    ck.check_true("features == config.BDT_FEATURES (as a set)", set(f) == set(config.BDT_FEATURES) and len(d["features"]) == len(config.BDT_FEATURES), "config.py", str(sorted(f)))
    frozen = load_json(BDT_FROZEN)
    ck.check_true("order == importance_sorted of ztautau_bdt.json", [r["name"] for r in d["features"]] == [k for k, _ in frozen["importance_sorted"]], "ztautau_bdt.json")
    ck.check("importance values == ztautau_bdt.json", [r["importance"] for r in d["features"]], [v for _, v in frozen["importance_sorted"]], 1e-12, "ztautau_bdt.json")
    for row in d["features"]:
        ck.check(f"{row['name']}: signal and fakes shapes sum to 1", [sum(row["signal"]), sum(row["fakes"])], [1.0, 1.0], 1e-6, "normalisation")
        nb = len(row["categories"]) if "categories" in row else len(row["edges"]) - 1
        ck.check_true(f"{row['name']}: {nb} bins, all arrays that long", all(len(row[k]) == nb for k in ("signal", "fakes", "signal_raw", "fakes_raw")), detail=str(nb))
    ck.check_true("signal dR median < fakes dR median (docs/09: signal at small dR)", f["dr_tt"]["median_signal"] < f["dr_tt"]["median_fakes"], D9,
                  f"signal {f['dr_tt']['median_signal']:.3f}, fakes {f['dr_tt']['median_fakes']:.3f}")
    ck.check_true("n_signal > 5000", d["n_signal"] > 5000, detail=str(d["n_signal"]))
    ck.check_true("every symbol is MathTex without words", all("\\" in r["symbol"] or "(" in r["symbol"] for r in d["features"]))
    ck.check_true("dm_labels for 0, 1, 10, 11", set(d["dm_labels"]) == {"0", "1", "10", "11"})
    ck.check_true("DM categories [0, 1, 10, 11]", f["t1_dm"]["categories"] == [0, 1, 10, 11] and f["t2_dm"]["categories"] == [0, 1, 10, 11], "config.TAU_DMS")
    ck.check("fakes sum of weights 13926 +- 50 (AR x FF x inclusive C, MC not subtracted)", d["sumw_fakes"], 13926, 50, "ztautau_bdt.json score fakes 13569 after MC subtraction (357 subtracted)", soft=True)
    ck.check("n_fakes == 120805 AR events", d["n_fakes"], 120805, 0, "docs/05 contamination table")
    neg = {r["name"]: [(i, round(v, 5)) for i, v in enumerate(r["signal"]) if v < 0] for r in d["features"] if min(r["signal"]) < 0}
    ck.check("no negative normalised signal bin (signed genWeight; informational, `signal_unit` is the non-negative shape)", len(neg), 0, 0,
             "aMC@NLO negative weights: " + str(neg), soft=True)
    ck.check_true("signal_unit shapes sum to 1 and are non-negative", all(abs(sum(r["signal_unit"]) - 1) < 1e-6 and min(r["signal_unit"]) >= 0 for r in d["features"]), "consistency")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_bdt_inputs.json")
    args = ap.parse_args()
    finalize(args, build, verify, "ztautau_bdt_inputs")


if __name__ == "__main__":
    main()
