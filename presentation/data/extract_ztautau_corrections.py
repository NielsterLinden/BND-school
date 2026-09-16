#!/usr/bin/env python
"""Freeze the corrections applied to the Z -> tautau simulation: the TauPOG tau_h ID scale factors and energy scales per decay
mode, a trigger-SF example curve, and the mean event weights on the fiducial signal in the signal region.

    source setup.sh && python presentation/data/extract_ztautau_corrections.py [--json PATH] [--check-only] [--no-means]

Reads (read-only): z-tautau/external/tau_pog_UL2016postVFP.json (TauPOG UL2016 postVFP, DeepTau2017v2p1), z-tautau/output/data/fakefactors.json
(AR contamination), z-tautau/output/data/yields.json (the DYtautau yield the weight chain must reproduce), and the four Drell-Yan
ntuples through ztautau.analysis (dy_norm stitching, pileup, L1 prefiring, corrections.id_sf / trigger_sf) for the means.
Writes presentation/data/ztautau_corrections.json. Anchors: docs/07-corrections-and-systematics.md, docs/05, docs/06, docs/09.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATASET_TAUTAU, TAUTAU_NTUPLES, VERSION_TAUTAU, ZTAUTAU, Checker, add_ztautau_path,  # noqa: E402
                             finalize, load_json, provenance, standard_args)

POG = ZTAUTAU / "external" / "tau_pog_UL2016postVFP.json"
FF_JSON = ZTAUTAU / "output" / "data" / "fakefactors.json"
YIELDS = ZTAUTAU / "output" / "data" / "yields.json"
DMS = ["0", "1", "10", "11"]
TRIG_EXAMPLE = {"pt": 50.0, "dm": "1"}
CURVE_PT = [40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 100.0]


def means_on_fiducial_signal():
    """Mean weight factors on the fiducial Z -> tautau (LHE tautau, gen_fid) signal-region events of the stitched DY samples,
    weighted by the generator weight (dy_norm x genWeight); plus the raw / corrected yields (corrected must equal RESULTS.md 3759)."""
    add_ztautau_path()
    from ztautau import analysis, corrections, samples
    keys = analysis.dy_stitch_keys()
    acc = {k: 0.0 for k in ("w0", "w_pu", "w_l1", "sf_id", "sf_trig", "all", "n")}
    w0_pos = 0.0
    per_sample = {}
    for key in keys:
        d, _ = analysis.load(key)
        r = analysis.regions(d, is_mc=True)
        comp = dict(analysis.mc_components(key))
        m = r["SR"] & comp[samples.SIGNAL]
        w0 = analysis.dy_norm(key, d)[m] * d["genWeight"][m].astype(np.float64)
        pu = analysis.pileup()(d["Pileup_nTrueInt"][m])
        l1 = d["L1PreFiringWeight_Nom"][m].astype(np.float64)
        sf_id = np.ones(int(m.sum()))
        sf_tr = np.ones(int(m.sum()))
        for i in (1, 2):
            dm, fl, ae = d[f"t{i}_dm"][m].astype(int), d[f"t{i}_genflav"][m].astype(int), np.abs(d[f"t{i}_eta"][m])
            sf_id = sf_id * corrections.id_sf(dm, fl, ae)
            sf_tr = sf_tr * corrections.trigger_sf(dm, d[f"t{i}_pt"][m], genflav=fl)
        allw = pu * l1 * sf_id * sf_tr
        full = analysis.weights(d, key)[m]
        assert np.allclose(w0 * allw, full), key
        per_sample[key] = {"n_events": int(m.sum()), "raw": float(w0.sum()), "corrected": float(full.sum())}
        acc["n"] += int(m.sum()); acc["w0"] += float(w0.sum()); w0_pos += float(np.abs(w0).sum())
        for name, arr in (("w_pu", pu), ("w_l1", l1), ("sf_id", sf_id), ("sf_trig", sf_tr), ("all", allw)):
            acc[name] += float((w0 * arr).sum())
    means = {"w_pu": acc["w_pu"] / acc["w0"], "w_l1": acc["w_l1"] / acc["w0"], "sf_id_both_legs": acc["sf_id"] / acc["w0"],
             "sf_trig_both_legs": acc["sf_trig"] / acc["w0"], "all": acc["all"] / acc["w0"]}
    return {**means, "n_events": acc["n"], "samples": keys, "per_sample": per_sample,
            "signal_yield": {"raw": acc["w0"], "corrected": acc["all"], "ratio": acc["all"] / acc["w0"]},
            "weighting": "generator weight (stitched dy_norm x genWeight, sign kept) over fiducial Z->tautau SR events of DY_NLO + DY_0J/1J/2J; "
                         "sf_trig = product of the two leg SFs (genuine legs), sf_id = product of the two Tight DM-binned ID SFs"}


def build(with_means=True):
    p = load_json(POG)
    ff = load_json(FF_JSON)
    y = load_json(YIELDS)
    tight = p["id_vsjet_dm"]["Tight"]
    trig = p["trigger_ditau_Tight"]
    grid = np.asarray(trig["pt"])

    def sf_at(dm, pt):
        i = int(np.clip(np.searchsorted(grid, pt, side="right") - 1, 0, len(grid) - 1))
        return {"pt": float(pt), "dm": int(dm), "sf": float(trig[str(dm)]["sf"][i]), "sf_err": float(trig[str(dm)]["sf_err"][i]),
                "eff_data": float(trig[str(dm)]["data"][i]), "eff_mc": float(trig[str(dm)]["mc"][i])}

    means = means_on_fiducial_signal() if with_means else None
    out = {
        "provenance": provenance("extract_ztautau_corrections.py", [POG, FF_JSON, YIELDS] + ([TAUTAU_NTUPLES / f"{k}.root" for k in (means["samples"] if means else [])]),
                                 dataset=DATASET_TAUTAU, version=VERSION_TAUTAU,
                                 pog_sources=p["sources"], pog_retrieved=p["retrieved"], tau_id=p["tau_id"], working_point="Tight (VSjet), VSe VVLoose, VSmu VLoose",
                                 anchors=["z-tautau/docs/07-corrections-and-systematics.md", "z-tautau/docs/05-fake-factors.md", "z-tautau/docs/06-cross-section.md",
                                          "z-tautau/docs/09-bdt.md", "z-tautau/output/RESULTS.md"]),
        "taupog": {
            "id_sf_tight_per_dm": {dm: {"value": tight[dm][0], "err": tight[dm][1]} for dm in DMS},
            "id_sf_medium_per_dm": {dm: {"value": p["id_vsjet_dm"]["Medium"][dm][0], "err": p["id_vsjet_dm"]["Medium"][dm][1]} for dm in DMS},
            "tes_per_dm": {dm: {"value": p["tes_dm"][dm][0], "err": p["tes_dm"][dm][1]} for dm in DMS},
            "vse_vvloose": p["vse"]["VVLoose"], "vsmu_vloose": p["vsmu"]["VLoose"],
            "trigger_sf_example": sf_at(TRIG_EXAMPLE["dm"], TRIG_EXAMPLE["pt"]),
            "trigger_sf_curve": {dm: {"pt": CURVE_PT, "sf": [sf_at(dm, pt)["sf"] for pt in CURVE_PT], "sf_err": [sf_at(dm, pt)["sf_err"] for pt in CURVE_PT],
                                      "eff_data": [sf_at(dm, pt)["eff_data"] for pt in CURVE_PT]} for dm in DMS},
            "id_pt_binned_check": p["id_vsjet_pt_check"],
            "note": "TauPOG UL2016 postVFP DeepTau2017v2p1 (external/tau_pog_UL2016postVFP.json); ID SF per DM for the Tight working point "
                    "(v3 nominal); trigger = di-tau trigger leg SF(pT, DM) = eff_data / eff_MC for Tight offline legs, applied to genuine legs only",
        },
        "means_on_fiducial_signal_in_sr": means,
        "constants": {
            "w_l1_mean_docs07": {"value": 0.990, "source": "docs/07 table: L1PreFiringWeight_Nom mean 0.990 for the signal"},
            "mc_subtraction_ar_fraction": {"value": ff["contamination"]["AR"]["fraction"], "source": "output/data/fakefactors.json contamination.AR (docs/05 table: 2.5 %)"},
            "signal_double_count_without_subtraction": {"value": 0.06, "source": "docs/05: ~6 % of the signal double counted without the MC subtraction"},
            "nonfiducial_fraction_of_selected_dy": {"value": 0.38, "source": "docs/06: selected Z/gamma*->tautau is 62 % fiducial; REVIEW.md 3.1 (38 %)"},
            "nonfiducial_fraction_from_yields": {"value": y["yields"]["DYtautau_nonfid"]["value"] / (y["yields"]["DYtautau"]["value"] + y["yields"]["DYtautau_nonfid"]["value"]),
                                                 "source": "output/data/yields.json: 2337 / (3759 + 2337)"},
            "fakes_fraction_medium": {"value": 0.80, "source": "docs/05, docs/09: 80 % fakes with DeepTau Medium (v1-v2.1)"},
            "fakes_fraction_tight": {"value": 0.64, "source": "docs/02, docs/05: 64 % fakes with DeepTau Tight (v3)"},
            "pileup_sigma_mb": {"value": 69.2, "source": "docs/07: minimum-bias cross section 69.2 mb, +-4.6 %"},
            "c_osss_bias_without_subtraction": {"value": 0.025, "source": "docs/05: C biased up by 2.5 % without the MC subtraction"},
        },
    }
    return out


def verify(d, ck: Checker):
    D7 = "z-tautau/docs/07-corrections-and-systematics.md"
    sf = d["taupog"]["id_sf_tight_per_dm"]
    ck.check("ID SF Tight per DM 0.902 / 0.892 / 0.938 / 0.812", [sf[dm]["value"] for dm in DMS], [0.902, 0.892, 0.938, 0.812], 5e-4, D7)
    ck.check("ID SF Tight err 0.126 / 0.053 / 0.148 / 0.149", [sf[dm]["err"] for dm in DMS], [0.126, 0.053, 0.148, 0.149], 5e-4, D7)
    tes = d["taupog"]["tes_per_dm"]
    ck.check("TES per DM 0.993 / 0.991 / 1.001 / 0.997", [tes[dm]["value"] for dm in DMS], [0.993, 0.991, 1.001, 0.997], 5e-4, D7)
    ck.check("TES err 0.009 / 0.007 / 0.007 / 0.016", [tes[dm]["err"] for dm in DMS], [0.009, 0.007, 0.007, 0.016], 5e-4, D7)
    med = d["taupog"]["id_sf_medium_per_dm"]
    ck.check("ID SF Medium (v2.1) 0.923 / 0.880 / 0.868 / 0.898", [med[dm]["value"] for dm in DMS], [0.923, 0.880, 0.868, 0.898], 5e-4, D7)
    ex = d["taupog"]["trigger_sf_example"]
    ck.check_true("trigger SF example in 0.5 .. 1.5", 0.5 < ex["sf"] < 1.5, D7 + " (0.94-0.97 at 50 GeV)", f"DM{ex['dm']} pT {ex['pt']}: {ex['sf']:.3f} +- {ex['sf_err']:.3f}")
    ck.check("trigger SF example ~ 0.94-0.97 at 50 GeV", ex["sf"], 0.955, 0.06, D7 + " (e.g. 0.94-0.97 at 50 GeV)", soft=True)
    c = d["constants"]
    ck.check("AR contamination 0.0248", c["mc_subtraction_ar_fraction"]["value"], 0.0248, 5e-4, "docs/05 (2.5 %)")
    ck.check("non-fiducial fraction from yields 0.38", c["nonfiducial_fraction_from_yields"]["value"], 0.38, 0.01, "docs/06 (62 % fiducial)")
    ck.check("w_l1 docs 0.990", c["w_l1_mean_docs07"]["value"], 0.990, 0, D7)
    m = d["means_on_fiducial_signal_in_sr"]
    if m is not None:
        ck.check("<w_L1> on the fiducial signal 0.990 +- 0.003", m["w_l1"], 0.990, 0.003, D7)
        ck.check("corrected fiducial signal yield == RESULTS.md 3759.14 +- 0.5", m["signal_yield"]["corrected"], 3759.14, 0.5, "output/RESULTS.md / yields.json DYtautau")
        ck.check("<all> == corrected / raw", m["all"], m["signal_yield"]["ratio"], 1e-9, "consistency", rel=True)
        ck.check_true("<SF_ID both legs> in 0.7 .. 0.95 (per-leg 0.81-0.94)", 0.7 < m["sf_id_both_legs"] < 0.95, D7, f"{m['sf_id_both_legs']:.4f}")
        ck.check_true("<w_PU> in 0.9 .. 1.1", 0.9 < m["w_pu"] < 1.1, "normalised pileup weights", f"{m['w_pu']:.4f}")
        ck.check_true("<SF_trig both legs> in 0.6 .. 1.1", 0.6 < m["sf_trig_both_legs"] < 1.1, D7 + " (turn-on: data efficiency 0.3-0.6 at 40 GeV)", f"{m['sf_trig_both_legs']:.4f}")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_corrections.json")
    ap.add_argument("--no-means", action="store_true", help="skip the ntuple pass (means_on_fiducial_signal_in_sr = null)")
    args = ap.parse_args()
    finalize(args, lambda: build(with_means=not args.no_means), verify, "ztautau_corrections")


if __name__ == "__main__":
    main()
