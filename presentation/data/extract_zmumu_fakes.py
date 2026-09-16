#!/usr/bin/env python
"""Freeze the Z -> mumu fake-factor estimate (FF maps, SR template, closure, same-sign data).

    source setup.sh && python presentation/data/extract_zmumu_fakes.py [--json PATH] [--check-only]

Reads (read-only): z-mumu/output/v2/fakes.json (scripts/v2_3_control.py, zmumu/fakes.py) and
z-mumu/output/v2/histograms.pkl (`Data|SS|mass_fit|nominal`). Writes presentation/data/zmumu_fakes.json.
Anchors: output/v2/RESULTS_v2.md:116-119, docs/13-fake-factor.md, handoff.md:36.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import V2, Checker, finalize, load_json, provenance, standard_args  # noqa: E402

FAKES_JSON = V2 / "fakes.json"
HIST_PKL = V2 / "histograms.pkl"

DOC_NOTE = ("docs/13-fake-factor.md quoted 3 825 / 2 077 / 1 922 from an earlier pass; the numbers of this JSON are the "
            "current output/v2/fakes.json (3 870 +- 80; closure 2 094 +- 20 vs 1 932 +- 79), to which the doc was corrected on 2026-09-15.")


def build():
    f = load_json(FAKES_JSON)
    with open(HIST_PKL, "rb") as fh:
        hall = pickle.load(fh)
    t = np.array(f["templates"]["nominal"], dtype=float)                 # 120 x 0.5 GeV, 60-120 GeV
    tv = np.array(f["templates"]["nominal_var"], dtype=float)            # (b_err * F)^2: fully correlated between bins
    t60 = t.reshape(60, 2).sum(1)
    tv60 = np.sqrt(tv).reshape(60, 2).sum(1) ** 2
    edges60 = np.linspace(60.0, 120.0, 61).tolist()
    edges120 = np.linspace(60.0, 120.0, 121).tolist()
    ss = np.array(hall["Data|SS|mass_fit|nominal"], dtype=float)
    assert np.all(ss == np.round(ss))
    out = {
        "provenance": provenance("extract_zmumu_fakes.py", [FAKES_JSON, HIST_PKL],
                                 method=("fake factor FF = N_tight / N_anti-tight (tightId, 0.20 < pfRelIso04 < 1.0) measured in same-sign "
                                         "tag+probe pairs (tag tight, pT > 26, trigger-matched; probe pT > 20; m > 12 GeV), prompt MC subtracted; "
                                         "applied to opposite-sign tight + anti-tight events (60 < m < 120) through a two-template fit "
                                         "D = a P + b F (zmumu/fakes.py:fit_application)"),
                                 histograms={"template_sr": "fakes.json templates.nominal (120 x 0.5 GeV) and its 2-bin sum (60 x 1 GeV)",
                                             "ss_data_mass": "histograms.pkl Data|SS|mass_fit|nominal (60 x 1 GeV, prescale-weighted data)"},
                                 anchors=["z-mumu/output/v2/RESULTS_v2.md:116-119", "z-mumu/handoff.md:36", "z-mumu/docs/13-fake-factor.md"]),
        "pt_edges": f["pt_edges"], "eta_edges": f["eta_edges"],
        "ff": {"nominal": f["maps"]["nominal"], "nominal_err": f["maps"]["nominal_err"], "region_c": f["maps"]["region_c"],
               "anti_alt": f["maps"]["anti_alt"], "shape": "[pt bin][|eta| bin] = (6, 4)",
               "note": "region_c (single-muon + jet, trigger-biased) is a documented cross-check only; anti_alt = 0.15 < iso < 0.60"},
        "template_sr": {"edges": edges60, "counts": t60.tolist(), "var": tv60.tolist(), "total": float(t.sum()),
                        "edges_halfgev": edges120, "counts_halfgev": t.tolist(), "var_halfgev": tv.tolist(),
                        "var_note": "template_var = (sigma_b F)^2 is fully correlated between bins; the 1 GeV var is (sum of sqrt)^2"},
        "yields": {"sr_fakes": f["yields"]["sr_fakes"], "sr_fakes_stat": f["yields"]["sr_fakes_stat"],
                   "sr_data_tight_tight": f["yields"]["sr_data_tight_tight"], "ss_data_tight_tight": f["yields"]["ss_data_tight_tight"],
                   "os_over_ss": f["yields"]["os_over_ss_fakes"], "method_rel_unc": f["yields"]["method_rel_unc"],
                   "method_worst": f["yields"]["method_worst"], "prompt_norm_fitted": f["yields"]["prompt_norm_fitted"],
                   "variants": f["yields"]["variants"], "method_deviations": f["yields"]["method_deviations"],
                   "fraction_of_sr_data": f["yields"]["sr_fakes"] / f["yields"]["sr_data_tight_tight"]},
        "closure": {"predicted": f["closure"]["predicted_ss"], "predicted_err": f["closure"]["predicted_ss_err"],
                    "observed": f["closure"]["observed_ss"], "observed_err": f["closure"]["observed_ss_err"],
                    "edges_halfgev": edges120, "predicted_hist": f["closure"]["predicted_hist"], "observed_hist": f["closure"]["observed_hist"],
                    "note": "same-sign application region: FF x (tight + anti-tight) minus prompt MC vs same-sign tight-tight data minus prompt MC"},
        "ss_data_mass": {"edges": edges60, "counts": [int(v) for v in ss], "total": int(ss.sum()),
                         "note": "same-sign tight-tight data (SS region of the analysis), no prompt subtraction"},
        "fit_nominal": {k: f["fits"]["nominal"][k] for k in ("a", "b", "chi2", "ndf", "n_fake", "n_fake_err", "n_prompt_fitted", "n_prompt_mc")},
        "doc_note": DOC_NOTE,
    }
    return out


def verify(d, ck: Checker):
    R = "z-mumu/output/v2/RESULTS_v2.md:116-119"
    ff = np.array(d["ff"]["nominal"])
    ck.check_true("FF map shape (6, 4)", ff.shape == (6, 4), detail=str(ff.shape))
    ck.check_true("FF err map shape (6, 4)", np.array(d["ff"]["nominal_err"]).shape == (6, 4))
    y = d["yields"]
    ck.check("SR fakes 3869.96", y["sr_fakes"], 3869.96, 5e-3, "fakes.json yields.sr_fakes")
    ck.check("SR fakes stat 80.45", y["sr_fakes_stat"], 80.45, 5e-3, "fakes.json")
    ck.check("SR fakes rounded 3870 +- 80", [y["sr_fakes"], y["sr_fakes_stat"]], [3870, 80], 0.5, R)
    ck.check("fakes fraction of SR data 0.037%", 100 * y["fraction_of_sr_data"], 0.037, 5e-4, "z-mumu/docs/13-fake-factor.md:38")
    ck.check("SR data tight-tight", y["sr_data_tight_tight"], 10378567, 0, "z-mumu/handoff.md:35")
    ck.check("SS data tight-tight 2302", y["ss_data_tight_tight"], 2302, 0, "RESULTS_v2.md yields table")
    ck.check("os_over_ss 1.848", y["os_over_ss"], 1.848, 5e-4, "fakes.json")
    ck.check("method_rel_unc 0.190", y["method_rel_unc"], 0.190, 5e-4, "RESULTS_v2.md:118 (19%)")
    c = d["closure"]
    ck.check("closure predicted 2093.99 +- 19.71", [c["predicted"], c["predicted_err"]], [2093.99, 19.71], 5e-3, R)
    ck.check("closure observed 1932.28 +- 78.91", [c["observed"], c["observed_err"]], [1932.28, 78.91], 5e-3, R)
    ck.check("closure rounded 2094 +- 20 / 1932 +- 79", [c["predicted"], c["predicted_err"], c["observed"], c["observed_err"]], [2094, 20, 1932, 79], 0.5, R)
    t = d["template_sr"]
    ck.check("template 60-bin sum == total", sum(t["counts"]), t["total"], 1e-9, rel=True)
    ck.check("template 120-bin sum == total", sum(t["counts_halfgev"]), t["total"], 1e-9, rel=True)
    ck.check("template total == sr_fakes", t["total"], y["sr_fakes"], 1e-9, rel=True)
    ck.check_true("template has no negative bins", min(t["counts_halfgev"]) >= 0, detail=f"min {min(t['counts_halfgev']):.3f}")
    ck.check("closure predicted_hist sum == predicted", sum(c["predicted_hist"]), c["predicted"], 1e-9, rel=True)
    ck.check("closure observed_hist sum == observed", sum(c["observed_hist"]), c["observed"], 1e-9, rel=True)
    ck.check("ss_data_mass sum == 2302", d["ss_data_mass"]["total"], 2302, 0, "histograms.pkl Data|SS")
    ck.check("ss_data_mass counts sum == total", sum(d["ss_data_mass"]["counts"]), d["ss_data_mass"]["total"], 0)
    barrel = ff[:4, 0]
    ck.check_true("barrel FF cells (pT < 50, |eta| < 0.9) in [0.08, 0.13]", bool(np.all((barrel > 0.08) & (barrel < 0.13))), detail=str(np.round(barrel, 4)))
    ck.check("barrel FF cells 0.093 / 0.103 / 0.087 / 0.118", barrel, [0.093, 0.103, 0.087, 0.118], 5e-4, "plan A.2 / fakes.json maps.nominal")
    ck.check("FF map range 0.087 .. 0.81", [ff.min(), ff.max()], [0.087, 0.806], 5e-3, "fakes.json")
    ck.check("fit_nominal b == os_over_ss", d["fit_nominal"]["b"], y["os_over_ss"], 1e-9, rel=True)
    ck.check_true("edges 60 x 1 GeV", len(t["edges"]) == 61 and t["edges"][0] == 60.0 and t["edges"][-1] == 120.0)
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "zmumu_fakes.json")
    args = ap.parse_args()
    finalize(args, build, verify, "zmumu_fakes")


if __name__ == "__main__":
    main()
