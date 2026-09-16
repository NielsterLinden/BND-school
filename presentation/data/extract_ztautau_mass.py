#!/usr/bin/env python
"""Freeze the di-tau mass shapes of the simulated genuine tau_h tau_h signal in the Z -> tautau signal region: m_vis, the
MET-likelihood mass m_tt and the collinear mass m_col (defined events only), plus the docs/04 performance table.

    source setup.sh && python presentation/data/extract_ztautau_mass.py [--json PATH] [--check-only] [--weighting auto|unit|weighted]

Reads (read-only, uproot through ztautau.analysis) the Drell-Yan ntuples of $BND_TAUTAU_CACHE/ntuples_v1 restricted to
70 < m_LHE < 110 GeV, both legs genuine tau_h (genPartFlav 5), SR cuts (OS, both Tight, pT > 40).

Three variants are stored, `signal` is the one the deck draws:
  * `signal_inclusive_unit`      inclusive aMC@NLO sample DY_NLO alone, unit weights (one entry per simulated event);
  * `signal_inclusive_weighted`  DY_NLO alone, sigma x L / sumw x genWeight (analysis.normalisation, sign kept);
  * `signal_stitched`            DY_NLO + DY_0J/1J/2J stitched per LHE_NpNLO bin (analysis.dy_norm) x genWeight -- the
                                 analysis weighting, kept for the record: LHE 0-jet events carry the largest stitched weight,
                                 which puts a shoulder at 100-120 GeV in m_tt and moves the m_vis mode to 80-90 GeV.
`--weighting auto` (default) takes the inclusive weighted shape if it passes the single-peak test below, else the inclusive
unit-weight shape (the aMC@NLO sample has 25 % negative-weight events, so with 5.4 k events the signed weighted histogram
is noisy: one bin rises by 6.6 %). The chosen variant is recorded in `signal_weighting`.
Writes presentation/data/ztautau_mass.json. Anchors: docs/04-ditau-mass.md performance table, config.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATASET_TAUTAU, TAUTAU_NTUPLES, VERSION_TAUTAU, ZTAUTAU, Checker, add_ztautau_path,  # noqa: E402
                             finalize, provenance, standard_args)

M_Z = 91.1876
EDGES = np.arange(0.0, 200.0 + 1e-9, 5.0)
M_LHE_LO, M_LHE_HI = 70.0, 110.0
INCLUSIVE = "DY_NLO"
SINGLE_PEAK_TOL = 0.03          # one bin may rise by less than this fraction between the peak and the end of the tested range
SINGLE_PEAK_END = {"m_tt": 130.0, "m_vis": 110.0}    # the bin starting here is the last one compared
# docs/04-ditau-mass.md, "Performance on Z->tautau simulation" (genuine tau_h tau_h in the SR, 70 < m_LHE < 110)
PERFORMANCE = {"m_vis": {"scale": 0.80, "core_res": 0.131, "width68": 0.183, "defined": 1.0, "source": "docs/04 table row 'm_vis'"},
               "m_col": {"scale": 1.05, "core_res": 0.16, "width68": None, "defined": 0.38, "source": "docs/04 table row 'm_col (defined events only)'"},
               "m_tt": {"scale": 0.99, "core_res": 0.109, "width68": 0.160, "defined": 1.0, "source": "docs/04 table row 'm_tt, posterior median, 1/m^2 prior (used)'"},
               "m_tt_mean_noprior": {"scale": 1.04, "core_res": 0.119, "width68": 0.191, "defined": 1.0, "source": "docs/04 table"},
               "m_tt_map": {"scale": 0.91, "core_res": 0.091, "width68": 0.141, "defined": 1.0, "source": "docs/04 table"}}


def wquantile(x, w, q):
    """Weighted quantile with signed weights: the first point where the normalised cumulative sum of the sorted weights
    reaches q (the cumulative sum is not monotone with negative aMC@NLO weights, so no searchsorted)."""
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    cw /= cw[-1]
    return float(x[o][min(int(np.argmax(cw >= q)), len(x) - 1)])


def stats(m, mlhe, w, unit=False):
    ww = np.ones_like(w) if unit else w
    r = m / mlhe
    med = wquantile(r, ww, 0.5)
    q25, q75 = wquantile(r, ww, 0.25), wquantile(r, ww, 0.75)
    q16, q84 = wquantile(r, ww, 0.16), wquantile(r, ww, 0.84)
    return {"median_over_mlhe": med, "core_res_iqr_half": 0.5 * (q75 - q25), "width68_half": 0.5 * (q84 - q16),
            "median_gev": wquantile(m, ww, 0.5), "median_over_mZ": wquantile(m, ww, 0.5) / M_Z}


def single_peak(counts, edges, end_gev, tol=SINGLE_PEAK_TOL):
    """(ok, detail): from the peak bin to the bin starting at `end_gev` the counts fall monotonically, except that one
    bin may rise by less than `tol` (relative to the previous bin)."""
    c, e = np.asarray(counts, dtype=float), np.asarray(edges, dtype=float)
    ip, iend = int(np.argmax(c)), int(np.searchsorted(e, end_gev))
    seg = c[ip:iend + 1]
    rises = [(float(e[ip + i + 1]), float(seg[i + 1] / seg[i] - 1.0)) for i in range(len(seg) - 1) if seg[i + 1] > seg[i]]
    ok = len(rises) == 0 or (len(rises) == 1 and rises[0][1] < tol)
    detail = f"peak bin {e[ip]:.0f}-{e[ip + 1]:.0f}, tested to {end_gev:.0f}: " + ("no rise" if not rises else "rises " + ", ".join(f"{r[1]:+.1%} at {r[0]:.0f}" for r in rises))
    return ok, detail


def load_parts(keys, weighting):
    """Concatenated arrays of the selected events of the DY ntuples with the requested per-event weights."""
    from ztautau import analysis
    parts = []
    for key in keys:
        d, _ = analysis.load(key)
        r = analysis.regions(d, is_mc=True)
        sel = r["SR"] & (d["t1_genflav"] == 5) & (d["t2_genflav"] == 5) & (d["gen_mll_lhe"] > M_LHE_LO) & (d["gen_mll_lhe"] < M_LHE_HI) \
            & (d["gen_lhe_flavour"] == 15)
        gw = d["genWeight"][sel].astype(np.float64)
        if weighting == "stitched":
            w = analysis.dy_norm(key, d)[sel] * gw
        elif weighting == "weighted":
            w = analysis.normalisation(key) * gw
        else:
            w = np.ones_like(gw)
        parts.append({k: d[k][sel].astype(np.float64) for k in ("m_vis", "m_tt", "m_col", "gen_mll_lhe", "gen_fid")} | {"w": w})
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def shapes(a):
    """The `signal` block of one weighting: histograms, statistics, fractions."""
    n = len(a["w"])
    col_ok = a["m_col"] > 0
    frac_col = float(a["w"][col_ok].sum() / a["w"].sum())
    hists = {}
    for name in ("m_vis", "m_tt", "m_col"):
        m = a[name] if name != "m_col" else a["m_col"][col_ok]
        w = a["w"] if name != "m_col" else a["w"][col_ok]
        h, _ = np.histogram(m, bins=EDGES, weights=w)
        h2, _ = np.histogram(m, bins=EDGES, weights=w ** 2)
        hu, _ = np.histogram(m, bins=EDGES)
        hists[name] = {"counts": h.tolist(), "var": h2.tolist(), "unit_counts": hu.tolist(), "overflow": float(w[m >= EDGES[-1]].sum()),
                       "total": float(w.sum())}
    args_of = lambda name: ((a[name], a["gen_mll_lhe"], a["w"]) if name != "m_col" else (a["m_col"][col_ok], a["gen_mll_lhe"][col_ok], a["w"][col_ok]))
    perf = {name: stats(*args_of(name)) for name in ("m_vis", "m_tt", "m_col")}
    perf_unit = {name: stats(*args_of(name), unit=True) for name in ("m_vis", "m_tt", "m_col")}
    return {"m_vis": hists["m_vis"]["counts"], "m_tt": hists["m_tt"]["counts"], "m_col": hists["m_col"]["counts"],
            "hists": hists, "n_events": int(n), "sumw": float(a["w"].sum()),
            "frac_col_defined": frac_col, "frac_col_defined_unit": float(col_ok.mean()),
            "median_over_mZ": {k: perf[k]["median_over_mZ"] for k in ("m_vis", "m_tt", "m_col")},
            "median_over_mlhe": {k: perf[k]["median_over_mlhe"] for k in ("m_vis", "m_tt", "m_col")},
            "core_res_iqr_half": {k: perf[k]["core_res_iqr_half"] for k in ("m_vis", "m_tt", "m_col")},
            "width68_half": {k: perf[k]["width68_half"] for k in ("m_vis", "m_tt", "m_col")},
            "median_gev": {k: perf[k]["median_gev"] for k in ("m_vis", "m_tt", "m_col")},
            "unit_weights": {"median_over_mlhe": {k: perf_unit[k]["median_over_mlhe"] for k in ("m_vis", "m_tt", "m_col")},
                             "core_res_iqr_half": {k: perf_unit[k]["core_res_iqr_half"] for k in ("m_vis", "m_tt", "m_col")}},
            "fid_fraction": float(a["w"][a["gen_fid"] > 0].sum() / a["w"].sum()),
            "peak_bin_start": {k: float(EDGES[int(np.argmax(hists[k]["counts"]))]) for k in ("m_vis", "m_tt", "m_col")},
            "single_peak": {k: dict(zip(("ok", "detail"), single_peak(hists[k]["counts"], EDGES, SINGLE_PEAK_END[k]))) for k in ("m_tt", "m_vis")}}


def build(weighting="auto"):
    add_ztautau_path()
    from ztautau import analysis, config
    stitch_keys = analysis.dy_stitch_keys()
    variants = {"inclusive_unit": shapes(load_parts([INCLUSIVE], "unit")),
                "inclusive_weighted": shapes(load_parts([INCLUSIVE], "weighted")),
                "stitched": shapes(load_parts(stitch_keys, "stitched"))}
    for k, v in variants.items():
        print(f"[ztautau_mass] {k:19s} n {v['n_events']:5d}  m_tt: {v['single_peak']['m_tt']['detail']}  |  m_vis: {v['single_peak']['m_vis']['detail']}")
    if weighting == "auto":
        w_ok = variants["inclusive_weighted"]["single_peak"]
        chosen = "inclusive_weighted" if (w_ok["m_tt"]["ok"] and w_ok["m_vis"]["ok"]) else "inclusive_unit"
    else:
        chosen = {"unit": "inclusive_unit", "weighted": "inclusive_weighted"}[weighting]
    print(f"[ztautau_mass] `signal` = {chosen}")
    weights_doc = {"inclusive_unit": "unit weights (one entry per simulated event of DY_NLO; the 25 % negative-weight aMC@NLO events count +1)",
                   "inclusive_weighted": "sigma x L / sumw(DY_NLO) x genWeight, sign kept (analysis.normalisation); no scale factors, no pileup, no prefiring",
                   "stitched": "stitched normalisation (sigma x L x f_j / sum_j, analysis.dy_norm) x genWeight, sign kept; no scale factors, no pileup, no prefiring"}
    out = {
        "provenance": provenance("extract_ztautau_mass.py", [TAUTAU_NTUPLES / f"{k}.root" for k in stitch_keys] + [ZTAUTAU / "docs" / "04-ditau-mass.md"],
                                 dataset=DATASET_TAUTAU, version=VERSION_TAUTAU,
                                 samples={"signal": [INCLUSIVE], "signal_inclusive_unit": [INCLUSIVE], "signal_inclusive_weighted": [INCLUSIVE], "signal_stitched": stitch_keys},
                                 selection=f"SR (OS, both DeepTau Tight, pT > 40, |eta| < 2.1), t1_genflav == 5 and t2_genflav == 5, LHE tautau, {M_LHE_LO} < m_LHE < {M_LHE_HI} GeV",
                                 weights={"signal": weights_doc[chosen], **{f"signal_{k}": v for k, v in weights_doc.items()}},
                                 why_inclusive_only=("the analysis weighting (signal_stitched) gives the LHE 0-jet events of the inclusive sample the largest "
                                                     "stitched weight (Z at rest, back-to-back tau_h with both pT > 40 -> m_vis >= 80 GeV), which puts a shoulder "
                                                     "at 100-120 GeV in m_tt and moves the m_vis mode to 80-90 GeV; the deck draws one clean Z -> tautau peak, so "
                                                     "`signal` is the inclusive aMC@NLO sample alone with a uniform normalisation. The signed genWeight version "
                                                     "(signal_inclusive_weighted) has 25 % negative-weight events and with 5.4 k events its m_tt still rises by "
                                                     "6.6 % in one bin, so the single-peak test selects the unit-weight histogram (signal_inclusive_unit)"),
                                 single_peak_test=(f"counts fall monotonically from the peak bin to the bin starting at {SINGLE_PEAK_END['m_tt']:.0f} GeV (m_tt) / "
                                                   f"{SINGLE_PEAK_END['m_vis']:.0f} GeV (m_vis); one bin may rise by < {SINGLE_PEAK_TOL:.0%}"),
                                 binning="40 bins of 5 GeV, 0-200 GeV, overflow NOT in the last bin (stored separately)",
                                 mass_definition={"m_vis": "invariant mass of the two tau_h", "m_col": "collinear approximation, undefined (-1) when the 2x2 system is singular or x outside (0,1]",
                                                  "m_tt": "posterior median of m_vis / sqrt(x1 x2) on a 60 x 60 grid, MET transfer function with the event MET covariance, "
                                                          "flat phase space in x on [m_vis,i^2/m_tau^2, 1], 1/m^2 prior (ztautau/mass.py likelihood_mass)"},
                                 anchors=["z-tautau/docs/04-ditau-mass.md performance table", "z-tautau/ztautau/config.py (MASS_GRID_N 60, MASS_PRIOR_POW 2)"]),
        "edges": EDGES.tolist(),
        "m_Z": M_Z,
        "signal_weighting": chosen,
        "signal": variants[chosen],
        "signal_inclusive_unit": variants["inclusive_unit"],
        "signal_inclusive_weighted": variants["inclusive_weighted"],
        "signal_stitched": variants["stitched"],
        "performance": PERFORMANCE,
        "grid_n": int(config.MASS_GRID_N), "prior_pow": float(config.MASS_PRIOR_POW), "m_tau": float(config.M_TAU),
    }
    return out


def verify(d, ck: Checker):
    D4 = "z-tautau/docs/04-ditau-mass.md"
    s = d["signal"]
    ck.check("41 edges 0..200 step 5", d["edges"], np.arange(0, 201, 5), 0, "brief")
    ck.check_true("40 bins in each shape", all(len(s[k]) == 40 for k in ("m_vis", "m_tt", "m_col")))
    ck.check_true("signal == the recorded variant", d["signal_weighting"] in ("inclusive_unit", "inclusive_weighted") and s["m_tt"] == d[f"signal_{d['signal_weighting']}"]["m_tt"],
                  "consistency", d["signal_weighting"])
    ck.check("performance constants (docs/04)", [d["performance"]["m_vis"]["scale"], d["performance"]["m_vis"]["core_res"], d["performance"]["m_col"]["scale"],
                                                 d["performance"]["m_col"]["core_res"], d["performance"]["m_col"]["defined"], d["performance"]["m_tt"]["scale"], d["performance"]["m_tt"]["core_res"]],
             [0.80, 0.131, 1.05, 0.16, 0.38, 0.99, 0.109], 0, D4)
    ck.check("median m_vis / m_LHE within 0.06 of 0.80", s["median_over_mlhe"]["m_vis"], 0.80, 0.06, D4, soft=True)
    ck.check("median m_tt / m_LHE within 0.06 of 0.99", s["median_over_mlhe"]["m_tt"], 0.99, 0.06, D4, soft=True)
    ck.check("median m_col / m_LHE within 0.06 of 1.05", s["median_over_mlhe"]["m_col"], 1.05, 0.06, D4, soft=True)
    ck.check("median m_vis / m_Z within 0.06 of 0.80", s["median_over_mZ"]["m_vis"], 0.80, 0.06, D4 + " (brief: median_over_mZ)", soft=True)
    ck.check("median m_tt / m_Z within 0.06 of 0.99", s["median_over_mZ"]["m_tt"], 0.99, 0.06, D4 + " (brief: median_over_mZ)", soft=True)
    ck.check("core resolution m_vis 13.1 %, m_tt 10.9 % (+- 1.5 %)", [s["core_res_iqr_half"]["m_vis"], s["core_res_iqr_half"]["m_tt"]], [0.131, 0.109], 0.015, D4, soft=True)
    ck.check("m_col defined in 38 % (+- 3 %)", s["frac_col_defined"], 0.38, 0.03, D4, soft=True)
    ck.check_true("m_tt peak bin near m_Z (80-100 GeV)", 80 <= d["edges"][int(np.argmax(s["m_tt"]))] < 100, D4, f"peak bin starts at {d['edges'][int(np.argmax(s['m_tt']))]}")
    for name in ("m_tt", "m_vis"):
        ok, detail = single_peak(s[name], d["edges"], SINGLE_PEAK_END[name])
        ck.check_true(f"{name} shape has a single peak (monotone fall to {SINGLE_PEAK_END[name]:.0f} GeV, one rise < {SINGLE_PEAK_TOL:.0%} allowed)", ok, "deck owner", detail)
    ck.check_true("weighted median m_vis below weighted median m_tt", s["median_over_mZ"]["m_vis"] < s["median_over_mZ"]["m_tt"], D4,
                  f"{s['median_over_mZ']['m_vis']:.3f} < {s['median_over_mZ']['m_tt']:.3f}")
    ck.check_true("m_vis peak bin at or below m_tt peak bin", int(np.argmax(s["m_vis"])) <= int(np.argmax(s["m_tt"])), D4,
                  f"m_vis peak {s['peak_bin_start']['m_vis']:.0f}, m_tt peak {s['peak_bin_start']['m_tt']:.0f}")
    ck.check_true("more than 5000 simulated events", s["n_events"] > 5000, detail=str(s["n_events"]))
    ck.check("m_vis, m_tt shapes have the same normalisation", sum(s["m_vis"]) + s["hists"]["m_vis"]["overflow"], sum(s["m_tt"]) + s["hists"]["m_tt"]["overflow"], 1e-9, "consistency", rel=True)
    ck.check("m_col normalisation == frac_col_defined x total", sum(s["m_col"]) + s["hists"]["m_col"]["overflow"], s["frac_col_defined"] * s["sumw"], 1e-9, "consistency", rel=True)
    st = d["signal_stitched"]
    ck.check_true("stitched variant kept (4 samples, > 5000 events)", st["n_events"] > 5000 and len(d["provenance"]["samples"]["signal_stitched"]) == 4, "record", str(st["n_events"]))
    ck.check("stitched m_tt / m_LHE median within 0.06 of 0.99", st["median_over_mlhe"]["m_tt"], 0.99, 0.06, D4 + " (the docs/04 weighting)", soft=True)
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_mass.json")
    ap.add_argument("--weighting", choices=["auto", "unit", "weighted"], default="auto", help="which inclusive-only variant becomes `signal` (default: auto by the single-peak test)")
    args = ap.parse_args()
    finalize(args, lambda: build(args.weighting), verify, "ztautau_mass")


if __name__ == "__main__":
    main()
