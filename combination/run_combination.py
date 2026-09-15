#!/usr/bin/env python
"""Combine the Z -> mu mu and Z -> tau_h tau_h cross sections and write every output.

    source ../setup.sh
    python run_combination.py                # checks, baseline, variations, plots, JSON
    python run_combination.py --check        # input validation and closure tests only
    python run_combination.py --no-plots     # numbers only (no matplotlib)

Writes output/combination_result.json and output/plots/*.{pdf,png}; `result.md` and the docs quote
those numbers. The physics is in comb/ -- this file only orchestrates and reports.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np

from comb import blue, inputs, likelihood, model, ratio

OUT = HERE / "output"
REFERENCE = "sigma(pp -> Z/gamma* -> ll, 60 < m < 120 GeV) per lepton flavour, assuming lepton universality"


# ------------------------------------------------------------------------------- checks
def check(verbose: bool = True) -> list[str]:
    """Validate the inputs against the channels' published numbers and run the closure tests."""
    notes = []

    def ok(label, got, want, tol, unit=""):
        good = abs(got - want) <= tol
        notes.append(f"{'pass' if good else 'FAIL'}  {label}: {got:.6g}{unit} (published {want:.6g}{unit})")
        if not good:
            raise AssertionError(notes[-1])

    mm_fit = inputs.load_mumu("shapefit")
    ok("z-mumu shape fit mu_Z", mm_fit.mu, 0.990297, 1e-6)
    ok("z-mumu shape fit sigma(60-120)", mm_fit.sigma, 1934.97, 0.05, " pb")
    mm = inputs.load_mumu("counting")
    ok("z-mumu counting mu_Z", mm.mu, 0.993539, 1e-6)
    ok("z-mumu counting sigma_fid", mm.extra["sigma_fid_pb"], 794.38, 0.02, " pb")
    ok("z-mumu counting luminosity impact", mm.groups["Luminosity"], 0.0120375, 1e-7)

    tt = inputs.load_tautau("nominal")
    ok("z-tautau mu_Z", tt.mu, 1.1596, 1e-4)
    ok("z-tautau sigma(60-120)", tt.sigma, 2255.29, 0.05, " pb")
    ok("z-tautau tau ID impact", tt.groups["Tau ID"], 0.125475, 1e-6)
    tt_mc = inputs.load_tautau("mcsub")
    ok("z-tautau mcsub mu_Z", tt_mc.mu, 1.20751, 1e-5)

    # every grouped-impact category must have a correlation assigned (model.build raises otherwise)
    spec = model.build({"mumu": mm, "tautau": tt})
    notes.append(f"pass  every category mapped: {len(spec.sources)} sources built from "
                 f"{len({k for c in spec.channels for k in c.groups})} categories")

    # closure 1: one channel through the combination machinery returns that channel
    for name, ch in (("mumu", mm), ("tautau", tt)):
        r1 = blue.single(spec, name)
        ok(f"closure: BLUE({name}) value", r1.value, ch.sigma, 1e-6, " pb")
        ok(f"closure: BLUE({name}) error", r1.error, ch.sigma_err_total, 1e-6, " pb")

    # closure 2: two identical measurements. Uncorrelated -> error / sqrt(2) and equal weights;
    # with every systematic correlated -> only the statistical part gains.
    def twin_spec(**kw):
        clone = mm.__class__(**{**mm.__dict__, "name": "tautau"})
        return model.build({"mumu": mm, "tautau": clone}, **kw)

    r2 = blue.combine(twin_spec(rho_override=0.0), iterate=False)
    ok("closure: identical inputs, rho=0 -> error/sqrt(2)", r2.error, mm.sigma_err_total / 2 ** 0.5, 1e-6, " pb")
    ok("closure: identical inputs, rho=0 -> equal weights", r2.weights["mumu"], 0.5, 1e-9)
    r3 = blue.combine(twin_spec(rho_override=1.0), iterate=False)
    expect = (mm.sigma_err_total ** 2 - mm.sigma_stat ** 2 / 2) ** 0.5
    ok("closure: identical inputs, systematics rho=1 -> only stat averages", r3.error, expect, 1e-6, " pb")

    # closure 3: the profile likelihood must reproduce BLUE when the errors are symmetric
    rb = blue.combine(spec, iterate=False)
    rl = likelihood.combine(spec, asymmetric=False, scan_points=3)
    ok("closure: likelihood vs BLUE, central value", rl.value, rb.value, 1e-3, " pb")
    ok("closure: likelihood vs BLUE, uncertainty", rl.error, rb.error, 2e-2 * rb.error, " pb")

    if verbose:
        print("\n".join(notes))
    return notes


# ------------------------------------------------------------------------------- variations
VARIATIONS = {
    "baseline": dict(),
    "tautau_mcsub": dict(tautau="mcsub",
                         why="z-tautau fake factor with genuine-tau subtraction (its open issue 1)"),
    "mumu_shapefit": dict(mumu="shapefit",
                          why="z-mumu 30-bin shape fit instead of the reviewed counting extraction"),
    "sigmodel_correlated": dict(build=dict(correlate_sigmodel=True),
                                why="generator comparison treated as correlated between channels"),
    "acc_scale_decorrelated": dict(build=dict(correlate_acc_scale=False),
                                   why="acceptance QCD-scale term uncorrelated between channels"),
    "rho_zero": dict(build=dict(rho_override=0.0), why="all systematics uncorrelated (naive)"),
    "rho_one": dict(build=dict(rho_override=1.0), why="all systematics fully correlated"),
}


def run_variation(mumu="counting", tautau="nominal", build=None, **_):
    chans = inputs.load_channels(mumu=mumu, tautau=tautau)
    spec = model.build(chans, **(build or {}))
    return chans, spec, blue.combine(spec)


# ------------------------------------------------------------------------------- reporting
def fmt(value, error, unit=" pb", digits=1):
    return f"{value:.{digits}f} +- {error:.{digits}f}{unit}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="input validation and closure tests only")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    print("=" * 78)
    print("  Input validation and closure tests")
    print("=" * 78)
    notes = check()
    if args.check:
        return

    chans, spec, res = run_variation()
    mm, tt = chans["mumu"], chans["tautau"]
    solo = {name: blue.single(spec, name) for name in spec.order}
    lik = likelihood.combine(spec, asymmetric=True)
    lik_sym = likelihood.combine(spec, asymmetric=False, scan_points=3)
    rat = ratio.compute(spec)

    print("\n" + "=" * 78)
    print(f"  {REFERENCE}")
    print("=" * 78)
    for c in spec.channels:
        print(f"  {c.name:>7s} ({c.variant:<9s}) mu_Z = {c.mu:6.4f}   sigma = "
              f"{fmt(c.sigma, c.sigma_err_total)}   prediction {c.sigma_pred:.1f} pb")
    print(f"\n  COMBINED  sigma = {fmt(res.value, res.error)}   ({100 * res.rel:.2f} %)")
    g = res.group_breakdown()
    print(f"            = {res.value:.0f} +- {g['statistical']:.1f} (stat) +- {g['other systematic']:.1f} "
          f"(syst) +- {g['acceptance']:.1f} (acc) +- {g['luminosity']:.1f} (lumi) pb")
    print(f"  weights   : " + ", ".join(f"{k} {v:+.4f}" for k, v in res.weights.items()))
    print(f"  chi2/ndf  : {res.chi2:.2f}/{res.ndf}   p = {res.pvalue:.3f}"
          f"   rho(mumu,tautau) = {res.correlation[0, 1]:.3f}")
    print(f"  gain over mu mu alone: {solo['mumu'].error:.2f} -> {res.error:.2f} pb "
          f"({100 * (1 - res.error / solo['mumu'].error):+.2f} %)")
    print(f"  profile likelihood (asymmetric): {lik.value:.1f} +{lik.error_up:.1f} -{lik.error_down:.1f} pb")
    print(f"  profile likelihood (symmetric) : {lik_sym.value:.1f} +- {lik_sym.error:.1f} pb  "
          f"(BLUE {res.value:.1f} +- {res.error:.1f})")

    print("\n  uncertainty breakdown of the combined value [pb]:")
    for name, value in res.breakdown.items():
        print(f"    {name:<38s} {value:7.2f}   ({100 * value / res.value:5.2f} %)")

    print(f"\n  lepton universality:  R = sigma(tautau)/sigma(mumu) = {rat.value:.3f} +- {rat.error:.3f}"
          f"   ({rat.z_from_unity:+.2f} sigma from 1, p = {rat.pvalue:.3f})")

    print("\n  variations:")
    variations = {}
    for key, cfg in VARIATIONS.items():
        why = cfg.pop("why", "the baseline")
        _, vspec, vres = run_variation(**cfg)
        cfg["why"] = why
        variations[key] = {"value": vres.value, "error": vres.error, "why": why,
                           "weights": vres.weights, "chi2": vres.chi2, "pvalue": vres.pvalue,
                           "shift": vres.value - res.value}
        print(f"    {key:<24s} {vres.value:7.1f} +- {vres.error:5.1f} pb "
              f"({vres.value - res.value:+6.1f})   {why}")

    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "what": REFERENCE,
        "combined": {"sigma_pb": res.value, "total_pb": res.error, "relative": res.rel,
                     "stat_pb": g["statistical"], "syst_pb": g["other systematic"],
                     "acc_pb": g["acceptance"], "lumi_pb": g["luminosity"],
                     "chi2": res.chi2, "ndf": res.ndf, "pvalue": res.pvalue,
                     "weights": res.weights, "iterations": res.iterations,
                     "correlation_mumu_tautau": float(res.correlation[0, 1]),
                     "breakdown_pb": res.breakdown},
        "channels": {c.name: {"variant": c.variant, "mu": c.mu, "mu_err_up": c.mu_err_up,
                              "mu_err_down": c.mu_err_down, "mu_stat": c.mu_stat,
                              "sigma_pb": c.sigma, "sigma_err_pb": c.sigma_err_total,
                              "sigma_pred_pb": c.sigma_pred, "gof_p": c.gof_p,
                              "groups": c.groups, "acceptance": c.acc,
                              "alone_pb": solo[c.name].value, "alone_err_pb": solo[c.name].error,
                              "provenance": c.provenance, "extra": c.extra} for c in spec.channels},
        "sources": [{"name": s.name, "kind": s.kind, "rho": s.rho, "sizes_pb": s.sizes}
                    for s in spec.sources],
        "likelihood": {"value_pb": lik.value, "err_up_pb": lik.error_up, "err_down_pb": lik.error_down,
                       "symmetric_pb": lik_sym.value, "symmetric_err_pb": lik_sym.error,
                       "valid": lik.valid},
        "ratio": {"value": rat.value, "error": rat.error, "z": rat.z_from_unity,
                  "pvalue": rat.pvalue, "breakdown": rat.breakdown},
        "variations": variations,
        "checks": notes,
    }
    (OUT / "combination_result.json").write_text(json.dumps(payload, indent=1, default=float))
    print(f"\n  wrote {(OUT / 'combination_result.json').relative_to(HERE)}")
    from comb import report
    report.write_result_md(json.loads((OUT / "combination_result.json").read_text()), HERE / "result.md")
    print(f"  wrote result.md")

    if not args.no_plots:
        from comb import plots
        made = plots.make_all(spec, res, solo, lik, rat, variations, OUT / "plots")
        print(f"  wrote {len(made)} figures to {(OUT / 'plots').relative_to(HERE)}")


if __name__ == "__main__":
    main()
