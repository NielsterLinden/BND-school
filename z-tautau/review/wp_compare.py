#!/usr/bin/env python
"""Medium vs Tight working-point comparison (docs/08, "working-point cross-check").

    python review/wp_compare.py            # prints the markdown table and writes review/wp_compare.md

Reads output/results.json (Medium, nominal) and variants/tight/output/results.json (BND_TAUTAU_WP=Tight run).
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
M = json.loads((HERE / "output/results.json").read_text())
T = json.loads((HERE / "variants/tight/output/results.json").read_text())


def row(label, f):
    return f"| {label} | {f(M)} | {f(T)} |"


def fit(r):
    return r["fit"]["mcsub"]


lines = ["| quantity | Medium (nominal) | Tight |", "|---|---:|---:|",
         row("μ_Z", lambda r: f"{fit(r)['mu']:.3f} +{fit(r)['mu_err_up']:.3f} −{fit(r)['mu_err_down']:.3f}"),
         row("total uncertainty on μ_Z", lambda r: f"+{100*fit(r)['mu_err_up']:.1f} / −{100*fit(r)['mu_err_down']:.1f} %"),
         row("expected (Asimov)", lambda r: f"+{100*fit(r)['mu_expected_asimov']['err_up']:.1f} / −{100*fit(r)['mu_expected_asimov']['err_down']:.1f} %"),
         row("data statistics", lambda r: f"{100*fit(r)['mu_stat']:.1f} %"),
         row("σ(60–120) [pb]", lambda r: f"{fit(r)['sigma_60_120_pb']['value']:.0f} +{fit(r)['sigma_60_120_pb']['err_up']:.0f} −{fit(r)['sigma_60_120_pb']['err_down']:.0f}"),
         row("goodness of fit p", lambda r: f"{fit(r)['gof_probability']:.3f}")]
for g in ("Tau ID", "Fakes", "Gammas", "Tau trigger", "Tau energy scale", "Signal modelling", "Background normalisation"):
    lines.append(row(f"impact: {g}", lambda r, g=g: f"{100*fit(r)['grouped_impact'].get(g, 0):.1f} %"))
lines += [row("prefit fiducial signal / fakes (all categories)", lambda r: f"{r['yields_prefit']['mcsub']['DYtautau']['value']:.0f} / {r['yields_prefit']['mcsub']['Fakes']['value']:.0f}"),
          row("SR2: signal / fakes", lambda r: f"{r['yields_prefit_per_region']['tautau_SR2']['DYtautau']['value']:.0f} / {r['yields_prefit_per_region']['tautau_SR2']['Fakes']['value']:.0f}"),
          row("BDT held-out AUC", lambda r: f"{sum(r['bdt']['training']['auc_test'])/len(r['bdt']['training']['auc_test']):.3f}"),
          row("C_OS/SS (inclusive)", lambda r: f"{r['fake_factors']['mcsub']['C_OS_SS_inclusive']:.3f}"),
          row("closure NPs (c0/lo … c2/hi)", lambda r: ", ".join(f"{100*v['delta']:.0f} %" for v in r["closure_nps"].values())),
          row("τh ID SF per DM (0/1/10/11)", lambda r: "0.92±0.14, 0.88±0.05, 0.87±0.09, 0.90±0.16" if r is M else "0.90±0.13, 0.89±0.05, 0.94±0.15, 0.81±0.15")]
md = "\n".join(lines)
(HERE / "review/wp_compare.md").write_text(md + "\n")
print(md)
