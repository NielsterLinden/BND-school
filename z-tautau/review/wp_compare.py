#!/usr/bin/env python
"""Medium (v2.1) vs Tight (v3 nominal) working-point comparison (docs/08).

    python review/wp_compare.py            # prints the markdown table and writes review/wp_compare.md

Reads output/results.json (v3, Tight) and the recorded v2.1 Medium numbers in slides/history.json.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
T = json.loads((HERE / "output/results.json").read_text())            # v3 nominal: Tight
H = json.loads((HERE / "slides/history.json").read_text())["v2.1"]      # the last Medium result (same chain)


def fit(r):
    return r["fit"]["mcsub"]


rows = [("μ_Z", f"{H['mu']:.3f} +{H['up']:.3f} −{H['down']:.3f}", f"{fit(T)['mu']:.3f} +{fit(T)['mu_err_up']:.3f} −{fit(T)['mu_err_down']:.3f}"),
        ("total uncertainty on μ_Z", f"+{100*H['up']:.1f} / −{100*H['down']:.1f} %", f"+{100*fit(T)['mu_err_up']:.1f} / −{100*fit(T)['mu_err_down']:.1f} %"),
        ("data statistics", f"{100*H['stat']:.1f} %", f"{100*fit(T)['mu_stat']:.1f} %"),
        ("σ(60–120) [pb]", f"{H['sigma60']} +{H['sigma_up']} −{H['sigma_down']}", f"{fit(T)['sigma_60_120_pb']['value']:.0f} +{fit(T)['sigma_60_120_pb']['err_up']:.0f} −{fit(T)['sigma_60_120_pb']['err_down']:.0f}"),
        ("goodness of fit p", f"{H['gof']:.2f}", f"{fit(T)['gof_probability']:.2f}")]
for g in ("Tau ID", "Fakes", "Gammas", "Tau trigger", "Tau energy scale", "Signal modelling", "Background normalisation"):
    rows.append((f"impact: {g}", f"{H['impacts'].get(g, 0):.1f} %", f"{100*fit(T)['grouped_impact'].get(g, 0):.1f} %"))
y = T["yields_prefit"]["mcsub"]; yr = T["yields_prefit_per_region"]["tautau_SR2"]
rows += [("prefit fiducial signal / fakes (all categories)", f"{H['yields']['DYtautau']} / {H['yields']['Fakes']}", f"{y['DYtautau']['value']:.0f} / {y['Fakes']['value']:.0f}"),
         ("SR2: signal / fakes", f"{H['yields']['SR2_DYtautau']} / {H['yields']['SR2_Fakes']}", f"{yr['DYtautau']['value']:.0f} / {yr['Fakes']['value']:.0f}"),
         ("BDT held-out AUC", f"{H['bdt_auc']:.3f}", f"{sum(T['bdt']['training']['auc_test'])/len(T['bdt']['training']['auc_test']):.3f}"),
         ("τh ID SF per DM (0/1/10/11)", H["sf"], "0.90±0.13, 0.89±0.05, 0.94±0.15, 0.81±0.15")]
lines = ["| quantity | Medium (v2.1) | Tight (v3, nominal) |", "|---|---:|---:|"] + [f"| {a} | {b} | {c} |" for a, b, c in rows]
md = "\n".join(lines)
(HERE / "review/wp_compare.md").write_text(md + "\n")
print(md)
