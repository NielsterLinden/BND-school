#!/usr/bin/env python
"""Figures made specifically for the deck (everything else is reused from the channels).

    source ../setup.sh && python figures.py

Writes figures/*.pdf and figures/*.png. Light background, to match beamer/metropolis and the
channels' own CMS-style plots. The combination figures come from ../combination/output/plots/ and
are not remade here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT = HERE / "figures"

C_MUMU, C_TAUTAU, C_COMB, C_PRED, C_GREY = "#2B6CB0", "#EB811B", "#23373B", "#14B03D", "#8C8C94"

plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "figure.facecolor": "white",
    "savefig.facecolor": "white", "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": ":",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  figures/{name}.pdf")


def mumu_fit_stability():
    """Why the combination uses the counting extraction: z-mumu/REVIEW.md section 4, table 1.

    The five fit configurations the review ran, their mu_Z, and the +-0.7 % lineshape term
    (half the spread) that the recommendation attaches to the counting value.
    """
    rows = [                                        # label, mu_Z, err, GoF p, is the chosen one
        ("30 x 2 GeV bins\n(channel headline)", 0.9903, 0.0140, "0.33", False),
        ("6 x 10 GeV bins", 0.9959, 0.0160, "0.0009", False),
        ("1 bin = counting", 0.993539, None, "n/a (1 bin)", True),
        ("30 bins, no SigModel", 1.00501, 0.0139, "0.0014", False),
        ("30 bins, no SigModel,\nno smoothing", 1.00634, 0.0139, "0.0003", False),
    ]
    chosen = next(r[1] for r in rows if r[4])
    fig, ax = plt.subplots(figsize=(8.4, 4.0))
    ax.axvspan(chosen * (1 - 0.007), chosen * (1 + 0.007), color=C_PRED, alpha=0.18, lw=0,
               label=r"assigned lineshape systematic $\pm$0.7 %")
    ax.axvline(chosen, color=C_COMB, lw=1.5, label=r"used here: counting, $\mu_Z$ = 0.9935")
    for i, (label, mu, err, gof, is_chosen) in enumerate(rows[::-1]):
        colour = C_COMB if is_chosen else C_MUMU
        ax.errorbar(mu, i, xerr=err, fmt="o" if is_chosen else "s", color=colour,
                    ms=8 if is_chosen else 6, capsize=4, lw=2.0)
        ax.text(1.0245, i, f"GoF p = {gof}", va="center", fontsize=9, color=C_GREY)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=9.5)
    ax.set_ylim(-0.6, len(rows) - 0.35)
    ax.set_xlim(0.972, 1.036)
    ax.set_xticks([0.98, 0.99, 1.00, 1.01])
    ax.set_xlabel(r"$\mu_Z$ (Z $\to\mu\mu$)")
    ax.legend(fontsize=9.5, loc="upper center", bbox_to_anchor=(0.45, -0.19), ncol=2)
    ax.grid(axis="y", visible=False)
    ax.set_title("the 30-bin fit is not robust: spread 0.990-1.006, only one has an acceptable GoF",
                 loc="left", fontsize=9.5, color=C_GREY)
    save(fig, "mumu_fit_stability")


def precision_budget():
    """One picture of why the combination gains nothing: the correlated floor."""
    payload = json.loads((REPO / "combination/output/combination_result.json").read_text())
    res, mm = payload["combined"], payload["channels"]["mumu"]
    srcs = {s["name"]: s for s in payload["sources"]}
    corr = sum(v ** 2 for k, v in res["breakdown_pb"].items() if srcs[k]["rho"] >= 0.999) ** 0.5
    uncorr = sum(v ** 2 for k, v in res["breakdown_pb"].items() if srcs[k]["rho"] < 0.999) ** 0.5

    fig, ax = plt.subplots(figsize=(7.8, 3.2))
    bars = [(r"$Z\to\mu\mu$ alone", mm["alone_err_pb"], C_MUMU),
            (r"$Z\to\tau_h\tau_h$ alone", payload["channels"]["tautau"]["alone_err_pb"], C_TAUTAU),
            ("combined", res["total_pb"], C_COMB)]
    y = np.arange(len(bars))
    ax.barh(y, [b[1] for b in bars], color=[b[2] for b in bars], height=0.55)
    for i, (label, value, _) in enumerate(bars):
        ax.text(value * 1.06, i, f"{value:.0f} pb", va="center", fontsize=10.5, color=C_GREY)
    ax.axvline(corr, color=C_PRED, lw=2.0, ls="--",
               label=f"correlated floor {corr:.0f} pb "
                     f"(luminosity {res['lumi_pb']:.0f} pb + shared theory)")
    ax.set_yticks(y)
    ax.set_yticklabels([b[0] for b in bars], fontsize=11)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(8, 900)
    ax.set_xticks([10, 20, 50, 100, 200, 500])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.legend(fontsize=9.5, loc="center right", bbox_to_anchor=(1.0, 0.62))
    ax.set_xlabel(r"total uncertainty on $\sigma(60<m<120)$  [pb]")
    ax.grid(axis="y", visible=False)
    ax.set_title(f"combining cannot go below the floor: the uncorrelated part is only {uncorr:.0f} pb",
                 loc="left", fontsize=10, color=C_GREY)
    save(fig, "precision_budget")


def main():
    if not (REPO / "combination/output/combination_result.json").exists():
        sys.exit("run ../combination/run_combination.py first")
    print("writing deck figures:")
    mumu_fit_stability()
    precision_budget()


if __name__ == "__main__":
    main()
