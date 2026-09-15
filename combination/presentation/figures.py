#!/usr/bin/env python
"""Figures made specifically for the deck (everything else is reused from the channels).

    source ../setup.sh && python figures.py

Writes figures/*.pdf and figures/*.png. Light background, to match beamer/metropolis and the
channels' own CMS-style plots. The combination figures come from ../output/plots/ and
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
REPO = HERE.parents[1]          # the repository root (this deck lives in combination/)
COMB = HERE.parent
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
    """Why the combination can now use the shape fit: z-mumu's own stability table.

    Read from `z-mumu/fit/results/stability.json`, the seven fit configurations the channel ran
    after repairing the `SigModel` template (`z-mumu/REVIEW.md` section 0). Configurations with an
    acceptable goodness of fit are shown filled, the ones the channel rejects are greyed.
    """
    rows = json.loads((REPO / "z-mumu/fit/results/stability.json").read_text())
    labels = {"zmumu": "12 x 5 GeV, two-sided SigModel\n(used here)",
              "stab_2gev": "30 x 2 GeV", "stab_10gev": "6 x 10 GeV",
              "stab_1bin": "1 bin = counting", "stab_nosig": "12 x 5 GeV, no SigModel",
              "stab_smooth": "12 x 5 GeV, smoothed",
              "stab_2gev_nosig": "30 x 2 GeV, no SigModel"}
    baseline = next(r for r in rows if r["tag"] == "zmumu")
    # the configurations the channel accepts: a real goodness of fit above 5 %
    good = [r for r in rows if r["gof_p"] > 0.05]
    lo, hi = min(r["mu"] for r in good), max(r["mu"] for r in good)

    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    ax.axvspan(lo, hi, color=C_PRED, alpha=0.18, lw=0,
               label=f"binnings with GoF p > 0.05: {lo:.4f}-{hi:.4f} "
                     rf"($\pm${100 * (hi - lo) / 2 / baseline['mu']:.1f} %)")
    ax.axvline(baseline["mu"], color=C_COMB, lw=1.5,
               label=rf"used here: $\mu_Z$ = {baseline['mu']:.4f} $\pm$ {baseline['err_up']:.4f}")
    for i, r in enumerate(rows[::-1]):
        chosen = r["tag"] == "zmumu"
        accepted = r["gof_p"] > 0.05
        colour = C_COMB if chosen else (C_MUMU if accepted else C_GREY)
        ax.errorbar(r["mu"], i, xerr=0.5 * (r["err_up"] + r["err_down"]),
                    fmt="o" if chosen else "s", color=colour, ms=8 if chosen else 6,
                    capsize=4, lw=2.0)
        gof = "n/a (1 bin)" if r["tag"] == "stab_1bin" else f"{r['gof_p']:.3g}"
        ax.text(1.0255, i, f"GoF p = {gof}", va="center", fontsize=9,
                color=C_GREY if not accepted else "#444")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([labels[r["tag"]] for r in rows[::-1]], fontsize=9.5)
    ax.set_ylim(-0.6, len(rows) - 0.35)
    ax.set_xlim(0.968, 1.042)
    ax.set_xticks([0.98, 0.99, 1.00, 1.01])
    ax.set_xlabel(r"$\mu_Z$ (Z $\to\mu\mu$)")
    ax.legend(fontsize=9.5, loc="upper center", bbox_to_anchor=(0.42, -0.17), ncol=1)
    ax.grid(axis="y", visible=False)
    ax.set_title("after the SigModel repair the fit is stable; grey = rejected by the channel "
                 "(GoF p < 0.05)", loc="left", fontsize=9.5, color=C_GREY)
    save(fig, "mumu_fit_stability")


def precision_budget():
    """One picture of why the combination gains nothing: the correlated floor."""
    payload = json.loads((COMB / "output/combination_result.json").read_text())
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
    if not (COMB / "output/combination_result.json").exists():
        sys.exit("run ../run_combination.py first")
    print("writing deck figures:")
    mumu_fit_stability()
    precision_budget()


if __name__ == "__main__":
    main()
