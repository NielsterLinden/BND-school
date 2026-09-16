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
import tarfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]          # the repository root (this deck lives in combination/)
COMB = HERE.parent
OUT = HERE / "figures"

C_MUMU, C_TAUTAU, C_EE = "#2B6CB0", "#EB811B", "#8E44AD"
C_COMB, C_PRED, C_GREY, C_WARN = "#23373B", "#14B03D", "#8C8C94", "#C0392B"

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
    """What adding a third channel of comparable precision does -- and where the floor is."""
    payload = json.loads((COMB / "output/combination_result.json").read_text())
    res = payload["combined"]
    srcs = {s["name"]: s for s in payload["sources"]}
    corr = sum(v ** 2 for k, v in res["breakdown_pb"].items() if srcs[k]["rho"] >= 0.999) ** 0.5
    uncorr = sum(v ** 2 for k, v in res["breakdown_pb"].items() if srcs[k]["rho"] < 0.999) ** 0.5
    ch = payload["channels"]
    no_ee = payload["variations"]["no_ee"]

    fig, ax = plt.subplots(figsize=(7.8, 3.9))
    bars = [(r"$Z\to\mu\mu$ alone", ch["mumu"]["alone_err_pb"], C_MUMU),
            (r"$Z\to ee$ alone", ch["ee"]["alone_err_pb"], C_EE),
            (r"$Z\to\tau_h\tau_h$ alone", ch["tautau"]["alone_err_pb"], C_TAUTAU),
            (r"$\mu\mu \oplus \tau\tau$", no_ee["error"], C_GREY),
            ("all three", res["total_pb"], C_COMB)]
    y = np.arange(len(bars))
    ax.barh(y, [b[1] for b in bars], color=[b[2] for b in bars], height=0.58)
    for i, (label, value, _) in enumerate(bars):
        ax.text(value * 1.06, i, f"{value:.0f} pb", va="center", fontsize=10.5, color=C_GREY)
    ax.axvline(corr, color=C_PRED, lw=2.0, ls="--")
    ax.text(corr * 0.94, -0.72, f"correlated floor {corr:.0f} pb "
            f"(luminosity {res['lumi_pb']:.0f} pb $\\oplus$ shared theory)",
            ha="right", va="center", fontsize=9.5, color=C_PRED)
    ax.set_yticks(y)
    ax.set_yticklabels([b[0] for b in bars], fontsize=10.5)
    ax.set_ylim(len(bars) - 0.4, -1.1)
    ax.set_xscale("log")
    ax.set_xlim(8, 900)
    ax.set_xticks([10, 20, 50, 100, 200, 500])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"total uncertainty on $\sigma(60<m<120)$  [pb]")
    ax.grid(axis="y", visible=False)
    gain = 100 * (1 - res["total_pb"] / no_ee["error"])
    ax.set_title(f"a second precise channel buys {gain:.0f} %; the uncorrelated part is "
                 f"only {uncorr:.0f} pb, so the floor stands",
                 loc="left", fontsize=10, color=C_GREY)
    save(fig, "precision_budget")


def zee_normalisation():
    """The z-ee signal-template normalisation problem, in one picture.

    Left: what each nuisance parameter does to the *signal normalisation* at the best fit,
    kappa = response(pull). Red are the two that scale the prediction itself (`fitting/
    CONVENTIONS.md` section 3 requires them to be renormalised to a constant yield), so their
    pull must not be absorbed into mu_hat; grey are the ones that scale L x A x C, where it must.
    Right: what that costs, against the mu mu measurement.
    """
    sys.path.insert(0, str(COMB / "tools"))
    from extract_zee import kappa                                   # noqa: E402

    zee = json.loads((COMB / "inputs/zee_fit_result.json").read_text())
    payload = json.loads((COMB / "output/combination_result.json").read_text())
    nps = {k: v for k, v in zee["nuisance_parameters"].items() if "norm_up" in v}
    theory = ("PDF", "QCDScale")
    items = sorted(((k, 100 * (kappa(v["pull"], v["norm_up"], v["norm_down"]) - 1), v["pull"])
                    for k, v in nps.items()), key=lambda t: abs(t[1]))

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 3.9),
                             gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.42})

    ax = axes[0]
    names = [k for k, _, _ in items]
    vals = [v for _, v, _ in items]
    cols = [C_WARN if k in theory else C_GREY for k in names]
    ax.barh(names, vals, color=cols, height=0.62)
    for i, (k, v, pull) in enumerate(items):
        ax.text(v + (0.5 if v > 0 else -0.5), i, f"{v:+.1f} %  ({pull:+.2f}$\\sigma$)",
                va="center", ha="left" if v > 0 else "right", fontsize=8.5,
                color=C_WARN if k in theory else "#555")
    ax.axvline(0, color="#999", lw=0.9)
    ax.set_xlim(-20, 12)
    ax.set_xlabel(r"effect on the $Z\to ee$ signal normalisation at the best fit  [%]")
    ax.set_title("red: scales the prediction (must not be absorbed by $\\mu$)\n"
                 r"grey: scales $L\times A\times C$ (must be)", loc="left", fontsize=9.5)
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    ee, mm = payload["channels"]["ee"], payload["channels"]["mumu"]
    rows = [("as published\n" + r"$\hat\mu\times\sigma^{\rm pred}$",
             ee["sigma_pb"], ee["sigma_err_pb"], C_EE),
            (r"counting" + "\n" + r"$(N_{\rm data}-N_{\rm bkg})/N_{\rm sig}$",
             ee["extra"]["sigma_counting_pb"], 0.0, C_GREY),
            ("normalisation fixed\n" + r"$\hat\mu\times\kappa\times\sigma^{\rm pred}$",
             payload["variations"]["ee_normfix"]["channels"]["ee"]["sigma_pb"],
             payload["variations"]["ee_normfix"]["channels"]["ee"]["sigma_err_pb"], C_WARN)]
    ax.axvspan(mm["sigma_pb"] - mm["sigma_err_pb"], mm["sigma_pb"] + mm["sigma_err_pb"],
               color=C_MUMU, alpha=0.18, lw=0,
               label=rf"$Z\to\mu\mu$: {mm['sigma_pb']:.0f} $\pm$ {mm['sigma_err_pb']:.0f} pb")
    ax.axvline(mm["sigma_pb"], color=C_MUMU, lw=1.4)
    for i, (lab, value, err, col) in enumerate(rows[::-1]):
        ax.errorbar(value, i, xerr=err if err else None, fmt="o", color=col, ms=8, capsize=4, lw=2.2)
        ax.text(value, i + 0.22, f"{value:.0f} pb", ha="center", fontsize=10, color=col)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=9)
    ax.set_ylim(-0.55, len(rows) - 0.3)
    ax.set_xlim(1770, 2140)
    ax.set_xlabel(r"$\sigma(Z\to ee,\ 60<m<120)$  [pb]")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="y", visible=False)
    ax.set_title(r"$\kappa$ = " f"{zee['theory_normalisation']['kappa']:.3f}"
                 f"  $\\Rightarrow$  {payload['variations']['ee_normfix']['shift']:+.0f} pb "
                 "on the combination", loc="left", fontsize=10)
    save(fig, "zee_normalisation")


#: plots the z-ee group made inside its TRExFitter job, pulled out of the committed tarball so the
#: deck can use them without a second copy living in git. name in the tarball -> name in figures/.
ZEE_PLOTS = {"Plots/SR_postFit.png": "zee_SR_postfit.png",
             "Plots/SR.png": "zee_SR_prefit.png",
             "Pulls/All/NuisPar.png": "zee_pulls.png",
             "Rankings/Ranking_mu_signal_Breakdown.png": "zee_ranking.png"}


def zee_channel_plots():
    """Extract the z-ee TRExFitter plots from `z-ee/Zee_fit.tar.gz` into figures/."""
    tar_path = REPO / "z-ee/Zee_fit.tar.gz"
    if not tar_path.exists():
        print(f"  !! {tar_path} missing, skipping the z-ee channel plots")
        return
    OUT.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        for inside, name in ZEE_PLOTS.items():
            member = tar.extractfile(f"Zee_fit/{inside}")
            if member is None:
                print(f"  !! Zee_fit/{inside} not in the tarball")
                continue
            (OUT / name).write_bytes(member.read())
            print(f"  figures/{name}")


def main():
    if not (COMB / "output/combination_result.json").exists():
        sys.exit("run ../run_combination.py first")
    print("writing deck figures:")
    mumu_fit_stability()
    precision_budget()
    zee_normalisation()
    zee_channel_plots()


if __name__ == "__main__":
    main()
