"""Figures for the combination: vector PDF (for the deck) plus PNG (for markdown).

Light background on purpose -- the deck is beamer/metropolis and the channels' own plots are
white, so a dark figure would show as a black rectangle on the slide.

    from comb import plots
    plots.make_all(spec, res, solo, lik, rat, variations, Path("output/plots"))
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# metropolis-compatible palette
C_MUMU, C_TAUTAU, C_COMB = "#2B6CB0", "#EB811B", "#23373B"
C_PRED, C_GREY, C_GREEN = "#14B03D", "#8C8C94", "#577D29"

#: published references for the comparison figure (see docs/04-results.md for the citations)
REFERENCES = [
    ("CMS 2024 (206 pb$^{-1}$)\n60 < m < 120 GeV", 1952.0, (4.0 ** 2 + 18.0 ** 2 + 45.0 ** 2) ** 0.5,
     "CMS-SMP-20-004, arXiv:2408.03744"),
    ("ATLAS 2016 (81 pb$^{-1}$)\n66 < m < 116 GeV", 1981.0, (7.0 ** 2 + 38.0 ** 2 + 42.0 ** 2) ** 0.5,
     "arXiv:1603.09222"),
]

plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": ":",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def _save(fig, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out.with_suffix(".pdf")


def _label(name: str) -> str:
    return {"mumu": r"$Z\to\mu\mu$", "tautau": r"$Z\to\tau_h\tau_h$"}.get(name, name)


# --------------------------------------------------------------------------------- figures
def forest(spec, res, solo, out: Path) -> Path:
    """The two channels and the combination on one axis, with the aMC@NLO reference."""
    fig, ax = plt.subplots(figsize=(7.6, 3.6))
    rows = [(c.name, c.sigma, c.sigma_err_total, c.sigma_stat) for c in spec.channels][::-1]
    rows.append(("combined", res.value, res.error, res.group_breakdown()["statistical"]))
    colours = {"mumu": C_MUMU, "tautau": C_TAUTAU, "combined": C_COMB}

    # the two channels normalise to slightly different aMC@NLO references (docs/01-inputs.md),
    # so the reference is drawn as the band they span rather than as one line
    lo = min(c.sigma_pred for c in spec.channels)
    hi = max(c.sigma_pred for c in spec.channels)
    ax.axvspan(lo, hi, color=C_PRED, alpha=0.20, lw=0,
               label=f"aMC@NLO reference {lo:.0f}-{hi:.0f} pb")

    for i, (name, value, err, stat) in enumerate(rows):
        c = colours[name]
        ax.errorbar(value, i, xerr=err, fmt="o", color=c, ms=7 if name == "combined" else 6,
                    capsize=4, lw=2.0, zorder=3)
        ax.errorbar(value, i, xerr=stat, fmt="none", color=c, capsize=8, lw=3.5, zorder=4)
        ax.text(value, i + 0.28, f"{value:.0f} $\\pm$ {err:.0f} pb", ha="center", va="bottom",
                fontsize=10, color=c, fontweight="bold" if name == "combined" else "normal")

    ax.axhspan(-0.5, 0.5, color=C_COMB, alpha=0.05, lw=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([_label(n) if n != "combined" else "combined" for n, *_ in rows], fontsize=12)
    ax.set_ylim(-0.6, len(rows) - 0.25)
    ax.set_xlabel(r"$\sigma(pp\to Z/\gamma^*\to\ell\ell,\ 60<m<120\,$GeV$)$  [pb]")
    ax.set_xlim(1550, 2720)
    ax.legend(loc="lower right", fontsize=9)
    ax.text(0.015, 0.06, f"$\\chi^2$/ndf = {res.chi2:.2f}/{res.ndf},  p = {res.pvalue:.2f}",
            transform=ax.transAxes, fontsize=9, color=C_GREY, va="bottom")
    ax.set_title("CMS Open Data 2016 (Run2016G+H), 16.4 fb$^{-1}$", loc="left", fontsize=10, color=C_GREY)
    return _save(fig, out)


def breakdown(spec, res, out: Path) -> Path:
    """Where the uncertainty of the combined number comes from, next to the channel inputs."""
    items = [(k, v) for k, v in res.breakdown.items() if v > 0.02][:14]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), gridspec_kw={"width_ratios": [1.35, 1]})

    ax = axes[0]
    names = [k for k, _ in items][::-1]
    values = [v for _, v in items][::-1]
    rho = {src.name: src.rho for src in spec.sources}
    colours = [C_GREEN if rho.get(n, 0.0) >= 0.999 else C_COMB for n in names]
    ax.barh(names, values, color=colours, height=0.65)
    for i, v in enumerate(values):
        ax.text(v + max(values) * 0.015, i, f"{v:.1f} ({100 * v / res.value:.2f} %)",
                va="center", fontsize=8.5, color="#444")
    ax.set_xlabel(r"contribution to $\delta\sigma_{\rm comb}$  [pb]")
    ax.set_xlim(0, max(values) * 1.38)
    ax.set_title(f"combined: {res.value:.0f} $\\pm$ {res.error:.0f} pb ({100 * res.rel:.2f} %)"
                 "   -- green: correlated between the channels", loc="left", fontsize=10.5)
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    cats = sorted({k for c in spec.channels for k in c.groups},
                  key=lambda k: -max(c.groups.get(k, 0) for c in spec.channels))
    y = np.arange(len(cats))
    w = 0.38
    for off, c, col in ((+w / 2, spec.channels[0], C_MUMU), (-w / 2, spec.channels[1], C_TAUTAU)):
        ax.barh(y + off, [100 * c.groups.get(k, 0) / c.mu for k in cats], height=w,
                color=col, label=_label(c.name))
    ax.set_yticks(y)
    ax.set_yticklabels(cats, fontsize=9)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.003, 40)
    ax.set_xlabel(r"impact on $\mu_Z$  [%]")
    ax.set_title("in-fit systematics, per channel", loc="left", fontsize=11)
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def sources(spec, res, out: Path) -> Path:
    """Every uncertainty source, its size in each channel, and whether it is correlated."""
    srcs = sorted(spec.sources, key=lambda s: -max(s.sizes.values()))
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    y = np.arange(len(srcs))
    w = 0.38
    ax.barh(y + w / 2, [s.sizes.get("mumu", 0) for s in srcs], height=w, color=C_MUMU, label=_label("mumu"))
    ax.barh(y - w / 2, [s.sizes.get("tautau", 0) for s in srcs], height=w, color=C_TAUTAU, label=_label("tautau"))
    def rho_tag(rho):
        return r"$\rho{=}1$" if rho >= 0.999 else (r"$\rho{=}0$" if rho <= 0.001 else fr"$\rho{{=}}{rho:.2f}$")

    ax.set_yticks(y)
    ax.set_yticklabels([f"{s.name}   {rho_tag(s.rho)}" for s in srcs], fontsize=8.5)
    for tick, s in zip(ax.get_yticklabels(), srcs):
        if s.rho >= 0.999:
            tick.set_color(C_GREEN)
            tick.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.15, 800)
    ax.set_xlabel(r"absolute uncertainty on $\sigma(60<m<120)$  [pb]")
    ax.set_title(f"correlated sources give $\\rho(\\mu\\mu,\\tau\\tau)$ = {res.correlation[0,1]:.2f}",
                 loc="left", fontsize=11)
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def likelihood_scan(res, lik, out: Path) -> Path:
    """-2 dlnL from the nuisance-parameter form against the BLUE parabola."""
    grid, curve = lik.scan
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(grid, curve, color=C_COMB, lw=2.2, label="profile likelihood (asymmetric)")
    ax.plot(grid, ((grid - res.value) / res.error) ** 2, color=C_MUMU, lw=1.6, ls="--",
            label="BLUE parabola")
    for level, text in ((1.0, r"$1\sigma$"), (4.0, r"$2\sigma$")):
        ax.axhline(level, color=C_GREY, lw=0.8, ls=":")
        ax.text(grid[0], level + 0.08, text, fontsize=9, color=C_GREY)
    ax.axvline(lik.value, color=C_COMB, lw=0.9, ls="-", alpha=0.5)
    ax.set_xlabel(r"$\sigma(pp\to Z/\gamma^*\to\ell\ell,\ 60<m<120\,$GeV$)$  [pb]")
    ax.set_ylabel(r"$-2\,\Delta\ln L$")
    ax.set_ylim(0, 9)
    ax.set_xlim(grid[0], grid[-1])
    ax.legend(fontsize=9.5, loc="upper center")
    ax.set_title(f"{lik.value:.0f} $^{{+{lik.error_up:.0f}}}_{{-{lik.error_down:.0f}}}$ pb   "
                 f"(BLUE {res.value:.0f} $\\pm$ {res.error:.0f} pb)", loc="left", fontsize=10.5)
    return _save(fig, out)


def universality(spec, rat, out: Path) -> Path:
    """R = sigma(tautau)/sigma(mumu): the value, and what cancels in it."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.6),
                             gridspec_kw={"width_ratios": [1, 1.3], "wspace": 0.55})

    ax = axes[0]
    ax.axvline(1.0, color=C_PRED, lw=1.6, ls="--", label="lepton universality")
    ax.axvspan(1 - 0.0, 1 + 0.0, color=C_PRED, alpha=0.1, lw=0)
    ax.errorbar(rat.value, 0, xerr=rat.error, fmt="o", color=C_TAUTAU, ms=8, capsize=5, lw=2.2)
    ax.text(rat.value, 0.10, f"$R$ = {rat.value:.2f} $\\pm$ {rat.error:.2f}", ha="center",
            fontsize=11.5, color=C_TAUTAU)
    ax.text(rat.value, -0.16, f"{rat.z_from_unity:+.1f}$\\sigma$ from 1  (p = {rat.pvalue:.2f})",
            ha="center", fontsize=9.5, color=C_GREY)
    ax.set_yticks([])
    ax.set_ylim(-0.45, 0.45)
    ax.set_xlim(0.55, 1.75)
    ax.set_xlabel(r"$R=\sigma(Z\to\tau\tau)\,/\,\sigma(Z\to\mu\mu)$")
    ax.legend(fontsize=9.5, loc="upper left")
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    items = [(k, v) for k, v in rat.breakdown.items() if v > 0.002][:10][::-1]
    cancels = {s.name: s.rho for s in spec.sources}
    ax.barh([k for k, _ in items], [v for _, v in items],
            color=[C_GREEN if cancels.get(k, 0) >= 0.999 else C_TAUTAU for k, _ in items], height=0.62)
    ax.set_xlabel(r"contribution to $\delta R$")
    ax.set_title(r"green: $\rho=1$ sources, which cancel in the ratio", loc="left", fontsize=9.5)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def comparison(spec, res, out: Path) -> Path:
    """This combination next to the published 13 TeV measurements and the reference prediction."""
    fig, ax = plt.subplots(figsize=(7.8, 3.8))
    rows = [(lab, v, e, C_GREY) for lab, v, e, _ in REFERENCES]
    rows.append(("this work (combined)\n$\\mu\\mu + \\tau_h\\tau_h$, 16.4 fb$^{-1}$",
                 res.value, res.error, C_COMB))
    pred = np.mean([c.sigma_pred for c in spec.channels])
    ax.axvline(pred, color=C_PRED, lw=1.4, ls="--", label=f"aMC@NLO reference {pred:.0f} pb")
    for i, (lab, value, err, colour) in enumerate(rows):
        ax.errorbar(value, i, xerr=err, fmt="o", color=colour, ms=7, capsize=4,
                    lw=2.4 if colour == C_COMB else 1.8)
        ax.text(value, i - 0.34, f"{value:.0f} $\\pm$ {err:.0f} pb", ha="center", fontsize=9.5,
                color=colour, fontweight="bold" if colour == C_COMB else "normal")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.5)
    ax.set_ylim(-0.75, len(rows) - 0.35)
    ax.set_xlim(1820, 2120)
    ax.set_xlabel(r"$\sigma\cdot\mathcal{B}(Z\to\ell\ell)$  [pb]   (uncertainties added in quadrature)")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def variations_plot(res, variations, out: Path) -> Path:
    """How far each variation of the combination moves the central value.

    Every variation carries essentially the baseline uncertainty (34.5-34.7 pb), so the plot shows
    the shifts on a zoomed axis; the baseline +-1 sigma band is marked at the edges for scale.
    """
    keys = [k for k in variations if k != "baseline"]
    fig, ax = plt.subplots(figsize=(7.8, 3.4))
    span = max(12.0, 1.35 * max(abs(variations[k]["shift"]) for k in keys))
    ax.axvline(res.value, color=C_COMB, lw=1.6, label=f"baseline {res.value:.0f} pb")
    for i, k in enumerate(keys[::-1]):
        v = variations[k]
        ax.plot(v["value"], i, "s", color=C_MUMU, ms=7)
        ax.annotate(f"{v['shift']:+.1f} pb", (v["value"], i), textcoords="offset points",
                    xytext=(0, 11), ha="center", fontsize=9, color=C_GREY)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([k.replace("_", " ") for k in keys[::-1]], fontsize=10)
    ax.set_ylim(-0.6, len(keys) - 0.25)
    ax.set_xlim(res.value - span, res.value + span)
    ax.set_xlabel(r"combined $\sigma(60<m<120)$  [pb]")
    top = ax.secondary_xaxis("top", functions=(lambda v: (v - res.value) / res.error,
                                               lambda t: res.value + t * res.error))
    top.set_xlabel(r"shift in units of the baseline uncertainty ($\pm$"
                   f"{res.error:.0f} pb)", fontsize=9.5)
    top.tick_params(labelsize=9)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def make_all(spec, res, solo, lik, rat, variations, outdir: Path) -> list[Path]:
    outdir = Path(outdir)
    return [forest(spec, res, solo, outdir / "forest"),
            breakdown(spec, res, outdir / "breakdown"),
            sources(spec, res, outdir / "sources"),
            likelihood_scan(res, lik, outdir / "likelihood_scan"),
            universality(spec, rat, outdir / "universality"),
            comparison(spec, res, outdir / "comparison"),
            variations_plot(res, variations, outdir / "variations")]
