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
C_MUMU, C_TAUTAU, C_EE, C_COMB = "#2B6CB0", "#EB811B", "#8E44AD", "#23373B"
C_PRED, C_GREY, C_GREEN = "#14B03D", "#8C8C94", "#577D29"

#: one colour per channel, plus the combination
COLOUR = {"mumu": C_MUMU, "tautau": C_TAUTAU, "ee": C_EE, "combined": C_COMB}

#: ratio sigma(60 < m < 120) / sigma(66 < m < 116) at LHE level in the same aMC@NLO DY sample
#: (z-mumu/output/v2/gensums.json, DY_NLO h_lhe_mll). Used to put ATLAS's narrower window on ours.
WINDOW_66_116_TO_60_120 = 1.01425

#: published 13 TeV measurements, all converted to sigma(60 < m_ll < 120 GeV) in pb.
#: (label, value, total uncertainty, citation). `scaled` marks the ones we moved windows.
REFERENCES = [
    ("CMS 2024 (206 pb$^{-1}$)\n60 < m < 120 GeV", 1952.0, (4.0 ** 2 + 18.0 ** 2 + 45.0 ** 2) ** 0.5,
     "CMS-SMP-20-004, arXiv:2408.03744"),
    ("ATLAS 2016 (81 pb$^{-1}$)\n66 < m < 116 GeV $\\to$ 60-120", 1981.0 * WINDOW_66_116_TO_60_120,
     (7.0 ** 2 + 38.0 ** 2 + 42.0 ** 2) ** 0.5 * WINDOW_66_116_TO_60_120, "arXiv:1603.09222"),
]

#: the same, per channel. Only these three exist: CMS SMP-20-004 fits ee and mu mu together and
#: publishes no per-channel cross section, so the CMS point is the combined one; ATLAS publishes
#: both lepton channels separately (its Table 2); and the only published inclusive Z -> tau tau
#: cross section at 13 TeV is CMS's, in exactly our 60-120 GeV window.
#: (channel, label, value pb, total uncertainty pb, citation)
CHANNEL_REFERENCES = [
    ("ee", "ATLAS $Z\\to ee$\n66-116 GeV $\\to$ 60-120", 1987.0 * WINDOW_66_116_TO_60_120,
     (11.0 ** 2 + 41.0 ** 2 + 42.0 ** 2) ** 0.5 * WINDOW_66_116_TO_60_120, "arXiv:1603.09222 Tab. 2"),
    ("mumu", "ATLAS $Z\\to\\mu\\mu$\n66-116 GeV $\\to$ 60-120", 1977.0 * WINDOW_66_116_TO_60_120,
     (9.0 ** 2 + 41.0 ** 2 + 42.0 ** 2) ** 0.5 * WINDOW_66_116_TO_60_120, "arXiv:1603.09222 Tab. 2"),
    ("tautau", "CMS $Z\\to\\tau\\tau$ (2.3 fb$^{-1}$)\n60-120 GeV", 1848.0,
     (12.0 ** 2 + 57.0 ** 2 + 35.0 ** 2) ** 0.5, "arXiv:1801.03535"),
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


#: TeX for the decay, without the surrounding $ -- so it can go inside a larger expression
TEX = {"mumu": r"Z\to\mu\mu", "tautau": r"Z\to\tau_h\tau_h", "ee": r"Z\to ee"}
#: the short form used in titles
SHORT_TEX = {"mumu": r"\mu\mu", "tautau": r"\tau\tau", "ee": r"ee"}


def _label(name: str) -> str:
    return f"${TEX[name]}$" if name in TEX else name


def _ratio_label(rat) -> str:
    return fr"$\sigma({TEX[rat.numerator]})\,/\,\sigma({TEX[rat.denominator]})$"


# --------------------------------------------------------------------------------- figures
def forest(spec, res, solo, out: Path) -> Path:
    """Every channel and the combination on one axis, with the aMC@NLO reference."""
    fig, ax = plt.subplots(figsize=(7.6, 1.5 + 0.75 * (len(spec.channels) + 1)))
    rows = [(c.name, c.sigma, c.sigma_err_total, c.sigma_stat) for c in spec.channels][::-1]
    rows.append(("combined", res.value, res.error, res.group_breakdown()["statistical"]))
    colours = COLOUR

    # the channels normalise to slightly different aMC@NLO references (docs/01-inputs.md),
    # so the reference is drawn as the band they span rather than as one line
    lo = min(c.sigma_pred for c in spec.channels)
    hi = max(c.sigma_pred for c in spec.channels)
    ax.axvspan(lo, hi, color=C_PRED, alpha=0.20, lw=0,
               label=f"aMC@NLO reference {lo:.0f}-{hi:.0f} pb")

    for i, (name, value, err, stat) in enumerate(rows):
        c = colours.get(name, C_GREY)
        ax.errorbar(value, i, xerr=err, fmt="o", color=c, ms=7 if name == "combined" else 6,
                    capsize=4, lw=2.0, zorder=3)
        ax.errorbar(value, i, xerr=stat, fmt="none", color=c, capsize=8, lw=3.5, zorder=4)
        ax.text(value, i + 0.28, f"{value:.0f} $\\pm$ {err:.0f} pb", ha="center", va="bottom",
                fontsize=10, color=c, fontweight="bold" if name == "combined" else "normal")

    ax.axhspan(len(rows) - 1.5, len(rows) - 0.5, color=C_COMB, alpha=0.06, lw=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([_label(n) if n != "combined" else "combined" for n, *_ in rows], fontsize=12)
    ax.set_ylim(-0.6, len(rows) - 0.25)
    ax.set_xlabel(r"$\sigma(pp\to Z/\gamma^*\to\ell\ell,\ 60<m<120\,$GeV$)$  [pb]")
    lo_x = min(v - e for _, v, e, _ in rows)
    hi_x = max(v + e for _, v, e, _ in rows)
    pad = 0.22 * (hi_x - lo_x)
    ax.set_xlim(lo_x - pad, hi_x + pad)
    ax.legend(loc="lower right", fontsize=9)
    ax.text(0.015, 0.94, f"$\\chi^2$/ndf = {res.chi2:.2f}/{res.ndf},  p = {res.pvalue:.3f}",
            transform=ax.transAxes, fontsize=9, color=C_GREY, va="top")
    ax.set_title("CMS Open Data 2016 (Run2016G+H), 16.4 fb$^{-1}$", loc="left", fontsize=10, color=C_GREY)
    return _save(fig, out)


def breakdown(spec, res, out: Path) -> Path:
    """Where the uncertainty of the combined number comes from, next to the channel inputs."""
    items = [(k, v) for k, v in res.breakdown.items() if v > 0.02][:16]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.0), gridspec_kw={"width_ratios": [1.35, 1]})

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
    # two lines: on one line this title runs into the right-hand axes
    ax.set_title(f"combined: {res.value:.0f} $\\pm$ {res.error:.0f} pb ({100 * res.rel:.2f} %)"
                 "\ngreen: correlated between the channels", loc="left", fontsize=10.5)
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    cats = sorted({k for c in spec.channels for k in c.groups},
                  key=lambda k: -max(c.groups.get(k, 0) for c in spec.channels))
    y = np.arange(len(cats))
    n = len(spec.channels)
    w = 0.78 / n
    for i, c in enumerate(spec.channels):
        off = (n - 1) / 2.0 * w - i * w
        ax.barh(y + off, [100 * c.groups.get(k, 0) / c.mu for k in cats], height=w,
                color=COLOUR.get(c.name, C_GREY), label=_label(c.name))
    ax.set_yticks(y)
    ax.set_yticklabels(cats, fontsize=9)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.003, 40)
    ax.set_xlabel(r"impact on $\mu_Z$  [%]")
    ax.set_title("in-fit systematics, per channel", loc="left", fontsize=11)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def sources(spec, res, out: Path) -> Path:
    """Every uncertainty source, its size in each channel, and whether it is correlated."""
    srcs = sorted(spec.sources, key=lambda s: -max(s.sizes.values()))
    fig, ax = plt.subplots(figsize=(8.4, 0.9 + 0.26 * len(srcs)))
    y = np.arange(len(srcs))
    n = len(spec.channels)
    w = 0.78 / n
    for i, c in enumerate(spec.channels):
        off = (n - 1) / 2.0 * w - i * w
        ax.barh(y + off, [s.sizes.get(c.name, 0) for s in srcs], height=w,
                color=COLOUR.get(c.name, C_GREY), label=_label(c.name))

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
    order = res.meta["order"]
    pairs = [(order[i], order[j], res.correlation[i, j])
             for i in range(len(order)) for j in range(i + 1, len(order))]
    tag = ",  ".join(fr"$\rho({SHORT_TEX.get(a, a)},{SHORT_TEX.get(b, b)})$ = {r:.2f}"
                     for a, b, r in pairs)
    ax.set_title("correlated sources give " + tag, loc="left", fontsize=10)
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


def universality(spec, ratios, out: Path) -> Path:
    """R = sigma(X)/sigma(mu mu) for every other channel: the values, and what cancels in them.

    `ratios` is a list of `comb.ratio.RatioResult`, one per numerator channel.
    """
    ratios = list(ratios)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 1.6 + 1.0 * len(ratios)),
                             gridspec_kw={"width_ratios": [1, 1.3], "wspace": 0.55})

    ax = axes[0]
    ax.axvline(1.0, color=C_PRED, lw=1.6, ls="--", label="lepton universality")
    for i, rat in enumerate(ratios[::-1]):
        col = COLOUR.get(rat.numerator, C_GREY)
        ax.errorbar(rat.value, i, xerr=rat.error, fmt="o", color=col, ms=8, capsize=5, lw=2.2)
        ax.text(rat.value, i + 0.16, f"{rat.value:.2f} $\\pm$ {rat.error:.2f}", ha="center",
                fontsize=10.5, color=col)
        ax.text(rat.value, i - 0.30, f"{rat.z_from_unity:+.1f}$\\sigma$  (p = {rat.pvalue:.3f})",
                ha="center", fontsize=9, color=C_GREY)
    ax.set_yticks(range(len(ratios)))
    ax.set_yticklabels([_ratio_label(r) for r in ratios[::-1]], fontsize=10)
    lo = min(r.value - 1.6 * r.error for r in ratios)
    hi = max(r.value + 1.6 * r.error for r in ratios)
    ax.set_xlim(min(0.80, lo), max(1.25, hi))
    ax.set_ylim(-0.6, len(ratios) - 0.35)
    ax.set_xlabel(r"$R = \sigma(Z\to X)\,/\,\sigma(Z\to\mu\mu)$")
    ax.legend(fontsize=9.5, loc="upper left")
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    rat = max(ratios, key=lambda r: abs(r.z_from_unity))   # the one that needs explaining
    items = [(k, v) for k, v in rat.breakdown.items() if v > 0.002][:10][::-1]
    cancels = {s.name: s.rho for s in spec.sources}
    ax.barh([k for k, _ in items], [v for _, v in items],
            color=[C_GREEN if cancels.get(k, 0) >= 0.999 else COLOUR.get(rat.numerator, C_GREY)
                   for k, _ in items], height=0.62)
    ax.set_xlabel("contribution to $\\delta R$  for " + _ratio_label(rat))
    ax.set_title(r"green: $\rho=1$ sources, which cancel in the ratio", loc="left", fontsize=9.5)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def comparison(spec, res, out: Path) -> Path:
    """This combination next to the published 13 TeV measurements and the reference prediction."""
    fig, ax = plt.subplots(figsize=(7.8, 3.8))
    rows = [(lab, v, e, C_GREY) for lab, v, e, _ in REFERENCES]
    what = " + ".join(SHORT_TEX.get(c.name, c.name) for c in spec.channels)
    rows.append((f"this work (combined)\n${what}$, 16.4 fb$^{{-1}}$", res.value, res.error, C_COMB))
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
    lo = min(v - e for _, v, e, _ in rows) if False else min(r[1] - r[2] for r in rows)
    hi = max(r[1] + r[2] for r in rows)
    pad = 0.18 * (hi - lo)
    ax.set_xlim(min(lo - pad, pred - 40), max(hi + pad, pred + 40))
    ax.set_xlabel(r"$\sigma\cdot\mathcal{B}(Z\to\ell\ell)$  [pb]   (uncertainties added in quadrature)")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="y", visible=False)
    return _save(fig, out)


def channel_comparison(spec, res, out: Path) -> Path:
    """Each channel next to the published measurement of the same decay.

    Everything is on our window, 60 < m < 120 GeV: the ATLAS points are scaled from 66-116 GeV by
    the LHE-level ratio of the same aMC@NLO sample (`WINDOW_66_116_TO_60_120` = 1.0143), which is
    a 1.4 % correction with a negligible model dependence at this precision.
    """
    by_channel = {}
    for name, label, value, err, cite in CHANNEL_REFERENCES:
        by_channel.setdefault(name, []).append((label, value, err))

    rows = []
    for c in spec.channels:
        rows.append(("this", _label(c.name), c.sigma, c.sigma_err_total, COLOUR.get(c.name, C_GREY)))
        for label, value, err in by_channel.get(c.name, []):
            rows.append(("pub", label, value, err, C_GREY))
    rows.append(("this", "combined", res.value, res.error, C_COMB))
    for label, value, err, _ in [(r[0], r[1], r[2], r[3]) for r in REFERENCES]:
        rows.append(("pub", label.replace("\n", " "), value, err, C_GREY))

    fig, ax = plt.subplots(figsize=(8.6, 0.9 + 0.44 * len(rows)))
    pred = np.mean([c.sigma_pred for c in spec.channels])
    ax.axvline(pred, color=C_PRED, lw=1.3, ls="--", label=f"aMC@NLO reference {pred:.0f} pb")
    for i, (kind, label, value, err, colour) in enumerate(rows[::-1]):
        ax.errorbar(value, i, xerr=err, fmt="o" if kind == "this" else "s", color=colour,
                    ms=7 if kind == "this" else 5.5, capsize=4,
                    lw=2.4 if kind == "this" else 1.6)
        ax.text(value, i + 0.20, f"{value:.0f} $\\pm$ {err:.0f}", ha="center", fontsize=8.5,
                color=colour, fontweight="bold" if kind == "this" else "normal")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([("this work: " if k == "this" else "") + l for k, l, *_ in rows[::-1]],
                       fontsize=8.5)
    for tick, (kind, *_) in zip(ax.get_yticklabels(), rows[::-1]):
        if kind == "this":
            tick.set_fontweight("bold")
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_xlabel(r"$\sigma\cdot\mathcal{B}(Z\to\ell\ell)$, $60<m<120$ GeV  [pb]"
                  "   (uncertainties added in quadrature)")
    lo = min(v - e for _, _, v, e, _ in rows)
    hi = max(v + e for _, _, v, e, _ in rows)
    pad = 0.10 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="y", visible=False)
    ax.set_title("ATLAS points scaled from 66-116 GeV by 1.0143; CMS SMP-20-004 fits $ee$ and "
                 r"$\mu\mu$ together", loc="left", fontsize=9, color=C_GREY)
    return _save(fig, out)


def variations_plot(res, variations, out: Path) -> Path:
    """How far each variation of the combination moves the central value.

    Every variation carries essentially the baseline uncertainty (34.5-34.7 pb), so the plot shows
    the shifts on a zoomed axis; the baseline +-1 sigma band is marked at the edges for scale.
    """
    keys = [k for k in variations if k != "baseline"]
    # one row per variation, plus room for the "+x.x pb" label that sits above the top marker
    fig, ax = plt.subplots(figsize=(7.8, 1.3 + 0.42 * len(keys)))
    span = max(12.0, 1.35 * max(abs(variations[k]["shift"]) for k in keys))
    ax.axvline(res.value, color=C_COMB, lw=1.6, label=f"baseline {res.value:.0f} pb")
    for i, k in enumerate(keys[::-1]):
        v = variations[k]
        ax.plot(v["value"], i, "s", color=C_MUMU, ms=7)
        side = -1 if v["shift"] > 0 else 1
        ax.annotate(f"{v['shift']:+.1f} pb", (v["value"], i), textcoords="offset points",
                    xytext=(-side * 10, 10), ha="right" if side > 0 else "left",
                    fontsize=9, color=C_GREY)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([k.replace("_", " ") for k in keys[::-1]], fontsize=10)
    ax.set_ylim(-0.7, len(keys) - 0.1)
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


def make_all(spec, res, solo, lik, ratios, variations, outdir: Path) -> list[Path]:
    outdir = Path(outdir)
    return [forest(spec, res, solo, outdir / "forest"),
            breakdown(spec, res, outdir / "breakdown"),
            sources(spec, res, outdir / "sources"),
            likelihood_scan(res, lik, outdir / "likelihood_scan"),
            universality(spec, ratios, outdir / "universality"),
            comparison(spec, res, outdir / "comparison"),
            channel_comparison(spec, res, outdir / "comparison_channels"),
            variations_plot(res, variations, outdir / "variations")]
