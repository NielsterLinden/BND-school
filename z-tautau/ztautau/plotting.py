"""CMS-style plots: stacked prediction vs data with a ratio panel, and simple line/efficiency plots."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402

from . import config  # noqa: E402

plt.style.use(hep.style.CMS)
LUMI_LABEL = f"{config.LUMI_PB / 1000:.1f}"
DARK_BG = "#222222"          # the slide-deck background (prompts/presentation_style.md): dark figures must match it
DARK_FG = "#E6E6E6"


def dark_rc():
    """rcParams for figures that go on the dark slide background."""
    return {"figure.facecolor": DARK_BG, "axes.facecolor": DARK_BG, "savefig.facecolor": DARK_BG,
            "axes.edgecolor": DARK_FG, "axes.labelcolor": DARK_FG, "xtick.color": DARK_FG, "ytick.color": DARK_FG,
            "text.color": DARK_FG, "legend.facecolor": DARK_BG, "legend.edgecolor": DARK_BG, "grid.color": "#555555"}


def save(fig, path, dark=False):
    """PNG + PDF (vector), on the dark background if `dark`."""
    path = str(path)
    stem = path[:-4] if path.endswith((".png", ".pdf")) else path
    fc = DARK_BG if dark else "white"
    fig.savefig(stem + ".png", bbox_inches="tight", dpi=120, facecolor=fc)
    fig.savefig(stem + ".pdf", bbox_inches="tight", facecolor=fc)
COLORS = {"DYtautau": "#ffcc66", "DYtautau_nonfid": "#c9963c", "DYlowmass": "#e6a23c", "Fakes": "#ff99cc", "DYee": "#4a90d9", "DYmumu": "#1f5fa8",
          "WJets": "#d95f5f", "TTbar": "#9999cc", "SingleTop": "#6b6bb3", "WW": "#66cc99", "WZ": "#3fa877",
          "ZZ": "#2d7a56", "Diboson": "#66cc99", "Top": "#9999cc", "DYll": "#4a90d9", "MC": "#ffcc66"}
LABELS = {"DYtautau": r"Z$\rightarrow\tau\tau$ (fid.)", "DYtautau_nonfid": r"Z/$\gamma^*\rightarrow\tau\tau$ (non-fid.)", "DYlowmass": r"Z/$\gamma^*\rightarrow\ell\ell$ ($m$<50)", "Fakes": r"jet$\rightarrow\tau_h$ (FF, data)",
          "DYee": r"Z$\rightarrow$ee", "DYmumu": r"Z$\rightarrow\mu\mu$", "WJets": "W+jets", "TTbar": r"t$\bar{t}$",
          "SingleTop": "single t", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ", "Diboson": "diboson", "Top": "top",
          "DYll": r"Z$\rightarrow$ee/$\mu\mu$", "MC": "simulation (genuine $\tau_1$)"}


def label(ax, extra: str | None = None, color="black"):
    """'CMS Open Data' + luminosity header, drawn by hand: mplhep 1.1's label helpers fail to render with
    the matplotlib of LCG_110."""
    ax.text(0.0, 1.01, r"$\bf{CMS}$ $\it{Open\ Data}$", transform=ax.transAxes, fontsize=20, ha="left", va="bottom", color=color)
    ax.text(1.0, 1.01, rf"{LUMI_LABEL} fb$^{{-1}}$ (13 TeV)", transform=ax.transAxes, fontsize=17, ha="right",
            va="bottom", color=color)
    if extra:
        ax.text(0.04, 0.95, extra, transform=ax.transAxes, fontsize=15, va="top", color=color)
    tag(ax, y=0.905 if extra else 0.95, color=color)


def tag(ax, y=0.95, x=0.04, color="black", fontsize=11.5, ha="left", va="top"):
    """Version + working-point stamp (config.PLOT_TAG) so every figure says which analysis it comes from."""
    muted = "#99999E" if color != "black" else "#555555"
    ax.text(x, y, config.PLOT_TAG, transform=ax.transAxes, fontsize=fontsize, va=va, ha=ha, color=muted)


def fig_tag(fig, dark=False, fontsize=10):
    """The same stamp in the top-right corner of a figure with several axes (or no stack_plot header)."""
    fig.text(0.995, 0.995, config.PLOT_TAG, ha="right", va="top", fontsize=fontsize, color="#99999E" if dark else "#555555")


def stack_plot(path, edges, data, stack, xlabel, title=None, logy=False, band=None, ratio_range=(0.5, 1.5),
               density=False, ylabel="Events", dark=False, overlay=None, figsize=(9, 9)):
    """data: (values, errors) or None; stack: list of (name, values, variances) drawn bottom-up;
    band: optional per-bin absolute uncertainty of the total prediction (else MC/fake stat);
    overlay: optional (values, errors, label) of a total prediction drawn as a line (e.g. post-fit) and used
    for the ratio; dark: slide-deck background."""
    with plt.rc_context(dark_rc() if dark else {}):
        return _stack_plot(path, edges, data, stack, xlabel, title, logy, band, ratio_range, density, ylabel, dark, overlay, figsize)


def _stack_plot(path, edges, data, stack, xlabel, title, logy, band, ratio_range, density, ylabel, dark, overlay, figsize):
    fg = DARK_FG if dark else "black"
    edges = np.asarray(edges)
    widths = np.diff(edges)
    scale = 1.0 / widths if density else np.ones(len(widths))
    fig, (ax, rax) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
                                  sharex=True)
    total = np.zeros(len(widths))
    total_var = np.zeros(len(widths))
    names, values, colors = [], [], []
    for name, v, var in stack:
        names.append(LABELS.get(name, name))
        values.append(np.asarray(v) * scale)
        colors.append(COLORS.get(name, None))
        total += v
        total_var += var
    if values:
        hep.histplot(values, edges, stack=True, histtype="fill", label=names, color=colors, ax=ax,
                     edgecolor=DARK_BG if dark else "black", linewidth=0.5)
    ref = total
    unc = np.sqrt(total_var) if band is None else np.asarray(band)
    if overlay is not None:
        ov, oe, olab = overlay
        ref, unc = np.asarray(ov), np.asarray(oe)
        hep.histplot(ref * scale, edges, histtype="step", color="#2BA4DD" if dark else "#1f5fa8", linewidth=2.2, label=olab, ax=ax)
    ax.fill_between(edges, np.append((ref - unc) * scale, 0), np.append((ref + unc) * scale, 0), step="post",
                    hatch="////", facecolor="none", edgecolor="#AAAAAA" if dark else "gray", linewidth=0, label="Unc.")
    if data is not None:
        dv, de = data
        centers = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(centers, dv * scale, yerr=de * scale, fmt="o", color=fg, label="Data", markersize=5)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(ref > 0, dv / ref, np.nan)
            re = np.where(ref > 0, de / ref, np.nan)
        rax.errorbar(centers, r, yerr=re, fmt="o", color=fg, markersize=5)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(ref > 0, unc / ref, 0)
    rax.fill_between(edges, np.append(1 - rel, 1), np.append(1 + rel, 1), step="post", hatch="////",
                     facecolor="none", edgecolor="#AAAAAA" if dark else "gray", linewidth=0)
    rax.axhline(1, color=fg, linewidth=0.8)
    rax.set_ylim(*ratio_range)
    rax.set_ylabel("Data / Pred.")
    rax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel + (" / GeV" if density else ""))
    ax.set_xlim(edges[0], edges[-1])
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(max(0.1, 0.5 * np.min(total[total > 0] * scale[total > 0])) if (total > 0).any() else 0.1,
                    50 * max(np.max(total * scale), 1))
    else:
        top = max(np.max((total + unc) * scale), np.max(data[0] * scale) if data is not None else 0)
        ax.set_ylim(0, 1.9 * top)
    handles, labels = ax.get_legend_handles_labels()
    # the region label sits on the first line inside the frame, the legend starts below it
    ax.legend(handles[::-1], labels[::-1], ncol=2, fontsize=12, loc="upper right", bbox_to_anchor=(1.0, 0.90))
    label(ax, title, color=fg)
    save(fig, path, dark)
    plt.close(fig)
