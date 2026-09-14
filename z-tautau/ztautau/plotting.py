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
COLORS = {"DYtautau": "#ffcc66", "DYlowmass": "#e6a23c", "Fakes": "#ff99cc", "DYee": "#4a90d9", "DYmumu": "#1f5fa8",
          "WJets": "#d95f5f", "TTbar": "#9999cc", "SingleTop": "#6b6bb3", "WW": "#66cc99", "WZ": "#3fa877",
          "ZZ": "#2d7a56", "Diboson": "#66cc99", "Top": "#9999cc", "DYll": "#4a90d9"}
LABELS = {"DYtautau": r"Z$\rightarrow\tau\tau$", "DYlowmass": r"Z/$\gamma^*\rightarrow\ell\ell$ ($m$<50)", "Fakes": r"jet$\rightarrow\tau_h$ (FF, data)",
          "DYee": r"Z$\rightarrow$ee", "DYmumu": r"Z$\rightarrow\mu\mu$", "WJets": "W+jets", "TTbar": r"t$\bar{t}$",
          "SingleTop": "single t", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ", "Diboson": "diboson", "Top": "top",
          "DYll": r"Z$\rightarrow$ee/$\mu\mu$"}


def label(ax, extra: str | None = None):
    """'CMS Open Data' + luminosity header, drawn by hand: mplhep 1.1's label helpers fail to render with
    the matplotlib of LCG_110."""
    ax.text(0.0, 1.01, r"$\bf{CMS}$ $\it{Open\ Data}$", transform=ax.transAxes, fontsize=20, ha="left", va="bottom")
    ax.text(1.0, 1.01, rf"{LUMI_LABEL} fb$^{{-1}}$ (13 TeV)", transform=ax.transAxes, fontsize=17, ha="right",
            va="bottom")
    if extra:
        ax.text(0.04, 0.95, extra, transform=ax.transAxes, fontsize=15, va="top")


def stack_plot(path, edges, data, stack, xlabel, title=None, logy=False, band=None, ratio_range=(0.5, 1.5),
               density=False, ylabel="Events"):
    """data: (values, errors) or None; stack: list of (name, values, variances) drawn bottom-up;
    band: optional per-bin absolute uncertainty of the total prediction (else MC/fake stat)."""
    edges = np.asarray(edges)
    widths = np.diff(edges)
    scale = 1.0 / widths if density else np.ones(len(widths))
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(9, 9), gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
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
                     edgecolor="black", linewidth=0.5)
    unc = np.sqrt(total_var) if band is None else np.asarray(band)
    ax.fill_between(edges, np.append((total - unc) * scale, 0), np.append((total + unc) * scale, 0), step="post",
                    hatch="////", facecolor="none", edgecolor="gray", linewidth=0, label="Unc.")
    if data is not None:
        dv, de = data
        centers = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(centers, dv * scale, yerr=de * scale, fmt="o", color="black", label="Data", markersize=5)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(total > 0, dv / total, np.nan)
            re = np.where(total > 0, de / total, np.nan)
        rax.errorbar(centers, r, yerr=re, fmt="o", color="black", markersize=5)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(total > 0, unc / total, 0)
    rax.fill_between(edges, np.append(1 - rel, 1), np.append(1 + rel, 1), step="post", hatch="////",
                     facecolor="none", edgecolor="gray", linewidth=0)
    rax.axhline(1, color="black", linewidth=0.8)
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
    label(ax, title)
    fig.savefig(path, bbox_inches="tight", dpi=110)
    if str(path).endswith(".png"):
        fig.savefig(str(path)[:-4] + ".pdf", bbox_inches="tight")
    plt.close(fig)
