"""Histogram booking, persistence and CMS-style plotting.

Every claim in this analysis is backed by a histogram, so the plotting helpers
are shared rather than rewritten per step. Histograms are `hist.Hist` objects
and are persisted with pickle into output/data/<step>.pkl.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")            # headless: we only ever write files

import hist
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

from . import config

plt.style.use(hep.style.CMS)


# --------------------------------------------------------------------------
# Booking
# --------------------------------------------------------------------------
def mass_hist(name: str = "mass", lo: float = None, hi: float = None, nbins: int = None):
    """Dimuon invariant-mass histogram with the analysis binning."""
    lo = config.MASS_LO if lo is None else lo
    hi = config.MASS_HI if hi is None else hi
    nbins = config.MASS_NBINS if nbins is None else nbins
    return hist.Hist(hist.axis.Regular(nbins, lo, hi, name=name,
                                       label=r"$m_{\mu\mu}$ [GeV]"),
                     storage=hist.storage.Double())


def regular(nbins, lo, hi, label, name="x"):
    return hist.Hist(hist.axis.Regular(nbins, lo, hi, name=name, label=label),
                     storage=hist.storage.Double())


def variable(edges, label, name="x"):
    return hist.Hist(hist.axis.Variable(list(edges), name=name, label=label),
                     storage=hist.storage.Double())


def hist2d(xedges, yedges, xlabel, ylabel, xname="x", yname="y"):
    """Two-dimensional histogram on variable binning (used for eff. maps)."""
    return hist.Hist(
        hist.axis.Variable(list(xedges), name=xname, label=xlabel),
        hist.axis.Variable(list(yedges), name=yname, label=ylabel),
        storage=hist.storage.Double(),
    )


def eff_hist_pair(edges, label, name="x"):
    """(pass, total) histogram pair sharing one variable axis."""
    return variable(edges, label, name), variable(edges, label, name)


# --------------------------------------------------------------------------
# Persistence
# --------------------------------------------------------------------------
def save(obj, filename: str) -> Path:
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = config.DATA_DIR / filename
    with open(path, "wb") as fh:
        pickle.dump(obj, fh)
    return path


def load(filename: str):
    path = config.DATA_DIR / filename
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing -- run the step that produces it first "
            f"(see README.md for the order)."
        )
    with open(path, "rb") as fh:
        return pickle.load(fh)


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
def _decorate(ax, lumi_fb: float = None):
    lumi_fb = config.LUMI_PB / 1000.0 if lumi_fb is None else lumi_fb
    hep.cms.label("Open Data", data=True, lumi=round(lumi_fb, 1), year=2016, ax=ax)


def _title(ax, text, y=1.012):
    """Caption above the axes, left-aligned.

    Not `ax.set_title`: mplhep puts the CMS label there, and the two collide.
    """
    if text:
        ax.text(0.0, y, text, transform=ax.transAxes, fontsize=13,
                style="italic", va="bottom", ha="left")


def save_fig(fig, filename: str) -> Path:
    config.PLOT_DIR.mkdir(parents=True, exist_ok=True)
    path = config.PLOT_DIR / filename
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_mass(hists: dict, filename: str, title: str = "", logy: bool = True,
              ylabel: str = "Events / 0.5 GeV"):
    """Overlay one or more mass histograms. `hists` maps label -> Hist."""
    fig, ax = plt.subplots(figsize=(9, 7))
    for label, h in hists.items():
        hep.histplot(h, ax=ax, label=label, histtype="step", linewidth=1.6)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=13)
    _decorate(ax)
    _title(ax, title)
    return save_fig(fig, filename)


def plot_efficiency(edges, eff, lo, hi, filename, xlabel, ylabel="Efficiency",
                    title="", ylim=(0.0, 1.05), reference=None):
    """Efficiency vs a binned variable, with asymmetric Clopper-Pearson bars."""
    edges = np.asarray(edges, dtype=float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    widths = 0.5 * np.diff(edges)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.errorbar(centres, eff, yerr=[lo, hi], xerr=widths, fmt="o",
                color="black", markersize=5, capsize=0, linewidth=1.3)
    if reference is not None:
        ax.axhline(reference, linestyle="--", color="crimson", linewidth=1.3,
                   label=f"average = {reference:.4f}")
        ax.legend(fontsize=13)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(*ylim)
    _decorate(ax)
    _title(ax, title)
    return save_fig(fig, filename)


def plot_cutflow(labels, counts, filename, title="Cutflow"):
    """Horizontal cutflow with absolute and relative survival printed."""
    counts = np.asarray(counts, dtype=float)
    fig, ax = plt.subplots(figsize=(10, 0.62 * len(labels) + 3))
    ypos = np.arange(len(labels))[::-1]
    ax.barh(ypos, counts, color="#3f7fb5", height=0.68)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xscale("log")
    ax.set_xlabel("Events")
    span = counts.max() / max(counts.min(), 1.0)
    for y, n, prev in zip(ypos, counts, np.concatenate([[counts[0]], counts[:-1]])):
        frac = 100.0 * n / counts[0] if counts[0] else 0.0
        rel = 100.0 * n / prev if prev else 0.0
        ax.text(n * (1.0 + 0.08 * np.log10(span)), y, f"{n:,.0f}  ({frac:.2f}% abs, {rel:.2f}% rel)",
                va="center", fontsize=10)
    ax.set_xlim(counts.min() * 0.5, counts.max() * 12)
    _decorate(ax)
    _title(ax, title)
    return save_fig(fig, filename)


def plot_stack_with_ratio(data_h, components: dict, filename, title="",
                          xlabel=r"$m_{\mu\mu}$ [GeV]", logy=True):
    """Data with stacked background estimates and a data/(bkg) ratio panel."""
    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(9, 9), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07},
    )
    labels = list(components)
    stack = [components[k] for k in labels]
    if stack:
        hep.histplot(stack, ax=ax, stack=True, histtype="fill",
                     label=labels, alpha=0.85, edgecolor="black", linewidth=0.5)
    hep.histplot(data_h, ax=ax, histtype="errorbar", color="black",
                 label="Data", markersize=4)
    ax.set_ylabel("Events / 0.5 GeV")
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=12, ncol=2)
    _decorate(ax)
    _title(ax, title)

    total = np.zeros(data_h.shape)
    for h in stack:
        total = total + h.values()
    data_vals = data_h.values()
    centres = data_h.axes[0].centers
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(total > 0, data_vals / total, np.nan)
        err = np.where(total > 0, np.sqrt(np.maximum(data_vals, 0)) / total, np.nan)
    rax.errorbar(centres, ratio, yerr=err, fmt="o", color="black", markersize=3)
    rax.axhline(1.0, linestyle="--", color="grey")
    rax.set_ylabel("Data / bkg")
    rax.set_xlabel(xlabel)
    rax.set_yscale("log")
    return save_fig(fig, filename)
