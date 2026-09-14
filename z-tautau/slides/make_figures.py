#!/usr/bin/env python
"""Figures made only for the slides (everything else is taken from output/plots/).

    python slides/make_figures.py

fig_mass_estimators.pdf   m_vis, m_col, m_tt divided by the generator ttau mass (Z -> tau_h tau_h simulation)
fig_ff_closure.pdf        same-sign closure vs N_jets and per era: first iteration (DM x pT FF) vs final binning
fig_ff_dm1.pdf            fake factors of the most common decay mode (DM1) per era and jet multiplicity
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from ztautau import analysis, config, fakes  # noqa: E402

plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})


def mass_estimators():
    d, _ = analysis.load("DY_NLO")
    r = analysis.regions(d, is_mc=True)
    w = analysis.weights(d, "DY_NLO")
    sel = r["SR"] & (d["gen_lhe_flavour"] == 15) & (d["t1_genflav"] == 5) & (d["t2_genflav"] == 5)
    sel &= (d["gen_mll_lhe"] > 70) & (d["gen_mll_lhe"] < 110)
    ref = d["gen_mll_lhe"][sel]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    bins = np.linspace(0.3, 2.0, 69)
    for var, label, color in (("m_vis", r"$m_\mathrm{vis}$", "#4a90d9"), ("m_col", r"$m_\mathrm{col}$ (38 % defined)", "#999999"),
                              ("m_tt", r"$m_{\tau\tau}$ (MET likelihood)", "#d9822b")):
        x = d[var][sel] / ref
        ok = d[var][sel] > 0
        med = np.median(x[ok])
        q25, q75 = np.percentile(x[ok], [25, 75])
        ax.hist(x[ok], bins=bins, weights=w[sel][ok], histtype="step", lw=2, color=color, density=True,
                label=f"{label}: median {med:.2f}, IQR/2/med {100 * (q75 - q25) / 2 / med:.0f}%")
    ax.axvline(1.0, color="black", lw=0.8, ls="--")
    ax.set_xlabel(r"reconstructed mass / generator $m_{\tau\tau}$")
    ax.set_ylabel("normalised")
    ax.set_title(r"Z$\rightarrow\tau_h\tau_h$ simulation, signal region, $70<m_{\tau\tau}^{gen}<110$ GeV", fontsize=11)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.45)
    ax.legend(fontsize=9.5, frameon=False, loc="upper right")
    fig.savefig(HERE / "fig_mass_estimators.pdf", bbox_inches="tight")
    plt.close(fig)


def ff_closure():
    data = analysis.load_data()
    reg = analysis.regions(data)
    # iteration 1: FF binned in DM x pT only (both eras, all jets together)
    pt_bins = np.asarray([40, 45, 50, 55, 60, 70, 80, 100, 150, 1000.0])
    ipt = np.clip(np.searchsorted(pt_bins, data["t1_pt"], side="right") - 1, 0, len(pt_bins) - 2)
    idm = np.searchsorted(np.asarray(config.TAU_DMS), data["t1_dm"])
    num = np.zeros((4, len(pt_bins) - 1))
    den = np.zeros_like(num)
    np.add.at(num, (idm[reg["SS_T"]], ipt[reg["SS_T"]]), 1)
    np.add.at(den, (idm[reg["SS_L"]], ipt[reg["SS_L"]]), 1)
    ff1 = num / np.maximum(den, 1)
    w1 = np.where(reg["SS_L"], ff1[np.clip(idm, 0, 3), ipt], 0.0)
    ffres = json.loads((config.DATA_DIR / "fakefactors.json").read_text())
    table = fakes.from_json(ffres["nominal"]["ff"])
    w3 = fakes.fake_weights(data, reg["SS_L"], table, 1.0)
    cats = [("0 jets", data["njets"] == 0), ("1 jet", data["njets"] == 1), ("2 jets", data["njets"] == 2),
            (r"$\geq$3 jets", data["njets"] >= 3), ("era G", data["era"] == 0), ("era H", data["era"] == 1)]
    r1, r3, err = [], [], []
    for _, m in cats:
        obs = (reg["SS_T"] & m).sum()
        r1.append(obs / w1[m].sum())
        r3.append(obs / w3[m].sum())
        err.append(1 / np.sqrt(obs))
    x = np.arange(len(cats))
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.axhline(1, color="black", lw=0.8)
    ax.axvline(3.5, color="#bbbbbb", lw=0.8)
    ax.errorbar(x - 0.12, r1, yerr=err, fmt="s", color="#999999", ms=7, label=r"iteration 1: FF(DM, $p_T$)")
    ax.errorbar(x + 0.12, r3, yerr=err, fmt="o", color="#d9822b", ms=7, label=r"final: FF(era, DM, $N_{jets}$, $p_T$)")
    ax.set_xticks(x, [c for c, _ in cats], fontsize=11)
    ax.set_ylabel("same-sign data / FF prediction")
    ax.set_ylim(0.65, 1.2)
    ax.legend(fontsize=10, frameon=False, loc="lower left")
    fig.savefig(HERE / "fig_ff_closure.pdf", bbox_inches="tight")
    plt.close(fig)
    print("closure iteration 1:", np.round(r1, 3), "final:", np.round(r3, 3))


def ff_dm1():
    t = json.loads((config.DATA_DIR / "fakefactors.json").read_text())["nominal"]["ff"]
    ff, err = np.asarray(t["ff"]), np.asarray(t["err"])
    edges = np.asarray(config.FF_PT_BINS[:-1] + [110.0])
    centers = 0.5 * (edges[1:] + edges[:-1])
    i = t["dms"].index(1)
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.7), sharey=True)
    for e, (ax, era) in enumerate(zip(axes, t["eras"])):
        for j, (lab, col) in enumerate(zip(["0 jets", "1 jet", r"$\geq$2 jets"], ["black", "#d9822b", "#4a90d9"])):
            ax.errorbar(centers + 1.2 * (j - 1), ff[e, i, j], yerr=err[e, i, j], xerr=np.diff(edges) / 2, fmt="o", ms=4,
                        color=col, label=lab, elinewidth=1)
        ax.set_title(f"Run2016{era}, decay mode 1", fontsize=11)
        ax.set_xlabel(r"$p_T(\tau_1)$ [GeV]  (last bin > 80)", fontsize=10)
        ax.set_ylim(0.1, 0.32)
    axes[0].set_ylabel("fake factor", fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=9, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.1))
    fig.savefig(HERE / "fig_ff_dm1.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    mass_estimators()
    ff_closure()
    ff_dm1()
