#!/usr/bin/env python
"""v4 figures for the slide deck (dark vector figures, slides/figs/v4_*.pdf + .png), docs/10-v4-plan.md.

    source ../setup.sh && python slides/make_figures_v4.py

Inputs: output/results.json, fit/fitinputs/ztautau.root (+ meta), fit/results/ztautau/ (post-fit
totals, ranking / pull plots copied as PNG), output/data/fakes_*.json, external/trigger_insitu_v4.json.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from ztautau import config, plotting  # noqa: E402
from ztautau.plotting import DARK_FG, dark_rc, save  # noqa: E402

FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)
BLUE, GREEN, ORANGE, PINK, MUTED, RED = "#2BA4DD", "#8AC63F", "#D9822B", "#E377AC", "#99999E", "#CC4C4C"
R = json.loads((config.OUTPUT_DIR_V4 / "results.json").read_text())
FIT = R["fit"]["combined"]
META = json.loads((config.FIT_DIR_V4 / "fitinputs" / f"{config.JOB_V4}.root.meta.json").read_text())
REGION_LABELS = {"tautau_SR0": r"$\tau_h\tau_h$, BDT < 0.55 (m > 110)", "tautau_SR1": r"$\tau_h\tau_h$, 0.55 < BDT < 0.90", "tautau_SR2": r"$\tau_h\tau_h$, BDT > 0.90",
                 "mutau_SR": r"$\mu\tau_h$", "etau_SR": r"$e\tau_h$", "emu_SR": r"$e\mu$", "emu_CRtt": r"$e\mu$ $t\bar{t}$ control region"}
for _ch, _lab in (("mutau", r"$\mu\tau_h$"), ("etau", r"$e\tau_h$")):
    for _dm in config.TAU_DMS:
        REGION_LABELS[f"{_ch}_SR_dm{_dm}"] = f"{_lab}, decay mode {_dm}"
STACK_GROUPS = [("DYll", ["DYee", "DYmumu"]), ("DYlowmass", ["DYlowmass"]), ("Diboson", ["WW", "WZ", "ZZ"]), ("Top", ["TTbar", "SingleTop"]),
                ("WJets", ["WJets"]), ("Fakes", ["Fakes"]), ("DYtautau_nonfid", ["DYtautau_out"]), ("DYtautau", ["DYtautau"])]
plt.rcParams.update({"font.size": 13})


def dark_fig(figsize=(7, 4.4)):
    plt.rcParams.update(dark_rc())
    fig, ax = plt.subplots(figsize=figsize)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    plotting.tag(ax, x=1.0, y=1.01, ha="right", va="bottom", color=DARK_FG, fontsize=10)
    return fig, ax


# ------------------------------------------------------------------ 1. cross-section summary
def summary():
    rows = [("combined: " + r"$\tau_h\tau_h + \mu\tau_h + e\tau_h + e\mu$", FIT["sigma_60_120_pb"], BLUE)]
    for ch, lab in (("mutau", r"$\mu\tau_h$ alone"), ("etau", r"$e\tau_h$ alone"), ("emu", r"$e\mu$ alone"), ("tautau", r"$\tau_h\tau_h$ alone")):
        if ch in R["fit"]:
            f = R["fit"][ch]
            rows.append((lab + (r" (POG $\tau_h$ ID SF)" if f.get("tau_id_fixed") and ch != "emu" else ""), f["sigma_60_120_pb"], MUTED))
    if R.get("v3_reference"):
        rows.append((r"v3: $\tau_h\tau_h$, POG $\tau_h$ ID SF", R["v3_reference"], PINK))
    fig, ax = dark_fig((8.4, 1.2 + 0.62 * len(rows)))
    pred = FIT["sigma_60_120_pb"]["prediction"]
    ax.axvspan(pred * 0.96, pred * 1.04, color=ORANGE, alpha=0.25)
    ax.axvline(pred, color=ORANGE, lw=1.5)
    ax.text(pred + 12, 0.35, f"NLO prediction {pred:.0f} pb ($\\pm$4%)", fontsize=9, color=ORANGE)
    for i, (lab, s, col) in enumerate(rows):
        yv = len(rows) - i
        ax.errorbar([s["value"]], [yv], xerr=[[s["err_down"]], [s["err_up"]]], fmt="o", color=col, capsize=4, lw=2, ms=6)
        if s.get("stat"):
            ax.errorbar([s["value"]], [yv], xerr=[[s["stat"]], [s["stat"]]], fmt="none", color=RED, lw=5)
        ax.text(1330, yv + 0.3, lab, fontsize=10.5, color=DARK_FG)
        ax.text(2870, yv + 0.3, f"{s['value']:.0f}  +{s['err_up']:.0f} $-${s['err_down']:.0f} pb", fontsize=10.5, color=DARK_FG, ha="right")
    ax.set_xlim(1300, 2900); ax.set_ylim(0.2, len(rows) + 0.95); ax.set_yticks([])
    ax.set_xlabel(r"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau,\ 60<m<120\ \mathrm{GeV})$ [pb]   (red: statistical)")
    save(fig, FIGS / "v4_summary", dark=True)


# ------------------------------------------------------------------ 2. fitted tau ID scale factors
def tau_id():
    sf = FIT.get("tau_id_sf") or {}
    fig, ax = dark_fig((6.4, 4.2))
    x = np.arange(len(sf))
    for i, (dm, v) in enumerate(sf.items()):
        ax.errorbar(i - 0.12, v["pog"][0], v["pog"][1], fmt="s", color=MUTED, ms=7, capsize=4, label="TauPOG (external)" if i == 0 else None)
        ax.errorbar(i + 0.12, v["value"], [[v["err_down"]], [v["err_up"]]], fmt="o", color=BLUE, ms=7, capsize=4, label="fitted in situ (v4)" if i == 0 else None)
    ax.set_xticks(x); ax.set_xticklabels([k for k in sf]); ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax.set_ylabel(r"DeepTau VSjet Tight scale factor"); ax.set_ylim(0.5, 1.4); ax.legend(loc="upper left", fontsize=10)
    ax.set_title(r"$\tau_h$ identification scale factors per decay mode", fontsize=12, pad=18)
    save(fig, FIGS / "v4_tauid", dark=True)
    es = FIT.get("tau_es") or {}
    if es:
        fig, ax = dark_fig((6.4, 3.6))
        for i, (dm, v) in enumerate(es.items()):
            ax.errorbar(i, v["pull"], v["constraint"], fmt="o", color=GREEN, ms=7, capsize=4)
            ax.text(i, v["pull"] + v["constraint"] + 0.12, f"{v['constraint'] * v['prior_pct']:.1f}%", ha="center", fontsize=10, color=DARK_FG)
        ax.axhspan(-1, 1, color=MUTED, alpha=0.15); ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xticks(range(len(es))); ax.set_xticklabels(list(es)); ax.set_ylim(-2.2, 2.4)
        ax.set_ylabel("pull (units of the 3% prior)"); ax.set_title(r"$\tau_h$ energy scale: post-fit pull and constraint (label: absolute)", fontsize=11)
        save(fig, FIGS / "v4_taues", dark=True)


# ------------------------------------------------------------------ 3. post-fit regions
def read_fitinputs():
    out = {}
    with uproot.open(config.FIT_DIR_V4 / "fitinputs" / f"{config.JOB_V4}.root") as f:
        for k in f.keys():
            name = k.split(";")[0]
            if name.count("__") == 1 and name != "meta_json":
                h = f[name]
                region, sample = name.split("__")
                base = sample.split("_tDM")[0]
                key = f"{region}__{base}"
                v, var = h.values(), h.variances()
                out[key] = (out[key][0] + v, out[key][1] + var) if key in out else (v, var)
    return out


def stack_from(hists, region):
    stack = []
    for g, members in STACK_GROUPS:
        parts = [hists[f"{region}__{m}"] for m in members if f"{region}__{m}" in hists]
        if parts:
            v = sum(p[0] for p in parts); var = sum(p[1] for p in parts)
            stack.append((g, np.maximum(v, 0), var))
    return stack


def regions():
    hists = read_fitinputs()
    for region in META["regions"]:
        edges = np.asarray(META["bins"][region])
        data = hists[f"{region}__Data"][0]
        stack = stack_from(hists, region)
        post = config.FIT_DIR_V4 / f"results/{config.JOB_V4}/Histograms/{region}_postFit.root"
        overlay = None
        if post.exists():
            with uproot.open(post) as f:
                tot = f["h_tot_postFit"]; overlay = (tot.values(), tot.errors(), "post-fit total")
        lab = REGION_LABELS.get(region, region)
        if region == "tautau_SR0":
            keep = edges[:-1] >= config.SIDEBAND_REGION_MTT_MIN
            drop = ~keep
            data = np.where(drop, 0, data)
            stack = [(g, np.where(drop, 0, v), np.where(drop, 0, var)) for g, v, var in stack]
            if overlay:
                overlay = (np.where(drop, 0, overlay[0]), np.where(drop, 0, overlay[1]), overlay[2])
        plotting.stack_plot(FIGS / f"v4_postfit_{region}.png", edges, (data, np.sqrt(data)), stack, r"$m_{\tau\tau}$ [GeV]",
                            title=f"post-fit total over the prefit stack, {lab}", dark=True, overlay=overlay, density=True)
        plotting.stack_plot(FIGS / f"v4_prefit_{region}_log.png", edges, (data, np.sqrt(data)), stack, r"$m_{\tau\tau}$ [GeV]",
                            title=f"prefit, {lab}", dark=True, logy=True)


# ------------------------------------------------------------------ 4. impacts, ranking (TRExFitter PNGs)
def impacts():
    gi = {k: v for k, v in FIT["grouped_impact"].items() if k != "FullSyst"}
    if FIT.get("mu_stat"):
        gi["Data statistics"] = FIT["mu_stat"]
    items = sorted(gi.items(), key=lambda kv: kv[1])
    fig, ax = dark_fig((7, 0.33 * len(items) + 1.4))
    ax.barh([k for k, _ in items], [100 * v for _, v in items], color=BLUE)
    for i, (_, val) in enumerate(items):
        ax.text(100 * val + 0.05, i, f"{100 * val:.2f}%", va="center", fontsize=9, color=DARK_FG)
    ax.set_xlabel(r"impact on $\mu_Z$ [%]"); ax.set_title("grouped uncertainties, four-channel fit", fontsize=12)
    save(fig, FIGS / "v4_impacts", dark=True)
    base = config.FIT_DIR_V4 / f"results/{config.JOB_V4}"
    for src, dst in (("Rankings/RankingSysts_mu_Z_systs.png", "v4_ranking.png"), ("NuisPar.png", "v4_pulls.png"), ("NormFactors.png", "v4_normfactors.png"),
                     ("Plots/Summary_postFit.png", "v4_summary_postfit.png"), ("CorrMatrix.png", "v4_corrmatrix.png")):
        cand = sorted(base.glob(f"**/{Path(src).name}"))
        if cand:
            shutil.copy(cand[0], FIGS / dst)


# ------------------------------------------------------------------ 5. fakes and triggers
def fakes():
    for ch in ("mutau", "etau"):
        p = config.DATA_DIR_V4 / f"fakes_{ch}.json"
        if not p.exists():
            continue
        fk = json.loads(p.read_text())
        pt = np.asarray(fk["tables"]["qcd"]["pt_bins"]); x = 0.5 * (pt[1:] + pt[:-1]); x[-1] = pt[-2] + 15
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
        plt.rcParams.update(dark_rc())
        for ax, dm in zip(axes, (1, 10)):
            i = list(config.TAU_DMS).index(dm)
            for proc, col, lab in (("qcd", BLUE, "multijet (SS, pure)"), ("w", RED, "W+jets (OS, $m_T>70$)"), ("w_ss", ORANGE, "W+jets (SS)"), ("tt", GREEN, r"$t\bar{t}$ (simulation)")):
                t = fk["tables"][proc]
                ax.errorbar(x, np.asarray(t["ff"])[i, 0], np.asarray(t["err"])[i, 0], fmt="o", color=col, ms=4, label=lab if dm == 1 else None)
            ax.set_title(f"decay mode {dm}, 0 jets", color=DARK_FG); ax.set_xlim(30, 120); ax.set_ylim(0, 0.25); ax.grid(alpha=0.2)
            ax.set_xlabel(r"$p_T(\tau_h)$ [GeV]")
        axes[0].set_ylabel("fake factor (Tight / VVVLoose-not-Tight)"); axes[0].legend(fontsize=9)
        for ax in axes:
            ax.set_facecolor(plotting.DARK_BG)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        plotting.tag(axes[1], x=1.0, y=1.01, ha="right", va="bottom", color=DARK_FG, fontsize=10)
        fig.patch.set_facecolor(plotting.DARK_BG)
        fig.suptitle(f"{ch}: fake factors per process; C(OS/SS) = {fk['osss']['C']:.2f} ± {fk['osss']['stat']:.2f}, same-sign closure {fk['closure_ss']['ratio']:.3f} ± {fk['closure_ss']['stat']:.3f}", color=DARK_FG, fontsize=11)
        fig.tight_layout()
        save(fig, FIGS / f"v4_ff_{ch}", dark=True)
    p = config.DATA_DIR_V4 / "fakes_emu.json"
    if p.exists():
        fk = json.loads(p.read_text())
        e = np.asarray(fk["osss"]["dr_edges"]); x = 0.5 * (e[1:] + e[:-1])
        fig, ax = dark_fig((6.4, 4))
        ax.errorbar(x, fk["SB1"]["ratio"], fk["SB1"]["stat"], fmt="o", color=BLUE, label="SB1: both $I_{rel}$ < 0.5, one > 0.15")
        ax.errorbar(x + 0.06, fk["SB2"]["ratio"], fk["SB2"]["stat"], fmt="s", color=ORANGE, label="SB2: one $I_{rel}$ > 0.3")
        ax.set_xlabel(r"$\Delta R(e,\mu)$"); ax.set_ylabel("multijet OS / SS"); ax.set_ylim(0, 3.5); ax.legend(fontsize=9)
        ax.set_title(f"e$\\mu$ multijet: {fk['ss_region']['data']} same-sign data, {fk['ss_region']['mc']:.0f} simulated -> {fk['sr_multijet']:.0f} in the SR", fontsize=11)
        save(fig, FIGS / "v4_emu_osss", dark=True)


def triggers():
    p = HERE.parent / "external" / "trigger_insitu_v4.json"
    if not p.exists():
        return
    t = json.loads(p.read_text())
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    plt.rcParams.update(dark_rc())
    for ax, (name, title) in zip(axes, (("ele27", "Ele27_WPTight (e$\\tau_h$)"), ("emu_e", "e$\\mu$ cross trigger, electron legs"), ("emu_mu", "e$\\mu$ cross trigger, muon legs"))):
        tb = t[name]; xe = np.asarray(tb["x_edges"]); x = 0.5 * (xe[1:] + xe[:-1]); x[-1] = xe[-2] * 1.3
        for j in range(len(tb["y_edges"]) - 1):
            ax.errorbar(x * (1 + 0.02 * j), np.asarray(tb["sf"])[:, j], np.asarray(tb["err"])[:, j], fmt="o", ms=4, color=[BLUE, GREEN, ORANGE, PINK][j],
                        label=f"|η| {tb['y_edges'][j]}-{tb['y_edges'][j + 1]}")
        ax.set_xscale("log"); ax.set_ylim(0.7, 1.2); ax.axhline(1, color=MUTED, lw=0.8, ls=":"); ax.set_title(title, color=DARK_FG, fontsize=11)
        ax.set_xlabel(f"{tb['x'].replace('_pt', '')} $p_T$ [GeV]"); ax.set_facecolor(plotting.DARK_BG); ax.grid(alpha=0.2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("SF = data / simulation (e$\\mu$ events, in situ)"); axes[0].legend(fontsize=8)
    plotting.tag(axes[2], x=1.0, y=1.01, ha="right", va="bottom", color=DARK_FG, fontsize=10)
    fig.patch.set_facecolor(plotting.DARK_BG); fig.tight_layout()
    save(fig, FIGS / "v4_trigger", dark=True)


if __name__ == "__main__":
    summary(); tau_id(); regions(); impacts(); fakes(); triggers()
    print("v4 figures ->", FIGS)
