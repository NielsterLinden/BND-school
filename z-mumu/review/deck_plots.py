#!/usr/bin/env python
"""Figures for the Z->mumu review (REVIEW.md) and the LaTeX deck, on the dark deck background.

    source ../../setup.sh && python deck_plots.py [--pass review_pass.pkl]

Reads the v2 outputs (output/v2/histograms.pkl, tnp/fakes/momentum JSON, fit results) and the
review skim pass (review/review_pass.pkl, produced by review/review_pass.py). Writes
review/figures/<name>.pdf and .png (vector PDF for the deck).
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

from zmumu import config, histograms as H, tnp as TNP, momentum as MOM, fakes as FK

OUT_V2 = ROOT / "output" / "v2"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

BG, FG, MUTED, BLUE, GREEN, RULE = "#222222", "#E6E6E6", "#99999E", "#2BA4DD", "#8AC63F", "#8C8C94"
plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "axes.edgecolor": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG, "text.color": FG,
    "legend.frameon": False, "axes.titlecolor": FG, "font.size": 15, "legend.fontsize": 12,
    "axes.labelsize": 17, "xtick.labelsize": 14, "ytick.labelsize": 14, "hatch.color": MUTED,
})
LUMI_FB = config.LUMI_PB_NORMTAG / 1000.0

STACK_ORDER = ["Fakes", "WW", "WZ", "ZZ", "SingleTop", "TTbar", "DYee", "DYtautau", "DYmumu"]
COLOURS = {"DYmumu": "#f2b134", "DYtautau": "#a04cb0", "DYee": "#4c9ee0", "TTbar": "#d62728", "SingleTop": "#e07b7b",
           "WW": "#2ca6a4", "WZ": "#4dc0be", "ZZ": "#7fd4d2", "Fakes": "#9a9a9a", "NonPrompt": "#c8c8c8"}
LABELS = {"DYmumu": r"Z/$\gamma^*\to\mu\mu$", "DYtautau": r"Z/$\gamma^*\to\tau\tau$", "DYee": r"Z/$\gamma^*\to ee$",
          "TTbar": r"t$\bar{t}$", "SingleTop": "tW", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ", "Fakes": "non-prompt (FF)",
          "NonPrompt": "MC non-prompt (removed)"}
XL = {"mass_fit": r"$m_{\mu\mu}$ [GeV]", "mass_fine": r"$m_{\mu\mu}$ [GeV]", "pt1": r"leading muon $p_T$ [GeV]",
      "pt2": r"subleading muon $p_T$ [GeV]", "eta1": r"leading muon $\eta$", "zpt": r"$p_T^{\mu\mu}$ [GeV]",
      "zy": r"$y^{\mu\mu}$", "npv": "good primary vertices", "met": r"$p_T^{miss}$ [GeV]",
      "njet": r"jets ($p_T$ > 30 GeV, not lepton-cleaned)", "njet_clean": r"jets ($p_T$ > 30 GeV, $\Delta R(\ell)$ > 0.4)",
      "iso1": "leading muon rel. iso.", "nfsr": "FSR photons", "pt_el": r"electron $p_T$ [GeV]", "mass": r"$m_{e\mu}$ [GeV]"}
MCS = ("DYmumu", "DYee", "DYtautau", "TTbar", "SingleTop", "WW", "WZ", "ZZ")


def label(ax, title="", lumi=True):
    ax.text(0.0, 1.005, r"$\bf{CMS}$ $\it{Open\ Data}$", transform=ax.transAxes, fontsize=15, va="bottom", ha="left", color=FG)
    if lumi:
        ax.text(1.0, 1.005, f"{LUMI_FB:.1f} fb$^{{-1}}$ (13 TeV, 2016)", transform=ax.transAxes, fontsize=12, va="bottom", ha="right", color=FG)
    if title:
        ax.text(0.0, 1.06, title, transform=ax.transAxes, fontsize=13, style="italic", va="bottom", ha="left", color=GREEN)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight", facecolor=BG, dpi=170 if ext == "png" else None)
    plt.close(fig)
    print("  wrote", name)


def stack(ax, rax, e, comps, labels, cols, data, ylabel, xlabel, logy=False, ratio_range=(0.8, 1.2), extra=None):
    c = 0.5 * (e[1:] + e[:-1])
    tot = np.sum(comps, axis=0); tot_var = extra["tot_var"] if extra and "tot_var" in extra else np.zeros_like(tot)
    unit = 1e3 if (not logy and max(tot.max(), data.max()) >= 2e4) else 1.0
    hep.histplot([x / unit for x in comps], bins=e, ax=ax, stack=True, histtype="fill", label=labels, color=cols, edgecolor=BG, linewidth=0.4)
    ax.fill_between(e, np.append(tot - np.sqrt(tot_var), 0) / unit, np.append(tot + np.sqrt(tot_var), 0) / unit, step="post",
                    facecolor="none", hatch="////", edgecolor=MUTED, linewidth=0, label="MC stat.")
    ax.errorbar(c, data / unit, yerr=np.sqrt(np.maximum(data, 0)) / unit, fmt="o", color=FG, markersize=3.5, label="Data")
    ax.set_ylabel(ylabel + (r" [$\times 10^3$]" if unit != 1 else ""))
    if logy:
        ax.set_yscale("log"); ax.set_ylim(max(tot.min() * 0.1, 0.5), tot.max() * 40)
    else:
        ax.set_ylim(0, max(tot.max(), data.max()) / unit * 1.4)
    ax.legend(ncol=2, loc="upper right", fontsize=11)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(tot > 0, data / tot, np.nan); er = np.where(tot > 0, np.sqrt(np.maximum(data, 0)) / tot, np.nan)
        band = np.where(tot > 0, np.sqrt(tot_var) / tot, 0)
    rax.fill_between(e, np.append(1 - band, 1), np.append(1 + band, 1), step="post", facecolor="none", hatch="////", edgecolor=MUTED, linewidth=0)
    rax.errorbar(c, r, yerr=er, fmt="o", color=FG, markersize=3)
    rax.axhline(1, ls="--", color=MUTED)
    rax.set_ylim(*ratio_range); rax.set_ylabel("Data / pred."); rax.set_xlabel(xlabel)


def datamc(hall, region, var, name, fakes_fine=None, logy=False, title="", ratio_range=(0.8, 1.2)):
    e = H.edges(region, var)
    data = hall[f"Data|{region}|{var}|nominal"]
    comps, labels, cols, tot_var = [], [], [], np.zeros(len(e) - 1)
    for smp in STACK_ORDER:
        k = f"{smp}|{region}|{var}|nominal"
        if smp == "Fakes":
            if region == "SR" and fakes_fine is not None and var in ("mass_fit", "mass_fine"):
                f = np.asarray(fakes_fine); v = np.clip(f.reshape(-1, len(f) // (len(e) - 1)).sum(1), 0, None)
            else:
                continue
        elif k not in hall or hall[k].sum() <= 0:
            continue
        else:
            v = hall[k]; tot_var += hall[k + "|w2"]
        comps.append(v); labels.append(LABELS[smp]); cols.append(COLOURS[smp])
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8.5, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    gev = " GeV" if var in ("mass_fit", "mass_fine", "pt1", "pt2", "zpt", "met", "pt_el") else ""
    stack(ax, rax, e, comps, labels, cols, data, f"Events / {e[1]-e[0]:g}{gev}", XL.get(var, var), logy, ratio_range, {"tot_var": tot_var})
    label(ax, title or {"SR": "signal region", "SS": "same-sign region", "CRemu": r"e$\mu$ region (validation)"}[region])
    save(fig, name)


# --------------------------------------------------------------------------------------- review figures
def emu_review(R):
    """e-mu OS: prompt MC + MC non-prompt vs data; SS: data vs MC."""
    edges = {"pt_el": np.linspace(20, 200, 37), "met": np.linspace(0, 150, 31), "njet_clean": np.arange(0, 9, 1.0), "mass": np.linspace(60, 120, 13), "npv": np.linspace(0, 50, 51)}
    mcs = ("WW", "WZ", "ZZ", "SingleTop", "TTbar", "DYee", "DYtautau", "DYmumu")
    npsrc = ("WJets", "DYmumu", "DYee", "DYtautau", "TTbar", "SingleTop", "WW", "WZ", "ZZ")
    g = lambda s, ch, cat, var: np.asarray(R.get(f"emu|{s}|{ch}|{cat}|{var}", np.zeros(len(edges[var]) - 1)))
    for charge, name in (("os", "emu_os_with_nonprompt"), ("ss", "emu_ss")):
        fig, axes = plt.subplots(2, 4, figsize=(24, 9.5), sharex="col", gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06, "wspace": 0.28})
        for k, var in enumerate(("pt_el", "met", "njet_clean", "mass")):
            e = edges[var]; ax, rax = axes[0, k], axes[1, k]
            comps, labels, cols = [], [], []
            for s in mcs:
                v = g(s, charge, "pp", var)
                if v.sum() > 0:
                    comps.append(v); labels.append(LABELS[s]); cols.append(COLOURS[s])
            npv = sum(g(s, charge, "np", var) for s in npsrc)
            comps.insert(0, npv); labels.insert(0, LABELS["NonPrompt"]); cols.insert(0, COLOURS["NonPrompt"])
            data = g("Data", charge, "all", var)
            xl = XL[var]
            stack(ax, rax, e, comps, labels, cols, data, f"Events / bin", xl, logy=False, ratio_range=(0.5, 1.5) if charge == "ss" else (0.8, 1.2))
            if k:
                ax.legend().remove()
            if var in ("njet_clean",):
                ax.set_xlim(0, 6)
            if var == "pt_el":
                ax.set_xlim(20, 120)
            label(ax, ("opposite-sign e$\\mu$: prompt MC + MC non-prompt" if charge == "os" else "same-sign e$\\mu$: prompt MC + MC non-prompt") if k == 0 else "", lumi=(k == 3))
        save(fig, name)
    # composition bar chart of the OS non-prompt (MC) and the SS-based estimate
    comp = {}
    for s in npsrc:
        v = g(s, "os", "np", "mass").sum()
        if v > 0:
            comp[s] = v
    fig, ax = plt.subplots(figsize=(9, 5.5))
    names = sorted(comp, key=comp.get, reverse=True)
    ax.barh(range(len(names)), [comp[n] for n in names], color=[COLOURS.get(n, "#bbbbbb") for n in names])
    ax.set_yticks(range(len(names))); ax.set_yticklabels([{"WJets": "W+jets (jet → e)"}.get(n, LABELS.get(n, n)) for n in names])
    ax.invert_yaxis(); ax.set_xlabel("MC non-prompt events in the OS e$\\mu$ region")
    for i, n in enumerate(names):
        ax.text(comp[n], i, f"  {comp[n]:.0f}", va="center", fontsize=12, color=FG)
    ax.set_xlim(0, max(comp.values()) * 1.25)
    label(ax, "what the prompt-prompt requirement removes from the e$\\mu$ MC")
    save(fig, "emu_nonprompt_composition")
    return comp


def pileup_review(R, pu):
    """SR npv: data vs MC with the nominal profile and a 2 % higher profile."""
    w_nom = np.array(pu["weights"]["nominal"]); mcp = np.array(pu["mc_profile"]); mcn = mcp / mcp.sum()
    Hh = R["sr|DYmumu|ntrue_npv"]; dn = R["sr|Data|npv"]
    npvc = np.arange(60) + 0.5; e = np.arange(0, 61, 1.0)

    def weights_for(f):
        src = np.array(pu["data"]["nominal"]["hist"]); newp = np.zeros(100)
        for i, v in enumerate(src):
            newp[min(int((i + 0.5) * f), 99)] += v
        d = newp / newp.sum(); w = np.where(mcn > 0, d / np.where(mcn > 0, mcn, 1), 0.0); w = w / np.sum(w * mcn)
        return np.clip(w, 0, 10)
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8.5, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    hep.histplot(dn / 1e3, bins=e, ax=ax, histtype="errorbar", color=FG, label="Data", markersize=4)
    for f, col, lab in ((None, MUTED, "MC, no pileup weight"), (1.0, BLUE, r"MC, nominal profile ($\langle\mu\rangle$ = 24.7)"), (1.03, GREEN, r"MC, profile $\times$ 1.03")):
        w = np.ones(100) if f is None else weights_for(f)
        m = (Hh * w[:, None]).sum(0); m = m * dn.sum() / m.sum()
        hep.histplot(m / 1e3, bins=e, ax=ax, histtype="step", color=col, label=lab, linewidth=1.8)
        rax.plot(npvc, dn / np.maximum(m, 1), "o", color=col, markersize=3)
    ax.set_ylabel(r"Events / 1 [$\times 10^3$]"); ax.legend(fontsize=11); ax.set_xlim(0, 45); ax.set_ylim(0, dn.max() / 1e3 * 1.4)
    rax.axhline(1, ls="--", color=MUTED); rax.set_ylim(0.7, 1.3); rax.set_ylabel("Data / MC"); rax.set_xlabel("good primary vertices")
    label(ax, "signal region: pileup reweighting check")
    save(fig, "sr_npv_pileup_check")


def met_review(R, pu):
    w_nom = np.array(pu["weights"]["nominal"])
    metc = np.arange(30) * 5 + 2.5; e = np.linspace(0, 150, 31)
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
    for ax, (name, Dk, Mk) in zip(axes, (("PF $p_T^{miss}$", "sr|Data|npv_jet_met", "sr|DYmumu|ntrue_jet_met"),
                                        ("Puppi $p_T^{miss}$", "sr|Data|npv_jet_puppimet", "sr|DYmumu|ntrue_jet_puppimet"))):
        D, M = R[Dk], R[Mk]
        for j, col, lab in ((0, BLUE, "0 jets"), (1, GREEN, "1 jet"), (2, "#e07b7b", r"$\geq$ 2 jets")):
            d = D[:, j, :].sum(0); m = (M[:, j, :] * w_nom[:, None]).sum(0); m = m * d.sum() / m.sum()
            ax.errorbar(metc, d / np.maximum(m, 1e-9), yerr=np.sqrt(np.maximum(d, 0)) / np.maximum(m, 1e-9), fmt="o", color=col, markersize=4, label=lab + f" (data {d.sum()/D.sum():.0%})")
        ax.axhline(1, ls="--", color=MUTED); ax.set_ylim(0.6, 1.8); ax.set_xlabel(name + " [GeV]"); ax.set_ylabel("Data / MC (shape)")
        ax.legend(title="lepton-cleaned jets, $p_T$ > 30 GeV", fontsize=11, title_fontsize=11)
        label(ax, f"signal region, DY MC: {name} by jet multiplicity")
    save(fig, "sr_met_by_jets")
    # PF MET by pileup, 0 jets
    Hh = R["sr|DYmumu|ntrue_npv"]; prof = Hh.sum(1); avg_npv = (Hh * (np.arange(60) + 0.5)[None, :]).sum(1) / np.maximum(prof, 1e-9)
    D, M = R["sr|Data|npv_jet_met"], R["sr|DYmumu|ntrue_jet_met"]
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    for (lo, hi), col in (((0, 15), BLUE), ((15, 22), GREEN), ((22, 60), "#e07b7b")):
        d = D[lo:hi, 0, :].sum(0); sel = (avg_npv >= lo) & (avg_npv < hi); m = (M[sel, 0, :] * w_nom[sel][:, None]).sum(0); m = m * d.sum() / m.sum()
        ax.errorbar(metc, d / np.maximum(m, 1e-9), yerr=np.sqrt(np.maximum(d, 0)) / np.maximum(m, 1e-9), fmt="o", color=col, markersize=4, label=f"{lo} $\\leq$ N$_{{PV}}$ < {hi}")
    ax.axhline(1, ls="--", color=MUTED); ax.set_ylim(0.6, 1.8); ax.set_xlabel(r"PF $p_T^{miss}$ [GeV]"); ax.set_ylabel("Data / MC (shape)"); ax.legend(fontsize=11)
    label(ax, "signal region, 0 jets: PF $p_T^{miss}$ by pileup")
    save(fig, "sr_met_by_pileup")


def sigmodel_review(hall):
    import uproot
    f = uproot.open(ROOT / "fit" / "fitinputs" / "zmumu.root")
    e = np.linspace(60, 120, 31); c = 0.5 * (e[1:] + e[:-1])
    n = f["mumu_SR__DYmumu"].values(); d = f["mumu_SR__Data"].values()
    tot = sum(f[k].values() for k in f.keys() if k.startswith("mumu_SR__") and k.count("__") == 1 and "Data" not in k)
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.errorbar(c, d / tot, yerr=np.sqrt(d) / tot, fmt="o", color=FG, markersize=4, label="data / pre-fit prediction")
    for k, col, lab in (("SigModelUp", GREEN, "SigModel template (powheg) / nominal"), ("QCDScaleDown", BLUE, "QCDScale down / nominal"), ("PS_FSRDown", "#e07b7b", "PS FSR down / nominal")):
        v = f[f"mumu_SR__DYmumu__{k}"].values()
        ax.step(e, np.append(v / n, (v / n)[-1]), where="post", color=col, label=lab, linewidth=1.8)
    try:
        pf = uproot.open(ROOT / "fit" / "results" / "zmumu" / "Histograms" / "mumu_SR_postFit.root")["h_tot_postFit"].values()
        ax.errorbar(c, d / pf, yerr=np.sqrt(d) / pf, fmt="s", color=BLUE, markersize=3.5, label="data / post-fit prediction", alpha=0.9)
    except Exception:
        pass
    ax.axhline(1, ls="--", color=MUTED); ax.set_ylim(0.85, 1.08); ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("ratio"); ax.legend(fontsize=11, loc="lower right")
    label(ax, "signal-region shape: what the fit uses to absorb the low-mass tail")
    save(fig, "fit_sigmodel_shape")


def fit_review(res, extra):
    """Pulls/constraints and post-fit ranking from the fit-result JSON; binning-variant table."""
    nps = [k for k in res["nps"] if not k.startswith("gamma")]
    order = ["Lumi", "L1Prefiring", "MuonReco", "MuonID", "MuonIso", "MuonTrigger", "MuonScale", "MuonRes", "Pileup", "PDF", "AlphaS", "QCDScale", "PS_ISR", "PS_FSR", "SigModel",
             "XS_TTbar", "XS_SingleTop", "XS_WW", "XS_WZ", "XS_ZZ", "XS_DYtautau", "FakeStat_mumu", "FakeMethod_mumu"]
    nps = [k for k in order if k in res["nps"]]
    fig, ax = plt.subplots(figsize=(9, 8.5))
    y = np.arange(len(nps))
    ax.axvspan(-1, 1, color="#3a3a3a", zorder=0); ax.axvline(0, color=MUTED, lw=1)
    for i, k in enumerate(nps):
        p, up, dn = res["nps"][k]
        ax.errorbar(p, i, xerr=[[dn], [up]], fmt="o", color=GREEN if up < 0.5 else FG, markersize=5, capsize=3)
    ax.set_yticks(y); ax.set_yticklabels(nps, fontsize=12); ax.invert_yaxis(); ax.set_xlim(-2.2, 2.2)
    ax.set_xlabel(r"$(\hat\theta - \theta_0)/\Delta\theta$"); label(ax, "nuisance-parameter pulls and constraints (green: constrained below 0.5)", lumi=False)
    save(fig, "fit_pulls")
    # ranking
    rk = [r for r in res["ranking"] if not r["name"].startswith("gamma")][:14]
    fig, ax = plt.subplots(figsize=(9, 7.5))
    for i, r in enumerate(rk):
        ax.barh(i, 100 * r["dpoi_up_pre"], color="#555555", height=0.8); ax.barh(i, 100 * r["dpoi_down_pre"], color="#555555", height=0.8)
        ax.barh(i, 100 * r["dpoi_up_post"], color=BLUE, height=0.5); ax.barh(i, 100 * r["dpoi_down_post"], color=GREEN, height=0.5)
    ax.set_yticks(range(len(rk))); ax.set_yticklabels([r["name"] for r in rk], fontsize=12); ax.invert_yaxis()
    ax.set_xlabel(r"$\Delta\mu_Z$ [%]  (grey: pre-fit, colour: post-fit)"); ax.set_xlim(-3, 3)
    label(ax, "ranking: impact on $\\mu_Z$", lumi=False)
    save(fig, "fit_ranking")


def summary_fig(res):
    rows = [("aMC@NLO prediction", 799.6, 4.9, BLUE), ("madgraph LO prediction", 825.4, 8.3, "#6fa8c9"),
            ("v1: data-only counting", 773.2, 11.9, MUTED), ("v1 revised (review)", 776.9, 14.8, MUTED), ("reviewer's analysis", 797.2, 14.1, MUTED),
            ("v2 counting $(N-B)/(CL)$", res["meta"]["counting"]["sigma_fid_pb"], res["sigma_fid_tot_pb"], GREEN),
            ("v2: TRExFitter fit (30 bins)", res["sigma_fid_pb"], res["sigma_fid_tot_pb"], FG)]
    fig, ax = plt.subplots(figsize=(10, 6))
    for y, (lab, v, err, col) in enumerate(rows[::-1]):
        ax.errorbar(v, y, xerr=err, fmt="o", color=col, capsize=4, markersize=7)
        ax.text(v, y + 0.25, f"{v:.1f} $\\pm$ {err:.1f} pb", ha="center", fontsize=11, color=col)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=12); ax.set_ylim(-0.6, len(rows) - 0.2)
    ax.set_xlabel(r"$\sigma_{fid}(pp \to Z/\gamma^* \to \mu\mu)$ [pb]"); label(ax)
    save(fig, "summary_sigma")


def tnp_figs(tres):
    pt, eta = np.array(tres["pt_edges"]), np.array(tres["eta_edges"])
    c, hw = 0.5 * (pt[1:] + pt[:-1]), 0.5 * np.diff(pt)
    for eff, nice in (("id", "Tight ID"), ("iso", "Isolation (< 0.15)")):
        fig, axes = plt.subplots(1, 4, figsize=(24, 5.8), sharey=True)
        for j, ax in enumerate(axes):
            for key, lab, col, mk in ((f"{eff}_data_nominal", "data (fit)", FG, "o"), (f"{eff}_mc_nominal", "DY MC (fit)", BLUE, "s"), (f"{eff}_mc_truth", "DY MC gen-matched", GREEN, "^")):
                e = np.array(tres["eff"][key])[:, j]; er = np.array(tres["eff_err"][key])[:, j]
                ax.errorbar(c, e, yerr=er, xerr=hw, fmt=mk, color=col, markersize=4, label=lab)
            ax.set_xscale("log"); ax.set_xticks([20, 30, 50, 100, 200]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.set_xlabel(r"probe $p_T$ [GeV]"); ax.set_ylim(0.75, 1.03); ax.set_title(f"{eta[j]:.1f} < |$\\eta$| < {eta[j+1]:.1f}", fontsize=13, color=FG)
            if j == 0:
                ax.set_ylabel(f"{nice} efficiency"); ax.legend(fontsize=11, loc="lower right")
        save(fig, f"tnp_eff_{eff}")
        sf = np.array(tres["sf"][eff]); err = np.array(tres["sf_err"][eff])
        fig, ax = plt.subplots(figsize=(9.5, 6.5))
        mesh = ax.pcolormesh(pt, eta, sf.T, cmap="coolwarm", vmin=0.92, vmax=1.08)
        for i in range(len(pt) - 1):
            for j in range(len(eta) - 1):
                ax.text(np.sqrt(pt[i] * pt[i + 1]), 0.5 * (eta[j] + eta[j + 1]), f"{sf[i,j]:.3f}\n±{err[i,j]:.3f}", ha="center", va="center", fontsize=7, color="black")
        cb = fig.colorbar(mesh, ax=ax); cb.set_label(f"{nice} scale factor", color=FG); cb.ax.yaxis.set_tick_params(color=FG); plt.setp(cb.ax.get_yticklabels(), color=FG)
        ax.set_xscale("log"); ax.set_xticks([20, 30, 50, 100, 200]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xlabel(r"$p_T$ [GeV]"); ax.set_ylabel(r"|$\eta$|"); label(ax, f"{nice} data/MC scale factor")
        save(fig, f"tnp_sf_{eff}_map")
    tpt = np.array(tres["trig_pt_edges"]); tc, thw = 0.5 * (tpt[1:] + tpt[:-1]), 0.5 * np.diff(tpt)
    fig, axes = plt.subplots(1, 4, figsize=(24, 5.8), sharey=True)
    for j, ax in enumerate(axes):
        for s, col, lab in (("data", FG, "data"), ("mc", BLUE, "DY MC")):
            e = np.array(tres["eff"][f"trig_{s}"])[:, j]; er = np.array(tres["eff_err"][f"trig_{s}"])[:, j]
            ax.errorbar(tc, e, yerr=er, xerr=thw, fmt="o", color=col, markersize=4, label=lab)
        ax.set_xscale("log"); ax.set_xticks([20, 30, 50, 100, 200]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_ylim(0, 1.05); ax.set_xlabel(r"probe $p_T$ [GeV]"); ax.set_title(f"{eta[j]:.1f} < |$\\eta$| < {eta[j+1]:.1f}", fontsize=13, color=FG)
        if j == 0:
            ax.set_ylabel("IsoMu24 || IsoTkMu24 per-muon efficiency"); ax.legend(fontsize=11, loc="lower right")
    save(fig, "tnp_trigger")


def fakes_figs(fres, data, mc):
    pt, eta = np.array(fres["pt_edges"]), np.array(fres["eta_edges"])
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
    ff = np.array(fres["maps"]["nominal"]); ax = axes[0]
    mesh = ax.pcolormesh(pt, eta, ff.T, cmap="viridis", vmin=0, vmax=0.8)
    for i in range(len(pt) - 1):
        for j in range(len(eta) - 1):
            ax.text(np.sqrt(pt[i] * pt[i + 1]), 0.5 * (eta[j] + eta[j + 1]), f"{ff[i,j]:.2f}", ha="center", va="center", fontsize=9, color="white")
    cb = fig.colorbar(mesh, ax=ax); cb.set_label("fake factor N(tight)/N(anti-tight)", color=FG); plt.setp(cb.ax.get_yticklabels(), color=FG)
    ax.set_xscale("log"); ax.set_xticks([20, 30, 50, 100, 200]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel(r"muon $p_T$ [GeV]"); ax.set_ylabel(r"|$\eta$|"); label(ax, "fake factor, same-sign tag + probe region", lumi=False)
    ax = axes[1]; c = 0.5 * (FK.MASS_EDGES[1:] + FK.MASS_EDGES[:-1])
    obs = np.array(fres["closure"]["observed_hist"]); pred = np.array(fres["closure"]["predicted_hist"])
    k = 4; e2 = FK.MASS_EDGES[::k]; c2 = 0.5 * (e2[1:] + e2[:-1])
    ax.errorbar(c2, obs.reshape(-1, k).sum(1), yerr=np.sqrt(np.maximum(data["tt_ss_w2"].reshape(-1, k).sum(1), 0)), fmt="o", color=FG, markersize=4, label="SS tight-tight, data − prompt MC")
    ax.step(e2, np.append(pred.reshape(-1, k).sum(1), pred.reshape(-1, k).sum(1)[-1]), where="post", color=GREEN, label="FF × SS (tight + anti-tight)", linewidth=2)
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("events / 2 GeV"); ax.legend(fontsize=11); ax.set_ylim(-20, 160)
    label(ax, f"same-sign closure: predicted {fres['closure']['predicted_ss']:.0f} ± {fres['closure']['predicted_ss_err']:.0f}, observed {fres['closure']['observed_ss']:.0f} ± {fres['closure']['observed_ss_err']:.0f}", lumi=False)
    save(fig, "ff_map_and_closure")


def momentum_fig(cal, data, mc):
    fig, axes = plt.subplots(1, 4, figsize=(24, 5.8))
    eta = cal["eta_edges"]
    for j, ax in enumerate(axes):
        hd, hm = data[f"m_eta{j}"], mc[f"m_eta{j}"]; scale = hd.sum() / max(hm.sum(), 1e-9)
        fd, fm = MOM.fit_peak(hd, data[f"m_eta{j}_w2"]), MOM.fit_peak(hm, mc[f"m_eta{j}_w2"])   # refit (model not stored in JSON)
        ax.errorbar(MOM.CENTRES, hd, yerr=np.sqrt(np.maximum(hd, 0)), fmt="o", color=FG, markersize=3, label="data")
        ax.plot(MOM.CENTRES, np.array(fd["model"]), color=FG, label=f"fit: peak {fd['mu']:.2f}, σ {fd['sigma']:.2f}")
        ax.step(MOM.CENTRES, hm * scale, where="mid", color=BLUE, label="DY MC (scaled)")
        ax.plot(MOM.CENTRES, np.array(fm["model"]) * scale, color=BLUE, ls="--", label=f"fit: peak {fm['mu']:.2f}, σ {fm['sigma']:.2f}")
        ax.set_title(f"{eta[j]:.1f} < |η| < {eta[j+1]:.1f}: κ = {cal['kappa'][j]:.5f}, smear {100*cal['smear'][j]:.2f}%", fontsize=11, color=FG)
        ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.legend(fontsize=9)
    save(fig, "momentum_zpeak")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="pass_pkl", default=str(HERE / "review_pass.pkl"))
    args = ap.parse_args()
    hall = pickle.load(open(OUT_V2 / "histograms.pkl", "rb"))
    fres = json.load(open(OUT_V2 / "fakes.json"))
    fine = np.array(fres["templates"]["nominal"])
    res = json.load(open(ROOT / "fit" / "results" / "zmumu_fit_result.json"))
    tres = json.load(open(OUT_V2 / "tnp" / "tnp_result.json"))
    pu = json.load(open(OUT_V2 / "tnp" / "pileup_weights.json"))
    cal = json.load(open(OUT_V2 / "momentum.json"))
    print("[deck] data/MC")
    datamc(hall, "SR", "mass_fit", "sr_mass", fine, ratio_range=(0.9, 1.1))
    datamc(hall, "SR", "mass_fit", "sr_mass_log", fine, logy=True, ratio_range=(0.9, 1.1))
    for var in ("met", "npv", "zpt", "pt1", "njet"):
        datamc(hall, "SR", var, f"sr_{var}", fine)
    for var in ("mass_fit", "met", "pt_el", "njet", "pt1"):
        datamc(hall, "CRemu", var, f"cr_{var.replace('_fit', '')}", None)
    datamc(hall, "SS", "mass_fit", "ss_mass", None, ratio_range=(0, 2))
    print("[deck] T&P, fakes, momentum")
    tnp_figs(tres)
    from zmumu import batch
    parts = {}
    for key in ("data_2016G", "data_2016H", "DY_NLO", "TTbar", "ST_tW_top", "ST_tW_antitop", "WW", "WZ_3LNu", "WZ_2Q2L", "ZZ_4L", "ZZ_2L2Nu", "ZZ_2Q2L", "WJets"):
        d = OUT_V2 / "control_parts" / key
        if d.exists():
            parts[key] = batch.load_parts(d, {})
    def merged(keys, sub):
        tot = {}
        for k in keys:
            if k in parts:
                batch.merge(tot, parts[k][sub])
        return tot
    dat_ff = merged(("data_2016G", "data_2016H"), "ff"); mc_ff = merged([k for k in parts if not k.startswith("data")], "ff")
    fakes_figs(fres, dat_ff, mc_ff)
    dat_m = merged(("data_2016G", "data_2016H"), "mom"); mc_m = merged(("DY_NLO",), "mom")
    momentum_fig(cal, dat_m, mc_m)
    print("[deck] fit")
    fit_review(res, None); summary_fig(res); sigmodel_review(hall)
    if Path(args.pass_pkl).exists():
        print("[deck] review pass")
        R = pickle.load(open(args.pass_pkl, "rb"))
        comp = emu_review(R); pileup_review(R, pu); met_review(R, pu)
        json.dump({"emu_nonprompt_mc_os": comp}, open(HERE / "review_numbers_plots.json", "w"), indent=1)
    print("done")


if __name__ == "__main__":
    main()
