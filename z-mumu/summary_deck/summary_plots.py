#!/usr/bin/env python
"""Dark-style figures for the Z->mumu summary deck (summary_deck/zmumu_summary.tex).

    source ../../fitting/setup.sh && python summary_plots.py

Reuses the review's figure code (review/deck_plots.py) on the current v2 outputs and adds the
pileup-matching, lineshape, stability and post-fit-ratio figures. Writes figures/<name>.pdf/.png.
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "review"))

import numpy as np
import deck_plots as D
from deck_plots import BG, FG, MUTED, BLUE, GREEN, label, save, plt
from zmumu import batch, histograms as H, pileup

D.FIG = HERE / "figures"
D.STACK_ORDER = ["Fakes", "WJets", "WW", "WZ", "ZZ", "SingleTop", "TTbar", "DYee", "DYtautau", "DYmumu"]
D.COLOURS["WJets"] = "#8c6d4f"; D.LABELS["WJets"] = r"W+jets (jet$\to$e)"
D.XL["njet"] = r"lepton-cleaned jets ($p_T$ > 30 GeV)"; D.XL["puppimet"] = r"Puppi $p_T^{miss}$ [GeV]"; D.XL["met"] = r"PF $p_T^{miss}$ [GeV]"
OUT = ROOT / "output" / "v2"


def pileup_match(pu_json, npv_pkl):
    pu = pileup.PileupWeights.from_json(pu_json); d = pickle.load(open(npv_pkl, "rb"))
    h2, data = d["mc_ntrue_npv"], d["data_npv"]; e = np.arange(0.0, 61.0)
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8.5, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    dn = data / data.sum(); c = 0.5 * (e[1:] + e[:-1])
    ax.errorbar(c, dn, yerr=np.sqrt(data) / data.sum(), fmt="o", color=FG, markersize=3.5, label="data, signal region")
    for name, lab, col, w in (("none", "MC, no pileup weight", "#f2b134", None), ("csv", "MC, luminosity-CSV profile", MUTED, pileup._weights(pu.data["csv_raw"]["hist"], pu.mc_norm)),
                              ("matched", r"MC, $N_{PV}$-matched profile", BLUE, pu.weights["nominal"])):
        pred = h2.sum(0) if w is None else (w[:, None] * h2).sum(0); pred = pred / pred.sum()
        ax.step(e[:-1], pred, where="post", color=col, label=lab, linewidth=1.8)
        with np.errstate(divide="ignore", invalid="ignore"):
            rax.step(e[:-1], np.where(pred > 0, dn / pred, np.nan), where="post", color=col, linewidth=1.8)
    ax.set_ylabel("fraction of events"); ax.legend(fontsize=12); ax.set_xlim(0, 45); ax.set_ylim(0, dn.max() * 1.45)
    m = pu.data["nominal"]; label(ax, f"pileup profile: scale {pu.scale:.3f}, bunch spread {pu.rel_smear:.2f}; mean {pu.data['csv_raw']['mean']:.1f} -> {m['mean']:.1f}")
    rax.axhline(1, ls="--", color=MUTED); rax.set_ylim(0.75, 1.25); rax.set_ylabel("data / MC"); rax.set_xlabel("good primary vertices")
    save(fig, "pileup_match")


def lineshape(ls):
    e = np.linspace(60, 120, 13); c = 0.5 * (e[1:] + e[:-1])
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.step(e, np.append(ls["lhe_ratio_5gev"], ls["lhe_ratio_5gev"][-1]), where="post", color="#e07b7b", linewidth=2, label="powheg / aMC@NLO, generator level (LHE mass)")
    ax.step(e, np.append(ls["reco_ratio_5gev"], ls["reco_ratio_5gev"][-1]), where="post", color=GREEN, linewidth=2, ls="--", label=r"powheg / aMC@NLO, reconstructed, $50<m_{LHE}<120$")
    r = np.array(ls["data_over_pred_prefit_5gev"])
    ax.errorbar(c, r, yerr=0.002, fmt="o", color=FG, markersize=4, label=r"data / prediction, pre-fit ($\mu_Z$ = 1)")
    ax.axhline(1, ls=":", color=MUTED); ax.set_xlim(60, 120); ax.set_ylim(0.9, 1.05); ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("ratio (shapes normalised 60-120 GeV)")
    ax.legend(fontsize=11, loc="lower right"); label(ax, "signal lineshape: two NLO generators and the data")
    save(fig, "lineshape")


def stability(stab):
    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    n = len(stab)
    for i, x in enumerate(stab):
        ok = (x["gof_p"] or 0) > 0.05 or "1 bin" in x["label"]
        col = GREEN if i == 0 else (FG if ok else "#CC4C4C")
        ax.errorbar(x["mu"], n - 1 - i, xerr=[[x["err_down"]], [x["err_up"]]], fmt="o", color=col, capsize=4, markersize=6)
        ax.text(1.028, n - 1 - i, "GoF p = " + ("--" if "1 bin" in x["label"] else f"{x['gof_p']:.2g}"), va="center", fontsize=11, color=col)
    ax.axvspan(stab[0]["mu"] - stab[0]["err_down"], stab[0]["mu"] + stab[0]["err_up"], color="#3a3a3a", zorder=0)
    ax.set_yticks(range(n)); ax.set_yticklabels([x["label"].replace("nominal: ", "").replace(", SigModel two-sided (LHE window), no smoothing", " (nominal)") for x in stab[::-1]], fontsize=13)
    ax.set_xlim(0.96, 1.045); ax.set_xlabel(r"$\mu_Z$"); label(ax, "stability against the fit configuration (red: fit does not describe the data)", lumi=False)
    save(fig, "stability")


def postfit_ratio():
    import uproot
    f = uproot.open(ROOT / "fit" / "fitinputs" / "zmumu.root")
    e = np.linspace(60, 120, 13); c = 0.5 * (e[1:] + e[:-1])
    d = f["mumu_SR__Data"].values()
    tot = sum(f[k].values() for k in f.keys() if k.startswith("mumu_SR__") and k.count("__") == 1 and "Data" not in k)
    pf = uproot.open(ROOT / "fit" / "results" / "zmumu" / "Histograms" / "mumu_SR_postFit.root")["h_tot_postFit"].values()
    n = f["mumu_SR__DYmumu"].values(); su = f["mumu_SR__DYmumu__SigModelUp"].values()
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.step(e, np.append(su / n, (su / n)[-1]), where="post", color=GREEN, linewidth=1.8, label=r"SigModel $+1\sigma$ template / nominal")
    ax.errorbar(c, d / tot, yerr=np.sqrt(d) / tot, fmt="o", color=FG, markersize=4, label=r"data / pre-fit prediction ($\mu_Z$ = 1)")
    ax.errorbar(c, d / pf, yerr=np.sqrt(d) / pf, fmt="s", color=BLUE, markersize=4, label="data / post-fit prediction")
    ax.axhline(1, ls="--", color=MUTED); ax.set_ylim(0.9, 1.05); ax.set_xlim(60, 120); ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("ratio")
    ax.legend(fontsize=11, loc="lower right"); label(ax, "signal region, 12 bins: before and after the fit")
    save(fig, "postfit_ratio")
    return float(tot.sum()), float(pf.sum())


def summary(res):
    rows = [("aMC@NLO prediction", 799.6, 4.9, BLUE), ("madgraph LO prediction", 825.4, 8.3, "#6fa8c9"),
            ("v1: data-only counting", 773.2, 11.9, MUTED), ("external reviewer, same volume", 797.2, 14.1, MUTED),
            (r"v2 counting $(N-B)/(CL)$", res["meta"]["counting"]["sigma_fid_pb"], res["sigma_fid_tot_pb"], FG),
            (r"v2: TRExFitter fit (12 $\times$ 5 GeV)", res["sigma_fid_pb"], res["sigma_fid_tot_pb"], GREEN)]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for y, (lab, v, err, col) in enumerate(rows[::-1]):
        ax.errorbar(v, y, xerr=err, fmt="o", color=col, capsize=4, markersize=7)
        ax.text(v, y + 0.25, f"{v:.1f} $\\pm$ {err:.1f} pb", ha="center", fontsize=11, color=col)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=12); ax.set_ylim(-0.6, len(rows) - 0.2)
    ax.set_xlabel(r"$\sigma_{fid}(pp \to Z/\gamma^* \to \mu\mu)$ [pb]"); label(ax)
    save(fig, "summary_sigma")


def main():
    hall = pickle.load(open(OUT / "histograms.pkl", "rb"))
    fres = json.load(open(OUT / "fakes.json")); fine = np.array(fres["templates"]["nominal"])
    res = json.load(open(ROOT / "fit" / "results" / "zmumu_fit_result.json"))
    tres = json.load(open(OUT / "tnp" / "tnp_result.json")); cal = json.load(open(OUT / "momentum.json"))
    r2 = json.load(open(OUT / "results_v2.json"))
    D.datamc(hall, "SR", "mass_fit", "sr_mass", fine, ratio_range=(0.9, 1.1))
    D.datamc(hall, "SR", "mass_fit", "sr_mass_log", fine, logy=True, ratio_range=(0.9, 1.1))
    for var in ("met", "puppimet", "npv", "zpt", "njet", "pt1"):
        D.datamc(hall, "SR", var, f"sr_{var}", fine)
    for var in ("mass_fit", "pt_el", "njet", "met"):
        D.datamc(hall, "CRemu", var, f"cr_{var.replace('_fit', '')}", None, title=r"e$\mu$ region (OS), all MC")
    D.datamc(hall, "SSemu", "mass_fit", "ssemu_mass", None, title=r"e$\mu$ region (SS), all MC", ratio_range=(0.5, 1.5))
    D.tnp_figs(tres)
    parts = {k: batch.load_parts(OUT / "control_parts" / k, {}) for k in ("data_2016G", "data_2016H", "DY_NLO", "TTbar", "ST_tW_top", "ST_tW_antitop", "WW", "WZ_3LNu", "WZ_2Q2L", "ZZ_4L", "ZZ_2L2Nu", "ZZ_2Q2L", "WJets") if (OUT / "control_parts" / k).exists()}
    def merged(keys, sub):
        tot = {}
        for k in keys:
            if k in parts: batch.merge(tot, parts[k][sub])
        return tot
    D.fakes_figs(fres, merged(("data_2016G", "data_2016H"), "ff"), merged([k for k in parts if not k.startswith("data")], "ff"))
    D.momentum_fig(cal, merged(("data_2016G", "data_2016H"), "mom"), merged(("DY_NLO",), "mom"))
    D.fit_review(res, None); summary(res)
    pileup_match(OUT / "tnp" / "pileup_weights.json", OUT / "tnp" / "pileup_npv.pkl")
    lineshape(r2["lineshape"]); stability(r2["stability"])
    pre, post = postfit_ratio()
    y = r2["yields"]
    print("SR prefit total", pre, "postfit", post, "data", y["SR"]["Data"])
    for reg in ("CRemu", "SSemu"):
        mc = sum(v for k, v in y[reg].items() if k not in ("Data", "DYmumu_powheg", "DYother")); print(reg, y[reg]["Data"], round(mc), round(y[reg]["Data"] / mc, 3), {k: round(v) for k, v in y[reg].items() if k not in ("DYmumu_powheg",)})
    print("fakes", fres["yields"]["sr_fakes"], fres["yields"]["sr_fakes_stat"], fres["closure"]["predicted_ss"], fres["closure"]["observed_ss"])
    print("tnp", tres["meta"], "id/iso SF", np.mean(tres["sf"]["id"]), np.mean(tres["sf"]["iso"]))
    print("mom", cal["kappa"], cal["smear"])


if __name__ == "__main__":
    main()
