#!/usr/bin/env python
"""Vector figures for the slide deck, all on the deck background (prompts/presentation_style.md).

    source ../setup.sh && python slides/make_figures.py        # -> slides/figs/*.pdf (+ .png fallbacks)

Every figure is drawn on DARK_BG = #222222 so it blends with the slides. Inputs: output/results.json,
output/data/{fakefactors,bdt,yields}.json, fit/fitinputs/ztautau.root, fit/results/ztautau/ (post-fit total
histograms and yield tables) and, for the two control plots of the signal-dominated category, the ntuples.
The deck itself is assembled by slides/build_deck.py (PyMuPDF, run in the betterplottingtool venv).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402
import yaml  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from ztautau import analysis, bdt, config, fakes, plotting, samples  # noqa: E402
from ztautau.plotting import COLORS, DARK_BG, DARK_FG, LABELS, dark_rc, save  # noqa: E402

FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)
NOM = "mcsub"
BLUE, GREEN, ORANGE, PINK, MUTED = "#2BA4DD", "#8AC63F", "#D9822B", "#E377AC", "#99999E"
R = json.loads((HERE.parent / "output/results.json").read_text())
FIT = R["fit"][NOM]
YIELDS = json.loads((config.DATA_DIR / "yields.json").read_text())
BDT = json.loads((config.DATA_DIR / "bdt.json").read_text())
FF = json.loads((config.DATA_DIR / "fakefactors.json").read_text())
EDGES = np.asarray(config.FIT_BINS)
STACK_GROUPS = [("DYll", ["DYee", "DYmumu"]), ("DYlowmass", ["DYlowmass"]), ("Diboson", ["WW", "WZ", "ZZ"]),
                ("Top", ["TTbar", "SingleTop"]), ("WJets", ["WJets"]), ("Fakes", ["Fakes"]),
                ("DYtautau_nonfid", ["DYtautau_nonfid"]), ("DYtautau", ["DYtautau"])]
TITLES = {"Z#rightarrow#tau#tau (fiducial)": "DYtautau", "Z/#gamma*#rightarrow#tau#tau (non-fid.)": "DYtautau_nonfid",
          "Z#rightarrowee": "DYee", "Z#rightarrow#mu#mu": "DYmumu", "Z/#gamma*#rightarrowll (m<50)": "DYlowmass",
          "W+jets": "WJets", "t#bar{t}": "TTbar", "single t": "SingleTop", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ",
          "jet#rightarrow#tau_{h} (FF)": "Fakes"}
plt.rcParams.update({"font.size": 13})


def dark_fig(figsize=(7, 4.4)):
    plt.rcParams.update(dark_rc())
    fig, ax = plt.subplots(figsize=figsize)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return fig, ax


def stack_from(hists, region):
    stack = []
    for g, members in STACK_GROUPS:
        parts = [hists[f"{region}__{m}"] for m in members if f"{region}__{m}" in hists]
        if parts:
            v = sum(p[0] for p in parts); var = sum(p[1] for p in parts)
            stack.append((g, np.maximum(v, 0), var))
    return stack


# ------------------------------------------------------------------ 1. method flow diagram
def flow():
    fig, ax = dark_fig((6.2, 7.2))
    ax.set_axis_off()
    boxes = [("NanoAOD:  160 GB data, 4 Drell-Yan + background samples", MUTED), ("skims (trigger, $\\geq$2 $\\tau_h$ candidates):  1.7 GB", MUTED),
             ("ntuples: pair, trigger match, vetoes, masses:  330 MB (laptop)", MUTED),
             ("fake factors (era $\\times$ DM $\\times$ $N_{jets}$ $\\times$ $p_T$), MC subtracted,\nclosure corrections $f(|\\eta_1|)\\,g(p_{T,2})$, $C_{OS/SS}$", PINK),
             ("k-fold BDT (5 folds by event number)  $\\rightarrow$  SR0 / SR1 / SR2", BLUE),
             ("templates: fiducial signal, non-fiducial DY as background,\n39 nuisance parameters + 42 $\\gamma$", MUTED),
             ("TRExFitter v1.8.0: fit of $m_{\\tau\\tau}$ in 3 categories  $\\rightarrow$  $\\mu_Z$, $\\sigma$", ORANGE)]
    n = len(boxes); h = 1.0 / n
    for i, (txt, col) in enumerate(boxes):
        y = 1 - (i + 1) * h + 0.012
        ax.add_patch(plt.Rectangle((0.02, y), 0.96, h - 0.024, transform=ax.transAxes, facecolor=col, alpha=0.22, edgecolor=col, lw=2))
        ax.text(0.5, y + (h - 0.024) / 2, txt, transform=ax.transAxes, ha="center", va="center", fontsize=11.5, color=DARK_FG)
        if i < n - 1:
            ax.annotate("", xy=(0.5, y - 0.003), xytext=(0.5, y + 0.012 - 0.0), xycoords="axes fraction",
                        arrowprops=dict(arrowstyle="-|>", color=DARK_FG, lw=1.5))
    save(fig, FIGS / "flow", dark=True); plt.close(fig)


# ------------------------------------------------------------------ 2. result summary and impacts
def result_summary():
    fig, ax = dark_fig((9, 5.4))
    pred = R["prediction"]["sigma_tautau_60_120_pb"]
    s = FIT["sigma_60_120_pb"]
    H = json.loads((HERE / "history.json").read_text())
    rows = [("v3 nominal: DeepTau Tight (this result)", s["value"], s["stat"], s["err_down"], s["err_up"], BLUE)]
    for key, col in (("v2.1", ORANGE), ("v2", MUTED), ("v1", MUTED)):
        h = H[key]; rows.append((h["label"], h["sigma60"], h["sigma_stat"], h["sigma_down"], h["sigma_up"], col))
    rows.append((r"Z$\rightarrow\mu\mu$ v2 (60$-$120)", 1935, 5, 30, 30, GREEN))
    ax.axvspan(pred * 0.96, pred * 1.04, color=ORANGE, alpha=0.25); ax.axvline(pred, color=ORANGE)
    for i, (lab, v, st, dn, up, col) in enumerate(rows):
        y = len(rows) - i
        ax.errorbar([v], [y], xerr=[[dn], [up]], fmt="o", color=col, capsize=4, lw=1.6, ms=7)
        ax.errorbar([v], [y], xerr=[[st], [st]], fmt="none", color="#CC4C4C", lw=5)
        ax.text(1330, y + 0.28, lab, fontsize=11, color=col)
        ax.text(3230, y + 0.28, f"{v:.0f} +{up:.0f} $-${dn:.0f} pb", fontsize=11, color=col, ha="right")
    ax.text(pred + 12, 0.45, f"NNLO {pred:.0f} pb", fontsize=10, color=ORANGE)
    ax.set_ylim(0.3, len(rows) + 0.9); ax.set_xlim(1300, 3250); ax.set_yticks([])
    ax.set_xlabel(r"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau)$, $60<m<120$ GeV [pb]   (red: statistical)")
    save(fig, FIGS / "result_summary", dark=True); plt.close(fig)

    gi = {k: v for k, v in FIT["grouped_impact"].items() if k != "FullSyst"}
    gi["Data statistics"] = FIT["mu_stat"]
    items = sorted(gi.items(), key=lambda kv: kv[1])
    fig, ax = dark_fig((7, 5))
    ax.barh([k for k, _ in items], [100 * v for _, v in items], color=[GREEN if k == "Data statistics" else BLUE for k, _ in items])
    for i, (_, val) in enumerate(items):
        ax.text(100 * val + 0.15, i, f"{100 * val:.1f}%", va="center", fontsize=10, color=DARK_FG)
    ax.set_xlabel(r"impact on $\mu_Z$ [%]"); ax.set_xlim(0, 100 * max(gi.values()) * 1.25)
    save(fig, FIGS / "impacts", dark=True); plt.close(fig)


def ranking_and_pulls():
    rk = FIT["ranking"]
    titles = {"TauID_DM": r"$\tau_h$ ID DM", "TauTrigger_DM": r"$\tau_h$ trigger DM", "TauES_DM": r"$\tau_h$ energy scale DM",
              "FakeOSSS_tautau": "FF OS/SS extrapolation", "FakeClosure_tautau_c": "FF non-closure SR", "XS_DYtautau_nonfid": r"$\sigma$(non-fid. DY)",
              "MCStatNorm_WJets_tautau": "W+jets MC stat.", "MET_Unclustered": "MET unclustered", "PS_FSR": "PS FSR", "PS_ISR": "PS ISR",
              "QCDScale": "QCD scales", "Pileup": "pileup", "Lumi": "luminosity", "PDF": "PDF", "L1Prefiring": "L1 prefiring",
              "TauFakeEle": r"e$\rightarrow\tau_h$", "TauFakeMu": r"$\mu\rightarrow\tau_h$", "XS_": r"$\sigma$("}

    def nice(n):
        if n.startswith("gamma_stat_tautau_"):
            return "MC stat. " + n[len("gamma_stat_tautau_"):].replace("_bin_", " bin ")
        for k, v in titles.items():
            if n.startswith(k):
                rest = n[len(k):]
                if k == "FakeClosure_tautau_c":
                    return v + rest[0] + (" m<110" if rest.endswith("lo") else " m>110")
                if k == "FakeOSSS_tautau" and rest.startswith("_c"):
                    return v + " SR" + rest[2:]
                if k == "XS_":
                    return v + rest + ")"
                return v + rest
        return n
    fig, ax = dark_fig((8, 6.5))
    n = len(rk); ys = np.arange(n)[::-1]
    for y, r in zip(ys, rk):
        ax.barh(y, r["impact_up"], color=BLUE, alpha=0.9, height=0.7)
        ax.barh(y, r["impact_down"], color=GREEN, alpha=0.9, height=0.7)
    ax2 = ax.twiny()
    ax2.errorbar([r["pull"] for r in rk], ys, xerr=[r["constraint"] for r in rk], fmt="o", color=DARK_FG, ms=5, capsize=3)
    ax2.set_xlim(-2.2, 2.2); ax2.axvline(-1, color=MUTED, ls="--", lw=0.8); ax2.axvline(1, color=MUTED, ls="--", lw=0.8)
    ax2.set_xlabel(r"$(\hat\theta-\theta_0)/\Delta\theta$ (markers)", color=DARK_FG)
    lim = 1.15 * max(abs(r["impact_up"]) for r in rk); ax.set_xlim(-lim, lim)
    ax.set_yticks(ys, [nice(r["name"]) for r in rk], fontsize=10.5)
    ax.set_xlabel(r"post-fit impact on $\mu_Z$ (bars: $+1\sigma$ blue, $-1\sigma$ green)")
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    save(fig, FIGS / "ranking", dark=True); plt.close(fig)

    pulls = [(k, v) for k, v in FIT["pulls"].items()]
    fig, ax = dark_fig((8, 9))
    ys = np.arange(len(pulls))[::-1]
    ax.axvspan(-1, 1, color=GREEN, alpha=0.15); ax.axvspan(-2, 2, color=ORANGE, alpha=0.08)
    ax.errorbar([v[0] for _, v in pulls], ys, xerr=[v[1] for _, v in pulls], fmt="o", color=DARK_FG, ms=4, capsize=2)
    ax.set_yticks(ys, [nice(k) for k, _ in pulls], fontsize=8.5); ax.set_xlim(-2.5, 2.5)
    ax.axvline(0, color=MUTED, lw=0.8); ax.set_xlabel(r"$(\hat\theta-\theta_0)/\Delta\theta$")
    save(fig, FIGS / "pulls", dark=True); plt.close(fig)


# ------------------------------------------------------------------ 3. pre- and post-fit m_tt per category
def read_fitinputs():
    out = {}
    with uproot.open(config.FIT_DIR / "fitinputs" / "ztautau.root") as f:
        for k in f.keys():
            name = k.split(";")[0]
            if name.count("__") == 1 and name != "meta_json":
                h = f[name]; out[name] = (h.values(), h.variances())
    return out


def category_plots():
    hists = read_fitinputs()
    post = yaml.safe_load(open(config.FIT_DIR / "results/ztautau/Tables/Table_postfit.yaml"))
    pre = yaml.safe_load(open(config.FIT_DIR / "results/ztautau/Tables/Table_prefit.yaml"))
    for k, region in enumerate(config.REGIONS):
        data = hists[f"{region}__Data"][0]
        stack = stack_from(hists, region)
        lab = config.REGION_LABELS[k]
        for logy in (False, True):
            plotting.stack_plot(FIGS / f"prefit_{region}{'_log' if logy else ''}.png", EDGES, (data, np.sqrt(data)), stack,
                                r"$m_{\tau\tau}$ [GeV]", title=f"prefit, {lab}", dark=True, logy=logy, density=not logy)
        # post-fit: prefit templates scaled to the post-fit yields per sample, exact post-fit total from TRExFitter
        sf = {}
        for entry_pre, entry_post in zip(pre[k]["Samples"], post[k]["Samples"]):
            if "Sample" in entry_pre and entry_pre["Sample"] in TITLES:
                s = TITLES[entry_pre["Sample"]]
                sf[s] = entry_post["Yield"] / entry_pre["Yield"] if entry_pre["Yield"] > 0 else 1.0
        scaled = {}
        for name, (v, var) in hists.items():
            rg, s = name.split("__")
            if rg == region and s in sf:
                scaled[name] = (v * sf[s], var * sf[s] ** 2)
        with uproot.open(config.FIT_DIR / f"results/ztautau/Histograms/{region}_postFit.root") as f:
            tot = f["h_tot_postFit"]; tv, te = tot.values(), tot.errors()
        plotting.stack_plot(FIGS / f"postfit_{region}.png", EDGES, (data, np.sqrt(data)), stack_from(scaled, region),
                            r"$m_{\tau\tau}$ [GeV]", title=f"post-fit, {lab}", dark=True, overlay=(tv, te, "post-fit total"), density=True)
    # all categories together, log: signal composition
    tot = {}
    for name, (v, var) in hists.items():
        rg, s = name.split("__")
        tot[f"all__{s}"] = (tot.get(f"all__{s}", (0, 0))[0] + v, tot.get(f"all__{s}", (0, 0))[1] + var)
    data = tot["all__Data"][0]
    plotting.stack_plot(FIGS / "prefit_all_log.png", EDGES, (data, np.sqrt(data)), stack_from(tot, "all"), r"$m_{\tau\tau}$ [GeV]",
                        title="signal region, all categories, prefit", dark=True, logy=True)


# ------------------------------------------------------------------ 4. BDT figures (from bdt.json)
def bdt_figures():
    t = BDT["training"]; e = np.asarray(BDT["sr_score"]["edges"]); c = 0.5 * (e[1:] + e[:-1])
    # importance
    items = list(t["importance"].items())[::-1]
    fig, ax = dark_fig((6.4, 4.6))
    ax.barh([k for k, _ in items], [v for _, v in items], color=BLUE)
    ax.set_xlabel("feature importance (gain, mean over folds)")
    save(fig, FIGS / "bdt_importance", dark=True); plt.close(fig)
    # SR score, log
    sr = BDT["sr_score"]; d = np.asarray(sr["data"], float)
    mc = {k: (np.asarray(v), np.zeros(len(v))) for k, v in sr["mc"].items()}
    mc["Fakes"] = (np.maximum(np.asarray(sr["fakes"]), 0), np.zeros(len(d)))
    stack = []
    for g, members in STACK_GROUPS:
        parts = [mc[m] for m in members if m in mc]
        if parts:
            stack.append((g, np.maximum(sum(p[0] for p in parts), 0), np.zeros(len(d))))
    plotting.stack_plot(FIGS / "bdt_sr_score.png", e, (d, np.sqrt(d)), stack, "BDT score", title="signal region, prefit", dark=True, logy=True)
    # SS closure in the score
    ss = BDT["ss_closure_score"]; o = np.asarray(ss["obs"], float); om = np.asarray(ss["obs_mc"]); p = np.asarray(ss["pred"]); pv = np.asarray(ss["pred_var"])
    plotting.stack_plot(FIGS / "bdt_closure_ss.png", e, (o, np.sqrt(o)), [("MC", om, np.zeros(len(o))), ("Fakes", np.maximum(p, 0), pv)],
                        "BDT score", title="same-sign closure of the FF in the score", dark=True, logy=True)
    # score shapes are not stored per event: draw the normalised SR components instead
    fig, ax = dark_fig((6.4, 4.2))
    for name, col, lab in (("DYtautau", ORANGE, r"Z$\rightarrow\tau\tau$, fiducial"), ("DYtautau_nonfid", BLUE, r"Z/$\gamma^*\rightarrow\tau\tau$, non-fiducial")):
        v = np.asarray(sr["mc"][name], float); ax.stairs(v / v.sum(), e, color=col, lw=2, label=lab)
    v = np.maximum(np.asarray(sr["fakes"], float), 0); ax.stairs(v / v.sum(), e, color=PINK, lw=2, label="fakes (AR x FF)")
    for x in config.BDT_CATEGORY_EDGES[1:-1]:
        ax.axvline(x, color=MUTED, ls="--", lw=0.8)
    ax.set_xlabel("BDT score (held-out fold)"); ax.set_ylabel("normalised"); ax.legend(fontsize=10, frameon=False)
    ax.set_title(f"held-out AUC {np.mean(t['auc_test']):.3f} (train {np.mean(t['auc_train']):.3f})", fontsize=11, color=DARK_FG)
    save(fig, FIGS / "bdt_shapes", dark=True); plt.close(fig)


# ------------------------------------------------------------------ 5. fake-factor figures (from fakefactors.json)
def ff_figures():
    tbl = FF[NOM]["ff"]; ff, err = np.asarray(tbl["ff"]), np.asarray(tbl["err"])
    edges = np.asarray(config.FF_PT_BINS[:-1] + [110.0]); centers = 0.5 * (edges[1:] + edges[:-1])
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    plt.rcParams.update(dark_rc())
    for e_, (ax, era) in enumerate(zip(axes, tbl["eras"])):
        ax.set_facecolor(DARK_BG)
        i = tbl["dms"].index(1)
        for j, (lab, col) in enumerate(zip(["0 jets", "1 jet", r"$\geq$2 jets"], [DARK_FG, ORANGE, BLUE])):
            ax.errorbar(centers + 1.2 * (j - 1), ff[e_, i, j], yerr=err[e_, i, j], xerr=np.diff(edges) / 2, fmt="o", ms=4, color=col, label=lab, elinewidth=1)
        ax.set_title(f"Run2016{era}, decay mode 1", fontsize=11, color=DARK_FG); ax.set_xlabel(r"$p_T(\tau_1)$ [GeV]  (last bin > 80)")
        ax.tick_params(colors=DARK_FG); [sp.set_color(DARK_FG) for sp in ax.spines.values()]
    axes[0].set_ylabel("fake factor"); axes[0].legend(fontsize=9, frameon=False)
    fig.patch.set_facecolor(DARK_BG); save(fig, FIGS / "ff_dm1", dark=True); plt.close(fig)
    # closure corrections
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, (var, xl) in zip(axes, (("eta", r"$|\eta(\tau_1)|$"), ("pt2", r"$p_T(\tau_2)$ [GeV] (last bin: > 80)"))):
        ax.set_facecolor(DARK_BG); ax.tick_params(colors=DARK_FG); [sp.set_color(DARK_FG) for sp in ax.spines.values()]
        for variant, col, lab in ((NOM, BLUE, "MC subtracted (nominal)"), ("nosub", MUTED, "no subtraction")):
            c = FF[variant]["ff"]["closure"][var]; ed = np.asarray(c["edges"]); ed[-1] = min(ed[-1], 120.0)
            ax.errorbar(0.5 * (ed[1:] + ed[:-1]), c["values"], yerr=c["err"], xerr=np.diff(ed) / 2, fmt="o", color=col, label=lab)
        ax.axhline(1, color=MUTED, lw=0.8); ax.set_xlabel(xl); ax.set_ylabel("same-sign obs / FF prediction"); ax.set_ylim(0.8, 1.2); ax.legend(fontsize=9, frameon=False)
    fig.patch.set_facecolor(DARK_BG); fig.tight_layout(); save(fig, FIGS / "closure_corrections", dark=True); plt.close(fig)
    # closure before / after in eta(tau1), pT(tau2), N_jets, m_tt
    for var, xl in (("t1_eta", r"$\eta(\tau_1)$"), ("t2_pt", r"$p_T(\tau_2)$ [GeV]"), ("njets", r"$N_{jets}$"), ("m_tt", r"$m_{\tau\tau}$ [GeV]")):
        c = FF[NOM]["closure"][var]; ed = np.asarray(c["edges"]); obs = np.asarray(c["obs"], float)
        for tag, pred in (("before", c["pred_nocorr"]), ("after", c["pred"])):
            stack = [("MC", np.asarray(c["obs_mc"]), np.zeros(len(obs))), ("Fakes", np.maximum(np.asarray(pred), 0), np.asarray(c["pred_var"]))]
            plotting.stack_plot(FIGS / f"closure_{var}_{tag}.png", ed, (obs, np.sqrt(obs)), stack, xl,
                                title=f"same-sign closure, {'before' if tag == 'before' else 'after'} the closure corrections", dark=True,
                                logy=(var == "m_tt"))


# ------------------------------------------------------------------ 6. from the ntuples: mass estimators, SR2 control plots
def ntuple_figures():
    d, _ = analysis.load(samples.DY_INCLUSIVE)
    r = analysis.regions(d, is_mc=True); w = analysis.weights(d, samples.DY_INCLUSIVE)
    sel = r["SR"] & (d["gen_lhe_flavour"] == 15) & (d["t1_genflav"] == 5) & (d["t2_genflav"] == 5) & (d["gen_mll_lhe"] > 70) & (d["gen_mll_lhe"] < 110)
    ref = d["gen_mll_lhe"][sel]
    fig, ax = dark_fig((6.4, 4.2))
    bins = np.linspace(0.3, 2.0, 69)
    for var, lab, col in (("m_vis", r"$m_\mathrm{vis}$", BLUE), ("m_col", r"$m_\mathrm{col}$ (38% defined)", MUTED), ("m_tt", r"$m_{\tau\tau}$ (MET likelihood)", ORANGE)):
        x = d[var][sel] / ref; ok = d[var][sel] > 0
        med = np.median(x[ok]); q25, q75 = np.percentile(x[ok], [25, 75])
        ax.hist(x[ok], bins=bins, weights=w[sel][ok], histtype="step", lw=2, color=col, density=True, label=f"{lab}: median {med:.2f}, IQR/2 {100 * (q75 - q25) / 2 / med:.0f}%")
    ax.axvline(1.0, color=DARK_FG, lw=0.8, ls="--"); ax.set_xlabel(r"reconstructed mass / generator $m_{\tau\tau}$"); ax.set_ylabel("normalised")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.45); ax.legend(fontsize=9, frameon=False)
    save(fig, FIGS / "mass_estimators", dark=True); plt.close(fig)
    # SR2 control plots: pT(tau1), pT(tau2), dR
    table = fakes.from_json(FF[NOM]["ff"]); C = np.asarray(FF[NOM]["osss"]["C"])
    data = analysis.load_data(); reg = analysis.regions(data); cat = analysis.categories(data, "data")
    wf = fakes.fake_weights(data, reg["AR"], table, C)
    k = len(config.REGIONS) - 1
    sr, ar = reg["SR"] & (cat == k), reg["AR"] & (cat == k)
    vars_ = {"t1_pt": (np.array([40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 120, 150, 200, 300]), r"$p_T(\tau_1)$ [GeV]"),
             "t2_pt": (np.array([40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 120, 150]), r"$p_T(\tau_2)$ [GeV]"),
             "dr_tt": (np.linspace(0.5, 4.0, 15), r"$\Delta R(\tau_1,\tau_2)$")}
    books = {v: {} for v in vars_}

    def add(book, name, x, w_, ed):
        xc = np.clip(x, ed[0], ed[-1] - 1e-6)
        v = np.histogram(xc, bins=ed, weights=w_)[0]; v2 = np.histogram(xc, bins=ed, weights=w_ ** 2)[0]
        book[name] = (book[name][0] + v, book[name][1] + v2) if name in book else (v, v2)
    for v, (ed, _) in vars_.items():
        add(books[v], "Data", data[v][sr], np.ones(sr.sum()), ed); add(books[v], "Fakes", data[v][ar], wf[ar], ed)
    for key in analysis.available_mc():
        dm, _ = analysis.load(key)
        if not len(dm):
            continue
        rm = analysis.regions(dm, is_mc=True); wm = analysis.weights(dm, key); cm_ = analysis.categories(dm, key)
        wsub = analysis.subtraction_weights(key, wm)
        for name, cm in analysis.mc_components(key):
            m = rm["SR"] & cm & (cm_ == k); ma = rm["AR"] & cm & (cm_ == k)
            wfa = fakes.fake_weights(dm, ma, table, C) * wsub
            for v, (ed, _) in vars_.items():
                add(books[v], name, dm[v][m], wm[m], ed)
                xc = np.clip(dm[v][ma], ed[0], ed[-1] - 1e-6)
                sub = np.histogram(xc, bins=ed, weights=wfa[ma])[0]
                books[v]["Fakes"] = (books[v]["Fakes"][0] - sub, books[v]["Fakes"][1])
    for v, (ed, xl) in vars_.items():
        b = books[v]; dv = b["Data"][0]
        stack = []
        for g, members in STACK_GROUPS:
            parts = [b[m] for m in members if m in b]
            if parts:
                stack.append((g, np.maximum(sum(p[0] for p in parts), 0), sum(p[1] for p in parts)))
        plotting.stack_plot(FIGS / f"SR2_{v}.png", ed, (dv, np.sqrt(dv)), stack, xl, title=f"prefit, {config.REGION_LABELS[k]}", dark=True,
                            density=(v != "dr_tt"))


if __name__ == "__main__":
    only = sys.argv[1:]
    for fn in (flow, result_summary, ranking_and_pulls, category_plots, bdt_figures, ff_figures, ntuple_figures):
        if not only or fn.__name__ in only:
            fn()
    print(f"figures -> {FIGS}: {len(list(FIGS.glob('*.pdf')))} PDFs")
