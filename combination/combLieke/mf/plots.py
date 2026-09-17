"""Figures of the combination, from output/result.json only. PDF (slides) and PNG (markdown).

The look is the one of the earlier combination figures (white background, channel colours below),
with the CMS Open Data header and as little text inside the axes as possible.
"""

from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .paths import OUTPUT, PLOTS

C_MUMU, C_TAUTAU, C_EE, C_COMB = "#2B6CB0", "#EB811B", "#8E44AD", "#23373B"
C_PRED, C_GREY = "#14B03D", "#8C8C94"
COLOUR = {"ee": C_EE, "mumu": C_MUMU, "tautau": C_TAUTAU, "combined": C_COMB}
LABEL = {"ee": r"$Z\to ee$", "mumu": r"$Z\to\mu\mu$", "tautau": r"$Z\to\tau\tau$", "combined": r"$Z\to\ell\ell$ combined"}
XLABEL = r"$\sigma(\mathrm{pp}\to Z/\gamma^{*}\to\ell\ell)$, $60 < m_{\ell\ell} < 120$ GeV  [pb]"
ORDER = ("ee", "mumu", "tautau")

plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "figure.facecolor": "white", "savefig.facecolor": "white", "axes.facecolor": "white",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": ":",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def _header(ax, right="16.4 fb$^{-1}$ (13 TeV)"):
    ax.text(0.0, 1.015, r"$\mathbf{CMS}$ $\mathit{Open\ Data}$", transform=ax.transAxes, ha="left", va="bottom", fontsize=13)
    if right:
        ax.text(1.0, 1.015, right, transform=ax.transAxes, ha="right", va="bottom", fontsize=11)


def _save(fig, name):
    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTS / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(PLOTS / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _prediction_band(ax, pred, centre=None, label=True):
    centre = pred["sigma_pb"] if centre is None else centre
    lo = centre * (1 - pred["err_down_pb"] / pred["sigma_pb"])
    hi = centre * (1 + pred["err_up_pb"] / pred["sigma_pb"])
    ax.axvspan(lo, hi, color=C_PRED, alpha=0.18, lw=0, zorder=0,
               label="aMC@NLO (scale $\\oplus$ PDF $\\oplus$ $\\alpha_s$)" if label else None)
    ax.axvline(centre, color=C_PRED, lw=1.3, zorder=1)


def _value_text(v, up, down, digits=0):
    if abs(up - down) < 0.5 * 10 ** (-digits):
        return f"{v:.{digits}f} $\\pm$ {0.5 * (up + down):.{digits}f}"
    return f"{v:.{digits}f} $^{{+{up:.{digits}f}}}_{{-{down:.{digits}f}}}$"


def _rows(ax, rows, pred_centre, pred, xlim, value_column=True):
    """rows: dicts label, value, up, down, colour, own (bool), stat (optional). Top to bottom."""
    n = len(rows)
    _prediction_band(ax, pred, pred_centre)
    for i, r in enumerate(rows):
        y = n - 1 - i
        own = r.get("own", False)
        ax.errorbar(r["value"], y, xerr=[[r["down"]], [r["up"]]], fmt="o" if own else "s", color=r["colour"],
                    ms=7 if own else 5.5, capsize=4, lw=2.2 if own else 1.5, zorder=3)
        if r.get("stat"):
            ax.errorbar(r["value"], y, xerr=r["stat"], fmt="none", color=r["colour"], lw=4.5, capsize=0, zorder=4)
        if value_column:
            ax.text(1.02, y, _value_text(r["value"], r["up"], r["down"]), transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=10, color=r["colour"] if own else "#444444")
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["label"] for r in rows][::-1])
    for tick, r in zip(ax.get_yticklabels(), rows[::-1]):
        if r.get("own"):
            tick.set_color(r["colour"])
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(*xlim)
    ax.grid(axis="y", visible=False)


# ------------------------------------------------------------------------------------ figures
def summary(res):
    """Each channel alone and the combination, against the aMC@NLO prediction."""
    rows = []
    for k in ORDER:
        s = res["channels"][k]["standalone"]
        rows.append(dict(label=LABEL[k], value=s["sigma_pb"], up=s["err_up_pb"], down=s["err_down_pb"], colour=COLOUR[k], own=True))
    c = res["combined"]
    rows.append(dict(label="Combined", value=c["sigma_pb"], up=c["err_up_pb"], down=c["err_down_pb"],
                     colour=C_COMB, own=True, stat=c.get("stat_pb")))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    _rows(ax, rows, None, res["prediction"], (1650, 2400))
    ax.axhline(0.5, color="#999999", lw=0.8)
    ax.set_xlabel(XLABEL)
    ax.legend(loc="lower left", fontsize=9)
    _header(ax)
    _save(fig, "summary")


def channels_vs_published(res):
    """Every channel next to the published 13 TeV measurement(s) of the same decay."""
    pub = res["published"]["channel"]
    xl = {"ee": (1750, 2300), "mumu": (1700, 2150), "tautau": (1050, 2450)}
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.6), gridspec_kw={"hspace": 0.55})
    for ax, k in zip(axes, ORDER):
        s = res["channels"][k]["standalone"]
        rows = [dict(label="This work", value=s["sigma_pb"], up=s["err_up_pb"], down=s["err_down_pb"], colour=COLOUR[k], own=True)]
        for p in pub[k]:
            detail = p["detail"] + (" (66–116 GeV)" if p["window_scale"] != 1.0 else "")
            rows.append(dict(label=f'{p["label"]}\n{detail}', value=p["value_60_120_pb"], up=p["err_60_120_pb"],
                             down=p["err_60_120_pb"], colour=C_GREY))
        pred_centre = res["prediction"]["sigma_pb"] * res["channels"][k]["sigma_reference_pb"] / res["poi_reference_pb"]
        _rows(ax, rows, pred_centre, res["prediction"], xl[k])
        ax.tick_params(axis="y", labelsize=9.5)
        ax.text(0.5, 1.015, LABEL[k], transform=ax.transAxes, ha="center", va="bottom", fontsize=12, color=COLOUR[k])
        ax.set_xlabel(r"$\sigma$, $60 < m < 120$ GeV  [pb]", fontsize=11)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=1, fontsize=9, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(bottom=0.09)
    _header(axes[0])
    _save(fig, "channels_vs_published")


def combined_vs_published(res):
    c = res["combined"]
    rows = [dict(label="This work", value=c["sigma_pb"], up=c["err_up_pb"], down=c["err_down_pb"], colour=C_COMB, own=True,
                 stat=c.get("stat_pb"))]
    for p in res["published"]["combined"]:
        detail = p["detail"] + (" (66–116 GeV)" if p["window_scale"] != 1.0 else "")
        rows.append(dict(label=f'{p["label"]}\n{detail}', value=p["value_60_120_pb"], up=p["err_60_120_pb"],
                         down=p["err_60_120_pb"], colour=C_GREY))
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    _rows(ax, rows, None, res["prediction"], (1800, 2150))
    ax.set_xlabel(XLABEL)
    ax.legend(loc="lower right", fontsize=9)
    _header(ax)
    _save(fig, "combined_vs_published")


def nll_scan(res):
    scan = res["combined"].get("nll_scan")
    if not scan:
        return
    ref = res["poi_reference_pb"]
    x = np.array([p["X"] for p in scan]) * ref
    y = 2.0 * np.array([p["minusdeltaNLL"] for p in scan])
    y -= y.min()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(x, y, color=C_COMB, lw=2)
    for lvl, lab in ((1, r"68%"), (4, r"95%")):
        ax.axhline(lvl, color=C_GREY, lw=0.9, ls="--")
        ax.text(x.max(), lvl, lab, ha="right", va="bottom", fontsize=9, color=C_GREY)
    c = res["combined"]
    ax.set_xlim(max(x.min(), c["sigma_pb"] - 3.2 * c["err_down_pb"]), min(x.max(), c["sigma_pb"] + 3.2 * c["err_up_pb"]))
    ax.set_ylim(0, 10)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(r"$-2\,\Delta\ln L$")
    _header(ax)
    _save(fig, "nll_scan")


#: free normalisation factors of the tautau channel: fitted values, not pulls of a constrained parameter
FREE_FACTORS = ("TauIDSF_", "mu_ttbar")


def _np_label(name):
    if name == "mu_ttbar":
        return r"$\mu_{t\bar{t}}$ ($\tau\tau$, free)"
    return name.replace("_eeShape", " (ee shape)").replace("_", " ")


def pulls(res):
    """Post-fit values and constraints of every nuisance parameter of the combined fit (no gammas)."""
    cat = _categories()
    nps = {k: v for k, v in res["combined"]["pulls"].items() if not k.startswith(("gamma", "xsref") + FREE_FACTORS)}
    names = sorted(nps, key=lambda n: (cat.get(n, "zz"), n))
    n = len(names)
    fig, ax = plt.subplots(figsize=(6.2, 0.2 * n + 1.2))
    y = np.arange(n)[::-1]
    ax.axvspan(-2, 2, color="#FFE066", alpha=0.45, lw=0)
    ax.axvspan(-1, 1, color="#7CCB7C", alpha=0.45, lw=0)
    ax.axvline(0, color="k", lw=0.8)
    vals = np.array([nps[k]["value"] for k in names])
    errs = np.array([nps[k]["err"] for k in names])
    ax.errorbar(vals, y, xerr=errs, fmt="o", color="k", ms=3.5, lw=1.2, capsize=0)
    ax.set_yticks(y)
    ax.set_yticklabels([_np_label(k) for k in names], fontsize=7.5)
    prev = None
    for yi, k in zip(y, names):
        if prev is not None and cat.get(k) != prev:
            ax.axhline(yi + 0.5, color="#888888", lw=0.5)
        prev = cat.get(k)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel(r"$(\hat\theta - \theta_0)/\Delta\theta$")
    ax.grid(axis="y", visible=False)
    _header(ax, right=None)
    _save(fig, "pulls")


def _categories():
    """NP name -> Category string, from the generated channel configs of the common likelihood."""
    from . import trexcfg
    out = {}
    for cfg in sorted((OUTPUT.parent / "work" / "common").glob("*.config")):
        if cfg.name == "multifit.config":
            continue
        for b in trexcfg.read(cfg):
            if b.kind == "Systematic":
                out[b.opts.get("NuisanceParameter", b.name)] = b.opts.get("Category", "")
    return out


def impacts(res, top=20):
    """Pulls (left) and post-fit impact on sigma (right) of the nuisance parameters that matter most (no gammas)."""
    ranking = res["combined"].get("ranking")
    if not ranking:
        return
    ref = res["poi_reference_pb"]
    rows = [r for r in ranking if not r["name"].startswith("gamma")]
    rows = sorted(rows, key=lambda r: -(abs(r["post_up"]) + abs(r["post_down"])))[:top][::-1]
    n = len(rows)
    y = np.arange(n)
    fig, (axp, axi) = plt.subplots(1, 2, figsize=(8.6, 0.33 * n + 1.3), sharey=True,
                                   gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.06})
    axp.axvspan(-2, 2, color="#FFE066", alpha=0.45, lw=0)
    axp.axvspan(-1, 1, color="#7CCB7C", alpha=0.45, lw=0)
    axp.axvline(0, color="k", lw=0.8)
    con = [i for i, r in enumerate(rows) if not r["name"].startswith(FREE_FACTORS)]
    axp.errorbar([rows[i]["pull"] for i in con], y[con], xerr=[[rows[i]["err_down"] for i in con], [rows[i]["err_up"] for i in con]],
                 fmt="o", color="k", ms=4, lw=1.2, capsize=0)
    for i, r in enumerate(rows):          # a free factor has no prior to be pulled from: print what was fitted
        if r["name"].startswith(FREE_FACTORS):
            v = res["combined"]["pulls"][r["name"]]
            axp.text(0, y[i], f"free: {v['value']:.3f} $\\pm$ {v['err']:.3f}", ha="center", va="center", fontsize=8.5, color="#333333")
    axp.set_xlim(-2.6, 2.6)
    axp.set_xlabel(r"$(\hat\theta - \theta_0)/\Delta\theta$")
    axp.set_yticks(y)
    axp.set_yticklabels([_np_label(r["name"]) for r in rows], fontsize=9)
    axp.set_ylim(-0.7, n - 0.3)
    axp.grid(axis="y", visible=False)

    up = np.array([r["post_up"] for r in rows]) * ref
    dn = np.array([r["post_down"] for r in rows]) * ref
    lim = max(np.max(np.abs(up)), np.max(np.abs(dn)))
    axi.barh(y, up, height=0.72, color=C_MUMU, alpha=0.85, label=r"$\hat\theta+\Delta\hat\theta$")
    axi.barh(y, dn, height=0.72, color="#5FA8D3", alpha=0.85, label=r"$\hat\theta-\Delta\hat\theta$")
    axi.axvline(0, color="k", lw=0.8)
    axi.set_xlim(-1.35 * lim, 1.35 * lim)
    axi.set_xlabel(r"$\Delta\sigma$  [pb]")
    axi.grid(axis="y", visible=False)
    axi.tick_params(axis="y", length=0)
    axi.legend(loc="lower right", fontsize=8.5)
    _header(axp, right=None)
    axi.text(1.0, 1.015, "16.4 fb$^{-1}$ (13 TeV)", transform=axi.transAxes, ha="right", va="bottom", fontsize=11)
    _save(fig, "impacts")


def breakdown(res):
    """Uncertainty of the combined cross section by group (stat-only refits with each group fixed)."""
    g = res["combined"].get("grouped_impacts")
    if not g:
        return
    rename = {"Gammas": "MC statistics"}
    items = [(rename.get(k, k), v["impact_pb"]) for k, v in g.items() if k not in ("FullSyst", "NormFactors") and v["impact_pb"] >= 0.1]
    c = res["combined"]
    total = 0.5 * (c["err_up_pb"] + c["err_down_pb"])
    items.append(("Data statistics", c.get("stat_pb", float("nan"))))
    items = sorted(items, key=lambda kv: kv[1])
    fig, ax = plt.subplots(figsize=(6.4, 0.34 * len(items) + 1.4))
    names = [k for k, _ in items]
    vals = [v for _, v in items]
    ax.barh(names, vals, color=C_COMB, height=0.65)
    ax.axvline(total, color=C_GREY, ls="--", lw=1)
    ax.text(total, len(items) - 0.4, " total", color=C_GREY, fontsize=9, va="center")
    ax.set_xlabel(r"uncertainty on $\sigma$  [pb]")
    ax.set_xlim(0, total * 1.12)
    ax.grid(axis="y", visible=False)
    _header(ax)
    _save(fig, "breakdown")


def tautau_channels(res):
    """What the tautau line is made of: z-tautau's own fits of its final states, and the line as it enters the combination."""
    ch = res["channels"]["tautau"]
    sub = ch.get("published_submeasurements") or {}
    if not sub:
        return
    labels = {"tau_h tau_h": r"$\tau_h\tau_h$ alone (POG $\tau_h$ ID SF)", "mu tau_h": r"$\mu\tau_h$ alone (POG $\tau_h$ ID SF)",
              "e tau_h": r"$e\tau_h$ alone (POG $\tau_h$ ID SF)", "e mu": r"$e\mu$ alone (no $\tau_h$)",
              "tau channels, SF free": r"$\tau_h$ channels, free $\tau_h$ ID SF", "four channels": "four channels (z-tautau result)",
              "four channels, no tau_h ID pT uncertainty": r"four channels, no $\tau_h$ ID $p_T$ uncertainty",
              "four channels, pT-split SF": r"four channels, $p_T$-split $\tau_h$ ID SF",
              "four channels, e mu trigger prior x2": r"four channels, $e\mu$ trigger prior $\times 2$"}
    rows = [dict(label=labels.get(k, k), value=v["sigma_pb"], up=v["err_up_pb"], down=v["err_down_pb"],
                 colour=C_TAUTAU if k == "four channels" else C_GREY) for k, v in sub.items()]
    s = ch["standalone"]
    rows.append(dict(label="in the combination\n(standalone fit, same model)", value=s["sigma_pb"], up=s["err_up_pb"], down=s["err_down_pb"],
                     colour=C_TAUTAU, own=True))
    fig, ax = plt.subplots(figsize=(7.2, 0.5 * len(rows) + 1.4))
    pred_centre = res["prediction"]["sigma_pb"] * ch["sigma_reference_pb"] / res["poi_reference_pb"]
    _rows(ax, rows, pred_centre, res["prediction"], (1600, 2650))
    ax.axhline(0.5, color="#999999", lw=0.8)
    ax.axhline(len(rows) - 4.5, color="#999999", lw=0.8, ls=":")
    ax.set_xlabel(r"$\sigma(\mathrm{pp}\to Z/\gamma^{*}\to\tau\tau)$, $60 < m_{\tau\tau} < 120$ GeV  [pb]")
    ax.legend(loc="lower right", fontsize=8.5)
    _header(ax)
    _save(fig, "tautau_channels")


def variations(res):
    """The combined cross section under the alternative likelihoods of config/channels.json."""
    v = res["variations"]
    if not v:
        return
    labels = {"without_ee": r"without $Z\to ee$", "without_tautau": r"without $Z\to\tau\tau$",
              "split_all_channels": "shape/norm. split in all channels",
              "ee_split_shared_only": "ee: split shared NPs only",
              "ee_electron_id_1p2": r"ee: electron ID norm. $\pm$1.2%",
              "tautau_without_tauidpt": r"$\tau\tau$ without the $\tau_h$ ID $p_T$ uncertainty",
              "tautau_ptsplit": r"$\tau\tau$: $p_T$-split $\tau_h$ ID SF model",
              "tautau_emu_only": r"$\tau\tau$: $e\mu$ channel only",
              "tautau_tauh_only": r"$\tau\tau$: $\tau_h$ channels only",
              "tautau_emutrig2x": r"$\tau\tau$: $e\mu$ trigger prior $\times 2$",
              "tautau_leptons_correlated": r"$\tau\tau$: $\mu$/e NPs shared with $\mu\mu$/ee",
              "tautau_theory_decorrelated": r"$\tau\tau$: theory NPs not shared",
              "ttbar_cr_for_all_channels": r"$t\bar{t}$: $e\mu$ control region for all channels",
              "ttbar_xs_constrained": r"$t\bar{t}$: 6% prior everywhere, no free $\mu_{t\bar{t}}$"}
    c = res["combined"]
    rows = [dict(label="Baseline", value=c["sigma_pb"], up=c["err_up_pb"], down=c["err_down_pb"], colour=C_COMB, own=True)]
    for k, x in v.items():
        rows.append(dict(label=labels.get(k, k), value=x["sigma_pb"], up=x["err_up_pb"], down=x["err_down_pb"], colour=C_GREY))
    fig, ax = plt.subplots(figsize=(7.2, 0.5 * len(rows) + 1.3))
    _rows(ax, rows, None, res["prediction"], (1780, 2060))
    ax.axvspan(c["sigma_pb"] - c["err_down_pb"], c["sigma_pb"] + c["err_up_pb"], color=C_COMB, alpha=0.07, lw=0)
    ax.set_xlabel(XLABEL)
    ax.legend(loc="upper left", fontsize=8.5)
    _header(ax)
    _save(fig, "variations")


def make_all():
    res = json.loads((OUTPUT / "result.json").read_text())
    for f in (summary, channels_vs_published, combined_vs_published, tautau_channels, breakdown, impacts, pulls, nll_scan, variations):
        f(res)
    print("figures in", PLOTS)
