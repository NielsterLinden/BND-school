#!/usr/bin/env python
"""Every figure of the slide deck, drawn dark and vector into slides/figs/ (PDF + PNG).

    source ../fitting/setup.sh && python slides/make_figures.py

Nothing here is typed by hand: the numbers come from output/results.json, the fit results in
fit/results/*_fit_result.json, the fit inputs and their post-fit histograms, output/data/fakes_*.json and
external/trigger_insitu_v4.json. Figures are drawn on `plotting.DARK_BG`, which is the deck background, so
an embedded plot has no visible edge (.claude/prompts/presentation_style.md).

The deck itself is slides/build_deck.py.
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from ztautau import config, plotting  # noqa: E402
from ztautau.plotting import DARK_BG, DARK_FG, dark_rc, save  # noqa: E402

FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)
BLUE, GREEN, ORANGE, PINK, MUTED, RED = "#2BA4DD", "#8AC63F", "#D9822B", "#E377AC", "#99999E", "#CC4C4C"

R = json.loads((config.OUTPUT_DIR_V4 / "results.json").read_text())
FIT = R["fit"]["combined"]
VAR = R.get("fit_variants", {})
FREE = R.get("fit_per_channel_free", {})
META = json.loads((config.FIT_DIR_V4 / "fitinputs" / f"{config.JOB_V4}.root.meta.json").read_text())
FITTED_REGIONS = [r for r in META["regions"] if META.get("region_sets", {}).get(r, "nominal") == "nominal"]

CH = {"tautau": r"$\tau_h\tau_h$", "mutau": r"$\mu\tau_h$", "etau": r"$e\tau_h$", "emu": r"$e\mu$"}
REGION_LABELS = {"tautau_SR0": r"$\tau_h\tau_h$, BDT < 0.55  ($m_{\tau\tau}$ > 110 GeV only)",
                 "tautau_SR1": r"$\tau_h\tau_h$, 0.55 < BDT < 0.90", "tautau_SR2": r"$\tau_h\tau_h$, BDT > 0.90",
                 "emu_SR": r"$e\mu$ signal region", "emu_CRtt": r"$e\mu$, $t\bar{t}$ control region"}
for _c, _l in (("mutau", r"$\mu\tau_h$"), ("etau", r"$e\tau_h$")):
    for _dm in config.TAU_DMS:
        REGION_LABELS[f"{_c}_SR_dm{_dm}"] = f"{_l}, decay mode {_dm}"
STACK_GROUPS = [("DYll", ["DYee", "DYmumu"]), ("DYlowmass", ["DYlowmass"]), ("Diboson", ["WW", "WZ", "ZZ"]),
                ("Top", ["TTbar", "SingleTop"]), ("WJets", ["WJets"]), ("Fakes", ["Fakes"]),
                ("DYtautau_nonfid", ["DYtautau_out"]), ("DYtautau", ["DYtautau"])]
plt.rcParams.update({"font.size": 13})


def dark_fig(figsize=(7, 4.4), n=1, **kw):
    plt.rcParams.update(dark_rc())
    fig, axes = plt.subplots(1, n, figsize=figsize, **kw)
    for ax in (axes if n > 1 else [axes]):
        ax.set_facecolor(DARK_BG)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.patch.set_facecolor(DARK_BG)
    plotting.tag(axes[-1] if n > 1 else axes, x=1.0, y=1.01, ha="right", va="bottom", color=DARK_FG, fontsize=10)
    return fig, axes


def _sig(f):
    return f["sigma_60_120_pb"]


# ------------------------------------------------------------------ 1. the result
def summary():
    """The measurement, the two sub-measurements it combines, and each channel on its own."""
    rows = [(r"**combined**: $\tau_h\tau_h + \mu\tau_h + e\tau_h + e\mu$".replace("**", ""), _sig(FIT), BLUE, 7)]
    if FREE.get("emu"):
        rows.append((r"$e\mu$ alone  (no $\tau_h$, no ID scale factor)", _sig(FREE["emu"]), GREEN, 6))
    if VAR.get("taulep"):
        rows.append((r"$\tau_h\tau_h + \mu\tau_h + e\tau_h$ alone  (SF free)", _sig(VAR["taulep"]), ORANGE, 6))
    for ch in config.CHANNELS:
        if ch in R["fit"]:
            rows.append((CH[ch] + r" alone, POG $\tau_h$ ID SF", _sig(R["fit"][ch]), MUTED, 5))
    fig, ax = dark_fig((9.2, 1.2 + 0.58 * len(rows)))
    pred = _sig(FIT)["prediction"]
    ax.axvspan(pred * 0.96, pred * 1.04, color=ORANGE, alpha=0.20)
    ax.axvline(pred, color=ORANGE, lw=1.5)
    ax.text(pred + 18, len(rows) + 0.45, f"NLO prediction {pred:.0f} pb ($\\pm$4%)", fontsize=9.5, color=ORANGE)
    for i, (lab, s, col, ms) in enumerate(rows):
        y = len(rows) - i
        ax.errorbar([s["value"]], [y], xerr=[[s["err_down"]], [s["err_up"]]], fmt="none", color=col, capsize=4, lw=2, zorder=2)
        if s.get("stat"):
            ax.errorbar([s["value"]], [y], xerr=[[s["stat"]], [s["stat"]]], fmt="none", color=RED, lw=5, zorder=3)
        ax.plot([s["value"]], [y], "o", color=col, ms=ms, zorder=4)
        ax.text(1330, y + 0.28, lab, fontsize=10.5, color=DARK_FG)
        ax.text(2860, y + 0.28, f"{s['value']:.0f}  +{s['err_up']:.0f} $-${s['err_down']:.0f} pb",
                fontsize=10.5, color=DARK_FG, ha="right")
    ax.set_xlim(1300, 2900); ax.set_ylim(0.2, len(rows) + 0.95); ax.set_yticks([])
    ax.set_xlabel(r"$\sigma(pp\rightarrow Z/\gamma^*\rightarrow\tau\tau,\ 60<m<120\ \mathrm{GeV})$ [pb]"
                  "      (red: statistical only)")
    save(fig, FIGS / "summary", dark=True)


def submeasurements():
    """Why the combined number is a compromise: mu_Z of the two halves and of the cross-check fits."""
    rows = []
    if FREE.get("emu"):
        rows.append((r"$e\mu$ alone", FREE["emu"], GREEN))
    if VAR.get("taulep"):
        rows.append((r"$\tau$ channels alone, SF free", VAR["taulep"], ORANGE))
    rows.append(("combined (the measurement)", FIT, BLUE))
    if VAR.get("emutrig2x"):
        rows.append((r"combined, $e\mu$ trigger prior doubled", VAR["emutrig2x"], MUTED))
    if VAR.get("ptsplit"):
        rows.append((r"combined, $\tau_h$ ID SF split at $p_T$ = 40 GeV", VAR["ptsplit"], PINK))
    fig, ax = dark_fig((8.6, 1.1 + 0.62 * len(rows)))
    ax.axvline(1.0, color=MUTED, lw=1.0, ls=":")
    for i, (lab, f, col) in enumerate(rows):
        y = len(rows) - i
        ax.errorbar([f["mu"]], [y], xerr=[[f["mu_err_down"]], [f["mu_err_up"]]], fmt="o", color=col,
                    capsize=4, lw=2.2, ms=7 if col == BLUE else 6)
        ax.text(0.845, y + 0.3, lab, fontsize=10.5, color=DARK_FG)
        ax.text(1.335, y + 0.3, f"{f['mu']:.3f} +{f['mu_err_up']:.3f} $-${f['mu_err_down']:.3f}",
                fontsize=10.5, color=DARK_FG, ha="right")
    t = R.get("tension_taulep_vs_emu")
    if t:
        ax.text(0.845, 0.45, f"the two sub-measurements are {t['n_sigma']:.1f}$\\sigma$ apart "
                             "(uncorrelated limit)", fontsize=10, color=RED)
    ax.set_xlim(0.84, 1.34); ax.set_ylim(0.15, len(rows) + 0.95); ax.set_yticks([])
    ax.set_xlabel(r"$\mu_Z$   (signal strength w.r.t. the NLO prediction)")
    save(fig, FIGS / "submeasurements", dark=True)


# ------------------------------------------------------------------ 2. what the fit measures in situ
def tau_id():
    sf = {k: v for k, v in (FIT.get("tau_id_sf") or {}).items() if not k.endswith("_lowpt")}
    fig, ax = dark_fig((6.8, 4.3))
    for i, (dm, v) in enumerate(sf.items()):
        ax.errorbar(i - 0.13, v["pog"][0], v["pog"][1], fmt="s", color=MUTED, ms=7, capsize=4,
                    label="TauPOG, external" if i == 0 else None)
        ax.errorbar(i + 0.13, v["value"], [[v["err_down"]], [v["err_up"]]], fmt="o", color=BLUE, ms=7, capsize=4,
                    label="this fit, in situ" if i == 0 else None)
    ax.set_xticks(range(len(sf))); ax.set_xticklabels([f"DM {k[2:]}" for k in sf])
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax.set_ylabel("DeepTau VSjet Tight scale factor"); ax.set_ylim(0.55, 1.35)
    ax.legend(loc="upper left", fontsize=10)
    ax.set_title(r"$\tau_h$ identification scale factor per decay mode", fontsize=12, pad=14, color=DARK_FG)
    save(fig, FIGS / "tauid", dark=True)

    es = FIT.get("tau_es") or {}
    if es:
        fig, ax = dark_fig((6.8, 3.7))
        for i, (dm, v) in enumerate(es.items()):
            ax.errorbar(i, v["pull"], v["constraint"], fmt="o", color=GREEN, ms=7, capsize=4)
            ax.text(i, v["pull"] + v["constraint"] + 0.14, f"{v['constraint'] * v['prior_pct']:.1f}%",
                    ha="center", fontsize=10, color=DARK_FG)
        ax.axhspan(-1, 1, color=MUTED, alpha=0.15); ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xticks(range(len(es))); ax.set_xticklabels([f"DM {k[2:]}" for k in es]); ax.set_ylim(-2.2, 2.6)
        ax.set_ylabel("pull, in units of the 3% prior")
        ax.set_title(r"$\tau_h$ energy scale: pull and post-fit constraint (label: absolute)",
                     fontsize=11.5, pad=12, color=DARK_FG)
        save(fig, FIGS / "taues", dark=True)


def ptsplit():
    """The one assumption the tau_h tau_h / (l tau_h)^2 lever rests on, measured."""
    v = VAR.get("ptsplit")
    if not v or not v.get("tau_id_sf"):
        return
    sf = v["tau_id_sf"]
    dms = [dm for dm in config.TAU_DMS if f"DM{dm}" in sf and f"DM{dm}_lowpt" in sf]
    fig, ax = dark_fig((7.4, 4.3))
    for i, dm in enumerate(dms):
        hi, lo = sf[f"DM{dm}"], sf[f"DM{dm}_lowpt"]
        ax.errorbar(i - 0.13, hi["value"], [[hi["err_down"]], [hi["err_up"]]], fmt="o", color=BLUE, ms=7, capsize=4,
                    label=r"$p_T(\tau_h)$ > 40 GeV" if i == 0 else None)
        ax.errorbar(i + 0.13, lo["value"], [[lo["err_down"]], [lo["err_up"]]], fmt="D", color=PINK, ms=7, capsize=4,
                    label=r"30 < $p_T(\tau_h)$ < 40 GeV" if i == 0 else None)
        ax.text(i, 1.28, f"{100 * (lo['value'] / hi['value'] - 1):+.0f}%", ha="center", fontsize=11, color=PINK)
    ax.set_xticks(range(len(dms))); ax.set_xticklabels([f"DM {d}" for d in dms])
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax.set_ylabel("DeepTau VSjet Tight scale factor"); ax.set_ylim(0.6, 1.4)
    ax.legend(loc="lower left", fontsize=10)
    ax.set_title(r"the scale factor is not flat in $p_T$:  $\mu_Z$ moves "
                 f"{v['mu'] - FIT['mu']:+.3f} to {v['mu']:.3f}", fontsize=11.5, pad=14, color=DARK_FG)
    save(fig, FIGS / "ptsplit", dark=True)


# ------------------------------------------------------------------ 3. uncertainties
def impacts():
    gi = {k: val for k, val in FIT["grouped_impact"].items() if k not in ("FullSyst", "Total")}
    if FIT.get("mu_stat"):
        gi["Data statistics"] = FIT["mu_stat"]
    items = sorted(gi.items(), key=lambda kv: kv[1])
    fig, ax = dark_fig((7.6, 0.32 * len(items) + 1.7))
    ax.barh([k for k, _ in items], [100 * v for _, v in items],
            color=[GREEN if k == "Data statistics" else BLUE for k, _ in items])
    for i, (_, val) in enumerate(items):
        ax.text(100 * val + 0.05, i, f"{100 * val:.2f}%", va="center", fontsize=9, color=DARK_FG)
    tot = 100 * 0.5 * (FIT["mu_err_up"] + FIT["mu_err_down"])
    ax.axvline(tot, color=RED, ls="--", lw=1.4)
    ax.text(tot + 0.05, len(items) - 1.2, f"MINOS total {tot:.1f}%", fontsize=10, color=RED)
    ax.set_xlabel(r"impact on $\mu_Z$ [%]")
    ax.set_title("the categories overlap: their quadrature sum exceeds the total by "
                 f"{1 / (FIT.get('grouped_impact_scale') or 1):.2f}", fontsize=11, pad=12, color=DARK_FG)
    save(fig, FIGS / "impacts", dark=True)


def ranking():
    rows = (FIT.get("ranking") or [])[:16]
    if not rows:
        return
    rows = sorted(rows, key=lambda r: max(abs(r["impact_up"]), abs(r["impact_down"])))
    fig, ax = dark_fig((8.2, 0.36 * len(rows) + 1.6))
    y = np.arange(len(rows))
    ax.barh(y, [100 * r["impact_up"] for r in rows], color=BLUE, height=0.62, label=r"+1$\sigma$")
    ax.barh(y, [100 * r["impact_down"] for r in rows], color=ORANGE, height=0.62, label=r"$-$1$\sigma$")
    ax.set_yticks(y); ax.set_yticklabels([r["name"] for r in rows], fontsize=9.5)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel(r"post-fit impact on $\mu_Z$ [%]"); ax.legend(fontsize=10, loc="lower right")
    ax.set_title("the parameters the data determine, and what they move", fontsize=11.5, pad=12, color=DARK_FG)
    save(fig, FIGS / "ranking", dark=True)

    pulls = {k: v for k, v in (FIT.get("pulls") or {}).items()
             if not k.startswith("gamma") and not k.startswith(("mu_", "TauIDSF"))}
    if not pulls:
        return
    names = sorted(pulls, key=lambda k: -abs(pulls[k][0]))[:24]
    names = names[::-1]
    fig, ax = dark_fig((8.2, 0.33 * len(names) + 1.6))
    y = np.arange(len(names))
    ax.axvspan(-2, 2, color=MUTED, alpha=0.10); ax.axvspan(-1, 1, color=MUTED, alpha=0.18)
    ax.errorbar([pulls[n][0] for n in names], y,
                xerr=[[abs(pulls[n][2]) for n in names], [abs(pulls[n][1]) for n in names]],
                fmt="o", color=BLUE, ms=5, capsize=3, lw=1.4)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=9)
    ax.set_xlim(-3, 3); ax.set_xlabel(r"pull  $(\hat{\theta}-\theta_0)/\Delta\theta$   (band: the prior)")
    ax.set_title("the 24 most-pulled nuisance parameters", fontsize=11.5, pad=12, color=DARK_FG)
    save(fig, FIGS / "pulls", dark=True)


# ------------------------------------------------------------------ 4. the fitted regions
# The post-fit content comes from TRExFitter's own per-bin yaml (Plots/<region>_postfit.yaml), which lists
# every sample, the total and its uncertainty. The `Histograms/<region>_postFit.root` file next to it holds
# `h_tot_postFit`, which is the post-fit total *without* the signal (the job sets PlotOptions: NOSIG, so
# TRExFitter's own plots draw the signal separately). Using it as the total, as the first version of these
# figures did, left Z -> tautau out of the stack and out of the data/prediction ratio.
TITLE_TO_SAMPLE = {v: k for k, v in META["titles"].items()}


def postfit_yaml(region):
    """(edges, data, [(sample, per-bin yield)], per-bin uncertainty of the total) of one region, post-fit."""
    import yaml
    p = config.FIT_DIR_V4 / f"results/{config.JOB_V4}/Plots/{region}_postfit.yaml"
    if not p.exists():
        return None
    d = yaml.safe_load(p.read_text())
    edges = np.asarray(d["Figure"][0]["BinEdges"], dtype=float)
    data = np.asarray(d["Data"][0]["Yield"], dtype=float)
    tot = d["Total"][0]
    unc = 0.5 * (np.abs(np.asarray(tot["UncertaintyUp"], dtype=float))
                 + np.abs(np.asarray(tot["UncertaintyDown"], dtype=float)))
    samples = {TITLE_TO_SAMPLE.get(s["Name"], s["Name"]): np.asarray(s["Yield"], dtype=float) for s in d["Samples"]}
    return edges, data, samples, unc


def regions():
    for region in FITTED_REGIONS:
        got = postfit_yaml(region)
        if got is None:
            print(f"  no post-fit yaml for {region}, skipped")
            continue
        edges, data, samples, unc = got
        stack = []
        for g, members in STACK_GROUPS:
            parts = [samples[m] for m in members if m in samples]
            if parts:
                v = np.maximum(sum(parts), 0)
                stack.append((g, v, np.zeros_like(v)))     # the band below is the full post-fit uncertainty
        total = sum(v for _, v, _ in stack)
        if abs(total.sum() - sum(s.sum() for s in samples.values())) > 1.0:
            missing = [k for k in samples if not any(k in m for _, m in STACK_GROUPS)]
            print(f"  WARNING {region}: samples not drawn: {missing}")
        plotting.stack_plot(FIGS / f"region_{region}.png", edges, (data, np.sqrt(np.maximum(data, 0))), stack,
                            r"$m_{\tau\tau}$ [GeV]", title=REGION_LABELS.get(region, region), dark=True,
                            band=unc, density=True, figsize=(5.6, 7.8))


# ------------------------------------------------------------------ 5. backgrounds measured from data
def fakes():
    for ch in ("mutau", "etau"):
        p = config.DATA_DIR_V4 / f"fakes_{ch}.json"
        if not p.exists():
            continue
        fk = json.loads(p.read_text())
        pt = np.asarray(fk["tables"]["qcd"]["pt_bins"]); x = 0.5 * (pt[1:] + pt[:-1]); x[-1] = pt[-2] + 15
        fig, axes = dark_fig((11, 4.3), n=2, sharey=True)
        for ax, dm in zip(axes, (1, 10)):
            i = list(config.TAU_DMS).index(dm)
            for proc, col, lab in (("qcd", BLUE, "multijet (same sign)"), ("w", RED, r"W+jets ($m_T$ > 70)"),
                                   ("w_ss", ORANGE, "W+jets (same sign)"), ("tt", GREEN, r"$t\bar{t}$ (simulation)")):
                if proc not in fk["tables"]:
                    continue
                t = fk["tables"][proc]
                ax.errorbar(x, np.asarray(t["ff"])[i, 0], np.asarray(t["err"])[i, 0], fmt="o", color=col, ms=4,
                            label=lab if dm == 1 else None)
            ax.set_title(f"decay mode {dm}, 0 jets", color=DARK_FG, fontsize=11)
            ax.set_xlim(30, 120); ax.set_ylim(0, 0.25); ax.grid(alpha=0.2)
            ax.set_xlabel(r"$p_T(\tau_h)$ [GeV]")
        axes[0].set_ylabel("fake factor  (Tight / VVVLoose-not-Tight)")
        axes[0].legend(fontsize=9)
        fig.suptitle(f"{CH[ch].strip('$')}: per-process fake factors — "
                     f"C(OS/SS) = {fk['osss']['C']:.2f} $\\pm$ {fk['osss']['stat']:.2f}, "
                     f"same-sign closure {fk['closure_ss']['ratio']:.3f} $\\pm$ {fk['closure_ss']['stat']:.3f}",
                     color=DARK_FG, fontsize=11.5)
        fig.tight_layout()
        save(fig, FIGS / f"fakefactors_{ch}", dark=True)

    p = config.DATA_DIR_V4 / "fakes_emu.json"
    if p.exists():
        fk = json.loads(p.read_text())
        e = np.asarray(fk["osss"]["dr_edges"]); x = 0.5 * (e[1:] + e[:-1])
        fig, ax = dark_fig((6.8, 4.2))
        ax.errorbar(x, fk["SB1"]["ratio"], fk["SB1"]["stat"], fmt="o", color=BLUE,
                    label=r"SB1: both $I_{rel}$ < 0.5, one > 0.15")
        ax.errorbar(x + 0.06, fk["SB2"]["ratio"], fk["SB2"]["stat"], fmt="s", color=ORANGE,
                    label=r"SB2: one $I_{rel}$ > 0.3")
        ax.set_xlabel(r"$\Delta R(e,\mu)$"); ax.set_ylabel("multijet OS / SS")
        ax.set_ylim(0, 3.5); ax.legend(fontsize=9)
        ax.set_title(f"$e\\mu$ multijet: {fk['ss_region']['data']} same-sign events in data, "
                     f"{fk['ss_region']['mc']:.0f} simulated $\\rightarrow$ {fk['sr_multijet']:.0f} in the SR",
                     fontsize=11, pad=12, color=DARK_FG)
        save(fig, FIGS / "emu_osss", dark=True)


def triggers():
    p = HERE.parent / "external" / "trigger_insitu_v4.json"
    if not p.exists():
        return
    t = json.loads(p.read_text())
    panels = (("ele27", r"Ele27_WPTight  ($e\tau_h$)"), ("emu_e", r"$e\mu$ cross trigger, electron leg"),
              ("emu_mu", r"$e\mu$ cross trigger, muon leg"))
    fig, axes = dark_fig((13, 4.2), n=3, sharey=True)
    for ax, (name, title) in zip(axes, panels):
        tb = t[name]; xe = np.asarray(tb["x_edges"]); x = 0.5 * (xe[1:] + xe[:-1]); x[-1] = xe[-2] * 1.3
        for j in range(len(tb["y_edges"]) - 1):
            ax.errorbar(x * (1 + 0.02 * j), np.asarray(tb["sf"])[:, j], np.asarray(tb["err"])[:, j], fmt="o", ms=4,
                        color=[BLUE, GREEN, ORANGE, PINK][j],
                        label=r"$|\eta|$ " + f"{tb['y_edges'][j]}–{tb['y_edges'][j + 1]}")
        ax.set_xscale("log"); ax.set_ylim(0.7, 1.2); ax.axhline(1, color=MUTED, lw=0.8, ls=":")
        ax.set_title(title, color=DARK_FG, fontsize=11); ax.grid(alpha=0.2)
        ax.set_xlabel(r"lepton $p_T$ [GeV]")
    axes[0].set_ylabel("scale factor, data / simulation")
    axes[0].legend(fontsize=8)
    axes[0].axvspan(config.EL_PT_MIN_ETAU, config.ELE27_PLATEAU_PT, color=RED, alpha=0.15)
    axes[0].text(30, 0.74, "turn-on:\nown NP", fontsize=8, color=RED)
    fig.tight_layout()
    save(fig, FIGS / "triggers", dark=True)


if __name__ == "__main__":
    summary(); submeasurements(); tau_id(); ptsplit(); impacts(); ranking(); regions(); fakes(); triggers()
    print(f"{len(list(FIGS.glob('*.pdf')))} figures -> {FIGS}")
