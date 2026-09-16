#!/usr/bin/env python
"""Uncertainty parity with a published measurement: budget, three estimate options, comparison, plots.

Shared by the three channels (procedure: fitting/UNCERTAINTY_PARITY.md; agent:
.claude/agents/uncertainty-parity-auditor.md). A channel writes one *parity file* (JSON) that maps
every uncertainty row of the published measurement onto its own analysis; this module turns it into

  * the total uncertainty of the channel's result for each of the three options used for the sources
    it could not measure ("estimate" rows): a = the published number, b = zero, c = our estimate;
  * the comparison with the published value (difference / sigma), without and with the theory rows
    that the two measurements share treated as correlated;
  * one figure per option (budget side by side + comparison), an overview of the three options, and a
    markdown table for the channel docs.

    python fitting/uncertainty_parity.py z-mumu/output/v2/cms_parity/parity_zmumu.json \\
        --plot-dir z-mumu/output/v2/plots --prefix cms_parity

Parity file (all uncertainties in percent of the respective central value):

    {"schema": "uncertainty-parity/1", "channel": "zmumu", "channel_label": "Z->mumu (this work)",
     "colour": "#2a78d6", "observable": "...", "unit": "pb",
     "reference": {"label": "CMS-SMP-20-004", "cite": "arXiv:2408.03744, Table 7", "value": 1952,
                   "stat_pct": 0.2, "lumi_pct": 2.3},
     "measurement": {"value": 1930.7, "stat_pct": 0.03, "syst_pct_measured": 1.04, "label": "..."},
     "rows": [{"id": "lumi", "source": "Luminosity", "group": "luminosity",
               "reference_pct": 2.3, "ours_pct": 1.2,
               "status": "similar|better|partial|missing|not-applicable|ours-only",
               "treatment": "measured|external|estimate|none",
               "ours_how": "...", "reference_how": "...", "rho_reference": 0.0,
               "components": {"name": pct, ...},
               "estimate": null | {"why_not_measured": "...", "why_not_cited": "...",
                                   "options": {"a": pct, "b": 0.0, "c": pct}, "derivation": "..."}}]}

`ours_pct` of an estimate row is ignored: the option decides. Rows are individual impacts and are shown as
such; post-fit impacts added in quadrature overshoot a fit's correlated total, so `syst_pct_measured` (if
given) is the channel's proper systematic uncertainty without luminosity and without the estimate rows,
and the total for an option is `syst_pct_measured` ⊕ its estimate rows (the row sum otherwise). `rho_reference` is the correlation of that
row between the two measurements (theory rows computed with the same generator and PDF set); rows that
the reference does not have (`reference_pct` null) contribute nothing to the correlated term.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

OPTIONS = {"a": "3a: published value", "b": "3b: set to zero", "c": "3c: our estimate"}
STATUS_TEXT = {"similar": "similar", "better": "ours better", "partial": "partial / eyeballed",
               "missing": "missing", "not-applicable": "n/a", "ours-only": "ours only", "measured-now": "measured now"}
INK, INK_2, GRID, REF, BG = "#1f1f1f", "#5c5c5c", "#e6e6e6", "#6b6b6b", "white"
THEMES = {"light": dict(INK="#1f1f1f", INK_2="#5c5c5c", GRID="#e6e6e6", REF="#6b6b6b", BG="white"),
          "dark": dict(INK="#E6E6E6", INK_2="#99999E", GRID="#3a3a40", REF="#8C8C94", BG="#222222")}


def set_theme(name):
    """'light' (docs, default) or 'dark' (the house slide style: #222222 ground, like the deck)."""
    globals().update(THEMES[name])


# ----------------------------------------------------------------------------- numbers
def load(path):
    with open(path) as fh:
        p = json.load(fh)
    if p.get("schema") != "uncertainty-parity/1":
        raise ValueError(f"{path}: not an uncertainty-parity/1 file")
    for r in p["rows"]:
        if r.get("treatment") == "estimate" and not r.get("estimate"):
            raise ValueError(f"row {r['id']}: treatment 'estimate' needs an 'estimate' block with three options")
    return p


def ours_pct(row, option):
    if row.get("treatment") == "estimate":
        return float(row["estimate"]["options"][option])
    return float(row.get("ours_pct") or 0.0)


def budget(p, option):
    """Dict with the total and its parts, in percent of our central value."""
    rows = {r["id"]: ours_pct(r, option) for r in p["rows"]}
    lumi = sum(v ** 2 for r, v in zip(p["rows"], rows.values()) if r.get("group") == "luminosity") ** 0.5
    row_sum = sum(v ** 2 for r, v in zip(p["rows"], rows.values()) if r.get("group") != "luminosity") ** 0.5
    measured = p["measurement"].get("syst_pct_measured")
    if measured is None:
        syst = row_sum
    else:
        est = sum(v ** 2 for r, v in zip(p["rows"], rows.values())
                  if r.get("group") != "luminosity" and r.get("treatment") == "estimate")
        syst = float(np.sqrt(measured ** 2 + est))
    stat = float(p["measurement"].get("stat_pct", 0.0))
    return {"rows": rows, "stat": stat, "syst": syst, "syst_row_sum": row_sum, "lumi": lumi,
            "total": float(np.sqrt(stat ** 2 + syst ** 2 + lumi ** 2))}


def reference_total_pct(p):
    ref = p["reference"]
    if "total_pct" in ref:
        return float(ref["total_pct"])
    syst = sum((r.get("reference_pct") or 0.0) ** 2 for r in p["rows"] if r.get("group") != "luminosity") ** 0.5
    return float(np.sqrt(ref.get("stat_pct", 0.0) ** 2 + syst ** 2 + ref.get("lumi_pct", 0.0) ** 2))


def compare(p, option, correlated=True):
    """Difference ours - reference and its uncertainty; theory rows with rho_reference > 0 enter with their
    correlation when `correlated`."""
    m, ref = p["measurement"]["value"], p["reference"]["value"]
    b = budget(p, option)
    s_o = m * b["total"] / 100
    s_r = ref * reference_total_pct(p) / 100
    cov = 0.0
    if correlated:
        for r in p["rows"]:
            rho, rp = float(r.get("rho_reference") or 0.0), r.get("reference_pct")
            if rho and rp:
                cov += rho * (m * ours_pct(r, option) / 100) * (ref * rp / 100)
    var = s_o ** 2 + s_r ** 2 - 2 * cov
    diff = m - ref
    return {"diff": diff, "sigma": float(np.sqrt(max(var, 0.0))), "pull": diff / np.sqrt(max(var, 1e-12)),
            "ours_total": s_o, "reference_total": s_r, "covariance": cov}


def summary(p):
    out = {}
    for o in OPTIONS:
        b = budget(p, o)
        out[o] = {"total_pct": b["total"], "syst_pct": b["syst"], "syst_row_sum_pct": b["syst_row_sum"],
                  "lumi_pct": b["lumi"], "stat_pct": b["stat"],
                  "total_abs": p["measurement"]["value"] * b["total"] / 100,
                  "uncorrelated": compare(p, o, correlated=False), "theory_correlated": compare(p, o, correlated=True)}
    out["reference_total_pct"] = reference_total_pct(p)
    return out


# ----------------------------------------------------------------------------- markdown
def markdown(p):
    unit = p.get("unit", "")
    lines = ["| source | " + p["reference"]["label"] + " [%] | ours [%] | status | how we do it | how they do it |",
             "|---|---:|---:|---|---|---|"]
    for r in p["rows"]:
        rp = "—" if r.get("reference_pct") is None else f"{r['reference_pct']:.2f}"
        if r.get("treatment") == "estimate":
            o = r["estimate"]["options"]
            op = f"a {o['a']:.2f} / b {o['b']:.2f} / **c {o['c']:.2f}**"
        else:
            op = f"{ours_pct(r, 'c'):.2f}"
        lines.append(f"| {r['source']} | {rp} | {op} | {STATUS_TEXT.get(r.get('status'), r.get('status'))} | "
                     f"{r.get('ours_how', '')} | {r.get('reference_how', '')} |")
    s = summary(p)
    lines += ["", f"| option | stat ⊕ syst [%] | total [%] | total [{unit}] | ours − {p['reference']['label']} [{unit}] | "
                  "pull (uncorrelated) | pull (theory correlated) |", "|---|---:|---:|---:|---:|---:|---:|"]
    for o, label in OPTIONS.items():
        d = s[o]
        lines.append(f"| {label} | {np.hypot(d['stat_pct'], d['syst_pct']):.3f} | {d['total_pct']:.2f} | {d['total_abs']:.1f} | "
                     f"{d['uncorrelated']['diff']:+.1f} | "
                     f"{d['uncorrelated']['pull']:+.2f}σ (±{d['uncorrelated']['sigma']:.1f}) | "
                     f"{d['theory_correlated']['pull']:+.2f}σ (±{d['theory_correlated']['sigma']:.1f}) |")
    return "\n".join(lines)


# ----------------------------------------------------------------------------- plots
def _style(ax):
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
    ax.tick_params(colors=INK_2)


def _default_style(fn):
    """Draw with matplotlib's defaults: channel code may have activated a global (mplhep) style."""
    import functools

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        import matplotlib.pyplot as plt
        with plt.style.context("default"):
            plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 13, "figure.facecolor": BG, "axes.facecolor": BG,
                                 "savefig.facecolor": BG, "text.color": INK, "axes.labelcolor": INK,
                                 "xtick.color": INK_2, "ytick.color": INK_2, "axes.edgecolor": INK_2,
                                 "legend.labelcolor": INK, "hatch.color": INK_2})
            return fn(*args, **kwargs)
    return wrapper


@_default_style
def plot_option(p, option, path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    colour = p.get("colour", "#2a78d6")
    rows = [r for r in p["rows"] if r.get("status") != "not-applicable" or (r.get("reference_pct") or 0) > 0]
    n = len(rows)
    fig = plt.figure(figsize=(17, 0.46 * n + 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.25, 1.0], wspace=0.08)
    ax = fig.add_subplot(gs[0])
    y = np.arange(n)[::-1]
    h = 0.36
    ref_vals = [r.get("reference_pct") or 0.0 for r in rows]
    our_vals = [ours_pct(r, option) for r in rows]
    ax.barh(y + h / 2 + 0.02, ref_vals, height=h, color=REF, label=p["reference"]["label"])
    for yi, v, r in zip(y, our_vals, rows):
        est = r.get("treatment") == "estimate"
        ax.barh(yi - h / 2 - 0.02, v, height=h, color=colour if not est else BG,
                edgecolor=colour, hatch="///" if est else None, linewidth=1.2)
    xmax = max(max(ref_vals), max(our_vals)) * 1.28 + 0.05
    for yi, rv, ov, r in zip(y, ref_vals, our_vals, rows):
        ax.text(max(rv, 0) + xmax * 0.008, yi + h / 2 + 0.02, "—" if r.get("reference_pct") is None else f"{rv:.2f}",
                va="center", fontsize=9, color=INK_2)
        tag = f"{ov:.2f}" + ("  (" + OPTIONS[option].split(":")[0] + ")" if r.get("treatment") == "estimate" else "")
        ax.text(max(ov, 0) + xmax * 0.008, yi - h / 2 - 0.02, tag, va="center", fontsize=9, color=INK)
        ax.text(xmax * 1.0, yi, STATUS_TEXT.get(r.get("status"), r.get("status")), va="center", ha="right",
                fontsize=9, color=INK_2, style="italic")
    ax.set_yticks(y)
    ax.set_yticklabels([r["source"] for r in rows], fontsize=10.5, color=INK)
    ax.set_xlim(0, xmax)
    ax.set_xlabel("uncertainty [% of the cross section]", color=INK)
    _style(ax)
    handles = [Patch(color=REF, label=p["reference"]["label"]), Patch(color=colour, label=p.get("channel_label", "ours")),
               Patch(facecolor=BG, edgecolor=colour, hatch="///", label="ours, not measurable: option " + option)]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.86, 0.0), fontsize=10, frameon=False)
    ax.set_title(f"Uncertainty budget — {OPTIONS[option]}", loc="left", fontsize=14, color=INK)

    # comparison panel
    ax2 = fig.add_subplot(gs[1])
    s = summary(p)[option]
    ref = p["reference"]
    m = p["measurement"]["value"]
    ref_tot = ref["value"] * reference_total_pct(p) / 100
    ref_nolumi = ref["value"] * np.sqrt(max(reference_total_pct(p) ** 2 - ref.get("lumi_pct", 0) ** 2, 0)) / 100
    b = budget(p, option)
    our_nolumi = m * np.hypot(b["stat"], b["syst"]) / 100
    pts = [(1.0, ref["value"], ref_nolumi, ref_tot, REF, ref["label"]),
           (0.0, m, our_nolumi, s["total_abs"], colour, p.get("channel_label", "ours"))]
    for yy, v, inner, outer, col, lab in pts:
        ax2.errorbar(v, yy, xerr=outer, fmt="none", ecolor=col, elinewidth=1.4, capsize=5)
        ax2.errorbar(v, yy, xerr=inner, fmt="o", color=col, ecolor=col, elinewidth=4.0, capsize=0, markersize=9,
                     markeredgecolor=BG, markeredgewidth=1.5)
        ax2.text(v, yy + 0.2, f"{lab}\n{v:.0f} ± {outer:.0f} {p.get('unit', '')}", ha="center", va="bottom", fontsize=10.5,
                 color=INK)
    lo = min(ref["value"] - 1.6 * ref_tot, m - 1.6 * s["total_abs"])
    hi = max(ref["value"] + 1.6 * ref_tot, m + 1.6 * s["total_abs"])
    ax2.set_xlim(lo, hi)
    ax2.set_ylim(-0.6, 1.9)
    ax2.set_yticks([])
    ax2.set_xlabel(p.get("observable", "") + f" [{p.get('unit', '')}]", color=INK)
    _style(ax2)
    ax2.spines["left"].set_visible(False)
    u, c = s["uncorrelated"], s["theory_correlated"]
    ax2.text(0.02, 0.02, f"ours − reference = {u['diff']:+.0f} {p.get('unit', '')}\n"
                         f"{u['pull']:+.2f}σ (uncorrelated), {c['pull']:+.2f}σ (shared theory correlated)\n"
                         f"thick bars: without luminosity; thin: total",
             transform=ax2.transAxes, fontsize=10, color=INK_2, va="bottom")
    ax2.set_title(f"total: ours {s['total_pct']:.2f} %, reference {reference_total_pct(p):.2f} %", loc="left",
                  fontsize=12, color=INK)
    fig.savefig(path, dpi=130, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    return path


@_default_style
def plot_overview(p, path):
    import matplotlib.pyplot as plt

    colour = p.get("colour", "#2a78d6")
    s = summary(p)
    ref = p["reference"]
    m = p["measurement"]["value"]
    fig, (ax, axt) = plt.subplots(1, 2, figsize=(15, 5.2), gridspec_kw={"width_ratios": [1.6, 1.0], "wspace": 0.05})
    ref_tot = ref["value"] * s["reference_total_pct"] / 100
    ax.axvspan(ref["value"] - ref_tot, ref["value"] + ref_tot, color=REF, alpha=0.18, lw=0)
    ax.axvline(ref["value"], color=REF, lw=1.5)
    ax.text(ref["value"], 3.55, f"{ref['label']}  {ref['value']:.0f} ± {ref_tot:.0f}", ha="center", color=INK, fontsize=11)
    labels = []
    for i, o in enumerate(("a", "b", "c")):
        yy = 2 - i
        d = s[o]
        nolumi = np.hypot(d["stat_pct"], d["syst_pct"])
        ax.errorbar(m, yy, xerr=d["total_abs"], fmt="none", ecolor=colour, elinewidth=1.4, capsize=6)
        ax.errorbar(m, yy, xerr=m * nolumi / 100, fmt="o", color=colour, ecolor=colour, elinewidth=4.0, capsize=0,
                    markersize=9, markeredgecolor=BG, markeredgewidth=1.5)
        labels.append(OPTIONS[o])
        axt.text(0.0, yy, f"± {d['total_abs']:.1f} {p.get('unit', '')} ({d['total_pct']:.2f} %)    "
                          f"{nolumi:.3f} %    "
                          f"{d['uncorrelated']['pull']:+.2f}σ / {d['theory_correlated']['pull']:+.2f}σ",
                 va="center", fontsize=11.5, color=INK)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels(labels, fontsize=11.5, color=INK)
    ax.set_ylim(-0.7, 3.9)
    ax.set_xlabel(p.get("observable", "") + f" [{p.get('unit', '')}]", color=INK)
    _style(ax)
    axt.set_ylim(-0.7, 3.9)
    axt.axis("off")
    axt.text(0.0, 3.2, "total                          without lumi    pull: uncorrelated / theory correlated",
             fontsize=10.5, color=INK_2)
    ax.text(0.01, 0.02, "thick bars: without luminosity; thin: total", transform=ax.transAxes, fontsize=10, color=INK_2)
    est = [r["source"] for r in p["rows"] if r.get("treatment") == "estimate"]
    fig.suptitle(f"{p.get('channel_label', '')} vs {ref['label']}: the options for the sources we cannot measure "
                 f"({', '.join(est)})", fontsize=13.5, color=INK, x=0.02, ha="left")
    fig.savefig(path, dpi=130, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parity", type=Path)
    ap.add_argument("--plot-dir", type=Path, required=True)
    ap.add_argument("--prefix", default="parity")
    args = ap.parse_args()
    p = load(args.parity)
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    for o in OPTIONS:
        print("wrote", plot_option(p, o, args.plot_dir / f"{args.prefix}_option_{o}.png"))
    print("wrote", plot_overview(p, args.plot_dir / f"{args.prefix}_options.png"))
    md = markdown(p)
    out_md = args.parity.with_suffix(".md")
    out_md.write_text(md + "\n")
    summ = summary(p)
    args.parity.with_name(args.parity.stem + "_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    print(md)


if __name__ == "__main__":
    main()
