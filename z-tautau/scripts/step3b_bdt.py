#!/usr/bin/env python
"""Step 3b -- train the k-fold BDT and validate the fake factor in its score (docs/09-bdt.md).

    python scripts/step3b_bdt.py [--no-train]

Signal: simulated Z -> tautau (fiducial part) in the SR of every available Drell-Yan sample; background:
application-region data weighted by the nominal fake factor (step 3). Five models (fold = event mod 5),
saved to $BND_TAUTAU_CACHE/bdt/. Every event downstream is scored by the model that never saw it.

Validation written to output/data/bdt.json and output/plots/step3b_*.png:
  * held-out AUC per fold, feature importance, score shapes (signal / non-fiducial DY / fakes);
  * same-sign closure of the FF in the score (SS_T vs FF x SS_L, MC subtracted): the FF must not depend
    on the score beyond the residual that the per-category closure nuisance parameters cover;
  * prefit signal-region score distribution (data vs prediction) and the category yields;
  * stat-only sensitivity of an m_tt fit inclusive vs in categories.
`fit/bdt_info.json` (committed) records the training summary.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import analysis, bdt, config, fakes, plotting, samples  # noqa: E402

NOMINAL = "mcsub" if config.FF_SUBTRACT_MC else "nosub"
SCORE_EDGES = np.linspace(0, 1, 21)


def h(x, w, edges):
    edges = np.asarray(edges, dtype=float)
    xc = np.clip(x, edges[0], edges[-1] - 1e-6)
    return np.histogram(xc, bins=edges, weights=w)[0], np.histogram(xc, bins=edges, weights=w ** 2)[0]


def sens(s, b):
    ok = (s + b) > 0
    return float(1.0 / np.sqrt(np.sum(s[ok] ** 2 / (s[ok] + b[ok]))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-train", action="store_true", help="only the validation, with the existing models")
    args = ap.parse_args()
    config.PLOT_DIR.mkdir(parents=True, exist_ok=True)
    ffres = json.loads((config.DATA_DIR / "fakefactors.json").read_text())
    table = fakes.from_json(ffres[NOMINAL]["ff"])
    C = np.asarray(ffres[NOMINAL]["osss"]["C"])
    data = analysis.load_data()
    reg = analysis.regions(data)
    ar = reg["AR"]
    wf = fakes.fake_weights(data, ar, table, C)

    mc = []
    for key in analysis.available_mc():
        d, _ = analysis.load(key)
        if len(d):
            mc.append((key, d, analysis.regions(d, is_mc=True), analysis.weights(d, key)))
    mc_sub = [(k, d, r, analysis.subtraction_weights(k, w)) for k, d, r, w in mc if k not in analysis.SUBTRACT_EXCLUDE]
    dy_keys = analysis.dy_stitch_keys()
    Xs, ws, fs = [], [], []
    for key, d, r, w in mc:
        if key not in dy_keys:
            continue
        m = r["SR"] & (d["gen_lhe_flavour"] == 15)
        if config.BDT_TRAIN_ON_FIDUCIAL:
            m &= d["gen_fid"].astype(bool)
        Xs.append(bdt.features(d)[m]); ws.append(w[m]); fs.append(bdt.folds(d)[m])
    Xs, ws, fs = np.concatenate(Xs), np.concatenate(ws), np.concatenate(fs)
    Xb, wb, fb = bdt.features(data)[ar], wf[ar], bdt.folds(data)[ar]
    print(f"signal: {len(Xs):,} events from {dy_keys} (sum w {ws.sum():.0f}); fakes: {len(Xb):,} AR events (sum w {wb.sum():.0f})")
    if not args.no_train:
        info = bdt.train(Xs, ws, fs, Xb, wb, fb)
        print(f"AUC per fold (held out): {np.round(info['auc_test'], 3).tolist()}  (train: {np.round(info['auc_train'], 3).tolist()})")
        print("importance: " + ", ".join(f"{k} {v:.2f}" for k, v in info["importance"].items()))
    else:
        info = json.loads((config.BDT_DIR / "bdt.json").read_text())

    # ---------------------------------------------------------------- scores of everything
    s_data = analysis.scores(data, "data")
    cat_data = bdt.category(s_data)
    cats = list(range(len(config.BDT_CATEGORY_EDGES) - 1))
    # signal region: data vs fakes (MC subtracted in the AR) + simulation
    sr = reg["SR"]
    obs, _ = h(s_data[sr], np.ones(sr.sum()), SCORE_EDGES)
    fake, fake2 = h(s_data[ar], wf[ar], SCORE_EDGES)
    stacks = {}
    edges = np.asarray(config.FIT_BINS)
    yields = {c: {"data": int((sr & (cat_data == c)).sum()), "fakes": float(wf[ar & (cat_data == c)].sum())} for c in cats}
    mtt_s = {c: np.zeros(len(edges) - 1) for c in cats}
    mtt_b = {c: h(data["m_tt"][ar & (cat_data == c)], wf[ar & (cat_data == c)], edges)[0] for c in cats}
    shapes = {"DYtautau": [], "DYtautau_nonfid": [], "Fakes": (s_data[ar], wf[ar])}
    for key, d, r, w in mc:
        s_mc = analysis.scores(d, key)
        c_mc = bdt.category(s_mc)
        if NOMINAL == "mcsub" and key not in analysis.SUBTRACT_EXCLUDE:
            m = r["AR"]
            wm = fakes.fake_weights(d, m, table, C) * analysis.subtraction_weights(key, w)
            f_, f2_ = h(s_mc[m], wm[m], SCORE_EDGES)
            fake -= f_
            for c in cats:
                yields[c]["fakes"] -= float(wm[m & (c_mc == c)].sum())
                mtt_b[c] -= h(d["m_tt"][m & (c_mc == c)], wm[m & (c_mc == c)], edges)[0]
        for name, cm in analysis.mc_components(key):
            m = r["SR"] & cm
            v, v2 = h(s_mc[m], w[m], SCORE_EDGES)
            stacks[name] = stacks.get(name, np.zeros((2, len(v)))) + np.stack([v, v2])
            for c in cats:
                yields[c][name] = yields[c].get(name, 0.0) + float(w[m & (c_mc == c)].sum())
                if name == samples.SIGNAL:
                    mtt_s[c] += h(d["m_tt"][m & (c_mc == c)], w[m & (c_mc == c)], edges)[0]
            if name in ("DYtautau", "DYtautau_nonfid"):
                shapes[name].append((s_mc[m], w[m]))
    groups = [("DYll", ["DYee", "DYmumu"]), ("DYlowmass", ["DYlowmass"]), ("Diboson", ["WW", "WZ", "ZZ"]),
              ("Top", ["TTbar", "SingleTop"]), ("WJets", ["WJets"]), ("Fakes", None), ("DYtautau_nonfid", ["DYtautau_nonfid"]),
              ("DYtautau", ["DYtautau"])]
    stack = []
    for g, members in groups:
        if members is None:
            stack.append(("Fakes", np.maximum(fake, 0), fake2))
        else:
            parts = [stacks[m] for m in members if m in stacks]
            if parts:
                tot = sum(parts)
                stack.append((g, np.maximum(tot[0], 0), tot[1]))
    for logy in (False, True):
        plotting.stack_plot(config.PLOT_DIR / f"step3b_bdt_SR_score{'_log' if logy else ''}.png", SCORE_EDGES,
                            (obs, np.sqrt(obs)), stack, "BDT score", title="signal region (prefit)", logy=logy)
    # same-sign closure in the score
    wss = fakes.fake_weights(data, reg["SS_L"], table, 1.0)
    o_ss, _ = h(s_data[reg["SS_T"]], np.ones(reg["SS_T"].sum()), SCORE_EDGES)
    p_ss, p2_ss = h(s_data[reg["SS_L"]], wss[reg["SS_L"]], SCORE_EDGES)
    o_ss = o_ss.astype(float)
    sig_ss = np.zeros(len(o_ss))
    if NOMINAL == "mcsub":
        for key, d, r, w in mc_sub:
            s_mc = analysis.scores(d, key)
            sig_ss += h(s_mc[r["SS_T"]], w[r["SS_T"]], SCORE_EDGES)[0]
            wm = fakes.fake_weights(d, r["SS_L"], table, 1.0) * w
            p_ss -= h(s_mc[r["SS_L"]], wm[r["SS_L"]], SCORE_EDGES)[0]
    plotting.stack_plot(config.PLOT_DIR / "step3b_bdt_closure_SS.png", SCORE_EDGES, (o_ss, np.sqrt(o_ss)),
                        [("MC", sig_ss, np.zeros_like(sig_ss)), ("Fakes", np.maximum(p_ss, 0), p2_ss)],
                        "BDT score", title="same-sign closure of the fake factor in the score", logy=True)
    # score shapes
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for name, col, lab in (("DYtautau", "#d9822b", r"Z$\rightarrow\tau\tau$, fiducial"),
                           ("DYtautau_nonfid", "#8a6200", r"Z/$\gamma^*\rightarrow\tau\tau$, non-fiducial"),
                           ("Fakes", "#e377ac", "fakes (AR x FF)")):
        if name == "Fakes":
            x, w = shapes["Fakes"]
        else:
            x = np.concatenate([a for a, _ in shapes[name]]); w = np.concatenate([b for _, b in shapes[name]])
        ax.hist(x, bins=SCORE_EDGES, weights=w / w.sum(), histtype="step", lw=2, color=col, label=lab)
    for e in config.BDT_CATEGORY_EDGES[1:-1]:
        ax.axvline(e, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("BDT score (held-out fold)"); ax.set_ylabel("normalised"); ax.legend(fontsize=9, frameon=False)
    ax.set_title(f"k-fold BDT: held-out AUC {np.mean(info['auc_test']):.3f}", fontsize=11)
    fig.savefig(config.PLOT_DIR / "step3b_bdt_score_shapes.png", bbox_inches="tight", dpi=120)
    fig.savefig(config.PLOT_DIR / "step3b_bdt_score_shapes.pdf", bbox_inches="tight")
    plt.close(fig)
    # importance plot
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    items = list(info["importance"].items())[::-1]
    ax.barh([k for k, _ in items], [v for _, v in items], color="#4a90d9")
    ax.set_xlabel("feature importance (gain, mean over folds)")
    fig.savefig(config.PLOT_DIR / "step3b_bdt_importance.png", bbox_inches="tight", dpi=120)
    fig.savefig(config.PLOT_DIR / "step3b_bdt_importance.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---------------------------------------------------------------- category table and sensitivity
    tot_s = sum(mtt_s.values()); tot_b = sum(mtt_b.values())
    incl = sens(tot_s, tot_b)
    cat_sens = 1.0 / np.sqrt(sum(1.0 / sens(mtt_s[c], mtt_b[c]) ** 2 for c in cats))
    win = (edges[:-1] >= 70) & (edges[:-1] < 110)
    print(f"{'cat':>4s} {'edges':>12s} {'data':>8s} {'fakes':>8s} {'Z->tautau':>9s} {'non-fid':>8s} {'S/B':>6s} {'S/B 70-110':>10s}")
    for c in cats:
        y = yields[c]
        print(f"{c:4d} {config.BDT_CATEGORY_EDGES[c]:5.2f}-{config.BDT_CATEGORY_EDGES[c+1]:5.2f} {y['data']:8d} {y['fakes']:8.0f} "
              f"{y.get('DYtautau', 0):9.0f} {y.get('DYtautau_nonfid', 0):8.0f} {y.get('DYtautau', 0) / max(y['fakes'], 1):6.2f} "
              f"{mtt_s[c][win].sum() / max(mtt_b[c][win].sum(), 1):10.2f}")
    print(f"stat-only d(mu)/mu, fakes fixed: inclusive {100 * incl:.2f}% -> categories {100 * cat_sens:.2f}%")
    out = {"training": info, "category_edges": config.BDT_CATEGORY_EDGES, "category_yields_prefit": yields,
           "stat_only_sensitivity": {"inclusive": incl, "categories": cat_sens},
           "ss_closure_score": {"edges": SCORE_EDGES.tolist(), "obs": o_ss.tolist(), "obs_mc": sig_ss.tolist(), "pred": p_ss.tolist(),
                                "pred_var": p2_ss.tolist()},
           "sr_score": {"edges": SCORE_EDGES.tolist(), "data": obs.tolist(), "fakes": fake.tolist(),
                        "mc": {k: v[0].tolist() for k, v in stacks.items()}}}
    (config.DATA_DIR / "bdt.json").write_text(json.dumps(out, indent=1, default=float))
    config.FIT_DIR.mkdir(parents=True, exist_ok=True)
    (config.FIT_DIR / "bdt_info.json").write_text(json.dumps({"training": info, "category_edges": config.BDT_CATEGORY_EDGES,
                                                                "category_yields_prefit": yields,
                                                                "stat_only_sensitivity": out["stat_only_sensitivity"]},
                                                               indent=1, default=float))
    print(f"-> {config.DATA_DIR / 'bdt.json'}, {config.FIT_DIR / 'bdt_info.json'}, models in {config.BDT_DIR}")


if __name__ == "__main__":
    main()
