#!/usr/bin/env python
"""Step 3 -- measure the fake factors, the closure corrections, the OS/SS correction and the closure
(docs/05-fake-factors.md).

    python scripts/step3_fakefactors.py

Reads the ntuples; writes output/data/fakefactors.json (variant "mcsub": simulated events with a genuine
leading tau subtracted from every FF region; --with-nosub adds the unsubtracted cross-check) with the
era x DM x N_jets x pT table, the factorised closure corrections in |eta(tau1)| and pT(tau2),
the OS/SS correction C, the same-sign closure in several variables and the genuine-tau contamination of
every region; and the plots output/plots/step3_*.png.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import analysis, config, fakes, plotting  # noqa: E402

VARS = {"m_tt": (config.FIT_BINS, r"$m_{\tau\tau}$ (MET likelihood) [GeV]"),
        "m_vis": (np.arange(0, 305, 15), r"$m_{vis}$ [GeV]"),
        "t1_pt": (np.array([40, 45, 50, 55, 60, 70, 80, 100, 150, 250]), r"$p_T(\tau_1)$ [GeV]"),
        "t2_pt": (np.array([40, 45, 50, 55, 60, 70, 80, 100, 150]), r"$p_T(\tau_2)$ [GeV]"),
        "t1_eta": (np.linspace(-2.1, 2.1, 15), r"$\eta(\tau_1)$"),
        "t2_eta": (np.linspace(-2.1, 2.1, 15), r"$\eta(\tau_2)$"),
        "njets": (np.arange(-0.5, 5.5, 1), r"$N_{jets}$"),
        "dr_tt": (np.linspace(0.5, 5.0, 19), r"$\Delta R(\tau_1, \tau_2)$"),
        "met": (np.arange(0, 155, 10), r"$p_T^{miss}$ [GeV]")}
NOMINAL = "mcsub" if config.FF_SUBTRACT_MC else "nosub"


def hist(x, w, edges):
    edges = np.asarray(edges, dtype=float)
    xc = np.clip(x, edges[0], edges[-1] - 1e-6)            # overflow into the last bin
    v, _ = np.histogram(xc, bins=edges, weights=w)
    v2, _ = np.histogram(xc, bins=edges, weights=w ** 2)
    return v, v2


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--with-nosub", action="store_true", help="also compute the unsubtracted variant")
    args = ap.parse_args()
    config.PLOT_DIR.mkdir(parents=True, exist_ok=True)
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    data = analysis.load_data()
    reg = analysis.regions(data)
    print({k: int(v.sum()) for k, v in reg.items()})

    # simulation with a genuine leading tau (the part the FF must not contain)
    mc = []
    for key in analysis.available_mc():
        d, _ = analysis.load(key)
        if not len(d):
            continue
        r = analysis.regions(d, is_mc=True)
        w = analysis.weights(d, key)
        mc.append((key, d, r, w))
    # subtracted from the FF regions (W+jets with uniform weights, analysis.subtraction_weights)
    mc_sub = [(k, d, r, analysis.subtraction_weights(k, w)) for k, d, r, w in mc if k not in analysis.SUBTRACT_EXCLUDE]
    contamination = {}
    for name in analysis.REGIONS:
        mcsum = sum(float(w[r[name]].sum()) for _, _, r, w in mc)
        contamination[name] = {"data": int(reg[name].sum()), "mc_genuine_tau1": mcsum,
                               "fraction": mcsum / max(int(reg[name].sum()), 1)}
        print(f"  {name:7s} data {int(reg[name].sum()):>9,d}   MC (genuine tau1) {mcsum:>10.1f}  "
              f"({100 * contamination[name]['fraction']:.2f}%)")

    results = {"contamination": contamination, "nominal": NOMINAL}
    variants = [("mcsub", True)] + ([("nosub", False)] if args.with_nosub else [])
    for variant, subtract in variants:
        sub_ff = [(d, r["SS_T"], r["SS_L"], w) for _, d, r, w in mc_sub] if subtract else []
        sub_reg = [(d, r, w) for _, d, r, w in mc_sub] if subtract else []
        table = fakes.measure(data, reg["SS_T"], reg["SS_L"], sub_ff)
        # closure plots *before* the corrections, for the record
        wraw = fakes.fake_weights(data, reg["SS_L"], {k: v for k, v in table.items() if k != "closure"}, 1.0)
        fakes.closure_corrections(data, reg, table, sub_reg)
        osss = fakes.osss_correction(data, reg, sub_reg)
        wpred = fakes.fake_weights(data, reg["SS_L"], table, 1.0)
        closure = {}
        for var, (edges, xl) in VARS.items():
            obs, _ = hist(data[var][reg["SS_T"]], np.ones(int(reg["SS_T"].sum())), edges)
            obs_mc = np.zeros_like(obs, dtype=float)
            preds = {}
            for tag, wp in (("", wpred), ("_nocorr", wraw)):
                pred, pred2 = hist(data[var][reg["SS_L"]], wp[reg["SS_L"]], edges)
                pred_mc = np.zeros_like(obs, dtype=float)
                if subtract:
                    for _, d, r, w in mc_sub:
                        if not tag:
                            obs_mc += hist(d[var][r["SS_T"]], w[r["SS_T"]], edges)[0]
                        tbl = table if not tag else {k: v for k, v in table.items() if k != "closure"}
                        wm = fakes.fake_weights(d, r["SS_L"], tbl, 1.0) * w
                        pred_mc += hist(d[var][r["SS_L"]], wm[r["SS_L"]], edges)[0]
                preds[tag] = (pred - pred_mc, pred2)
            closure[var] = {"edges": list(map(float, edges)), "obs": obs.tolist(), "pred": preds[""][0].tolist(),
                            "pred_var": preds[""][1].tolist(), "obs_mc": obs_mc.tolist(),
                            "pred_nocorr": preds["_nocorr"][0].tolist()}
            if variant == NOMINAL or var in ("m_tt", "t1_eta"):
                for tag in ("", "_nocorr") if var in ("t1_eta", "t2_pt", "m_tt") else ("",):
                    stack = [("Fakes", preds[tag][0], preds[tag][1])]
                    if subtract:
                        stack = [("DYtautau", obs_mc, np.zeros_like(obs_mc))] + stack
                    plotting.stack_plot(config.PLOT_DIR / f"step3_closure_SS_{var}_{variant}{tag}.png", edges,
                                        (obs, np.sqrt(obs)), stack, xl,
                                        title=f"same-sign closure ({'MC subtracted' if subtract else 'no MC subtraction'}"
                                              f"{', before closure corrections' if tag else ''})")
        # OS tau2-anti-isolated sideband: prediction with C
        tai = osss["table_ai"]
        m = reg["OSAI_L"]
        obs, _ = hist(data["m_tt"][reg["OSAI_T"]], np.ones(int(reg["OSAI_T"].sum())), config.FIT_BINS)
        w_ai = fakes.fake_weights(data, m, tai, osss["C"])
        pred, pred2 = hist(data["m_tt"][m], w_ai[m], config.FIT_BINS)
        stack = [("Fakes", pred, pred2)]
        if subtract:
            sig = np.zeros(len(obs))
            for _, d, r, w in mc_sub:
                sig += hist(d["m_tt"][r["OSAI_T"]], w[r["OSAI_T"]], config.FIT_BINS)[0]
                wm = fakes.fake_weights(d, r["OSAI_L"], tai, osss["C"]) * w
                stack[0] = ("Fakes", stack[0][1] - hist(d["m_tt"][r["OSAI_L"]], wm[r["OSAI_L"]], config.FIT_BINS)[0], pred2)
            stack = [("DYtautau", sig, np.zeros_like(sig))] + stack
        plotting.stack_plot(config.PLOT_DIR / f"step3_OSAI_m_tt_{variant}.png", config.FIT_BINS, (obs, np.sqrt(obs)),
                            stack, VARS["m_tt"][1], title=r"OS, $\tau_2$ anti-isolated" + f" (C = {osss['inclusive']['C']:.3f})")
        # residual non-closure in the fit variable (same sign, inclusive; the per-category version is in step 4)
        c = closure["m_tt"]
        o, p = np.asarray(c["obs"], float) - np.asarray(c["obs_mc"]), np.asarray(c["pred"], float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(p > 0, o / p, 1.0)
            rerr = np.where((p > 0) & (o > 0), ratio * np.sqrt(1 / np.maximum(np.asarray(c["obs"]), 1) + np.asarray(c["pred_var"]) / p ** 2), 1.0)
        results[variant] = {"ff": fakes.to_json(table), "osss": {k: v for k, v in osss.items() if k != "table_ai"},
                            "ff_ai": fakes.to_json(tai), "closure": closure,
                            "nonclosure_m_tt": {"ratio": ratio.tolist(), "err": rerr.tolist()}}
        print(f"[{variant}] C_OS/SS inclusive = {osss['inclusive']['C']:.4f} +- {osss['inclusive']['stat']:.4f}")
        for var, corr in table["closure"].items():
            print(f"   closure correction {var}: " + " ".join(f"{v:.3f}" for v in corr["values"]))
        for e, era in enumerate(table["eras"]):
            print(f"   era {era}: C per N_jets " + ", ".join(f"{c:.4f} +- {x:.4f}" for c, x in zip(osss["C"][e], osss["stat"][e])))
            for i, dm in enumerate(fakes.DMS):
                for j, nj in enumerate(fakes.NJ_BINS):
                    print(f"     DM{dm:>2d} nj{nj} FF: " + " ".join(f"{v:.3f}" for v in table["ff"][e, i, j])
                          + "   N(num): " + " ".join(f"{v:.0f}" for v in table["num"][e, i, j]))
        for e, era in enumerate(table["eras"]):
            me = data["era"] == e
            wpe = fakes.fake_weights(data, reg["SS_L"] & me, table, 1.0)
            print(f"   era {era}: SS closure obs/pred = {(reg['SS_T'] & me).sum() / wpe.sum():.4f}")

    # FF plot: rows = eras, columns = DMs, one series per jet multiplicity; open markers = no subtraction
    eras = results[NOMINAL]["ff"]["eras"]
    fig, axes = plt.subplots(len(eras), 4, figsize=(22, 6.2 * len(eras)), squeeze=False)
    fig.subplots_adjust(hspace=0.38, wspace=0.22)
    edges = np.asarray(config.FF_PT_BINS[:-1] + [120.0])
    centers = 0.5 * (edges[1:] + edges[:-1])
    colors = ["black", "tab:red", "tab:blue"]
    for e, era in enumerate(eras):
        for i, dm in enumerate(fakes.DMS):
            ax = axes[e, i]
            for j, nj in enumerate(fakes.NJ_BINS):
                lab = f"{nj} jets" if j < len(fakes.NJ_BINS) - 1 else rf"$\geq${nj} jets"
                for variant, mfc, off in [v for v in (("mcsub", colors[j], -1.0), ("nosub", "none", 1.0)) if v[0] in results]:
                    t = results[variant]["ff"]
                    ax.errorbar(centers + off + 0.6 * (j - 1), np.asarray(t["ff"])[e, i, j],
                                yerr=np.asarray(t["err"])[e, i, j], fmt="o", color=colors[j], mfc=mfc,
                                label=lab if variant == "mcsub" else None)
            ax.set_title(f"Run2016{era}, decay mode {dm}", fontsize=16)
            ax.set_xlim(40, 120)
            ax.set_ylim(0, None)
            if i == 0:
                ax.set_ylabel("FF = N(Medium) / N(VVVLoose & !Medium)", fontsize=13)
            ax.set_xlabel(r"$p_T(\tau_1)$ [GeV] (last bin: > 80)", fontsize=13)
            if e == 0 and i == 0:
                ax.legend(fontsize=12, title="filled: MC subtracted (nominal), open: no subtraction", title_fontsize=11)
    fig.savefig(config.PLOT_DIR / "step3_fakefactors.png", bbox_inches="tight", dpi=90)
    fig.savefig(config.PLOT_DIR / "step3_fakefactors.pdf", bbox_inches="tight")
    plt.close(fig)
    # closure-correction plot
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, (var, xl) in zip(axes, (("eta", r"$|\eta(\tau_1)|$"), ("pt2", r"$p_T(\tau_2)$ [GeV] (last bin: > 80)"))):
        for variant, mfc in [v for v in (("mcsub", "black"), ("nosub", "none")) if v[0] in results]:
            c = results[variant]["ff"]["closure"][var]
            ed = np.asarray(c["edges"]); ed[-1] = min(ed[-1], 120.0)
            ax.errorbar(0.5 * (ed[1:] + ed[:-1]), c["values"], yerr=c["err"], xerr=np.diff(ed) / 2, fmt="o", color="black",
                        mfc=mfc, label="MC subtracted (nominal)" if mfc == "black" else "no subtraction")
        ax.axhline(1, color="grey", lw=0.8)
        ax.set_xlabel(xl); ax.set_ylabel("same-sign obs / FF prediction"); ax.set_ylim(0.8, 1.2)
        ax.legend(fontsize=9)
    fig.suptitle("factorised closure corrections of the fake factor")
    fig.tight_layout()
    fig.savefig(config.PLOT_DIR / "step3_closure_corrections.png", bbox_inches="tight", dpi=110)
    fig.savefig(config.PLOT_DIR / "step3_closure_corrections.pdf", bbox_inches="tight")
    plt.close(fig)
    (config.DATA_DIR / "fakefactors.json").write_text(json.dumps(results, indent=1))
    print(f"-> {config.DATA_DIR / 'fakefactors.json'}")


if __name__ == "__main__":
    main()
