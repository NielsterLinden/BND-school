#!/usr/bin/env python
"""Step 3 (v4) -- fake factors of the mu tau_h and e tau_h channels and the e mu multijet estimate
(docs/10-v4-plan.md section 5).

    python scripts/step3d_fakes_lepton.py [--channels mutau etau emu] [--no-plots]

Per lepton channel: FF_qcd (same-sign DR), FF_w (m_T > 70 DR, no b jet), FF_tt (ttbar/single-top simulation),
the AR fractions R_p(N_jets, m_T), C_OS/SS from the anti-isolated-lepton sidebands, r_W(DM) from the W+jets
simulation, and the same-sign validation (obs/pred -> FakeClosure_<ch>). e mu: the OS/SS ratio of the multijet
background from the isolation sidebands SB1 (nominal) and SB2 (systematic) in bins of dR(e, mu).
Output: output_v4/data/fakes_<channel>.json, output_v4/plots/step3_<channel>_*.png.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import analysis_v4 as an, config, fakes_v4 as fk, plotting, samples  # noqa: E402

W_LIKE = {"WJets", "DYee", "DYmumu", "DYtautau", "DYtautau_out", "DYlowmass", "WW", "WZ", "ZZ"}     # jet fakes use the W FF
TT_LIKE = {"TTbar", "SingleTop"}
DR_BINS = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 6.0])


def load_channel(channel):
    data = an.load_data(channel)
    reg_d = an.regions(data, channel)
    mc = []
    for key in an.available_mc(channel):
        d, _ = an.load(key, channel)
        d = dict(d)
        w = an.weights(d, key, channel)
        r_all = an.regions(d, channel, is_mc=False)
        r_gen = an.regions(d, channel, is_mc=True)
        for name, cm in an.mc_components(key, channel):
            mc.append((key, name, d, cm, w, r_all, r_gen))
    return data, reg_d, mc


def ltau(channel, plots: bool):
    data, reg, mc = load_channel(channel)
    jet = lambda d: d["t_genflav"] == 0        # noqa: E731
    # genuine / lepton-faked tau_h simulation, subtracted from every data region
    # the multijet FF is the FF of *pure* multijet events: genuine / lepton-faked tau_h of every sample and the
    # jet fakes of W+jets, Drell-Yan (W-like) and top are subtracted from the same-sign DR (the same-sign region
    # with an isolated lepton is about half W+jets; without this the W part was double counted, 16% non-closure)
    sub_T = [(d, r_gen["QCD_T"] & cm, r_gen["QCD_L"] & cm, w) for _, _, d, cm, w, _, r_gen in mc]
    sub_T += [(d, r_all["QCD_T"] & cm & jet(d), r_all["QCD_L"] & cm & jet(d), w) for _, name, d, cm, w, r_all, _ in mc if name in W_LIKE | TT_LIKE]
    sub_W = [(d, r_gen["W_T"] & cm, r_gen["W_L"] & cm, w) for _, _, d, cm, w, _, r_gen in mc]
    sub_WSS = [(d, r_gen["WSS_T"] & cm, r_gen["WSS_L"] & cm, w) for _, _, d, cm, w, _, r_gen in mc]
    tables = {"qcd": fk.measure(data, reg["QCD_T"], reg["QCD_L"], sub_T),
              "w": fk.measure(data, reg["W_T"], reg["W_L"], sub_W),
              "w_ss": fk.measure(data, reg["WSS_T"], reg["WSS_L"], sub_WSS)}
    tt_parts = [(d, r_all["SR"] & cm & jet(d), r_all["AR"] & cm & jet(d), w) for _, name, d, cm, w, r_all, _ in mc if name in TT_LIKE]
    tables["tt"] = fk.coarsen_tt(fk.measure_mc(tt_parts))
    # fractions in the AR
    parts = []
    for _, name, d, cm, w, r_all, _ in mc:
        proc = "w" if name in W_LIKE else ("tt" if name in TT_LIKE else "other")
        parts.append((d, r_all["AR"] & cm, jet(d), w, proc))
    fr = fk.fractions(data, reg["AR"], parts)
    # corrections
    osss_sub = [(d, {k: v & cm for k, v in r_gen.items()}, w) for _, _, d, cm, w, _, r_gen in mc]
    osss_sub += [(d, {k: v & cm & jet(d) for k, v in r_all.items()}, w) for _, name, d, cm, w, r_all, _ in mc if name in W_LIKE | TT_LIKE]
    osss = fk.osss_correction(data, reg, osss_sub)
    rw = fk.w_mt_correction([(d, {k: v & cm for k, v in r_all.items()}, jet(d), w) for _, name, d, cm, w, r_all, _ in mc if name == "WJets"])
    # same-sign validation: SS_SR (data - genuine MC) vs SS_AR x w (multijet FF without C, fractions of the SS AR)
    fr_ss = fk.fractions(data, reg["SS_AR"], [(d, r_all["SS_AR"] & cm, jet(d), w, ("w" if n in W_LIKE else "tt" if n in TT_LIKE else "other")) for _, n, d, cm, w, r_all, _ in mc])
    osss_one = {**osss, "C": 1.0}
    tables_ss = {**tables, "w": tables["w_ss"]}          # same-sign W FF for the same-sign prediction
    wss, wss2 = fk.fake_weights(data, reg["SS_AR"], tables_ss, fr_ss, osss_one, rw, with_err=True)
    edges = np.asarray(config.FIT_BINS_LTAU)
    obs = np.histogram(np.clip(data["m_tt"][reg["SS_SR"]], 0, edges[-1] - 1e-6), bins=edges)[0].astype(float)
    pred = np.histogram(np.clip(data["m_tt"][reg["SS_AR"]], 0, edges[-1] - 1e-6), bins=edges, weights=wss[reg["SS_AR"]])[0]
    pred2 = np.histogram(np.clip(data["m_tt"][reg["SS_AR"]], 0, edges[-1] - 1e-6), bins=edges, weights=wss2[reg["SS_AR"]] ** 2)[0]
    obs_mc, pred_mc = np.zeros(len(edges) - 1), np.zeros(len(edges) - 1)
    for _, _, d, cm, w, _, r_gen in mc:
        m = r_gen["SS_SR"] & cm
        obs_mc += np.histogram(np.clip(d["m_tt"][m], 0, edges[-1] - 1e-6), bins=edges, weights=w[m])[0]
        m = r_gen["SS_AR"] & cm
        wm = fk.fake_weights(d, m, tables_ss, fr_ss, osss_one, rw) * w
        pred_mc += np.histogram(np.clip(d["m_tt"][m], 0, edges[-1] - 1e-6), bins=edges, weights=wm[m])[0]
    o, p = (obs - obs_mc).sum(), (pred - pred_mc).sum()
    ratio = o / p if p > 0 else 1.0
    stat = ratio * np.sqrt(1 / max(o, 1) + pred2.sum() / max(p, 1e-9) ** 2)
    closure = {"ratio": float(ratio), "stat": float(stat), "delta": float(np.hypot(ratio - 1, stat)), "obs": float(o), "pred": float(p),
               "obs_hist": (obs - obs_mc).tolist(), "pred_hist": (pred - pred_mc).tolist(), "edges": edges.tolist()}
    # AR composition and the SR estimate
    wf, wf2 = fk.fake_weights(data, reg["AR"], tables, fr, osss, rw, with_err=True)
    n_fake = float(wf[reg["AR"]].sum())
    n_sub = sum(float((fk.fake_weights(d, r_gen["AR"] & cm, tables, fr, osss, rw) * w)[r_gen["AR"] & cm].sum()) for _, _, d, cm, w, _, r_gen in mc)
    out = {"channel": channel, "tables": fk.to_json(tables), "fractions": fr, "osss": osss, "w_mt": rw, "closure_ss": closure,
           "sr_fakes": n_fake - n_sub, "sr_fakes_before_subtraction": n_fake, "sr_subtracted_mc": n_sub,
           "n_data": {k: int(v.sum()) for k, v in reg.items()}}
    print(f"[{channel}] FF_qcd (DM1, 0 jet): {np.round(tables['qcd']['ff'][1, 0], 3)}  FF_w: {np.round(tables['w']['ff'][1, 0], 3)}  "
          f"FF_w(SS): {np.round(tables['w_ss']['ff'][1, 0], 3)}  FF_tt 1-prong: {np.round(tables['tt']['ff'][1, 0], 3)}")
    print(f"[{channel}] C_OS/SS = {osss['C']:.3f} +- {osss['stat']:.3f};  r_W = " + ", ".join(f"{g}: {v['r']:.3f} +- {v['rel_stat']:.2f}" for g, v in rw.items()))
    print(f"[{channel}] AR fractions (0 jet, m_T bins): qcd {np.round(fr['fractions']['qcd'][0], 2)} w {np.round(fr['fractions']['w'][0], 2)} tt {np.round(fr['fractions']['tt'][0], 2)}")
    print(f"[{channel}] SS closure obs/pred = {ratio:.3f} +- {stat:.3f} -> NP {100 * closure['delta']:.1f}%;  SR fakes {n_fake - n_sub:.0f} (subtracted MC {n_sub:.0f}); data SR {reg['SR'].sum()}")
    (config.DATA_DIR_V4 / f"fakes_{channel}.json").write_text(json.dumps(fk.to_json(out), indent=1, default=float))
    if plots:
        plot_ltau(channel, out)
    return out


def plot_ltau(channel, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pt = np.asarray(out["tables"]["qcd"]["pt_bins"])
    x = 0.5 * (pt[1:] + pt[:-1]); x[-1] = pt[-2] + 15
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for ax, dm in zip(axes.flat, config.TAU_DMS):
        i = list(config.TAU_DMS).index(dm)
        for proc, col in (("qcd", "C0"), ("w", "C3"), ("w_ss", "C1"), ("tt", "C2")):
            t = out["tables"][proc]
            for j, nj in enumerate(fk.NJ_BINS):
                ax.errorbar(x + 1.5 * j, np.asarray(t["ff"])[i, j], np.asarray(t["err"])[i, j], fmt="o" if j == 0 else ("s" if j == 1 else "^"),
                            color=col, ms=4, alpha=0.8, label=f"{proc} N_jets{'=' if j < 2 else '>='}{nj}" if dm == 0 else None)
        ax.set_title(f"DM {dm}"); ax.set_ylim(0, 1.0); ax.set_xlim(30, 120); ax.grid(alpha=0.3)
        ax.set_xlabel(r"$p_T(\tau_h)$ [GeV]"); ax.set_ylabel("FF")
    axes[0, 0].legend(fontsize=7, ncol=3)
    plotting.tag(axes[0, 0])
    fig.suptitle(f"{channel}: fake factors (tight / loose-not-tight) per process")
    fig.tight_layout(); plotting.fig_tag(fig)
    fig.savefig(config.PLOT_DIR_V4 / f"step3_{channel}_ff.png", dpi=110); plt.close(fig)
    fr = out["fractions"]
    mt = np.asarray(fr["mt_bins"]); xm = 0.5 * (mt[1:] + mt[:-1]); xm[-1] = mt[-2] + 30
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for j, ax in enumerate(axes):
        bottom = np.zeros(len(xm))
        for proc, col in (("qcd", "C0"), ("w", "C3"), ("tt", "C2")):
            v = np.asarray(fr["fractions"][proc])[j]
            ax.bar(xm, v, width=np.diff(mt) * 0.9, bottom=bottom, color=col, label=proc, alpha=0.8)
            bottom += v
        ax.set_title(f"N_jets {'=' if j < 2 else '>='} {fk.NJ_BINS[j]}"); ax.set_xlabel(r"$m_T(\ell, MET)$ [GeV]"); ax.set_xlim(0, 180)
    axes[0].set_ylabel("fraction of jet fakes in the AR"); axes[0].legend()
    fig.suptitle(f"{channel}: AR composition"); fig.tight_layout(); plotting.fig_tag(fig)
    fig.savefig(config.PLOT_DIR_V4 / f"step3_{channel}_fractions.png", dpi=110); plt.close(fig)
    c = out["closure_ss"]
    e = np.asarray(c["edges"])
    plotting.stack_plot(config.PLOT_DIR_V4 / f"step3_{channel}_closure_SS.png", e, (np.asarray(c["obs_hist"]), np.sqrt(np.abs(c["obs_hist"]))),
                        [("Fakes", np.asarray(c["pred_hist"]), np.zeros(len(e) - 1))], r"$m_{\tau\tau}$ [GeV]",
                        title=f"{channel}: same-sign validation (data - genuine MC vs FF prediction)", density=True)


def emu(plots: bool):
    data, reg, mc = load_channel("emu")
    edges = DR_BINS
    out = {"channel": "emu"}
    for sb in ("SB1", "SB2"):
        res = {}
        for tag in ("OS", "SS"):
            m = reg[f"{sb}_{tag}"]
            h = np.histogram(np.clip(data["dr"][m], 0, edges[-1] - 1e-6), bins=edges)[0].astype(float)
            hmc = np.zeros(len(edges) - 1)
            for _, _, d, cm, w, r_all, _ in mc:
                mm = r_all[f"{sb}_{tag}"] & cm
                hmc += np.histogram(np.clip(d["dr"][mm], 0, edges[-1] - 1e-6), bins=edges, weights=w[mm])[0]
            res[tag] = {"data": h.tolist(), "mc": hmc.tolist(), "net": np.maximum(h - hmc, 0).tolist()}
        os_, ss = np.asarray(res["OS"]["net"]), np.asarray(res["SS"]["net"])
        ratio = np.where(ss > 0, os_ / np.maximum(ss, 1e-9), 1.0)
        stat = ratio * np.sqrt(1 / np.maximum(np.asarray(res["OS"]["data"]), 1) + 1 / np.maximum(np.asarray(res["SS"]["data"]), 1))
        res["ratio"] = ratio.tolist(); res["stat"] = stat.tolist()
        res["inclusive"] = float(os_.sum() / max(ss.sum(), 1e-9))
        out[sb] = res
    r1, r2 = np.asarray(out["SB1"]["ratio"]), np.asarray(out["SB2"]["ratio"])
    out["osss"] = {"dr_edges": edges.tolist(), "ratio": r1.tolist(), "stat": out["SB1"]["stat"],
                   "syst": np.abs(r2 - r1).tolist(), "rel_unc": np.hypot(np.asarray(out["SB1"]["stat"]) / np.maximum(r1, 1e-9), np.abs(r2 - r1) / np.maximum(r1, 1e-9)).tolist()}
    # SS template in the signal region selection: data - MC
    n_ss = int(reg["SS"].sum())
    n_ss_mc = sum(float(w[r_all["SS"] & cm].sum()) for _, _, d, cm, w, r_all, _ in mc)
    out["ss_region"] = {"data": n_ss, "mc": n_ss_mc, "net": n_ss - n_ss_mc}
    idx = np.clip(np.searchsorted(edges, data["dr"][reg["SS"]], side="right") - 1, 0, len(edges) - 2)
    out["sr_multijet"] = float((r1[idx]).sum() - 0.0) - float(sum((w * r1[np.clip(np.searchsorted(edges, d["dr"], side="right") - 1, 0, len(edges) - 2)])[r_all["SS"] & cm].sum() for _, _, d, cm, w, r_all, _ in mc))
    print(f"[emu] OS/SS per dR bin SB1: {np.round(r1, 2)}  SB2: {np.round(r2, 2)}  inclusive {out['SB1']['inclusive']:.2f} / {out['SB2']['inclusive']:.2f}")
    print(f"[emu] SS region: data {n_ss}, MC {n_ss_mc:.0f} -> multijet in the SR {out['sr_multijet']:.0f}; data SR {reg['SR'].sum()}, CRtt {reg['CRtt'].sum()}")
    (config.DATA_DIR_V4 / "fakes_emu.json").write_text(json.dumps(fk.to_json(out), indent=1, default=float))
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 5))
        x = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(x, r1, np.asarray(out["SB1"]["stat"]), fmt="o", label="SB1: both iso < 0.6, one > 0.15")
        ax.errorbar(x + 0.05, r2, np.asarray(out["SB2"]["stat"]), fmt="s", label="SB2: one iso > 0.6")
        ax.set_xlabel(r"$\Delta R(e, \mu)$"); ax.set_ylabel("OS / SS (multijet)"); ax.set_ylim(0, 3); ax.grid(alpha=0.3); ax.legend()
        plotting.label(ax, "e#mu multijet OS/SS")
        fig.tight_layout(); plotting.fig_tag(fig)
        fig.savefig(config.PLOT_DIR_V4 / "step3_emu_osss.png", dpi=110); plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--channels", nargs="*", default=["mutau", "etau", "emu"])
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    config.DATA_DIR_V4.mkdir(parents=True, exist_ok=True)
    config.PLOT_DIR_V4.mkdir(parents=True, exist_ok=True)
    for ch in args.channels:
        if ch == "emu":
            emu(not args.no_plots)
        else:
            ltau(ch, not args.no_plots)


if __name__ == "__main__":
    main()
