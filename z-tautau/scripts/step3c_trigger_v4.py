#!/usr/bin/env python
"""Step 3c (v4) -- in-situ trigger efficiencies and the b-tag efficiency map (docs/10-v4-plan.md section 7).

    python scripts/step3c_trigger_v4.py [--no-plots]

Trigger efficiencies from e mu events (ttbar and Z -> tautau -> e mu), which are unbiased for the lepton that
did not fire the selecting trigger:
  * `ele27`  : SingleMuon-triggered e mu events (emu_mu ntuples: muon pT > 26 matched to IsoMu24), efficiency of
               "HLT_Ele27_WPTight_Gsf fired and the electron matched" vs (pT, |eta_SC|) of the electron, in data and
               simulation; SF = ratio.
  * `emu_e`  : the same events, efficiency of the e mu cross-trigger OR vs the electron pT (electron legs + DZ).
  * `emu_mu` : SingleElectron-triggered e mu events (emu_el ntuples), efficiency of the cross-trigger OR vs the muon pT.
  * `iso24`  : cross-check of the POG IsoMu24 scale factor from the emu_el events.
The e mu trigger SF is SF_e(pT_e) x SF_mu(pT_mu) (the two measurements each contain the DZ filter once; the
double counting is at the per-mille level for the OR of two paths and is covered by the 2% added in
analysis_v4._insitu_sf). Output: external/trigger_insitu_v4.json (committed) and output_v4/plots/step3c_*.png.

b-tag efficiencies eps(flavour, pT) of the DeepJet medium working point from the ttbar and Drell-Yan simulation
(jets in the b-tag acceptance of the e mu selection): output_v4/data/btag_eff.json, used by analysis_v4.btag_weight.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import analysis_v4 as an, config, plotting  # noqa: E402

PT_EDGES_E = np.array([13.0, 20.0, 24.0, 29.0, 35.0, 45.0, 60.0, 100.0, 500.0])
PT_EDGES_MU = np.array([10.0, 15.0, 20.0, 24.0, 26.0, 30.0, 40.0, 60.0, 100.0, 500.0])
ETA_EDGES_E = np.array([0.0, 0.8, 1.479, 2.1, 2.5])
ETA_EDGES_MU = np.array([0.0, 0.9, 1.2, 2.1, 2.4])
BTAG_PT_EDGES = np.array([20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 1000.0])


def _pairs(d, channel):
    """Events of the side sample with a clean e mu pair (both isolated, OS, no extra lepton)."""
    r = an.regions(d, "emu")             # e mu regions apply to the side samples as well (same fields)
    base = (d["n_veto_mu"] == 0) & (d["n_veto_el"] == 0) & d["os"] & (d["mu_iso"] < config.EMU_ISO) & (d["el_iso"] < config.EMU_ISO)
    if channel == "emu_mu":
        base &= (d["mu_pt"] > config.MU_PT_MIN_MUTAU) & d["mu_match_iso24"] & (d["el_pt"] > 13)
    else:
        base &= (d["el_pt"] > config.EL_PT_MIN_ETAU) & d["el_match_ele27"] & (d["mu_pt"] > 10)
    return base


def _eff(x, y, passed, w, xe, ye):
    num = np.histogram2d(np.clip(x, xe[0], xe[-1] - 1e-3), np.clip(np.abs(y), ye[0], ye[-1] - 1e-3), bins=[xe, ye], weights=w * passed)[0]
    den = np.histogram2d(np.clip(x, xe[0], xe[-1] - 1e-3), np.clip(np.abs(y), ye[0], ye[-1] - 1e-3), bins=[xe, ye], weights=w)[0]
    den2 = np.histogram2d(np.clip(x, xe[0], xe[-1] - 1e-3), np.clip(np.abs(y), ye[0], ye[-1] - 1e-3), bins=[xe, ye], weights=w ** 2)[0]
    eff = np.where(den > 0, num / np.maximum(den, 1e-9), 0.0)
    n_eff = np.where(den2 > 0, den ** 2 / np.maximum(den2, 1e-9), 0.0)          # effective entries
    err = np.sqrt(np.maximum(eff * (1 - eff), 1e-6) / np.maximum(n_eff, 1.0))
    return eff, err, den


def measure(name, channel, x_field, y_field, passed_fn, xe, ye):
    data = an.load_data(channel)
    sel = _pairs(data, channel)
    ed, eed, nd = _eff(data[x_field][sel], data[y_field][sel], passed_fn(data)[sel].astype(float), np.ones(sel.sum()), xe, ye)
    num = np.zeros_like(ed); den = np.zeros_like(ed); den2 = np.zeros_like(ed)
    for key in an.available_mc(channel):
        d, _ = an.load(key, channel)
        d = dict(d)
        s = _pairs(d, channel)
        if not s.any():
            continue
        w = an.weights(d, key, "emu")
        x, y, p, ww = d[x_field][s], np.abs(d[y_field][s]), passed_fn(d)[s].astype(float), w[s]
        num += np.histogram2d(np.clip(x, xe[0], xe[-1] - 1e-3), np.clip(y, ye[0], ye[-1] - 1e-3), bins=[xe, ye], weights=ww * p)[0]
        den += np.histogram2d(np.clip(x, xe[0], xe[-1] - 1e-3), np.clip(y, ye[0], ye[-1] - 1e-3), bins=[xe, ye], weights=ww)[0]
        den2 += np.histogram2d(np.clip(x, xe[0], xe[-1] - 1e-3), np.clip(y, ye[0], ye[-1] - 1e-3), bins=[xe, ye], weights=ww ** 2)[0]
    em = np.where(den > 0, num / np.maximum(den, 1e-9), 0.0)
    eem = np.sqrt(np.maximum(em * (1 - em), 1e-6) / np.maximum(np.where(den2 > 0, den ** 2 / np.maximum(den2, 1e-9), 0), 1.0))
    sf = np.where((em > 0) & (ed > 0), ed / np.maximum(em, 1e-9), 1.0)
    err = np.where((em > 0) & (ed > 0), sf * np.hypot(eed / np.maximum(ed, 1e-9), eem / np.maximum(em, 1e-9)), 0.05)
    return {"x_edges": xe.tolist(), "y_edges": ye.tolist(), "abs_y": True, "x": x_field, "y": y_field,
            "eff_data": ed.tolist(), "eff_data_err": eed.tolist(), "eff_mc": em.tolist(), "eff_mc_err": eem.tolist(),
            "sf": sf.tolist(), "err": err.tolist(), "n_data": nd.tolist(), "n_mc": den.tolist()}


def btag_efficiencies():
    num = {"b": np.zeros(len(BTAG_PT_EDGES) - 1), "c": np.zeros(len(BTAG_PT_EDGES) - 1), "light": np.zeros(len(BTAG_PT_EDGES) - 1)}
    den = {k: np.zeros_like(v) for k, v in num.items()}
    for key in ("TTTo2L2Nu", "DY_NLO", "TTToSemiLeptonic", "WW"):
        if not an.has_ntuple(key, "emu"):
            continue
        d, _ = an.load(key, "emu")
        if not len(d):
            continue
        w = an.weights(dict(d), key, "emu")
        for i in range(1, 5):
            pt, fl, tag = d[f"bj{i}_pt"], d[f"bj{i}_flav"].astype(int), d[f"bj{i}_tag"].astype(bool)
            has = pt > 0
            ib = np.clip(np.searchsorted(BTAG_PT_EDGES, pt, side="right") - 1, 0, len(BTAG_PT_EDGES) - 2)
            for f, name in ((5, "b"), (4, "c"), (0, "light")):
                sel = has & (fl == f)
                np.add.at(den[name], ib[sel], w[sel])
                np.add.at(num[name], ib[sel & tag], w[sel & tag])
    out = {"pt_edges": BTAG_PT_EDGES.tolist()}
    for k in num:
        out[k] = np.where(den[k] > 0, num[k] / np.maximum(den[k], 1e-9), 0.5).tolist()
        out[k + "_n"] = den[k].tolist()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    config.DATA_DIR_V4.mkdir(parents=True, exist_ok=True)
    config.PLOT_DIR_V4.mkdir(parents=True, exist_ok=True)
    emu_or = lambda d: np.any([d[t] for t in config.EMU_TRIGGERS], axis=0)    # noqa: E731
    out = {}
    out["ele27"] = measure("ele27", "emu_mu", "el_pt", "el_sceta", lambda d: d["HLT_Ele27_WPTight_Gsf"] & d["el_match_ele27"], PT_EDGES_E, ETA_EDGES_E)
    out["emu_e"] = measure("emu_e", "emu_mu", "el_pt", "el_sceta", emu_or, PT_EDGES_E, ETA_EDGES_E)
    out["emu_mu"] = measure("emu_mu", "emu_el", "mu_pt", "mu_eta", emu_or, PT_EDGES_MU, ETA_EDGES_MU)
    out["iso24"] = measure("iso24", "emu_el", "mu_pt", "mu_eta", lambda d: (d["HLT_IsoMu24"] | d["HLT_IsoTkMu24"]) & d["mu_match_iso24"], PT_EDGES_MU, ETA_EDGES_MU)
    for name, t in out.items():
        sf, e, n = np.asarray(t["sf"]), np.asarray(t["err"]), np.asarray(t["n_data"])
        ok = n > 50
        rng = f"{sf[ok].min():.3f}-{sf[ok].max():.3f}, median err {np.median(e[ok]):.3f}" if ok.any() else "no bin with > 50 pairs"
        print(f"[trigger] {name}: eff_data(plateau, |eta| bin 0) {np.round(np.asarray(t['eff_data'])[-3:, 0], 3)}, SF range {rng}, data pairs {n.sum():.0f}")
    an.TRIG_INSITU.write_text(json.dumps(out, indent=1))
    print(f"-> {an.TRIG_INSITU}")
    beff = btag_efficiencies()
    an.BTAG_EFF.write_text(json.dumps(beff, indent=1))
    print("[btag] efficiencies (M): b", np.round(beff["b"], 3), "c", np.round(beff["c"], 3), "light", np.round(beff["light"], 4))
    if args.no_plots:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for name, t in out.items():
        xe = np.asarray(t["x_edges"]); x = 0.5 * (xe[1:] + xe[:-1]); x[-1] = xe[-2] * 1.3
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        for j in range(len(t["y_edges"]) - 1):
            axes[0].errorbar(x, np.asarray(t["eff_data"])[:, j], np.asarray(t["eff_data_err"])[:, j], fmt="o", ms=4, color=f"C{j}", label=f"data |eta| {t['y_edges'][j]}-{t['y_edges'][j + 1]}")
            axes[0].errorbar(x, np.asarray(t["eff_mc"])[:, j], np.asarray(t["eff_mc_err"])[:, j], fmt="s", ms=3, mfc="none", color=f"C{j}")
            axes[1].errorbar(x, np.asarray(t["sf"])[:, j], np.asarray(t["err"])[:, j], fmt="o", ms=4, color=f"C{j}")
        axes[0].set_ylabel("efficiency (filled: data, open: simulation)"); axes[1].set_ylabel("SF = data / simulation")
        for ax in axes:
            ax.set_xlabel(f"{t['x']} [GeV]"); ax.set_xscale("log"); ax.grid(alpha=0.3)
        axes[0].set_ylim(0, 1.1); axes[1].set_ylim(0.7, 1.3); axes[0].legend(fontsize=7)
        plotting.label(axes[0], f"in-situ trigger efficiency: {name}")
        fig.tight_layout(); plotting.fig_tag(fig)
        fig.savefig(config.PLOT_DIR_V4 / f"step3c_trigger_{name}.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
