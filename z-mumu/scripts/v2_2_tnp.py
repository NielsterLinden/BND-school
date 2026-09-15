#!/usr/bin/env python
"""v2 step 2 -- muon efficiencies and scale factors from tag-and-probe with fits.

    python scripts/v2_2_tnp.py                 # fill (data + DY NLO skims, parallel), fit, plot
    python scripts/v2_2_tnp.py --summarise-only

Inputs: the v2 skims ($BND_SKIM_DIR). Outputs: output/v2/tnp/tnp_parts/*.pkl (per file),
output/v2/tnp/tnp_result.json (efficiency and scale-factor maps read by zmumu.tnp.ScaleFactors),
output/v2/tnp/tnp_fits.pkl, output/v2/plots/tnp_*.png. See docs/12-tag-and-probe-fits.md.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import numpy as np

from zmumu import batch, config, hists, pileup, skim, tnp

OUT = config.OUTPUT_DIR / "v2" / "tnp"
config.PLOT_DIR = config.OUTPUT_DIR / "v2" / "plots"   # v2 plots live next to the v2 results
BRANCHES = ["run", "PV_npvsGood", "HLT_IsoMu24", "HLT_IsoTkMu24", "Muon_*", "TrigObj_*", "Flag_*",
            "genWeight", "Pileup_nTrueInt", "skim_cat"]


# ----------------------------------------------------------------------------- per file
def process_file(source, is_mc, pu_json):
    import uproot
    out = tnp.blank()
    out["n_files"] = 1
    pu = pileup.PileupWeights.from_json(pu_json) if (is_mc and pu_json) else None
    with uproot.open(source) as f:
        if "Events" not in f:
            return out
        for ev in f["Events"].iterate(filter_name=BRANCHES, step_size="200 MB", library="ak"):
            w = None
            if is_mc:
                w = np.sign(np.asarray(ev.genWeight, dtype=float))
                if pu is not None:
                    w = w * pu(np.asarray(ev.Pileup_nTrueInt))
            tnp.fill(ev, out, weight=w, is_mc=is_mc)
    return out


# ----------------------------------------------------------------------------- summary
def _sum_eras(parts):
    total = tnp.blank()
    for p in parts:
        batch.merge(total, p)
    return total


def measure(data, mc, mc_alt=None):
    """All fits and counts -> result dict with maps and per-cell details."""
    npt, neta = len(tnp.PT_EDGES) - 1, len(tnp.ETA_EDGES) - 1
    res = {"pt_edges": tnp.PT_EDGES.tolist(), "eta_edges": tnp.ETA_EDGES.tolist(),
           "trig_pt_edges": tnp.TRIG_PT_EDGES.tolist(), "sf": {}, "sf_err": {}, "sf_stat": {}, "sf_syst": {},
           "eff": {}, "eff_err": {}, "details": {}, "meta": {"n_pairs_data": int(data["n_pairs"]),
                                                             "n_pairs_mc": float(mc["n_pairs"])}}
    variants = [("nominal", dict(bkg="expo", rng=(60, 120), tag="nom", tmpl="gen")),
                ("bkg_cmsshape", dict(bkg="cmsshape", rng=(60, 120), tag="nom", tmpl="gen")),
                ("range_70_110", dict(bkg="expo", rng=(70, 110), tag="nom", tmpl="gen")),
                ("tag_alt", dict(bkg="expo", rng=(60, 120), tag="alt", tmpl="gen"))]
    if mc_alt is not None:
        variants.append(("tmpl_powheg", dict(bkg="expo", rng=(60, 120), tag="nom", tmpl="alt")))
    for eff in ("id", "iso", "antiiso"):
        maps = {name: {"data": np.zeros((npt, neta)), "data_err": np.zeros((npt, neta)),
                       "mc": np.zeros((npt, neta)), "mc_err": np.zeros((npt, neta))} for name, _ in variants}
        truth = np.zeros((npt, neta)); truth_err = np.zeros((npt, neta))
        count_d = np.zeros((npt, neta)); count_m = np.zeros((npt, neta))
        details = {}
        t0 = time.time()
        for i in range(npt):
            for j in range(neta):
                key = f"{eff}_nom_os"
                gp, gf = mc[f"{key}_pass_gen"][i, j], mc[f"{key}_fail_gen"][i, j]
                truth[i, j], truth_err[i, j] = tnp.count_eff(gp, gf)
                count_d[i, j], _ = tnp.count_eff(data[f"{key}_pass"][i, j], data[f"{key}_fail"][i, j],
                                                 data[f"{eff}_nom_ss_pass"][i, j], data[f"{eff}_nom_ss_fail"][i, j])
                count_m[i, j], _ = tnp.count_eff(mc[f"{key}_pass"][i, j], mc[f"{key}_fail"][i, j],
                                                 mc[f"{eff}_nom_ss_pass"][i, j], mc[f"{eff}_nom_ss_fail"][i, j])
                for name, v in variants:
                    k = f"{eff}_{v['tag']}_os"
                    src = mc_alt if v["tmpl"] == "alt" else mc
                    tp, tf = src[f"{k}_pass_gen"][i, j], src[f"{k}_fail_gen"][i, j]
                    for sample, h in (("data", data), ("mc", mc)):
                        r = tnp.fit_cell(h[f"{k}_pass"][i, j], h[f"{k}_fail"][i, j], tp, tf,
                                         pass_var=h[f"{k}_pass_w2"][i, j] if sample == "mc" else None,
                                         fail_var=h[f"{k}_fail_w2"][i, j] if sample == "mc" else None,
                                         bkg=v["bkg"], mass_range=v["rng"])
                        maps[name][sample][i, j] = r["eps"]
                        maps[name][f"{sample}_err"][i, j] = r["err"]
                        if name == "nominal":
                            details[(sample, i, j)] = {kk: r[kk] for kk in ("eps", "err", "N", "valid", "chi2", "ndf",
                                                                           "params", "model_pass", "model_fail", "x",
                                                                           "bkg_frac_fail")}
        nom = maps["nominal"]
        with np.errstate(divide="ignore", invalid="ignore"):
            sf = np.where(nom["mc"] > 0, nom["data"] / np.where(nom["mc"] > 0, nom["mc"], 1), 1.0)
            stat = sf * np.sqrt((nom["data_err"] / np.maximum(nom["data"], 1e-9)) ** 2
                                + (nom["mc_err"] / np.maximum(nom["mc"], 1e-9)) ** 2)
            syst2 = np.zeros_like(sf)
            for name, _ in variants[1:]:
                m = maps[name]
                sf_v = np.where(m["mc"] > 0, m["data"] / np.where(m["mc"] > 0, m["mc"], 1), 1.0)
                syst2 = np.maximum(syst2, (sf_v - sf) ** 2)
            closure = np.abs(nom["mc"] - truth) / np.maximum(truth, 1e-9) * sf
            syst = np.sqrt(syst2 + closure ** 2)
        res["sf"][eff] = sf.tolist(); res["sf_stat"][eff] = stat.tolist(); res["sf_syst"][eff] = syst.tolist()
        res["sf_err"][eff] = np.sqrt(stat ** 2 + syst ** 2).tolist()
        for name, _ in variants:
            res["eff"][f"{eff}_data_{name}"] = maps[name]["data"].tolist()
            res["eff_err"][f"{eff}_data_{name}"] = maps[name]["data_err"].tolist()
            res["eff"][f"{eff}_mc_{name}"] = maps[name]["mc"].tolist()
            res["eff_err"][f"{eff}_mc_{name}"] = maps[name]["mc_err"].tolist()
        res["eff"][f"{eff}_mc_truth"] = truth.tolist(); res["eff_err"][f"{eff}_mc_truth"] = truth_err.tolist()
        res["eff"][f"{eff}_data_count"] = count_d.tolist(); res["eff"][f"{eff}_mc_count"] = count_m.tolist()
        res["details"][eff] = details
        print(f"  {eff}: {npt*neta*len(variants)*2} fits in {time.time()-t0:.0f}s; "
              f"<SF> = {np.mean(sf):.4f}, <stat> {np.mean(stat):.4f}, <syst> {np.mean(syst):.4f}, "
              f"<MC fit - truth> {np.mean(nom['mc'] - truth):+.4f}", flush=True)

    # trigger: same-sign subtracted counting per cell, data and MC; L1 veto variant as syst (data)
    for sample, h in (("data", data), ("mc", mc)):
        e, err = {}, {}
        for variant in ("nominal", "l1veto"):
            p = h[f"trig_{variant}_os_pass"] - h[f"trig_{variant}_ss_pass"]
            t = h[f"trig_{variant}_os_tot"] - h[f"trig_{variant}_ss_tot"]
            with np.errstate(divide="ignore", invalid="ignore"):
                e[variant] = np.where(t > 0, p / np.where(t > 0, t, 1), 0.0)
                neff = np.where(t > 0, t ** 2 / np.maximum(h[f"trig_{variant}_os_tot_w2"], 1e-9), 0.0)
                err[variant] = np.sqrt(np.where(neff > 0, e[variant] * (1 - e[variant]) / np.maximum(neff, 1), 0.0))
        tot_err = err["nominal"]
        if sample == "data":
            tot_err = np.sqrt(err["nominal"] ** 2 + (e["l1veto"] - e["nominal"]) ** 2)
        res["eff"][f"trig_{sample}"] = e["nominal"].tolist()
        res["eff_err"][f"trig_{sample}"] = tot_err.tolist()
        res["eff"][f"trig_{sample}_l1veto"] = e["l1veto"].tolist()
        plateau = slice(list(tnp.TRIG_PT_EDGES).index(26.0), None)
        p = (h["trig_nominal_os_pass"] - h["trig_nominal_ss_pass"])[plateau].sum()
        t = (h["trig_nominal_os_tot"] - h["trig_nominal_ss_tot"])[plateau].sum()
        res["meta"][f"trig_plateau_{sample}"] = float(p / t) if t > 0 else 0.0
        print(f"  trigger {sample}: per-muon plateau efficiency (pT > 26) {p/t:.4f}")

    # isolation efficiency vs pileup
    for sample, h in (("data", data), ("mc", mc)):
        p = h["isonpv_os_pass"] - h["isonpv_ss_pass"]
        t = h["isonpv_os_tot"] - h["isonpv_ss_tot"]
        res["eff"][f"isonpv_{sample}"] = (p / np.maximum(t, 1e-9)).tolist()
    d = np.array(res["eff"]["isonpv_data"]); m = np.array(res["eff"]["isonpv_mc"])
    sf_npv = d / np.maximum(m, 1e-9)
    res["meta"]["iso_sf_vs_npv"] = sf_npv.tolist()
    res["meta"]["iso_sf_npv_spread"] = float(0.5 * (sf_npv.max() - sf_npv.min()))
    print(f"  iso SF in nPV bins {sf_npv.round(4)} -> half spread {res['meta']['iso_sf_npv_spread']:.4f}")
    return res


def plot(res, data, mc):
    import matplotlib.pyplot as plt
    pt, eta = np.array(res["pt_edges"]), np.array(res["eta_edges"])
    c, hw = 0.5 * (pt[1:] + pt[:-1]), 0.5 * np.diff(pt)
    for eff, nice in (("id", "Tight ID"), ("iso", "Isolation (< 0.15)"), ("antiiso", "Anti-isolation (0.20-1.0)")):
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))
        for j, ax in enumerate(axes.flat):
            for key, lab, col, mk in ((f"{eff}_data_nominal", "data (fit)", "black", "o"),
                                      (f"{eff}_mc_nominal", "DY MC (fit)", "#1f77b4", "s"),
                                      (f"{eff}_mc_truth", "DY MC gen-matched (truth)", "#d62728", "^")):
                e = np.array(res["eff"][key])[:, j]; er = np.array(res["eff_err"][key])[:, j]
                ax.errorbar(c, e, yerr=er, xerr=hw, fmt=mk, color=col, markersize=4, label=lab)
            e = np.array(res["eff"][f"{eff}_data_count"])[:, j]
            ax.plot(c, e, linestyle="none", marker="x", color="grey", label="data (SS-subtracted count)")
            hists.log_pt_axis(ax); ax.set_xlabel(r"probe $p_T$ [GeV]"); ax.set_ylabel("efficiency")
            lo = 0.0 if eff == "antiiso" else 0.75
            ax.set_ylim(lo, 1.05 if eff != "antiiso" else 0.3)
            ax.set_title(f"{eta[j]:.1f} < |$\\eta$| < {eta[j+1]:.1f}", fontsize=13)
            ax.legend(fontsize=9, loc="lower right")
        fig.suptitle(f"{nice} efficiency, tag-and-probe with pass/fail fits", fontsize=15)
        fig.tight_layout()
        hists.save_fig(fig, f"tnp_eff_{eff}_vs_pt.png")
        # SF map
        sf = np.array(res["sf"][eff]); err = np.array(res["sf_err"][eff])
        fig, ax = plt.subplots(figsize=(10, 7))
        mesh = ax.pcolormesh(pt, eta, sf.T, cmap="coolwarm", vmin=0.9, vmax=1.1)
        for i in range(len(pt) - 1):
            for j in range(len(eta) - 1):
                ax.text(np.sqrt(pt[i] * pt[i + 1]), 0.5 * (eta[j] + eta[j + 1]), f"{sf[i,j]:.3f}\n±{err[i,j]:.3f}",
                        ha="center", va="center", fontsize=7)
        fig.colorbar(mesh, ax=ax, label=f"{nice} scale factor")
        hists.log_pt_axis(ax); ax.set_xlabel(r"$p_T$ [GeV]"); ax.set_ylabel(r"|$\eta$|")
        hists._decorate(ax); hists._title(ax, f"{nice} data/MC scale factor")
        hists.save_fig(fig, f"tnp_sf_{eff}_map.png")
        # example fits: lowest pT bin in every eta bin, data and MC
        det = res["details"][eff]
        fig, axes = plt.subplots(4, 4, figsize=(20, 16))
        for j in range(4):
            for col, (sample, h) in enumerate((("data", data), ("mc", mc))):
                d = det[(sample, 0, j)]
                for pf, ax in (("pass", axes[j, 2 * col]), ("fail", axes[j, 2 * col + 1])):
                    obs = h[f"{eff}_nom_os_{pf}"][0, j]
                    ax.errorbar(tnp.MASS_CENTRES, obs, yerr=np.sqrt(np.maximum(obs, 0)), fmt="o", color="black", markersize=2)
                    ax.plot(d["x"], d[f"model_{pf}"], color="crimson")
                    ax.set_title(f"{sample} {pf}, 20-25 GeV, {eta[j]:.1f}<|η|<{eta[j+1]:.1f}: ε={d['eps']:.4f}±{d['err']:.4f}, χ²/ndf={d['chi2']:.0f}/{d['ndf']}", fontsize=8)
                    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
        fig.suptitle(f"{nice}: pass/fail fits in the lowest pT bin", fontsize=15)
        fig.tight_layout()
        hists.save_fig(fig, f"tnp_fits_{eff}_lowpt.png")
    # trigger
    tpt = np.array(res["trig_pt_edges"]); tc, thw = 0.5 * (tpt[1:] + tpt[:-1]), 0.5 * np.diff(tpt)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    for j, ax in enumerate(axes.flat):
        for sample, col in (("data", "black"), ("mc", "#1f77b4")):
            e = np.array(res["eff"][f"trig_{sample}"])[:, j]; er = np.array(res["eff_err"][f"trig_{sample}"])[:, j]
            ax.errorbar(tc, e, yerr=er, xerr=thw, fmt="o", color=col, markersize=4, label=sample)
        hists.log_pt_axis(ax); ax.set_ylim(0, 1.05); ax.set_xlabel(r"probe $p_T$ [GeV]")
        ax.set_ylabel("IsoMu24 || IsoTkMu24 per-muon efficiency")
        ax.set_title(f"{eta[j]:.1f} < |$\\eta$| < {eta[j+1]:.1f}", fontsize=13); ax.legend(fontsize=10)
    fig.suptitle("Trigger efficiency, trigger-object tag-and-probe (same-sign subtracted)", fontsize=15)
    fig.tight_layout(); hists.save_fig(fig, "tnp_trig_eff_vs_pt.png")
    # iso vs npv
    fig, ax = plt.subplots(figsize=(9, 7))
    npv = tnp.NPV_EDGES; nc = 0.5 * (npv[1:] + npv[:-1]); nc[-1] = 40
    ax.plot(nc, res["eff"]["isonpv_data"], "o", color="black", label="data")
    ax.plot(nc, res["eff"]["isonpv_mc"], "s", color="#1f77b4", label="DY MC (PU reweighted)")
    ax.plot(nc, res["meta"]["iso_sf_vs_npv"], "^", color="crimson", label="scale factor")
    ax.set_xlabel("good primary vertices"); ax.set_ylabel(r"isolation efficiency ($p_T$ > 26 GeV)")
    ax.set_ylim(0.85, 1.05); ax.legend(fontsize=12); hists._decorate(ax)
    hists._title(ax, "Isolation efficiency vs pileup")
    hists.save_fig(fig, "tnp_iso_vs_npv.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--plots-only", action="store_true")
    ap.add_argument("--skim-dir", default=None)
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--is-mc", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--pu-json", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:
        batch.write_part(process_file(args.one_file, args.is_mc, args.pu_json), args.out)
        return

    OUT.mkdir(parents=True, exist_ok=True)
    pu_json = OUT / "pileup_weights.json"
    if args.plots_only:
        saved = pickle.load(open(OUT / "tnp_fits.pkl", "rb"))
        plot(saved["result"], saved["data"], saved["mc"])
        return
    if not args.summarise_only:
        if pu_json.exists():        # N_PV-matched profile from scripts/v2_2_pileup.py
            pu = pileup.PileupWeights.from_json(pu_json)
        else:
            print("[tnp] no pileup_weights.json: run scripts/v2_2_pileup.py first (using the raw CSV profile)")
            pu = pileup.PileupWeights(pileup.mc_profile(("DY_NLO",), args.skim_dir))
            pu.to_json(pu_json)
        print(f"[tnp] pileup: data <mu> {pu.data['nominal']['mean']:.2f} (scale {pu.scale:.3f}, smear {pu.rel_smear:.2f}), MC <mu> "
              f"{(pu.mc * np.arange(len(pu.mc))).sum() / pu.mc.sum():.2f}, L = {pu.data['nominal']['lumi_pb']:.3f} pb^-1")
        for key, is_mc in (("data_2016G", False), ("data_2016H", False), ("DY_NLO", True), ("DY_powheg", True)):
            files = skim.skim_files(key, args.skim_dir)[: args.max_files]
            if not files:
                print(f"[tnp] {key}: no skim files, skipped")
                continue
            tasks = [(str(f), f"{key}__{f.stem}") for f in files]
            extra = ["--is-mc", "--pu-json", str(pu_json)] if is_mc else []
            print(f"[tnp] {key}: {len(tasks)} skim files", flush=True)
            batch.run_files(tasks, Path(__file__), OUT / "tnp_parts" / key, args.workers, extra_args=extra,
                            stall_s=1800, retries=1, log=lambda s: print(s, flush=True))
    totals = {}
    for key in ("data_2016G", "data_2016H", "DY_NLO", "DY_powheg"):
        d = OUT / "tnp_parts" / key
        if d.exists() and any(d.glob("*.pkl")):
            totals[key] = batch.load_parts(d, tnp.blank())
    data = _sum_eras([totals[k] for k in ("data_2016G", "data_2016H") if k in totals])
    mc = totals["DY_NLO"]
    print(f"[tnp] data pairs {data['n_pairs']:,}  MC pairs (weighted) {mc['n_pairs']:,.0f}")
    res = measure(data, mc, totals.get("DY_powheg"))
    with open(OUT / "tnp_fits.pkl", "wb") as fh:
        pickle.dump({"result": res, "data": data, "mc": mc}, fh)
    slim = {k: v for k, v in res.items() if k != "details"}
    with open(OUT / "tnp_result.json", "w") as fh:
        json.dump(slim, fh, indent=1)
    plot(res, data, mc)
    print(f"[tnp] wrote {OUT / 'tnp_result.json'}")


if __name__ == "__main__":
    main()
