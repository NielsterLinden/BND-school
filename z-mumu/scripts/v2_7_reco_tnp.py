#!/usr/bin/env python
"""v2 step 7 -- muon reconstruction efficiency by tag-and-probe on the unskimmed NanoAOD.

    python scripts/v2_7_reco_tnp.py                          # data G+H (all files) + DY NLO (--mc-files)
    python scripts/v2_7_reco_tnp.py --max-files 2 --mc-files 1 --workers 3    # smoke test
    python scripts/v2_7_reco_tnp.py --summarise-only

Replaces the assigned SF_reco = 1 +- 0.4 %/muon (docs/11) by a measurement, the gap found when
comparing with CMS-SMP-20-004 (docs/16). The skims cannot be used: a failing probe is by
definition not a loose muon, and skim category A needs two of them. See zmumu/recoeff.py.

Outputs: output/v2/tnp/reco_parts/<sample>/*.pkl, output/v2/tnp/reco_result.json,
output/v2/tnp/reco_fits.pkl, output/v2/plots/reco_*.png.
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

from zmumu import batch, config, hists, io, pileup, recoeff, samples, tnp

OUT = config.OUTPUT_DIR / "v2" / "tnp"
config.PLOT_DIR = config.OUTPUT_DIR / "v2" / "plots"
# Background shapes: the failing probes' background is sculpted into a broad peak under the Z by the
# kinematic cuts (a tag plus a random isolated track), so the nominal takes the shape from the same-sign
# pairs of the same cell ("ss"); CMSShape (erfc x exp) is the parametric alternative. A plain exponential
# cannot describe the peak: it is fitted and reported ("expo", chi2 in the log) but not used.
VARIANTS = [("nominal", dict(bkg="ss", rng=(60, 120), tag="nom")),
            ("bkg_cmsshape", dict(bkg="cmsshape", rng=(60, 120), tag="nom")),
            ("range_70_110", dict(bkg="ss", rng=(70, 110), tag="nom")),
            ("tag_alt", dict(bkg="ss", rng=(60, 120), tag="alt"))]
REPORTED = [("bkg_expo", dict(bkg="expo", rng=(60, 120), tag="nom"))]
# stand-alone pT is poorly measured: failing (stand-alone-only) probes migrate to high pT, so the tracking
# efficiency is binned in |eta| only; isolated-track probes have tracker pT and keep the pT binning
PT_EDGES_EFF = {"trk": np.array([20.0, 200.0]), "mutrk": recoeff.PT_EDGES, "mutrkT": recoeff.PT_EDGES}
MUTRK_NOMINAL = "mutrkT"      # the purer isolated-track probe; "mutrk" is the loose cross-check


# ----------------------------------------------------------------------------- per file
def process_file(source, is_mc, pu_json):
    import uproot
    out = recoeff.blank()
    out["n_files"] = 1
    pu = pileup.PileupWeights.from_json(pu_json) if (is_mc and pu_json) else None
    grl = None if is_mc else io.load_grl()
    wanted = recoeff.BRANCHES_MC if is_mc else recoeff.BRANCHES_DATA
    with uproot.open(source) as f:
        tree = f["Events"]
        present = [b for b in wanted if b in tree.keys()]
        for ev in tree.iterate(present, step_size="300 MB", library="ak"):
            w = None
            if is_mc:
                w = np.sign(np.asarray(ev.genWeight, dtype=float)) * pu(np.asarray(ev.Pileup_nTrueInt))
            recoeff.fill(ev, out, weight=w, is_mc=is_mc, grl=grl)
    return out


# ----------------------------------------------------------------------------- summary
def _rebin_pt(h, eff):
    """(npt, neta, nmass) in recoeff.PT_EDGES -> the binning of `eff`."""
    if eff == "trk":
        return h.sum(axis=0, keepdims=True)
    return h


def _fit(v, data_or_mc, k, i, j, tp, tf, is_mc, data_ss_fallback):
    h = data_or_mc
    pass_h, fail_h = h[f"{k}_pass"][i, j], h[f"{k}_fail"][i, j]
    pv = h[f"{k}_pass_w2"][i, j] if is_mc else None
    fv = h[f"{k}_fail_w2"][i, j] if is_mc else None
    if v["bkg"] == "ss":
        kss = k.replace("_os", "_ss")
        fb_f = data_ss_fallback[f"{kss}_fail"]
        bf = recoeff.ss_template(h[f"{kss}_fail"][i, j], fallback=fb_f)
        bp = recoeff.ss_template(h[f"{kss}_pass"][i, j], fallback=h[f"{kss}_fail"][i, j] if h[f"{kss}_fail"][i, j].sum() >= 30 else fb_f)
        wide = k.startswith("trk_")          # stand-alone-only failing probes: momentum from the muon system
        return recoeff.fit_cell_template(pass_h, fail_h, tp, tf, bp, bf, pass_var=pv, fail_var=fv, mass_range=v["rng"],
                                         fail_sigma_max=15.0 if wide else 3.0, fail_shift_max=8.0 if wide else 2.0)
    return tnp.fit_cell(pass_h, fail_h, tp, tf, pass_var=pv, fail_var=fv, bkg=v["bkg"], mass_range=v["rng"])


def measure(data_raw, mc_raw):
    res = {"eta_edges": recoeff.ETA_EDGES.tolist(), "pt_edges": {e: v.tolist() for e, v in PT_EDGES_EFF.items()},
           "sf": {}, "sf_stat": {}, "sf_syst": {}, "sf_err": {}, "eff": {}, "eff_err": {}, "chi2": {}, "details": {},
           "meta": {"n_events_data": int(data_raw["n_events"]), "n_events_mc": int(mc_raw["n_events"]),
                    "n_probes_data": data_raw["n_probes"], "n_probes_mc": mc_raw["n_probes"],
                    "variants": [n for n, _ in VARIANTS], "reported_only": [n for n, _ in REPORTED]}}
    neta = len(recoeff.ETA_EDGES) - 1
    for eff in recoeff.EFFS:
        npt = len(PT_EDGES_EFF[eff]) - 1
        data = {k: _rebin_pt(v, eff) for k, v in data_raw.items() if k.startswith(eff + "_")}
        mc = {k: _rebin_pt(v, eff) for k, v in mc_raw.items() if k.startswith(eff + "_")}
        # eta-inclusive same-sign shapes for cells with too few same-sign entries
        fallback_d = {k: v.sum(axis=(0, 1)) for k, v in data.items() if "_ss_" in k}
        allv = VARIANTS + REPORTED
        maps = {name: {k: np.zeros((npt, neta)) for k in ("data", "data_err", "mc", "mc_err")} for name, _ in allv}
        chi2 = {name: [] for name, _ in allv}
        count_d, count_m, gen_count = (np.zeros((npt, neta)) for _ in range(3))
        details = {}
        t0 = time.time()
        for i in range(npt):
            for j in range(neta):
                key = f"{eff}_nom_os"
                gen_count[i, j], _ = tnp.count_eff(mc[f"{key}_pass_gen"][i, j], mc[f"{key}_fail_gen"][i, j])
                count_d[i, j], _ = tnp.count_eff(data[f"{key}_pass"][i, j], data[f"{key}_fail"][i, j],
                                                 data[f"{eff}_nom_ss_pass"][i, j], data[f"{eff}_nom_ss_fail"][i, j])
                count_m[i, j], _ = tnp.count_eff(mc[f"{key}_pass"][i, j], mc[f"{key}_fail"][i, j],
                                                 mc[f"{eff}_nom_ss_pass"][i, j], mc[f"{eff}_nom_ss_fail"][i, j])
                for name, v in allv:
                    k = f"{eff}_{v['tag']}_os"
                    tp, tf = mc[f"{k}_pass_gen"][i, j], mc[f"{k}_fail_gen"][i, j]
                    if tf.sum() < 20:        # too few generator-matched failing probes: smear the pass shape
                        tf = tp
                    for sample, h, is_mc in (("data", data, False), ("mc", mc, True)):
                        r = _fit(v, h, k, i, j, tp, tf, is_mc, fallback_d)
                        maps[name][sample][i, j] = r["eps"]
                        maps[name][f"{sample}_err"][i, j] = r["err"]
                        if sample == "data":
                            chi2[name].append(r["chi2"] / max(r["ndf"], 1))
                        if name == "nominal":
                            details[(sample, i, j)] = {kk: r[kk] for kk in ("eps", "err", "N", "valid", "chi2", "ndf",
                                                                           "model_pass", "model_fail", "x", "bkg_frac_fail")}
        nom = maps["nominal"]
        sf = np.where(nom["mc"] > 0, nom["data"] / np.maximum(nom["mc"], 1e-9), 1.0)
        stat = sf * np.sqrt((nom["data_err"] / np.maximum(nom["data"], 1e-9)) ** 2
                            + (nom["mc_err"] / np.maximum(nom["mc"], 1e-9)) ** 2)
        syst2 = np.zeros_like(sf)
        for name, _ in VARIANTS[1:]:
            m = maps[name]
            sf_v = np.where(m["mc"] > 0, m["data"] / np.maximum(m["mc"], 1e-9), 1.0)
            syst2 = np.maximum(syst2, (sf_v - sf) ** 2)
        # closure: the fit on simulation against the generator-matched count of the same probes
        closure = np.abs(nom["mc"] - gen_count) / np.maximum(gen_count, 1e-9) * sf
        syst = np.sqrt(syst2 + closure ** 2)
        if eff == "trk":           # broadcast the |eta| binning onto the common pT x |eta| grid
            reps = len(recoeff.PT_EDGES) - 1
            sf, stat, syst = (np.repeat(x, reps, axis=0) for x in (sf, stat, syst))
        res["sf"][eff], res["sf_stat"][eff], res["sf_syst"][eff] = sf.tolist(), stat.tolist(), syst.tolist()
        res["sf_err"][eff] = np.sqrt(stat ** 2 + syst ** 2).tolist()
        for name, _ in allv:
            for sample in ("data", "mc"):
                res["eff"][f"{eff}_{sample}_{name}"] = maps[name][sample].tolist()
                res["eff_err"][f"{eff}_{sample}_{name}"] = maps[name][f"{sample}_err"].tolist()
        res["eff"][f"{eff}_mc_genmatched_count"] = gen_count.tolist()
        res["eff"][f"{eff}_data_count"] = count_d.tolist()
        res["eff"][f"{eff}_mc_count"] = count_m.tolist()
        truth_p, truth_t = _rebin_pt(mc_raw[f"truth_{eff}_pass"][:, :, None], eff)[..., 0], _rebin_pt(mc_raw[f"truth_{eff}_tot"][:, :, None], eff)[..., 0]
        res["eff"][f"{eff}_mc_truth"] = np.where(truth_t > 0, truth_p / np.maximum(truth_t, 1e-9), 1.0).tolist()
        res["chi2"][eff] = {name: float(np.median(c)) for name, c in chi2.items()}
        res["details"][eff] = details
        print(f"  {eff}: {npt*neta*len(allv)*2} fits in {time.time()-t0:.0f}s; <SF> {sf.mean():.5f} "
              f"<stat> {stat.mean():.5f} <syst> {syst.mean():.5f}; data eff {nom['data'].mean():.5f}, "
              f"MC fit {nom['mc'].mean():.5f}, MC gen-matched {gen_count.mean():.5f}; median data chi2/ndf "
              + ", ".join(f"{n} {c:.2f}" for n, c in res["chi2"][eff].items()), flush=True)
        for name, _ in allv[1:]:
            m = maps[name]
            sf_v = np.where(m["mc"] > 0, m["data"] / np.maximum(m["mc"], 1e-9), 1.0)
            print(f"      variant {name:14s}: <SF> {sf_v.mean():.5f}", flush=True)
        res["meta"][f"{eff}_variant_mean_sf"] = {
            name: float(np.mean(np.where(maps[name]["mc"] > 0, maps[name]["data"] / np.maximum(maps[name]["mc"], 1e-9), 1.0)))
            for name, _ in allv}

    # per-muon reconstruction SF and its per-event effect, weighted by the Z-muon population
    # (passing mutrk probes of the data, 70-110 GeV): the SR muons are these muons
    w = data_raw["mutrk_nom_os_pass"][:, :, tnp.COUNT_WINDOW].sum(axis=2)
    w = w / w.sum()
    mt = MUTRK_NOMINAL
    sf_mu = np.array(res["sf"]["trk"]) * np.array(res["sf"][mt])
    err_mu = sf_mu * np.hypot(np.array(res["sf_err"]["trk"]) / np.array(res["sf"]["trk"]),
                              np.array(res["sf_err"][mt]) / np.array(res["sf"][mt]))
    res["sf"]["reco"] = sf_mu.tolist()
    res["sf_err"]["reco"] = err_mu.tolist()
    mean_sf = float((w * sf_mu).sum())
    mean_err = float((w * err_mu).sum())               # coherent over cells (like MuonID/Iso/Trigger)
    res["meta"].update(population=w.tolist(), sf_reco_per_muon=mean_sf, sf_reco_per_muon_err=mean_err,
                       sf_reco_per_event=mean_sf ** 2, sf_reco_per_event_err=2 * mean_err * mean_sf,
                       sf_trk_per_muon=float((w * np.array(res["sf"]["trk"])).sum()),
                       sf_trk_per_muon_err=float((w * np.array(res["sf_err"]["trk"])).sum()),
                       mutrk_nominal=mt,
                       sf_mutrk_per_muon=float((w * np.array(res["sf"][mt])).sum()),
                       sf_mutrk_per_muon_err=float((w * np.array(res["sf_err"][mt])).sum()),
                       sf_mutrk_loose_per_muon=float((w * np.array(res["sf"]["mutrk"])).sum()),
                       sf_mutrk_loose_per_muon_err=float((w * np.array(res["sf_err"]["mutrk"])).sum()),
                       sf_spread_over_cells=float(np.sqrt((w * (sf_mu - mean_sf) ** 2).sum())))
    print(f"  per muon: SF_reco = {mean_sf:.5f} +- {mean_err:.5f} (tracking {res['meta']['sf_trk_per_muon']:.5f} +- "
          f"{res['meta']['sf_trk_per_muon_err']:.5f}, muon|track {res['meta']['sf_mutrk_per_muon']:.5f} +- "
          f"{res['meta']['sf_mutrk_per_muon_err']:.5f});  per event {mean_sf**2:.5f} +- {2*mean_err*mean_sf:.5f}"
          f"  (assigned before: 1 +- 0.004 per muon, 0.008 per event)")
    return res


def plot(res, data, mc):
    import matplotlib.pyplot as plt
    eta = np.array(res["eta_edges"])
    nice = {"trk": "Tracking (stand-alone probe)", "mutrk": "Muon | track (loose isolated-track probe)",
            "mutrkT": "Muon | track (tight isolated-track probe, nominal)"}
    fig, axes = plt.subplots(1, len(recoeff.EFFS), figsize=(8 * len(recoeff.EFFS), 6.5))
    for ax, eff in zip(axes, recoeff.EFFS):
        pt = np.array(res["pt_edges"][eff])
        c, hw = 0.5 * (pt[1:] + pt[:-1]), 0.5 * np.diff(pt)
        for j in range(len(eta) - 1):
            off = (j - 1.5) * 0.012 * c
            col = ("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")[j]
            e_d = np.array(res["eff"][f"{eff}_data_nominal"])[:, j]
            er_d = np.array(res["eff_err"][f"{eff}_data_nominal"])[:, j]
            e_m = np.array(res["eff"][f"{eff}_mc_nominal"])[:, j]
            ax.errorbar(c + off, e_d, yerr=er_d, xerr=None if eff != "trk" else hw, fmt="o", color=col, markersize=6,
                        label=f"data {eta[j]:.1f}<|η|<{eta[j+1]:.1f}")
            ax.plot(c + off, e_m, "s", markerfacecolor="white", markeredgecolor=col, markersize=6)
        if eff != "trk":
            hists.log_pt_axis(ax)
        ax.set_ylim(0.985, 1.002)
        ax.set_xlabel(r"probe $p_T$ [GeV]" if eff != "trk" else r"probe $p_T$ > 20 GeV (inclusive)")
        ax.set_ylabel("efficiency")
        ax.set_title(nice[eff] + "\nfilled: data, open: DY MC", loc="left", fontsize=13)
        ax.legend(fontsize=9, frameon=False, loc="lower left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle("Muon reconstruction efficiency, tag-and-probe on the unskimmed NanoAOD (same-sign background shapes)",
                 fontsize=14, x=0.01, ha="left")
    fig.tight_layout()
    hists.save_fig(fig, "reco_eff_vs_pt.png")

    pt = recoeff.PT_EDGES
    fig, axes = plt.subplots(1, 3, figsize=(21, 6.5))
    for ax, eff, title in zip(axes, ("trk", MUTRK_NOMINAL, "reco"),
                              ("tracking SF", "muon-given-track SF", "reconstruction SF (product)")):
        sf, err = np.array(res["sf"][eff]), np.array(res["sf_err"][eff])
        mesh = ax.pcolormesh(pt, eta, sf.T, cmap="RdBu", vmin=0.99, vmax=1.01)
        for i in range(len(pt) - 1):
            for j in range(len(eta) - 1):
                ax.text(np.sqrt(pt[i] * pt[i + 1]), 0.5 * (eta[j] + eta[j + 1]), f"{sf[i,j]:.4f}\n±{err[i,j]:.4f}",
                        ha="center", va="center", fontsize=8)
        fig.colorbar(mesh, ax=ax)
        hists.log_pt_axis(ax)
        ax.set_xlabel(r"$p_T$ [GeV]")
        ax.set_ylabel(r"|$\eta$|")
        ax.set_title(title, fontsize=14)
    m = res["meta"]
    fig.suptitle(f"Muon reconstruction scale factor: per muon {m['sf_reco_per_muon']:.4f} ± {m['sf_reco_per_muon_err']:.4f}"
                 f" (was 1 ± 0.004 assigned)", fontsize=15)
    fig.tight_layout()
    hists.save_fig(fig, "reco_sf_maps.png")

    for eff in recoeff.EFFS:
        det = res["details"][eff]
        i = 0 if eff == "trk" else 1
        pt_e = res["pt_edges"][eff]
        fig, axes = plt.subplots(4, 4, figsize=(20, 16))
        for j in range(4):
            for col, (sample, h) in enumerate((("data", data), ("mc", mc))):
                d = det[(sample, i, j)]
                for pf, ax in (("pass", axes[j, 2 * col]), ("fail", axes[j, 2 * col + 1])):
                    obs = h[f"{eff}_nom_os_{pf}"]
                    obs = (obs.sum(axis=0) if eff == "trk" else obs[i])[j]
                    ax.errorbar(tnp.MASS_CENTRES, obs, yerr=np.sqrt(np.maximum(obs, 0)), fmt="o", color="black", markersize=2)
                    ax.plot(d["x"], d[f"model_{pf}"], color="crimson")
                    ax.set_title(f"{sample} {pf}, {pt_e[i]:.0f}-{pt_e[i+1]:.0f} GeV, {eta[j]:.1f}<|η|<{eta[j+1]:.1f}: "
                                 f"ε={d['eps']:.5f}±{d['err']:.5f}, χ²/ndf={d['chi2']:.0f}/{d['ndf']}", fontsize=8)
                    ax.set_xlabel(r"$m$(tag, probe) [GeV]")
        fig.suptitle(f"{nice[eff]}: pass/fail fits with same-sign background shapes", fontsize=15)
        fig.tight_layout()
        hists.save_fig(fig, f"reco_fits_{eff}.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=None, help="data files per era")
    ap.add_argument("--mc-files", type=int, default=10, help="DY NLO parent files (of 41)")
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--plots-only", action="store_true")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--is-mc", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--pu-json", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:
        batch.write_part(process_file(args.one_file, args.is_mc, args.pu_json), args.out)
        return
    parts = OUT / "reco_parts"
    if args.plots_only:
        saved = pickle.load(open(OUT / "reco_fits.pkl", "rb"))
        plot(saved["result"], saved["data"], saved["mc"])
        return
    pu_json = OUT / "pileup_weights.json"
    if not args.summarise_only:
        for key, n in (("DY_NLO", args.mc_files), ("data_2016G", args.max_files), ("data_2016H", args.max_files)):
            srcs = samples.sources(key)[:n]
            alt = {Path(s["name"]).stem: s["fallback"] for s in srcs}
            tasks = [(s["primary"], Path(s["name"]).stem) for s in srcs]
            extra = ["--is-mc", "--pu-json", str(pu_json)] if samples.SAMPLES[key]["is_mc"] else []
            print(f"[reco] {key}: {len(tasks)} files", flush=True)
            batch.run_files(tasks, Path(__file__), parts / key, args.workers, extra_args=extra, stall_s=3600,
                            retries=1, fallback=lambda src, k: alt.get(k), log=lambda s: print(s, flush=True))
    data = recoeff.blank()
    for key in ("data_2016G", "data_2016H"):
        if (parts / key).exists():
            data = batch.load_parts(parts / key, data)
    mc = batch.load_parts(parts / "DY_NLO", recoeff.blank())
    print(f"[reco] data: {data['n_files']} files, {data['n_events']:,} events, probes {data['n_probes']}; "
          f"MC: {mc['n_files']} files, probes {mc['n_probes']}")
    res = measure(data, mc)
    with open(OUT / "reco_fits.pkl", "wb") as fh:
        pickle.dump({"result": res, "data": data, "mc": mc}, fh)
    with open(OUT / "reco_result.json", "w") as fh:
        json.dump({k: v for k, v in res.items() if k != "details"}, fh, indent=1)
    plot(res, data, mc)
    print(f"[reco] wrote {OUT / 'reco_result.json'}")


if __name__ == "__main__":
    main()
