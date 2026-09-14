#!/usr/bin/env python
"""v2 step 3 -- control regions: muon momentum calibration and the fake-factor estimate.

One pass over the data and MC skims fills (a) the Z-peak mass spectra per |eta| category
(data, DY) for the momentum scale/resolution calibration and (b) the fake-factor histograms
(measurement regions, application regions, tight-tight yields) with the prompt-prompt MC
contribution of every sample. The summary fits the peaks, derives the fake factor maps and
the fake template of the signal region, and makes the control plots.

    python scripts/v2_3_control.py [--workers 8] [--summarise-only]

Outputs: output/v2/momentum.json, output/v2/fakes.json, output/v2/fakes_templates.pkl,
output/v2/plots/momentum_*.png, output/v2/plots/ff_*.png. Docs: 13-fake-factor.md, 14-fit-and-systematics.md.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import numpy as np

from zmumu import batch, config, fakes, hists, momentum, pileup, regions, samples, skim, weights

OUT = config.OUTPUT_DIR / "v2"
config.PLOT_DIR = OUT / "plots"          # v2 plots live next to the v2 results
BRANCHES = ["run", "event", "PV_npvsGood", "HLT_IsoMu24", "HLT_IsoTkMu24", "Muon_*", "TrigObj_*", "Flag_*",
            "FsrPhoton_*", "genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "skim_cat", "skim_prescale", "gen_lhe_flavour"]
MC_KEYS = [k for k in samples.MC_KEYS if k != "DY_powheg"]


def blank():
    out = {"n_files": 0, "mom": momentum.blank(), "ff": fakes.blank()}
    return out


def process_file(source, key, is_mc, norm, pu_json):
    import uproot
    out = blank()
    out["n_files"] = 1
    pu = pileup.PileupWeights.from_json(pu_json) if is_mc else None
    with uproot.open(source) as f:
        if "Events" not in f:
            return out
        for ev in f["Events"].iterate(filter_name=BRANCHES, step_size="200 MB", library="ak"):
            if is_mc:
                w = (np.asarray(ev.genWeight, dtype=float) * norm * pu(np.asarray(ev.Pileup_nTrueInt))
                     * np.asarray(ev.L1PreFiringWeight_Nom, dtype=float))
            else:
                w = np.asarray(ev.skim_prescale, dtype=float)
            if key in ("data_2016G", "data_2016H", "DY_NLO"):
                keep = regions.fired(ev) & regions.event_clean(ev)
                e = ev[keep]
                if len(e):
                    masks = regions.muon_masks(e)
                    sr = regions.dimuon_regions(e, masks, regions.trigger_match(e), np.ones(len(e), bool))["SR"]
                    if sr is not None:
                        ww = w[keep][sr["idx"]]
                        if is_mc:      # DY -> mumu only for the peak
                            ww = np.where(np.asarray(e.gen_lhe_flavour)[sr["idx"]] == 13, ww, 0.0)
                        momentum.fill(out["mom"], sr, ww)
            fakes.fill(ev, out["ff"], w, is_mc, prompt_only=True)
    return out


def plot_momentum(cal, data, mc):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    eta = cal["eta_edges"]
    for j, ax in enumerate(axes.flat):
        fd, fm = cal["fits"][j]["data"], cal["fits"][j]["mc"]
        hd, hm = data[f"m_eta{j}"], mc[f"m_eta{j}"]
        scale = hd.sum() / max(hm.sum(), 1e-9)
        ax.errorbar(momentum.CENTRES, hd, yerr=np.sqrt(np.maximum(hd, 0)), fmt="o", color="black", markersize=3, label="data")
        ax.plot(momentum.CENTRES, np.array(fd["model"]), color="black", label=f"fit: peak {fd['mu']:.2f}, σ {fd['sigma']:.2f}")
        ax.step(momentum.CENTRES, hm * scale, where="mid", color="#1f77b4", label="DY MC (scaled)")
        ax.plot(momentum.CENTRES, np.array(fm["model"]) * scale, color="#1f77b4", linestyle="--",
                label=f"fit: peak {fm['mu']:.2f}, σ {fm['sigma']:.2f}")
        ax.set_title(f"both muons in {eta[j]:.1f} < |η| < {eta[j+1]:.1f}: κ = {cal['kappa'][j]:.5f}, smear {100*cal['smear'][j]:.2f}%", fontsize=11)
        ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.legend(fontsize=9)
    fig.suptitle("Z peak position and width, data vs simulation (momentum calibration)", fontsize=15)
    fig.tight_layout(); hists.save_fig(fig, "momentum_zpeak_fits.png")


def plot_fakes(res, data, mc):
    import matplotlib.pyplot as plt
    pt, eta = np.array(res["pt_edges"]), np.array(res["eta_edges"])
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    for ax, key, title in ((axes[0], "nominal", "FF, same-sign tag+probe region"), (axes[1], "region_c", "FF, single-muon + jet region")):
        ff = np.array(res["maps"][key])
        mesh = ax.pcolormesh(pt, eta, ff.T, cmap="viridis", vmin=0, vmax=max(ff.max(), 0.05))
        for i in range(len(pt) - 1):
            for j in range(len(eta) - 1):
                ax.text(np.sqrt(pt[i] * pt[i + 1]), 0.5 * (eta[j] + eta[j + 1]), f"{ff[i,j]:.3f}", ha="center", va="center", fontsize=8, color="white")
        fig.colorbar(mesh, ax=ax, label="fake factor N(tight)/N(anti-tight)")
        hists.log_pt_axis(ax); ax.set_xlabel(r"muon $p_T$ [GeV]"); ax.set_ylabel(r"|$\eta$|"); ax.set_title(title, fontsize=12)
    fig.tight_layout(); hists.save_fig(fig, "ff_maps.png")
    c = 0.5 * (fakes.MASS_EDGES[1:] + fakes.MASS_EDGES[:-1])
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    ax = axes[0]
    ax.errorbar(c, res["closure"]["observed_hist"], yerr=np.sqrt(np.maximum(data["tt_ss_w2"], 0)), fmt="o", color="black", markersize=3, label="SS tight-tight, data − prompt MC")
    ax.step(c, res["closure"]["predicted_hist"], where="mid", color="crimson", label="FF × SS (tight + anti-tight)")
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("events / 0.5 GeV"); ax.legend(fontsize=10)
    ax.set_title(f"same-sign closure: predicted {res['closure']['predicted_ss']:.0f} ± {res['closure']['predicted_ss_err']:.0f}, observed {res['closure']['observed_ss']:.0f} ± {res['closure']['observed_ss_err']:.0f}", fontsize=10)
    ax = axes[1]
    for key, col in (("nominal", "black"), ("region_c", "#1f77b4"), ("anti_alt", "#2ca02c"), ("mcsub_up", "grey"), ("mcsub_down", "grey")):
        ax.step(c, res["templates"][key], where="mid", color=col, label=f"{key}: {np.sum(res['templates'][key]):.0f}")
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("fake events / 0.5 GeV"); ax.legend(fontsize=9)
    ax.set_title("non-prompt template in the signal region", fontsize=12)
    fig.tight_layout(); hists.save_fig(fig, "ff_closure_and_template.png")
    # application region composition
    fig, ax = plt.subplots(figsize=(9, 7))
    d = data["app_os_anti"].sum(axis=(0, 1)); m = mc["app_os_anti"].sum(axis=(0, 1))
    ax.errorbar(c, d, yerr=np.sqrt(np.maximum(d, 0)), fmt="o", color="black", markersize=3, label="data: tight + anti-tight, OS")
    ax.step(c, m, where="mid", color="#1f77b4", label="prompt-prompt MC")
    ax.step(c, d - m, where="mid", color="crimson", label="difference (non-prompt)")
    ax.set_yscale("log"); ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]"); ax.set_ylabel("events / 0.5 GeV"); ax.legend(fontsize=11)
    hists._decorate(ax); hists._title(ax, "Fake-factor application region")
    hists.save_fig(fig, "ff_application_region.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--skim-dir", default=None)
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--sample", help=argparse.SUPPRESS)
    ap.add_argument("--norm", type=float, default=1.0, help=argparse.SUPPRESS)
    ap.add_argument("--pu-json", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:
        is_mc = samples.SAMPLES[args.sample]["is_mc"]
        batch.write_part(process_file(args.one_file, args.sample, is_mc, args.norm, args.pu_json), args.out)
        return

    OUT.mkdir(parents=True, exist_ok=True)
    pu_json = OUT / "tnp" / "pileup_weights.json"
    keys = args.samples or (samples.DATA_KEYS + MC_KEYS)
    if not args.summarise_only:
        if not pu_json.exists():
            pu = pileup.PileupWeights(pileup.mc_profile(("DY_NLO",), args.skim_dir)); pu_json.parent.mkdir(parents=True, exist_ok=True); pu.to_json(pu_json)
        for key in keys:
            files = skim.skim_files(key, args.skim_dir)[: args.max_files]
            if not files:
                print(f"[control] {key}: no skim files, skipped"); continue
            norm = 1.0
            if samples.SAMPLES[key]["is_mc"]:
                w = weights.Weighter(key, None, base=args.skim_dir); norm = w.norm
            tasks = [(str(f), f"{key}__{f.stem}") for f in files]
            print(f"[control] {key}: {len(tasks)} files, norm {norm:.4g}", flush=True)
            batch.run_files(tasks, Path(__file__), OUT / "control_parts" / key, args.workers,
                            extra_args=["--sample", key, "--norm", str(norm), "--pu-json", str(pu_json)],
                            stall_s=1800, retries=1, log=lambda s: print(s, flush=True))
    totals = {}
    for key in samples.DATA_KEYS + MC_KEYS:
        d = OUT / "control_parts" / key
        if d.exists() and any(d.glob("*.pkl")):
            totals[key] = batch.load_parts(d, blank())
    data = blank()
    for k in samples.DATA_KEYS:
        if k in totals:
            batch.merge(data, totals[k])
    mc_all = blank()
    for k in MC_KEYS:
        if k in totals:
            batch.merge(mc_all, totals[k])
    print("[control] momentum calibration (data vs DY NLO, Z -> mumu)")
    cal = momentum.calibrate(data["mom"], totals["DY_NLO"]["mom"])
    with open(OUT / "momentum.json", "w") as fh:
        json.dump({k: v for k, v in cal.items() if k != "fits"} | {"fits": {str(j): {s: {kk: vv for kk, vv in f.items() if kk != "model"} for s, f in v.items()} for j, v in cal["fits"].items()}}, fh, indent=1)
    plot_momentum(cal, data["mom"], totals["DY_NLO"]["mom"])
    antiiso_sf = None
    tnp_json = OUT / "tnp" / "tnp_result.json"
    if tnp_json.exists():
        t = json.load(open(tnp_json))
        # map the T&P anti-iso SF (T&P binning) onto the FF binning by cell centre
        sf_t = np.array(t["sf"]["antiiso"]); pte, ete = np.array(t["pt_edges"]), np.array(t["eta_edges"])
        antiiso_sf = np.zeros((len(fakes.PT_EDGES) - 1, len(fakes.ETA_EDGES) - 1))
        for i in range(antiiso_sf.shape[0]):
            for j in range(antiiso_sf.shape[1]):
                pc = np.sqrt(fakes.PT_EDGES[i] * fakes.PT_EDGES[i + 1]); ec = 0.5 * (fakes.ETA_EDGES[j] + fakes.ETA_EDGES[j + 1])
                antiiso_sf[i, j] = sf_t[min(np.digitize(pc, pte) - 1, len(pte) - 2), min(np.digitize(ec, ete) - 1, len(ete) - 2)]
        print(f"[control] anti-isolation SF from T&P applied to the prompt subtraction: mean {antiiso_sf.mean():.3f}")
    res = fakes.summarise(data["ff"], mc_all["ff"], antiiso_sf)
    print(f"[control] fakes in SR: {res['yields']['sr_fakes']:.0f} ± {res['yields']['sr_fakes_stat']:.0f} (stat), "
          f"method unc {100*res['yields']['method_rel_unc']:.0f}% (worst: {res['yields']['method_worst']}); "
          f"variants {{{', '.join(f'{k}: {v:.0f}' for k, v in res['yields']['variants'].items())}}}")
    print(f"[control] SS closure: predicted {res['closure']['predicted_ss']:.0f} ± {res['closure']['predicted_ss_err']:.0f}, "
          f"observed {res['closure']['observed_ss']:.0f} ± {res['closure']['observed_ss_err']:.0f}; data SS tight-tight {res['yields']['ss_data_tight_tight']:.0f}")
    with open(OUT / "fakes.json", "w") as fh:
        json.dump(res, fh, indent=1)
    with open(OUT / "fakes_templates.pkl", "wb") as fh:
        pickle.dump({"result": res, "data": data["ff"], "mc": mc_all["ff"]}, fh)
    plot_fakes(res, data["ff"], mc_all["ff"])


if __name__ == "__main__":
    main()
