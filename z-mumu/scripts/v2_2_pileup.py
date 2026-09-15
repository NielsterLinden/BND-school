#!/usr/bin/env python
"""v2 step 2a -- pileup profile and weights (N_PV-matched, docs/11-mc-weights.md).

    python scripts/v2_2_pileup.py [--workers 8] [--max-files N]

Fills, in the Z -> mu mu signal region, the (true pileup x N_PV) distribution of the DY
simulation (weighted with everything except the pileup weight) and the N_PV distribution of
the data, then chooses the scale and bunch-to-bunch smearing of the luminosity-CSV profile
that make the reweighted simulation reproduce the data (`pileup.match_npv`).
Writes output/v2/tnp/pileup_weights.json (read by the T&P, control and histogram steps)
and the plots pu_profile_data_mc.png, pu_npv_matching.png.
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

from zmumu import batch, config, hists, pileup, regions, samples, skim, weights

OUT = config.OUTPUT_DIR / "v2" / "tnp"
config.PLOT_DIR = config.OUTPUT_DIR / "v2" / "plots"
BRANCHES = ["PV_npvsGood", "HLT_IsoMu24", "HLT_IsoTkMu24", "Muon_*", "TrigObj_*", "Flag_*", "FsrPhoton_*",
            "genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "skim_prescale"]
NPV_EDGES = np.arange(0.0, 61.0)
NTRUE_EDGES = np.arange(0.0, pileup.NBINS + 1.0)


def process_file(source, is_mc):
    import awkward as ak
    import uproot
    out = {"n_files": 1, "mc_ntrue_npv": np.zeros((pileup.NBINS, len(NPV_EDGES) - 1)), "data_npv": np.zeros(len(NPV_EDGES) - 1)}
    with uproot.open(source) as f:
        if "Events" not in f:
            return out
        for ev in f["Events"].iterate(filter_name=BRANCHES, step_size="150 MB", library="ak"):
            ev = ev[regions.fired(ev) & regions.event_clean(ev)]
            if len(ev) == 0:
                continue
            masks = regions.muon_masks(ev)
            matched = regions.trigger_match(ev)
            sr = regions.dimuon_regions(ev, masks, matched, np.ones(len(ev), dtype=bool))["SR"]
            if sr is None or len(sr["idx"]) == 0:
                continue
            idx = sr["idx"]
            npv = np.asarray(ev.PV_npvsGood, dtype=float)[idx]
            if is_mc:
                w = (np.asarray(ev.genWeight, dtype=float) * np.asarray(ev.L1PreFiringWeight_Nom, dtype=float))[idx]
                pp = regions.is_prompt(sr["flav1"]) & regions.is_prompt(sr["flav2"])
                ntrue = np.asarray(ev.Pileup_nTrueInt, dtype=float)[idx]
                out["mc_ntrue_npv"] += np.histogram2d(ntrue[pp], npv[pp], bins=[NTRUE_EDGES, NPV_EDGES], weights=w[pp])[0]
            else:
                out["data_npv"] += np.histogram(npv, bins=NPV_EDGES, weights=np.asarray(ev.skim_prescale, dtype=float)[idx])[0]
    return out


def plot(pu, match, mc_h2, data_npv):
    import matplotlib.pyplot as plt
    x = np.arange(pileup.NBINS) + 0.5
    fig, ax = plt.subplots(figsize=(9, 6))
    for name, label, style in (("csv_raw", "data, luminosity CSV (69.2 mb)", dict(color="grey", linestyle="--")),
                               ("nominal", f"data, N$_{{PV}}$-matched (scale {pu.scale:.3f}, smear {pu.rel_smear:.2f})", dict(color="black")),
                               ("up", "sigma_mb +4.6%", dict(color="tab:red", linewidth=0.8)), ("down", "sigma_mb -4.6%", dict(color="tab:blue", linewidth=0.8))):
        h = pu.data[name]["hist"]
        ax.step(x, h / h.sum(), where="mid", label=label, **style)
    ax.step(x, pu.mc_norm, where="mid", color="tab:orange", label="MC (Pileup_nTrueInt, all generated)")
    ax.set_xlim(0, 70); ax.set_xlabel("true number of interactions"); ax.set_ylabel("fraction")
    ax.legend(fontsize=10); hists._decorate(ax, lumi_fb=config.LUMI_PB_NORMTAG / 1000); hists._title(ax, "Pileup profiles")
    hists.save_fig(fig, "pu_profile_data_mc.png")
    c = 0.5 * (NPV_EDGES[1:] + NPV_EDGES[:-1])
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})
    d = data_npv / data_npv.sum()
    ax.errorbar(c, d, yerr=np.sqrt(data_npv) / data_npv.sum(), fmt="o", color="black", markersize=3, label="data (SR)")
    for name, label, colour in (("unweighted", "MC, no pileup weight", "tab:orange"), ("csv", "MC, CSV profile", "grey"), ("matched", "MC, N$_{PV}$-matched profile", "tab:blue")):
        if name == "unweighted":
            pred = mc_h2.sum(axis=0)
        elif name == "csv":
            pred = (pileup._weights(pu.data["csv_raw"]["hist"], pu.mc_norm)[:, None] * mc_h2).sum(axis=0)
        else:
            pred = (pu.weights["nominal"][:, None] * mc_h2).sum(axis=0)
        pred = pred / pred.sum()
        ax.step(NPV_EDGES[:-1], pred, where="post", color=colour, label=label)
        with np.errstate(divide="ignore", invalid="ignore"):
            rax.step(NPV_EDGES[:-1], np.where(pred > 0, d / pred, np.nan), where="post", color=colour)
    ax.set_ylabel("fraction of events"); ax.legend(loc="upper right", fontsize=11); ax.set_xlim(0, 50); ax.set_ylim(0, 1.45 * d.max())
    hists._decorate(ax, lumi_fb=config.LUMI_PB_NORMTAG / 1000)
    hists._title(ax, f"Good primary vertices, Z -> mumu SR;  chi2/ndf {match['chi2_unmatched']:.0f} -> {match['chi2']:.0f} / {match['ndf']}")
    rax.axhline(1, color="grey", linestyle="--"); rax.set_ylim(0.7, 1.3); rax.set_ylabel("data / MC"); rax.set_xlabel("good primary vertices")
    hists.save_fig(fig, "pu_npv_matching.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--skim-dir", default=None)
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--is-mc", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.one_file:
        batch.write_part(process_file(args.one_file, args.is_mc), args.out)
        return
    OUT.mkdir(parents=True, exist_ok=True)
    if not args.summarise_only:
        for key in ("data_2016G", "data_2016H", "DY_NLO"):
            files = skim.skim_files(key, args.skim_dir)[: args.max_files]
            if not files:
                print(f"[pileup] {key}: no skim files, skipped"); continue
            tasks = [(str(f), f"{key}__{f.stem}") for f in files]
            print(f"[pileup] {key}: {len(tasks)} skim files", flush=True)
            batch.run_files(tasks, Path(__file__), OUT / "pileup_parts" / key, args.workers,
                            extra_args=["--is-mc"] if samples.SAMPLES[key]["is_mc"] else [], stall_s=1800, retries=1,
                            log=lambda s: print(s, flush=True))
    total = {}
    for key in ("data_2016G", "data_2016H", "DY_NLO"):
        d = OUT / "pileup_parts" / key
        if d.exists():
            batch.merge(total, batch.load_parts(d, {}))
    mc_hist = pileup.mc_profile(("DY_NLO",), args.skim_dir)
    match = pileup.match_npv(total["mc_ntrue_npv"], total["data_npv"], mc_hist)
    pu = pileup.PileupWeights(mc_hist, scale=match["scale"], rel_smear=match["rel_smear"])
    pu.to_json(OUT / "pileup_weights.json", extra={"npv_match": match})
    with open(OUT / "pileup_npv.pkl", "wb") as fh:
        pickle.dump(total, fh)
    npv_mean = lambda h: float(np.sum(h * (np.arange(len(h)) + 0.5)) / h.sum())
    mc_matched = (pu.weights["nominal"][:, None] * total["mc_ntrue_npv"]).sum(axis=0)
    mc_csv = (pileup._weights(pu.data["csv_raw"]["hist"], pu.mc_norm)[:, None] * total["mc_ntrue_npv"]).sum(axis=0)
    print(f"[pileup] N_PV match: scale {match['scale']:.3f}, rel. smear {match['rel_smear']:.2f}; chi2/ndf "
          f"{match['chi2_unmatched']:.0f} -> {match['chi2']:.0f} / {match['ndf']}")
    print(f"[pileup] profile <mu>: CSV {pu.data['csv_raw']['mean']:.2f} (rms {pu.data['csv_raw']['rms']:.2f}) -> matched "
          f"{pu.data['nominal']['mean']:.2f} (rms {pu.data['nominal']['rms']:.2f}); MC {npv_mean(pu.mc):.2f}")
    print(f"[pileup] <N_PV> in the SR: data {npv_mean(total['data_npv']):.2f}, MC unweighted {npv_mean(total['mc_ntrue_npv'].sum(0)):.2f}, "
          f"CSV profile {npv_mean(mc_csv):.2f}, matched {npv_mean(mc_matched):.2f}")
    plot(pu, match, total["mc_ntrue_npv"], total["data_npv"])
    print(f"[pileup] wrote {OUT / 'pileup_weights.json'}")


if __name__ == "__main__":
    main()
