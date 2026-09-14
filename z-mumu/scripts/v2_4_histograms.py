#!/usr/bin/env python
"""v2 step 4 -- fill every histogram (all samples, regions, variables, variations).

    python scripts/v2_4_histograms.py [--workers 8] [--samples KEY ...] [--merge-only]

Needs: the skims, output/v2/tnp/tnp_result.json (scale factors), output/v2/momentum.json.
Writes: output/v2/hist_parts/<sample>/*.pkl and output/v2/histograms.pkl (merged, keyed by
"<fit sample>|<region>|<variable>|<variation>"), plus the data/MC comparison plots.
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

from zmumu import batch, config, histograms as H, momentum, pileup, samples, skim, tnp, weights

OUT = config.OUTPUT_DIR / "v2"
BRANCHES = ["run", "event", "PV_npvsGood", "MET_pt", "HLT_IsoMu24", "HLT_IsoTkMu24", "Muon_*", "TrigObj_*",
            "Flag_*", "FsrPhoton_*", "Jet_pt", "Jet_eta", "Jet_jetId", "Electron_*", "genWeight", "Pileup_nTrueInt",
            "L1PreFiringWeight_*", "skim_cat", "skim_prescale", "gen_lhe_flavour", "LHEPdfWeight", "LHEScaleWeight",
            "PSWeight"]


class _FixedWeighter(weights.Weighter):
    """Weighter whose GenSums come from a JSON (avoids re-reading every skim file per worker)."""

    def __init__(self, key, gensums, pu):
        self.key = key
        self.sample = samples.SAMPLES[key]
        self.gensums = gensums
        self.pileup = pu
        self.sf = None
        sumw = float(gensums.get("sumw", 0.0)) if gensums else 0.0
        xsec = self.sample["xsec_pb"]
        self.norm = (xsec * weights.LUMI_PB / sumw) if (xsec and sumw) else 1.0


def process_file(source, key, aux_dir):
    import uproot
    s = samples.SAMPLES[key]
    is_mc = s["is_mc"]
    out = {"n_files": 1}
    aux = Path(aux_dir)
    pu = pileup.PileupWeights.from_json(aux / "tnp" / "pileup_weights.json") if is_mc else None
    gens = json.load(open(aux / "gensums.json")).get(key) if is_mc else None
    gens = {k: (np.asarray(v) if isinstance(v, list) else v) for k, v in gens.items()} if gens else None
    weighter = _FixedWeighter(key, gens, pu) if is_mc else None
    sf = tnp.ScaleFactors(aux / "tnp" / "tnp_result.json") if (is_mc and (aux / "tnp" / "tnp_result.json").exists()) else None
    calib = momentum.MomentumCalibration(aux / "momentum.json") if (is_mc and (aux / "momentum.json").exists()) else None
    fit_sample = s["fit_sample"]
    split = bool(s.get("split_lhe"))
    theory = key == "DY_NLO"
    variations = key != "DY_powheg"
    with uproot.open(source) as f:
        if "Events" not in f:
            return out
        for ev in f["Events"].iterate(filter_name=BRANCHES, step_size="150 MB", library="ak"):
            H.fill_chunk(ev, out, key, is_mc, weighter, sf, calib, split_flavour=split, theory=theory,
                         fit_sample=fit_sample, variations=variations)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument("--skim-dir", default=None)
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--sample", help=argparse.SUPPRESS)
    ap.add_argument("--aux", default=str(OUT), help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:
        batch.write_part(process_file(args.one_file, args.sample, args.aux), args.out)
        return

    keys = args.samples or list(samples.SAMPLES)
    if not args.merge_only:
        gens = {}
        for key in keys:
            if samples.SAMPLES[key]["is_mc"]:
                g = skim.load_gensums(key, args.skim_dir)
                gens[key] = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in g.items()}
        with open(OUT / "gensums.json", "w") as fh:
            json.dump(gens, fh)
        for key in keys:
            files = skim.skim_files(key, args.skim_dir)[: args.max_files]
            if not files:
                print(f"[hist] {key}: no skim files, skipped"); continue
            tasks = [(str(f), f"{key}__{f.stem}") for f in files]
            print(f"[hist] {key}: {len(tasks)} files", flush=True)
            batch.run_files(tasks, Path(__file__), OUT / "hist_parts" / key, args.workers,
                            extra_args=["--sample", key, "--aux", str(OUT)], stall_s=2400, retries=1,
                            log=lambda s: print(s, flush=True))
    total = {}
    for key in list(samples.SAMPLES):
        d = OUT / "hist_parts" / key
        if d.exists() and any(d.glob("*.pkl")):
            part = batch.load_parts(d, {})
            part.pop("n_files", None)
            batch.merge(total, part)
    with open(OUT / "histograms.pkl", "wb") as fh:
        pickle.dump(total, fh)
    yields = {}
    for k, v in total.items():
        smp, region, var, variation = k.split("|")[:4]
        if var == "mass_fit" and variation == "nominal" and not k.endswith("|w2"):
            yields[(region, smp)] = float(np.sum(v))
    for region in ("SR", "SS", "CRemu"):
        print(f"  {region}: " + ", ".join(f"{s} {yields[(r, s)]:,.1f}" for (r, s) in sorted(yields) if r == region))
    print(f"[hist] wrote {OUT / 'histograms.pkl'} ({len(total)} histograms)")


if __name__ == "__main__":
    main()
