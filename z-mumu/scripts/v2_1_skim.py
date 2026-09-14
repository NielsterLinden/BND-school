#!/usr/bin/env python
"""v2 step 1 -- skim the NanoAOD parents (data from dCache/EOS, MC from dCache or EOS).

    python scripts/v2_1_skim.py                          # every sample in zmumu/samples.py
    python scripts/v2_1_skim.py --samples DY_NLO data_2016G --max-files 2
    python scripts/v2_1_skim.py --manifest-only          # rebuild manifest.json

Each file runs in its own subprocess (zmumu/batch.py): resumable, a stalled dCache read is
killed after --stall-s seconds and the file is retried from EOS. Output:
$BND_SKIM_DIR/<sample>/<file>.root (+ .json provenance), and $BND_SKIM_DIR/manifest.json.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from zmumu import batch, samples, skim


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", nargs="*", default=None, help="sample keys (default: all)")
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--prefer", choices=["dcache", "eos"], default="dcache")
    ap.add_argument("--stall-s", type=float, default=2400)
    ap.add_argument("--out-dir", type=Path, default=None, help="default: $BND_SKIM_DIR")
    ap.add_argument("--manifest-only", action="store_true")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--sample", help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:                      # worker mode
        sample = dict(samples.SAMPLES[args.sample], key=args.sample)
        skim.process_file(args.one_file, args.out, sample)
        return

    base = skim.skim_dir(args.out_dir)
    keys = args.samples or list(samples.SAMPLES)
    if not args.manifest_only:
        for key in keys:
            srcs = samples.sources(key, prefer=args.prefer, max_files=args.max_files)
            tasks = [(s["primary"], Path(s["name"]).stem) for s in srcs]
            alt = {Path(s["name"]).stem: s["fallback"] for s in srcs}
            print(f"[skim] {key}: {len(tasks)} files -> {base / key}", flush=True)
            t0 = time.time()
            notes = batch.run_files(tasks, Path(__file__), base / key, args.workers,
                                    extra_args=["--sample", key], stall_s=args.stall_s, retries=2,
                                    fallback=lambda src, k: alt.get(k), suffix=".root",
                                    log=lambda s: print(s, flush=True))
            for n in notes:
                print("        " + n)
            print(f"[skim] {key} done in {time.time()-t0:.0f}s", flush=True)
    manifest = skim.build_manifest(keys, base)
    print(f"\n{'sample':16s} {'files':>5s} {'n_in':>12s} {'n_out':>10s} {'GB':>6s}")
    for key, e in manifest["samples"].items():
        print(f"{key:16s} {e['n_files']:>5d} {e['n_in']:>12,d} {e['n_out']:>10,d} {e['bytes']/1e9:>6.2f}")


if __name__ == "__main__":
    main()
