#!/usr/bin/env python
"""Step 1 -- skim the Tau data and the simulation (docs/03-skims.md).

    python scripts/step1_skim.py                              # every sample in ztautau/samples.py
    python scripts/step1_skim.py --samples DY_NLO --max-files 1
    python scripts/step1_skim.py --manifest-only
    python scripts/step1_skim.py --v4 --samples data_mu_2016G   # v4 four-channel skim -> skims_v4/ (docs/10-v4-plan.md)

Each parent file runs in its own subprocess (ztautau/batch.py): resumable; a stalled dCache read is
killed after --stall-s seconds and retried from EOS. Output: $BND_TAUTAU_CACHE/skims_v1/<sample>/
<file>.root with a <file>.root.json provenance sidecar, and skims_v1/manifest.json.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ztautau import batch, config, samples, skim  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--workers", type=int, default=config.N_WORKERS)
    ap.add_argument("--prefer", choices=["dcache", "eos"], default="dcache")
    ap.add_argument("--stall-s", type=float, default=3000)
    ap.add_argument("--manifest-only", action="store_true")
    ap.add_argument("--v4", action="store_true", help="four-channel preselection (skim.v4_categories) into skims_v4/")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--sample", help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:                                   # worker mode
        skim.process_file(args.one_file, args.out, samples.SAMPLES[args.sample], v4=args.v4)
        return

    skim_dir = config.SKIM_DIR_V4 if args.v4 else config.SKIM_DIR
    keys = args.samples or (samples.V4_SKIM_KEYS if args.v4 else [k for k in samples.SAMPLES if not samples.SAMPLES[k]["v4_only"]])
    if not args.manifest_only:
        for key in keys:
            srcs = samples.sources(key, prefer=args.prefer, max_files=args.max_files)
            tasks = [(s["primary"], s["stem"]) for s in srcs]
            alt = {s["stem"]: s["fallback"] for s in srcs}
            print(f"[skim] {key}: {len(tasks)} files -> {skim_dir / key}", flush=True)
            t0 = time.time()
            notes = batch.run_files(tasks, Path(__file__), skim_dir / key, args.workers,
                                    extra_args=["--sample", key] + (["--v4"] if args.v4 else []), stall_s=args.stall_s, retries=2,
                                    fallback=lambda src, k: alt.get(k), suffix=".root",
                                    log=lambda s: print(s, flush=True))
            for n in notes:
                print("        " + n)
            print(f"[skim] {key} done in {time.time() - t0:.0f}s", flush=True)
    manifest = skim.build_manifest([k for k in keys if skim.skim_files(k, skim_dir)], skim_dir, v4=args.v4)
    print(f"\n{'sample':18s} {'files':>9s} {'n_in':>12s} {'n_out':>10s} {'GB':>6s}")
    for key, e in manifest["samples"].items():
        print(f"{key:18s} {e['n_files']:>4d}/{e['n_files_parent']:<4d} {e['n_in']:>12,d} {e['n_out']:>10,d} "
              f"{e['bytes'] / 1e9:>6.2f}")


if __name__ == "__main__":
    main()
