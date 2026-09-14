#!/usr/bin/env python
"""Laptop bundle: one ROOT file per sample from the per-file v2 skims, plus a manifest.

    source ../setup.sh && python scripts/v2_bundle.py [--out DIR] [--samples KEY ...]

Uses `hadd` from the LCG view (no python ROOT import). The `GenSums` and `Runs` trees are
concatenated as well, so the per-sample sums are the sums over the entries of those trees.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from zmumu import samples, skim

DEFAULT_OUT = Path("/project/atlas/users/nterlind/BND-school-skims-lite")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--skim-dir", default=None)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = {"skim_version": skim.SKIM_VERSION, "samples": {}}
    for key in args.samples or list(samples.SAMPLES):
        files = skim.skim_files(key, args.skim_dir)
        if not files:
            continue
        target = args.out / f"{key}.root"
        if not target.exists():
            print(f"[bundle] {key}: hadd {len(files)} files", flush=True)
            subprocess.run(["hadd", "-f", "-v", "0", str(target), *map(str, files)], check=True)
        n_out = sum(json.load(open(str(f) + ".json"))["n_out"] for f in files)
        manifest["samples"][key] = {"file": target.name, "bytes": target.stat().st_size, "n_events": n_out,
                                    "n_input_files": len(files), "recid": samples.SAMPLES[key]["recid"],
                                    "xsec_pb": samples.SAMPLES[key]["xsec_pb"], "sha256": sha256(target)}
        print(f"[bundle] {key}: {n_out:,} events, {target.stat().st_size/1e9:.2f} GB")
    with open(args.out / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=1)
    total = sum(v["bytes"] for v in manifest["samples"].values()) / 1e9
    print(f"[bundle] {len(manifest['samples'])} samples, {total:.1f} GB in {args.out}")


if __name__ == "__main__":
    main()
