#!/usr/bin/env python
"""v2 step 0 -- file lists for every sample from the CERN Open Data API.

Writes filelists/<key>_<recid>.sizes.tsv (xrootd URL, size, adler32) for all samples in
zmumu/samples.py and reports how many of the files have a dCache copy.

Run:  python scripts/v2_0_filelists.py [--only KEY ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zmumu import samples


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--refresh", action="store_true", help="re-download even if the list exists")
    args = ap.parse_args()
    keys = args.only or list(samples.SAMPLES)
    print(f"{'sample':18s} {'recid':>6s} {'files':>5s} {'size GB':>8s} {'on dCache':>9s}")
    for key in keys:
        path = samples.filelist_path(key)
        if args.refresh or not path.exists():
            samples.write_filelist(key)
        src = samples.sources(key)
        size = sum(s["size"] for s in src) / 1e9
        local = sum(1 for s in src if s["local"])
        print(f"{key:18s} {samples.SAMPLES[key]['recid']:>6d} {len(src):>5d} {size:>8.1f} {local:>9d}")


if __name__ == "__main__":
    main()
