#!/usr/bin/env python3
"""
Set negative bin contents in ROOT histograms to a small positive floor, in-place, for every
ROOT file in a given folder.

Only bin *contents* are touched; bin variances/errors are left untouched, so the reported
error bar on a floored bin will still reflect its original (larger) statistical uncertainty.

Usage:
    python fix_negative_bins.py /path/to/histogram/folder
    python fix_negative_bins.py /path/to/histogram/folder --floor 1e-6 --pattern "*.root"
    python fix_negative_bins.py /path/to/histogram/folder --dry-run
    python fix_negative_bins.py /path/to/histogram/folder --backup
"""
import argparse
import glob
import os
import shutil
import sys

import uproot


def fix_histogram(h, floor):
    """h: a boost_histogram.Histogram (from uproot's .to_boost()).
    Returns (h, n_fixed) -- h is modified in place."""
    view = h.view(flow=True)
    # Weight storage -> structured array with .value/.variance; plain storage -> bare ndarray.
    values = view.value if hasattr(view, "value") else view
    neg_mask = values <= 0
    n_fixed = int(neg_mask.sum())
    if n_fixed:
        values[neg_mask] = floor
    return h, n_fixed


def process_file(path, floor, dry_run):
    with uproot.open(path) as fin:
        classnames = fin.classnames()
        hist_keys = sorted({k.split(";")[0] for k, cls in classnames.items() if cls.startswith("TH1")})

        if not hist_keys:
            return 0, 0

        fixed_hists = {}
        total_fixed = 0
        for k in hist_keys:
            h_boost = fin[k].to_boost()
            h_boost, n_fixed = fix_histogram(h_boost, floor)
            fixed_hists[k] = h_boost
            total_fixed += n_fixed

    if total_fixed and not dry_run:
        with uproot.recreate(path) as fout:
            for k, h in fixed_hists.items():
                fout[k] = h

    return len(hist_keys), total_fixed


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--folder", help="Folder containing ROOT histogram files")
    parser.add_argument("--floor", type=float, default=1e-6, help="Replacement value for negative bins (default: 1e-6)")
    parser.add_argument("--pattern", default="*.root", help="Glob pattern for files to process (default: *.root)")
    parser.add_argument("--dry-run", action="store_true", help="Report negative bins without modifying any file")
    parser.add_argument("--backup", action="store_true", help="Write a .bak copy of each file before modifying it")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.folder, args.pattern)))
    if not files:
        print(f"No files matching {args.pattern} in {args.folder}", file=sys.stderr)
        sys.exit(1)

    grand_total = 0
    for path in files:
        # Peek first so we only ever write a .bak when something will actually change.
        _, n_fixed_peek = process_file(path, args.floor, dry_run=True)

        if n_fixed_peek and args.backup and not args.dry_run:
            shutil.copy2(path, path + ".bak")

        n_hists, n_fixed = process_file(path, args.floor, args.dry_run)
        if n_fixed:
            tag = "[dry-run] would fix" if args.dry_run else "fixed"
            print(f"{os.path.basename(path)}: {tag} {n_fixed} negative bin(s) across {n_hists} histogram(s)")
        grand_total += n_fixed

    print(f"\nTotal negative bins {'found' if args.dry_run else 'fixed'}: {grand_total}")


if __name__ == "__main__":
    main()