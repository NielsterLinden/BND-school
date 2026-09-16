#!/usr/bin/env python3
"""
Merge the separate per-sample ROOT files (ttbar.root, Wjets.root, Diboson.root, DY_ee.root,
DY_tautau.root, Data.root, ...) already written by the notebook into a single ROOT file,
without recreating any histogram. Each histogram is renamed "<sample>_<variation>" so nothing
collides once everything lives in one file:

    h_mass                -> <sample>_nominal
    h_mass_PUUp            -> <sample>_PUUp
    h_mass_PDFUp            -> <sample>_PDFUp
    ...

<sample> is taken from each input file's name (without .root). Adjust --nominal-key /
--prefix if your histograms don't use the "h_mass" / "h_mass_<variation>" naming.

Usage:
    python merge_to_single_root.py /path/to/histogram/folder
    python merge_to_single_root.py /path/to/folder --output zee.root
    python merge_to_single_root.py /path/to/folder --pattern "*.root" --exclude Data.root
    python merge_to_single_root.py /path/to/folder --dry-run
"""
import argparse
import glob
import os
import sys

import uproot


def variation_name(key, nominal_key, prefix):
    if key == nominal_key:
        return None
    if key.startswith(prefix):
        return key[len(prefix):]
    return key  # fallback: histogram doesn't match the expected naming, keep as-is


def collect_histograms(files, nominal_key, prefix):
    """Returns (dict of {output_name: boost_histogram}, list of (source_file, key, output_name)
    collisions)."""
    merged = {}
    sources = {}
    collisions = []

    for path in files:
        sample = os.path.splitext(os.path.basename(path))[0]
        with uproot.open(path) as fin:
            classnames = fin.classnames()
            keys = sorted({k.split(";")[0] for k, cls in classnames.items() if cls.startswith("TH1")})
            for key in keys:
                var = variation_name(key, nominal_key, prefix)
                out_name = f"{sample}_{var}" if var is not None else f"{sample}"
                if out_name in merged:
                    collisions.append((path, key, out_name, sources[out_name]))
                    continue
                merged[out_name] = fin[key].to_boost()
                sources[out_name] = path

    return merged, collisions


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--folder", help="Folder containing the per-sample ROOT files")
    parser.add_argument("--output", default=None,
                         help="Output file path (default: <folder>/zee.root)")
    parser.add_argument("--pattern", default="*.root", help="Glob pattern for input files (default: *.root)")
    parser.add_argument("--exclude", default="", help="Comma-separated filenames to skip (e.g. an old zee.root)")
    parser.add_argument("--nominal-key", default="h_mass", help="Histogram name treated as nominal (default: h_mass)")
    parser.add_argument("--prefix", default="h_mass_", help="Prefix stripped to get the variation name (default: h_mass_)")
    parser.add_argument("--dry-run", action="store_true", help="Report what would be written without creating the output file")
    args = parser.parse_args()

    output_path = args.output or os.path.join(args.folder, "zee.root")
    exclude = {e.strip() for e in args.exclude.split(",") if e.strip()}

    files = sorted(glob.glob(os.path.join(args.folder, args.pattern)))
    files = [f for f in files if os.path.basename(f) not in exclude
             and os.path.abspath(f) != os.path.abspath(output_path)]
    if not files:
        print(f"No input files matching {args.pattern} in {args.folder} (after exclusions)", file=sys.stderr)
        sys.exit(1)

    merged, collisions = collect_histograms(files, args.nominal_key, args.prefix)

    print(f"Found {len(files)} input file(s):")
    for f in files:
        print(f"  {os.path.basename(f)}")

    if collisions:
        print(f"\n{len(collisions)} naming collision(s) -- these histograms were NOT included:")
        for path, key, out_name, first_source in collisions:
            print(f"  {os.path.basename(path)}:{key} -> '{out_name}' already came from {os.path.basename(first_source)}")

    print(f"\n{'Would write' if args.dry_run else 'Writing'} {len(merged)} histogram(s) to {output_path}")
    for name in sorted(merged):
        print(f"  {name}")

    if not args.dry_run:
        with uproot.recreate(output_path) as fout:
            for name, h in merged.items():
                fout[name] = h
        print(f"\nDone: {output_path}")


if __name__ == "__main__":
    main()