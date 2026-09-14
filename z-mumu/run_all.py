#!/usr/bin/env python
"""Run the complete Z -> mu mu cross-section measurement.

    .venv/bin/python run_all.py                 # full dataset, all 152 files
    .venv/bin/python run_all.py --max-files 4   # quick test on 4 files
    .venv/bin/python run_all.py --from 3        # resume from step 3

The steps must run in order; each reads the pickles written by the previous
ones from output/data/. Steps 1 and 2 are the expensive ones (they each read
the full 31.5 GB skim); steps 3-6 take seconds because they work on the
histograms those two produced.

Step 1  selection ...... three analysis regions, cutflow, mass spectra
Step 2  efficiency ..... tag-and-probe ID / isolation / trigger
Step 3  backgrounds .... same-sign and e-mu data-driven estimates
Step 4  signal ......... background-subtracted count, plus a fit cross-check
Step 5  cross section .. combine, propagate uncertainties
Step 6  report ......... results.json, RESULTS.md, summary plot
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYTHON = HERE / ".venv" / "bin" / "python"

STEPS = [
    (1, "selection", "step1_selection.py", True),
    (2, "efficiency", "step2_efficiency.py", True),
    (3, "backgrounds", "step3_backgrounds.py", False),
    (4, "signal extraction", "step4_signal.py", False),
    (5, "cross section", "step5_crosssection.py", False),
    (6, "report", "step6_report.py", False),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-files", type=int, default=None,
                    help="limit steps 1 and 2 to the first N skim files")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--from", dest="start", type=int, default=1,
                    help="first step to run (default 1)")
    ap.add_argument("--to", dest="stop", type=int, default=6)
    args = ap.parse_args()

    if not PYTHON.exists():
        sys.exit(f"venv interpreter not found at {PYTHON}\n"
                 "Create it with:  python3 -m venv .venv && "
                 ".venv/bin/pip install -r requirements.txt")

    t_start = time.time()
    for number, name, script, reads_data in STEPS:
        if not (args.start <= number <= args.stop):
            continue
        cmd = [str(PYTHON), str(HERE / "scripts" / script)]
        if reads_data:
            if args.max_files:
                cmd += ["--max-files", str(args.max_files)]
            if args.workers:
                cmd += ["--workers", str(args.workers)]

        print(f"\n{'='*74}\n  STEP {number}: {name}\n{'='*74}")
        t0 = time.time()
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(f"\nstep {number} ({name}) failed with code {result.returncode}")
        print(f"  -- step {number} done in {time.time()-t0:.0f}s")

    print(f"\nAll requested steps finished in {time.time()-t_start:.0f}s.")
    print("Results: output/RESULTS.md, output/results.json, output/plots/")


if __name__ == "__main__":
    main()
