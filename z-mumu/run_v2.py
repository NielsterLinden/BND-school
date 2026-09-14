#!/usr/bin/env python
"""Run the v2 (MC-based, TRExFitter) Z -> mu mu measurement end to end.

    source ../setup.sh
    python run_v2.py                     # everything (skims take ~1.5 h the first time)
    python run_v2.py --from 2            # after the skims exist
    python run_v2.py --from 2 --max-files 2 --workers 4     # quick smoke test

Steps: 0 filelists, 1 skim, 2 tag-and-probe, 3 control (momentum + fakes), 4 histograms,
5 fit, 6 report. Each step is also runnable on its own (scripts/v2_<n>_*.py).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = [(0, "file lists", "v2_0_filelists.py", []), (1, "skims", "v2_1_skim.py", ["--workers", "--max-files"]),
         (2, "tag-and-probe", "v2_2_tnp.py", ["--workers", "--max-files"]),
         (3, "control regions", "v2_3_control.py", ["--workers", "--max-files"]),
         (4, "histograms", "v2_4_histograms.py", ["--workers", "--max-files"]),
         (5, "fit", "v2_5_fit.py", []), (6, "report", "v2_6_report.py", [])]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=0)
    ap.add_argument("--to", dest="stop", type=int, default=6)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--max-files", type=int, default=None)
    args = ap.parse_args()
    t_start = time.time()
    for number, name, script, opts in STEPS:
        if not (args.start <= number <= args.stop):
            continue
        cmd = [sys.executable, str(HERE / "scripts" / script)]
        if "--workers" in opts and args.workers:
            cmd += ["--workers", str(args.workers)]
        if "--max-files" in opts and args.max_files:
            cmd += ["--max-files", str(args.max_files)]
        print(f"\n{'='*74}\n  v2 STEP {number}: {name}\n{'='*74}", flush=True)
        t0 = time.time()
        if subprocess.run(cmd).returncode != 0:
            sys.exit(f"step {number} ({name}) failed")
        print(f"  -- step {number} done in {time.time()-t0:.0f}s", flush=True)
    print(f"\nAll requested steps finished in {time.time()-t_start:.0f}s. See output/v2/RESULTS_v2.md")


if __name__ == "__main__":
    main()
