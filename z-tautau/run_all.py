#!/usr/bin/env python
"""Run the Z -> tau_h tau_h measurement end to end.

    source ../setup.sh                 # LCG_110 (python 3.13, uproot, awkward, hist, ROOT) + TRExFitter
    python run_all.py                  # steps 0-6
    python run_all.py --from 3         # from the fake factors on (ntuples exist)
    python run_all.py --only 5         # just the fit

Steps (each is a standalone script in scripts/, documented in docs/):
    0 external inputs and file lists        (network; seconds)
    1 skims of data and simulation          (dCache/EOS; ~1 h on 12 cores, resumable)
    2 flat ntuples                          (~5 min)
    3 fake factors, closure corrections, OS/SS correction, closure; then the k-fold BDT (step3b)
    4 histograms per BDT category, systematic variations, control plots, fit inputs   (both fake-factor variants)
    5 TRExFitter fits                       (both fake-factor variants)
    6 report: output/results.json, output/RESULTS.md, summary plots
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = {
    0: [["scripts/step0_external.py"]],
    1: [["scripts/step1_skim.py"]],
    2: [["scripts/step2_ntuples.py"]],
    3: [["scripts/step3_fakefactors.py"], ["scripts/step3b_bdt.py"]],
    4: [["scripts/step4_histograms.py", "--ff-variant", "mcsub"],
        ["scripts/step4_histograms.py", "--ff-variant", "nosub", "--no-plots"]],
    5: [["scripts/step5_fit.py", "--ff-variant", "mcsub"],
        ["scripts/step5_fit.py", "--ff-variant", "nosub", "--skip-ranking"]],
    6: [["scripts/step6_report.py"]],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=0)
    ap.add_argument("--to", dest="stop", type=int, default=6)
    ap.add_argument("--only", type=int, default=None)
    args = ap.parse_args()
    steps = [args.only] if args.only is not None else range(args.start, args.stop + 1)
    for step in steps:
        for cmd in STEPS[step]:
            t0 = time.time()
            print(f"==== step {step}: {' '.join(cmd)}", flush=True)
            subprocess.run([sys.executable, *cmd], cwd=HERE, check=True)
            print(f"==== step {step} done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
