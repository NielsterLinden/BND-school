#!/usr/bin/env python
"""Build the tau_h tau_h base of the measurement: skims, ntuples, fake factors, BDT and the templates
that the four-channel fit reads (`fit/fitinputs/tautau_base.root`).

    source ../fitting/setup.sh                 # LCG_110 (python 3.13, uproot, awkward, hist, ROOT) + TRExFitter
    python run_tautau_base.py          # steps 0-4
    python run_tautau_base.py --from 3 # from the fake factors on (ntuples exist, ~25 min)

Steps (each is a standalone script in scripts/, documented in docs/):
    0 external inputs and file lists        (network; seconds)
    1 skims of data and simulation          (dCache/EOS; ~1 h on 12 cores, resumable)
    2 flat ntuples                          (~5 min)
    3 fake factors, closure and OS/SS corrections; then the k-fold BDT (step3b)
    4 tau_h tau_h templates per BDT category: Data, the fake estimate and its nuisance parameters
      -> fit/fitinputs/tautau_base.root, from which step 4 of run_all.py copies them

This chain has no fit of its own any more: since v4 the measurement is the four-channel fit of
`run_all.py`, which rebuilds the tau_h tau_h simulation templates itself (v4 signal definition, free
tau_h ID scale factors, 3% tau energy-scale prior) and keeps only Data and Fakes from here.
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
    2: [["scripts/step2_ntuples_tautau.py"]],
    3: [["scripts/step3_fakefactors.py"], ["scripts/step3b_bdt.py"]],
    4: [["scripts/step4a_tautau_base.py"]],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=0)
    ap.add_argument("--to", dest="stop", type=int, default=4)
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
