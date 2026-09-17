#!/usr/bin/env python
"""Run the v4 four-channel measurement (docs/10-v4-plan.md) end to end.

    source ../setup.sh
    python run_v4.py                  # steps 1-6 (skims: hours; from the ntuples on: ~1.5 h)
    python run_v4.py --from 3         # from the fake factors / trigger efficiencies on
    python run_v4.py --only 5         # just the fit

Steps (each a standalone script in scripts/, documented in docs/10-v4-plan.md):
    1 v4 skims (SingleMuon, SingleElectron, MuonEG, all simulation)            step1_skim.py --v4
    2 flat ntuples of the lepton channels                                        step2_ntuples_v4.py
    3 in-situ trigger efficiencies + b-tag efficiencies; lepton-channel fakes    step3c_trigger_v4.py, step3_fakes_v4.py
    4 templates of all channels (tau_h tau_h from the v3 ntuples / fit inputs)   step4_histograms_v4.py, step4b_export_channels_v4.py
    5 TRExFitter fits: combined (with ranking), per channel (workspaces for the MultiFit, tau ID free)
      and per channel with the POG tau ID SFs fixed (cross-checks)               step5_fit_v4.py
    6 report                                                                     step6_report_v4.py
The tau_h tau_h inputs (steps 3-4 of run_all.py: fake factors, BDT, fit/fitinputs/ztautau.root) must exist.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = {
    1: [["scripts/step1_skim.py", "--v4"]],
    2: [["scripts/step2_ntuples_v4.py"]],
    3: [["scripts/step3c_trigger_v4.py"], ["scripts/step3_fakes_v4.py"]],
    4: [["scripts/step4_histograms_v4.py"], ["scripts/step4b_export_channels_v4.py"]],
    5: [["scripts/step5_fit_v4.py"]] + [["scripts/step5_fit_v4.py", "--channels", ch, "--job", f"ztautau_v4_{ch}", "--skip-ranking"] for ch in ("tautau", "mutau", "etau", "emu")]
       + [["scripts/step5_fit_v4.py", "--channels", ch, "--job", f"ztautau_v4_{ch}_fixedid", "--skip-ranking", "--fix-tauid"] for ch in ("tautau", "mutau", "etau", "emu")],
    6: [["scripts/step6_report_v4.py"]],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=1)
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
