#!/usr/bin/env python
"""Run the Z -> tautau measurement (four channels: tau_h tau_h, mu tau_h, e tau_h, e mu) end to end.

    source ../setup.sh
    python run_all.py                  # steps 1-6 (skims: hours; from the ntuples on: ~4 h with the ranking)
    python run_all.py --from 3         # from the fake factors / trigger efficiencies on
    python run_all.py --only 5         # just the fits

Steps (each a standalone script in scripts/, documented in docs/10-v4-plan.md):
    1 skims (SingleMuon, SingleElectron, MuonEG, all simulation)                 step1_skim.py --v4
    2 flat ntuples of the lepton channels                                        step2_ntuples_lepton.py
    3 in-situ trigger efficiencies + b-tag efficiencies; lepton-channel fakes     step3c_trigger.py, step3d_fakes_lepton.py
    4 templates of all channels (tau_h tau_h from fit/fitinputs/tautau_base.root) step4_histograms.py, step4b_export_channels.py
    5 TRExFitter fits: the measurement (with ranking), the cross-checks of REVIEW_v4.md, the per-channel
      workspaces and their MultiFit                                              step5_fit.py, step5b_multifit.py
    6 report, and the result block injected into README/handoff/docs                step6_report.py, update_docs.py

The tau_h tau_h base templates must exist: `python run_tautau_base.py --from 3` (fake factors, BDT and
fit/fitinputs/tautau_base.root, ~25 min from the ntuples).

The step-5 jobs, and why each exists (the first two run before the measurement, which reads them):
    ztautau_flatsf       the measurement's model without `TauIDpT_tautau`: one tau_h ID scale factor per decay
                         mode for every pT and no uncertainty on that assumption (the result quoted until
                         17 Sep 2026, 1981 +73 -70 pb)
    ztautau_ptsplit      the l tau_h regions split at pT(tau_h) = 40 GeV with their own scale factors below
                         it: the test of the assumption the lever rests on (finding 4). Its relative
                         difference to ztautau_flatsf is the size of `TauIDpT_tautau` (step5_fit.pt_model)
    ztautau              the measurement: four channels, tau_h ID scale factors free, `TauIDpT_tautau`,
                         ranking and impacts
    ztautau_<ch>         one channel, scale factors free: the workspaces the MultiFit (and a combination) reads
    ztautau_<ch>_fixedid one channel with the TauPOG scale factors fixed: the per-channel cross-check numbers
    ztautau_taulep       tau_h tau_h + mu tau_h + e tau_h without e mu, scale factors free: what the
                         tau_h tau_h / (l tau_h)^2 lever alone says about mu_Z (REVIEW_v4.md finding 2)
    ztautau_emutrig2x    the e mu trigger variation doubled (the 2% of the paper applied once per leg
                         instead of once per event, i.e. the v4 treatment before the review): how much of
                         the answer the trigger prior sets (finding 3)
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHANNELS = ("tautau", "mutau", "etau", "emu")
FIT = "scripts/step5_fit.py"
STEPS = {
    1: [["scripts/step1_skim.py", "--v4"]],
    2: [["scripts/step2_ntuples_lepton.py"]],
    3: [["scripts/step3c_trigger.py"], ["scripts/step3d_fakes_lepton.py"]],
    4: [["scripts/step4_histograms.py"], ["scripts/step4b_export_channels.py"]],
    5: [[FIT, "--job", "ztautau_flatsf", "--no-tauid-pt", "--skip-ranking", "--skip-asimov"],
        [FIT, "--job", "ztautau_ptsplit", "--region-set", "ptsplit", "--skip-ranking", "--skip-asimov"],
        [FIT]]
       + [[FIT, "--channels", ch, "--job", f"ztautau_{ch}", "--skip-ranking", "--skip-asimov"] for ch in CHANNELS]
       + [[FIT, "--channels", ch, "--job", f"ztautau_{ch}_fixedid", "--skip-ranking", "--skip-asimov", "--fix-tauid"] for ch in CHANNELS]
       + [[FIT, "--channels", "tautau", "mutau", "etau", "--job", "ztautau_taulep", "--skip-ranking", "--skip-asimov"],
          [FIT, "--job", "ztautau_emutrig2x", "--scale-syst", "EmuTrigger=2.0", "--skip-ranking", "--skip-asimov"],
          ["scripts/step5b_multifit.py"]],
    6: [["scripts/step6_report.py"], ["scripts/update_docs.py"]],
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
