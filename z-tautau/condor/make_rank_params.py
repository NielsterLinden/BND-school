#!/usr/bin/env python
"""Write condor/params_rank.txt: one line per Condor job of the parallel nuisance-parameter ranking.

    python condor/make_rank_params.py [<job>]

TRExFitter's `r` action ranks every floating parameter one after the other (two fits each), which is the
slowest step of the chain. `trex-fitter r <config> Ranking=<name>` does a single parameter and writes
`Fits/NPRanking_<name>_<poi>.txt`; `Ranking=plot` merges whatever files it finds. So the ranking fans out
over Condor with one job per parameter.

`Ranking=<name>` matches by **substring**, so a job for `TauIDSF_DM1` also covers `TauIDSF_DM10` and
`TauIDSF_DM11`. This script therefore emits the *minimal* set of names (those that contain no other
parameter name), which covers every parameter exactly once -- emitting all names would rank some of them
twice and duplicate them in the merged file.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHANNEL = HERE.parent
sys.path.insert(0, str(CHANNEL.parent))
sys.path.insert(0, str(CHANNEL))

from fitting import run_trex  # noqa: E402
from ztautau import config  # noqa: E402

SKIP = {"mu_Z"}          # the POI is not ranked against itself


def main():
    job = sys.argv[1] if len(sys.argv) > 1 else config.JOB_V4
    fit_txt = config.FIT_DIR_V4 / "results" / job / "Fits" / f"{job}.txt"
    names = [n for n in run_trex.parse_fit_txt(fit_txt)["nps"] if n not in SKIP]
    minimal = [n for n in names if not any(m != n and m in n for m in names)]
    # a parameter that contains two different minimal names would be ranked twice: warn (none so far)
    for n in names:
        hits = [m for m in minimal if m in n]
        if len(hits) > 1:
            print(f"WARNING: {n} is covered by {hits} and will appear more than once in the merged ranking")
    out = HERE / "params_rank.txt"
    out.write_text("\n".join(minimal) + "\n")
    print(f"{len(names)} floating parameters -> {len(minimal)} Condor jobs -> {out}")


if __name__ == "__main__":
    main()
