#!/usr/bin/env python
"""Were the di-tau trigger paths active and unprescaled in every certified lumisection? (docs/02)

    python review/trigger_rate.py

The selected events (any charge / isolation, QCD dominated) per recorded luminosity must be flat run by
run within an era; certified lumisections without triggered events would reveal an inactive path.
"""
import collections, csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ztautau import analysis, config, io

data = analysis.load_data(); grl = io.load_grl()
lumi = {}
with open(config.LUMIBYLS_CSV) as fh:
    for row in csv.reader(fh):
        if not row or row[0].startswith("#"):
            continue
        run = int(row[0].split(":")[0]); ls = int(row[1].split(":")[0])
        blocks = grl.get(run)
        if run < config.RUN_MIN or run > config.RUN_MAX or blocks is None or not any(lo <= ls <= hi for lo, hi in blocks):
            continue
        lumi[(run, ls)] = lumi.get((run, ls), 0.0) + float(row[6]) * 1e3
L_run = collections.Counter()
for (r, _), v in lumi.items():
    L_run[r] += v
cnt = collections.Counter(data["run"].tolist())
runs = sorted(L_run)
print(f"certified luminosity {sum(L_run.values()):.1f} pb-1 in {len(runs)} runs")
for e, name in ((0, "Run2016G"), (1, "Run2016H")):
    rs = [r for r in runs if (r >= config.ERA_H_FIRST_RUN) == bool(e)]
    n = sum(cnt.get(r, 0) for r in rs); L = sum(L_run[r] for r in rs); mean = n / L
    print(f"{name}: {n:,} events / {L:.1f} pb-1 = {mean:.1f} per pb-1, {len(rs)} runs")
    for r in rs:
        rate = cnt.get(r, 0) / L_run[r]
        if L_run[r] > 5 and abs(rate / mean - 1) > 0.3:
            print(f"   run {r}: {cnt.get(r, 0):6d} events, {L_run[r]:6.1f} pb-1, {rate / mean:.2f} x mean")
ls_ev = collections.Counter(zip(data["run"].tolist(), data["lumi"].tolist()))
big = [k for k, v in lumi.items() if v > 0.03]
zero = [k for k in big if ls_ev.get(k, 0) == 0]
print(f"lumisections with > 0.03 pb-1: {len(big)}, without a selected event: {len(zero)} "
      f"(Poisson expectation ~{len(big) * np.exp(-6.6):.0f}), luminosity {sum(lumi[k] for k in zero):.2f} pb-1")
print(f"selected events outside certified lumisections: {sum(1 for k in set(zip(data['run'].tolist(), data['lumi'].tolist())) if k not in lumi)}")
