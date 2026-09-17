#!/usr/bin/env python
"""Size of every systematic uncertainty on each channel's signal yield, as it enters the combination.

    source ../../../setup.sh
    python checks/systematics.py          # after `python run.py prepare`

For each Systematic of the generated channel configs (work/common/<channel>.config) on the signal sample:
the relative change of the signal yield summed over the fitted bins (DropBins respected), for the up and
down variation. OVERALL systematics are their configured values; HISTO ones are read from the fit inputs;
shape-only parts (the ee "<NP>_eeShape" blocks) are 0 by construction and skipped. This is the quantity
ATLAS tabulates as dC/C (arXiv:1603.09222, Table 1); the post-fit impacts on the combined cross section are
in output/result.json and output/plots/breakdown.pdf.

Writes checks/systematics.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import uproot

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from mf import trexcfg  # noqa: E402
from mf.paths import WORK, manifest  # noqa: E402


def _dropped(region):
    return {int(i) - 1 for i in trexcfg.listopt(region.opts.get("DropBins", ""))}


def channel_table(key: str, spec: dict) -> dict:
    blocks = trexcfg.read(WORK / "common" / f"{key}.config")
    job = trexcfg.first(blocks, "Job")
    regions = [b for b in blocks if b.kind == "Region" and b.opts.get("Type", "SIGNAL").upper() != "VALIDATION"]
    samples = {b.name: b for b in blocks if b.kind == "Sample"}
    signal = spec["signal"][0]
    s = samples[signal]
    f = uproot.open(Path(job.opts["HistoPath"]) / f'{s.opts.get("HistoFile", job.opts.get("HistoFile"))}.root')

    def integral(region, suffix=""):
        base = s.opts.get("HistoName", region.opts.get("HistoName", region.name)) + s.opts.get("HistoNameSuff", "")
        v = f[base + suffix].values()
        keep = [i for i in range(len(v)) if i not in _dropped(region)]
        return float(v[keep].sum())

    nominal = sum(integral(r) for r in regions)
    out = {}
    for b in blocks:
        if b.kind != "Systematic" or "DropNorm" in b.opts:
            continue
        smp = trexcfg.listopt(b.opts.get("Samples", "all"))
        if signal not in smp and "all" not in smp:
            continue
        name = b.opts.get("NuisanceParameter", b.name)
        sel = [r for r in regions if not b.opts.get("Regions") or r.name in trexcfg.listopt(b.opts["Regions"])]
        if not sel:                      # acts only in a validation region, not in the fit
            continue
        if b.opts["Type"].upper() == "OVERALL":
            up, down = float(b.opts["OverallUp"]), float(b.opts["OverallDown"])
        else:
            su, sd = b.opts.get("HistoNameSufUp"), b.opts.get("HistoNameSufDown")
            up = sum(integral(r, su) for r in sel) / sum(integral(r) for r in sel) - 1
            down = sum(integral(r, sd) for r in sel) / sum(integral(r) for r in sel) - 1 if sd else -up
        out[name] = {"category": b.opts.get("Category", ""), "up_pct": 100 * up, "down_pct": 100 * down,
                     "sym_pct": 50 * (abs(up) + abs(down))}
    return {"signal": signal, "nominal_yield": nominal, "systematics": out}


def main():
    m = manifest()
    table = {k: channel_table(k, spec) for k, spec in m["channels"].items()}
    (HERE / "checks" / "systematics.json").write_text(json.dumps(table, indent=1) + "\n")
    for k, t in table.items():
        print(f"== {k} ({t['signal']}, {t['nominal_yield']:.0f} events)")
        for name, v in sorted(t["systematics"].items(), key=lambda kv: -kv[1]["sym_pct"]):
            print(f"   {name:28s} {v['category']:26s} {v['up_pct']:+7.3f} {v['down_pct']:+7.3f}")


if __name__ == "__main__":
    main()
