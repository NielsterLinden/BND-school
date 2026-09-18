#!/usr/bin/env python
"""Size of every systematic uncertainty on each channel's signal yield, as it enters the combination.

    source ../../../fitting/setup.sh
    python checks/systematics.py          # after `python run.py prepare`

For each Systematic of the generated channel configs (work/common/<channel>.config) on the signal sample:
the relative change of the signal yield summed over the fitted bins (DropBins respected) of every fitted
region and every signal template, for the up and down variation. OVERALL systematics are their configured values; HISTO ones are read from the fit inputs;
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
    signals = trexcfg.signal_samples(spec)          # one template (ee, mumu) or one per tau_h decay mode (tautau)
    files = {}

    def histograms(sample):              # ee names its file per sample, mumu and tautau per job
        name = samples[sample].opts.get("HistoFile", job.opts.get("HistoFile"))
        if name not in files:
            f = uproot.open(Path(job.opts["HistoPath"]) / f"{name}.root")
            files[name] = (f, {k.split(";")[0] for k in f.keys()})
        return files[name]

    def in_region(block, region):
        where = trexcfg.listopt(block.opts.get("Regions", "all"))
        return "all" in where or region.name in where

    def integral(sample, region, suffix=""):
        """Signal yield of one template in the fitted bins of one region; None if that histogram does not exist."""
        smp = samples[sample].opts
        name = smp.get("HistoName", region.opts.get("HistoName", region.name)) + smp.get("HistoNameSuff", "") + suffix
        f, keys = histograms(sample)
        if name not in keys:
            return None
        v = f[name].values()
        return float(v[[i for i in range(len(v)) if i not in _dropped(region)]].sum())

    cells = {(smp, r.name): integral(smp, r) for smp in signals for r in regions if in_region(samples[smp], r)}
    cells = {k: v for k, v in cells.items() if v is not None}
    nominal = sum(cells.values())
    by_name = {r.name: r for r in regions}
    out = {}
    for b in blocks:
        if b.kind != "Systematic" or "DropNorm" in b.opts:
            continue
        on = trexcfg.listopt(b.opts.get("Samples", "all"))
        hit = [(smp, r) for (smp, r) in cells if ("all" in on or smp in on) and in_region(b, by_name[r])]
        if not hit:                      # not on the signal, or only in a validation region
            continue
        if b.opts["Type"].upper() == "OVERALL":
            frac = sum(cells[c] for c in hit) / nominal
            up, down = float(b.opts["OverallUp"]) * frac, float(b.opts["OverallDown"]) * frac
        else:
            def varied(suffix):
                got = {c: integral(c[0], by_name[c[1]], suffix) for c in hit}
                return sum(cells[c] if got[c] is None else got[c] for c in hit) + sum(v for c, v in cells.items() if c not in hit)
            up = varied(b.opts["HistoNameSufUp"]) / nominal - 1
            down = varied(b.opts["HistoNameSufDown"]) / nominal - 1 if b.opts.get("HistoNameSufDown") else -up
        name = b.opts.get("NuisanceParameter", b.name)
        out[name] = {"category": b.opts.get("Category", ""), "up_pct": 100 * up, "down_pct": 100 * down,
                     "sym_pct": 50 * (abs(up) + abs(down))}
    return {"signal": signals[0] if len(signals) == 1 else f"{len(signals)} templates ({signals[0]}, ...)",
            "nominal_yield": nominal, "systematics": out}


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
