#!/usr/bin/env python
"""Step 4b (v4) -- one fit-input file per channel, for the per-channel workspaces and the MultiFit combination
(fitting/CONVENTIONS.md sections 1 and 5; docs/10-v4-plan.md section 13).

    python scripts/step4b_export_channels.py

Reads fit_v4/fitinputs/ztautau_v4.root (+ .meta.json, step 4) and writes, for every channel,
fit_v4/fitinputs/ztautau_v4_<channel>.root with exactly the histograms of that channel's regions and a
.meta.json restricted to them (regions, bins, templates, systematic registry). The histograms are copied
bit for bit, so a channel workspace built from its own file and the combined workspace are the same model.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import hist
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fitting import trexhist  # noqa: E402
from ztautau import config  # noqa: E402


def main():
    src = config.FIT_DIR_V4 / "fitinputs" / f"{config.JOB_V4}.root"
    meta = json.loads(Path(str(src) + ".meta.json").read_text())
    with uproot.open(src) as f:
        names = [k.split(";")[0] for k in f.keys() if k.split(";")[0] != "meta_json"]
        hists = {n: (f[n].values(), f[n].variances()) for n in names}
    for ch in config.CHANNELS:
        regions = [r for r in meta["regions"] if r.split("_")[0] == ch]
        if not regions:
            print(f"[export] {ch}: no regions, skipped")
            continue
        out = {}
        for n, (v, var) in hists.items():
            region = n.split("__")[0]
            if region in regions:
                h = hist.Hist(hist.axis.Variable(np.asarray(meta["bins"][region], dtype=float)), storage=hist.storage.Weight())
                h.view().value[...] = v
                h.view().variance[...] = var
                out[n] = h
        samples = {t: {**i, "regions": [r for r in i["regions"] if r in regions]} for t, i in meta["samples"].items()
                   if any(r in regions for r in i["regions"])}
        systs = {}
        for name, s in meta["systs"].items():
            regs = [r for r in s["regions"] if r in regions]
            smp = [x for x in s["samples"] if x in samples or x == "Fakes"]
            if regs and smp:
                systs[name] = {**s, "regions": regs, "samples": smp}
        m = {**meta, "channel": ch, "regions": regions, "bins": {r: meta["bins"][r] for r in regions}, "samples": samples, "systs": systs,
             "source": src.name}
        path = config.FIT_DIR_V4 / "fitinputs" / f"{config.JOB_V4}_{ch}.root"
        rep = trexhist.write_fitinputs(path, out, meta=m)
        print(f"[export] {ch}: {len(regions)} regions, {len(samples)} templates, {len(systs)} systematics, {rep['n_hists']} histograms -> {path.name}")


if __name__ == "__main__":
    main()
