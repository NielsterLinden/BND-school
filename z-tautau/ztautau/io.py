"""Certification mask, trigger OR and MET filters (tolerant of branches missing from a file).

The golden JSON must be applied before counting: only certified lumisections are in the luminosity,
so N/L is only meaningful for the same set. The 2016 HLT menu changed during data taking, so trigger
branches differ between files (Run2016G has `HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg`, Run2016H
`HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg`): a missing trigger counts as not fired.
"""

from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np

from . import config


def load_grl(path: Path | None = None) -> dict[int, np.ndarray]:
    with open(path or config.GRL_PATH) as fh:
        raw = json.load(fh)
    return {int(run): np.asarray(blocks, dtype=np.int64) for run, blocks in raw.items()}


def lumi_mask(run, lumi, grl) -> np.ndarray:
    """True where (run, luminosityBlock) is certified; loops over the few runs of a chunk."""
    run = np.asarray(ak.to_numpy(run))
    lumi = np.asarray(ak.to_numpy(lumi))
    keep = np.zeros(run.shape, dtype=bool)
    for r in np.unique(run):
        blocks = grl.get(int(r))
        if blocks is None:
            continue
        sel = run == r
        ls = lumi[sel]
        ok = np.zeros(ls.shape, dtype=bool)
        for lo, hi in blocks:
            ok |= (ls >= lo) & (ls <= hi)
        keep[sel] = ok
    return keep


def trigger_or(events, names) -> np.ndarray:
    fields = set(ak.fields(events))
    fired = np.zeros(len(events), dtype=bool)
    for name in names:
        if name in fields:
            fired |= ak.to_numpy(events[name]).astype(bool)
    return fired


def met_filters(events) -> np.ndarray:
    fields = set(ak.fields(events))
    ok = np.ones(len(events), dtype=bool)
    for name in config.MET_FILTERS:
        if name in fields:
            ok &= ak.to_numpy(events[name]).astype(bool)
    return ok
