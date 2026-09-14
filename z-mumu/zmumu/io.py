"""Reading the skimmed NanoAOD: file discovery, chunking, and the lumi mask.

Two things deserve comment.

**The golden JSON.** CMS certifies data per (run, luminosityBlock). Only
certified lumisections are in the luminosity that `brilcalc` reports, so the
numerator and denominator of a cross section must use the same set. Applying
the mask is therefore not optional bookkeeping -- it is what makes
N_events / L meaningful.

**Branch sets differ between files.** The 2016 HLT menu changed during
data-taking, so 39 of the 328 trigger branches are absent from some files. Any
code that asks for a branch by name must tolerate its absence; `read_chunks`
handles this by intersecting the request with what the file actually has and
reporting missing branches as all-False.
"""

from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from . import config


def find_files(skim_dir: Path = None, eras=None) -> list[Path]:
    """Return every skimmed ROOT file, sorted, as absolute paths."""
    skim_dir = Path(skim_dir or config.SKIM_DIR)
    eras = eras or config.ERAS
    files: list[Path] = []
    for era in eras:
        found = sorted((skim_dir / era).glob("*.root"))
        if not found:
            raise FileNotFoundError(f"no ROOT files under {skim_dir / era}")
        files.extend(found)
    return files


def era_of(path: Path) -> str:
    """Era label ('Run2016G__30530') for a skim file."""
    return Path(path).parent.name


# --------------------------------------------------------------------------
# Golden JSON
# --------------------------------------------------------------------------
def load_grl(path: Path = None) -> dict[int, np.ndarray]:
    """Load the certification JSON into {run: array of [first, last] blocks}."""
    path = Path(path or config.GRL_PATH)
    with open(path) as fh:
        raw = json.load(fh)
    return {int(run): np.asarray(blocks, dtype=np.int64) for run, blocks in raw.items()}


def lumi_mask(run, luminosity_block, grl: dict[int, np.ndarray]) -> np.ndarray:
    """Boolean mask: True where (run, lumi) is in the certified list.

    Loops over the runs actually present in the chunk (a few at most), so this
    stays fast even though the GRL has ~400 runs.
    """
    run = np.asarray(ak.to_numpy(run))
    lumi = np.asarray(ak.to_numpy(luminosity_block))
    keep = np.zeros(run.shape, dtype=bool)

    for r in np.unique(run):
        blocks = grl.get(int(r))
        if blocks is None:
            continue                      # run not certified at all
        in_run = run == r
        lumis_here = lumi[in_run]
        ok = np.zeros(lumis_here.shape, dtype=bool)
        for lo, hi in blocks:
            ok |= (lumis_here >= lo) & (lumis_here <= hi)
        keep[in_run] = ok
    return keep


# --------------------------------------------------------------------------
# Chunked reading
# --------------------------------------------------------------------------
def read_chunks(path: Path, branches: list[str], chunk_size: str = None):
    """Yield awkward chunks of `branches` from one skim file.

    Branches missing from this particular file (trigger menu changes) are
    filled with False rather than raising, and the set of missing names is
    yielded alongside each chunk so callers can record it.
    """
    chunk_size = chunk_size or config.CHUNK_SIZE
    tree = uproot.open(f"{path}:Events")
    available = set(tree.keys())
    present = [b for b in branches if b in available]
    missing = [b for b in branches if b not in available]

    for chunk in tree.iterate(present, step_size=chunk_size):
        if missing:
            n = len(chunk)
            filler = {name: np.zeros(n, dtype=bool) for name in missing}
            chunk = ak.with_field(chunk, ak.Array(np.zeros(n, dtype=bool)), missing[0])
            for name, arr in filler.items():
                chunk = ak.with_field(chunk, ak.Array(arr), name)
        yield chunk, missing


def trigger_or(events, names: list[str]) -> ak.Array:
    """Logical OR of trigger branches, tolerating ones absent from this file."""
    fields = set(ak.fields(events))
    usable = [n for n in names if n in fields]
    if not usable:
        return ak.Array(np.zeros(len(events), dtype=bool))
    fired = events[usable[0]]
    for name in usable[1:]:
        fired = fired | events[name]
    return fired


def pass_met_filters(events) -> ak.Array:
    """AND of the recommended 2016 MET filter flags present in this file."""
    fields = set(ak.fields(events))
    usable = [f for f in config.MET_FILTERS if f in fields]
    ok = ak.Array(np.ones(len(events), dtype=bool))
    for name in usable:
        ok = ok & events[name]
    return ok
