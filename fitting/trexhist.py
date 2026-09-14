"""Write and check TRExFitter histogram inputs in the shared convention.

See CONVENTIONS.md. Histogram names are ``<region>__<sample>`` and
``<region>__<sample>__<syst>Up|Down``; every histogram is a TH1D with the sum of
weights squared stored, which needs ``hist.storage.Weight()`` before uproot writes it.

    from fitting import trexhist
    trexhist.write_fitinputs("fit/fitinputs/zmumu.root", {"mumu_SR__DYmumu": h, ...},
                             meta={"lumi_pb": 16393.381, ...})
    trexhist.check_fitinputs("fit/fitinputs/zmumu.root", ["mumu_SR__DYmumu", ...])

No ROOT import anywhere: uproot writes the TH1D objects.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import hist
import numpy as np
import uproot

NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")


def hname(region: str, sample: str, syst: str | None = None, direction: str | None = None) -> str:
    """Histogram name for (region, sample[, systematic, 'Up'/'Down'])."""
    name = f"{region}__{sample}"
    if syst:
        if direction not in ("Up", "Down"):
            raise ValueError("direction must be 'Up' or 'Down' when syst is given")
        name += f"__{syst}{direction}"
    if not NAME_RE.match(name):
        raise ValueError(f"illegal histogram name {name!r} (letters, digits, underscores only)")
    return name


def split_name(name: str):
    """Inverse of hname: (region, sample, syst or None, direction or None)."""
    parts = name.split("__")
    if len(parts) == 2:
        return parts[0], parts[1], None, None
    if len(parts) == 3:
        for d in ("Up", "Down"):
            if parts[2].endswith(d):
                return parts[0], parts[1], parts[2][: -len(d)], d
    raise ValueError(f"cannot parse histogram name {name!r}")


def with_weight_storage(h: hist.Hist) -> hist.Hist:
    """Return a copy of `h` with Weight storage (sum of w^2 kept), whatever it had."""
    if h.storage_type is hist.storage.Weight:
        return h.copy()
    out = hist.Hist(*h.axes, storage=hist.storage.Weight())
    values = h.values(flow=True)
    variances = h.variances(flow=True)
    if variances is None:            # unweighted storage: Poisson, w = 1
        variances = values
    view = out.view(flow=True)
    view.value[...] = values
    view.variance[...] = variances
    return out


def write_fitinputs(path, hists: dict, meta: dict | None = None, clip_negative: bool = True) -> dict:
    """Write `hists` ({name: hist.Hist}) as TH1D to `path`; returns a report dict.

    Checks: names follow the convention, all histograms of a region share the binning,
    nominal templates have no negative bins (clipped to 0 when `clip_negative`, reported),
    no NaN/inf. `meta` (a JSON-serialisable dict) is stored as a TObjString `meta_json`
    and next to the file as `<path>.meta.json`.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    report = {"n_hists": 0, "clipped": {}, "regions": {}}
    edges_by_region = {}
    prepared = {}
    for name, h in hists.items():
        region, sample, syst, direction = split_name(name)
        if h.ndim != 1:
            raise ValueError(f"{name}: only 1D histograms are fit inputs")
        h = with_weight_storage(h)
        edges = h.axes[0].edges
        ref = edges_by_region.setdefault(region, edges)
        if not np.array_equal(ref, edges):
            raise ValueError(f"{name}: binning differs from other histograms of region {region}")
        vals = h.values()
        if not np.all(np.isfinite(vals)) or not np.all(np.isfinite(h.variances())):
            raise ValueError(f"{name}: non-finite bin content")
        neg = vals < 0
        if neg.any():
            if not clip_negative:
                raise ValueError(f"{name}: {int(neg.sum())} negative bins")
            report["clipped"][name] = float(vals[neg].sum())
            view = h.view()
            view.value[neg] = 0.0
        prepared[name] = h
        report["regions"].setdefault(region, []).append(name)
    with uproot.recreate(path) as f:
        for name, h in prepared.items():
            f[name] = h
        if meta is not None:
            f["meta_json"] = json.dumps(meta, indent=1, sort_keys=True, default=_json_default)
    if meta is not None:
        with open(str(path) + ".meta.json", "w") as fh:
            json.dump(meta, fh, indent=1, sort_keys=True, default=_json_default)
    report["n_hists"] = len(prepared)
    return report


def _json_default(obj):
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not JSON serialisable: {type(obj)}")


def check_fitinputs(path, names) -> dict:
    """Assert every name in `names` exists in `path` as a TH1 with fSumw2 filled."""
    problems = []
    found = {}
    with uproot.open(path) as f:
        keys = {k.split(";")[0] for k in f.keys()}
        for name in names:
            if name not in keys:
                problems.append(f"missing: {name}")
                continue
            obj = f[name]
            sumw2 = obj.member("fSumw2")
            if sumw2 is None or len(sumw2) == 0:
                problems.append(f"no Sumw2: {name}")
            vals = obj.values()
            if (vals < 0).any():
                problems.append(f"negative bins: {name}")
            found[name] = float(vals.sum())
    if problems:
        raise RuntimeError("fit-input check failed:\n  " + "\n  ".join(problems))
    return found


def read_meta(path) -> dict:
    with uproot.open(path) as f:
        return json.loads(str(f["meta_json"]))
