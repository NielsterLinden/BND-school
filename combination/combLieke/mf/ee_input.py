"""The Z -> ee fit inputs, read from the TRExFitter job z-ee published as z-ee/Zee_fit.tar.gz.

The tarball holds TRExFitter's own copy of every input histogram (`Histograms/Zee_fit_histos.root`,
`<region>/<sample>/nominal/<region>_<sample>_orig` and `<region>/<sample>/<syst>/..._<syst>_Up_orig`),
i.e. the histograms exactly as TRExFitter read them from `datasets/z-ee/zee.root` before any
smoothing or symmetrisation. They are written back under the names `z-ee/fit.config` asks for
(`<HistoName><HistoNameSufUp|Down>`) so the z-ee config can be used unchanged. Where the delivered
`datasets/z-ee/zee.root` is present, every histogram is compared bin by bin with it.

The archive is read in memory; nothing in it is extracted to disk or executed.
"""

from __future__ import annotations

import io
import tarfile
from pathlib import Path

import hist
import numpy as np
import uproot

from . import trexcfg
from .paths import repo_path


def _members(tar_path: Path) -> dict[str, bytes]:
    wanted = ("Histograms/Zee_fit_histos.root", "reference_cross_section.txt", "Fits/Zee_fit.txt")
    out = {}
    with tarfile.open(tar_path, "r:gz") as tar:
        for m in tar.getmembers():
            if m.isfile() and m.name.endswith(wanted):
                out[m.name.split("/", 1)[1]] = tar.extractfile(m).read()
    return out


def _to_hist(h) -> hist.Hist:
    edges = h.axis().edges()
    out = hist.Hist(hist.axis.Variable(edges), storage=hist.storage.Weight())
    out.view(flow=True).value = h.values(flow=True)
    var = h.variances(flow=True)
    out.view(flow=True).variance = h.values(flow=True) if var is None else var
    return out


def extract(spec: dict, out_dir: Path) -> dict:
    """Write <out_dir>/<HistoFile>.root from the tarball; return provenance and checks."""
    tar_path = repo_path(spec["tarball"])
    members = _members(tar_path)
    blocks = trexcfg.read(repo_path(spec["config"]))
    job = trexcfg.first(blocks, "Job")
    region = trexcfg.first(blocks, "Region")
    samples = [b for b in blocks if b.kind == "Sample"]
    systs = [b for b in blocks if b.kind == "Systematic" and b.opts.get("Type", "").upper() == "HISTO"]
    histo_files = {s.opts.get("HistoFile", job.opts.get("HistoFile")) for s in samples}
    if len(histo_files) != 1:
        raise ValueError(f"z-ee: expected one HistoFile, found {histo_files}")
    histo_file = histo_files.pop()

    trex = uproot.open(io.BytesIO(members["Histograms/Zee_fit_histos.root"]))
    reg = region.name
    written = {}
    for s in samples:
        base = s.opts.get("HistoName", region.opts.get("HistoName"))
        written[base] = _to_hist(trex[f"{reg}/{s.name}/nominal/{reg}_{s.name}_orig"])
        if s.opts.get("Type", "").upper() == "DATA":
            continue
        for syst in systs:
            if s.name not in trexcfg.listopt(syst.opts.get("Samples", "all")) and syst.opts.get("Samples", "all") != "all":
                continue
            for d, key in (("Up", "HistoNameSufUp"), ("Down", "HistoNameSufDown")):
                name = f"{base}{syst.opts[key]}"
                written[name] = _to_hist(trex[f"{reg}/{s.name}/{syst.name}/{reg}_{s.name}_{syst.name}_{d}_orig"])

    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{histo_file}.root"
    with uproot.recreate(target) as f:
        for name, h in written.items():
            f[name] = h

    fit_txt = members["Fits/Zee_fit.txt"].decode()
    poi = job.opts["POI"]
    mu = next(line.split() for line in fit_txt.splitlines() if line.startswith(poi + " "))
    info = {"tarball": str(spec["tarball"]), "histograms": len(written), "written": str(target),
            "sigma_reference_pb": float(members["reference_cross_section.txt"].decode().split()[0]),
            "published_fit": [float(mu[1]), float(mu[2].lstrip("+")), float(mu[3].lstrip("-"))]}
    delivered = repo_path(spec["check_against"]) if spec.get("check_against") else None
    if delivered and delivered.exists():
        ref = uproot.open(delivered)
        worst = 0.0
        for name, h in written.items():
            r = ref[name]
            if not np.array_equal(r.axis().edges(), h.axes[0].edges):
                raise ValueError(f"z-ee: binning of {name} differs between the tarball and {delivered}")
            scale = np.maximum(np.abs(r.values()), 1e-12)
            worst = max(worst, float(np.max(np.abs(r.values() - h.values()) / scale)))
            if r.variances() is not None and name != "Data":
                worst = max(worst, float(np.max(np.abs(r.variances() - h.variances())
                                                / np.maximum(np.abs(r.variances()), 1e-24))))
        info["check_against"] = str(spec["check_against"])
        info["max_relative_difference"] = worst
        if worst > 1e-9:
            raise ValueError(f"z-ee: tarball histograms differ from {delivered} (max rel. diff {worst:.3g})")
    return info
