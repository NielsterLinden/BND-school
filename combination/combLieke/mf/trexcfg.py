"""Read the channels' TRExFitter configs, adapt them to the combined likelihood, write the MultiFit.

The reader handles the plain `Key: value` HIST configs of this repository (block lines, indented
options, `%`/`#` comments outside quotes). `adapt_channel` is the whole combination model for one
channel; every change it makes to a channel config is listed in its docstring and nothing else is
touched (samples, binning, the channels' own DropBins, smoothing, symmetrisation, MC-statistics treatment).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from .paths import repo_path

KINDS = ("Job", "Fit", "Region", "Sample", "NormFactor", "ShapeFactor", "Systematic", "MultiFit", "Options")
#: options whose value TRExFitter shows to a human; quoted on output (ROOT '#' markup inside)
QUOTED = ("Label", "Title", "Category", "SubCategory", "VariableTitle", "ShortLabel", "LegendLabel")


@dataclass
class Block:
    kind: str
    name: str
    opts: dict = field(default_factory=dict)


def _strip_comment(line: str) -> str:
    return re.split(r"[%#](?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", line)[0].strip()


def read(path: Path) -> list[Block]:
    blocks: list[Block] = []
    for number, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("%", "#")):
            continue
        line = _strip_comment(line)
        key, sep, value = line.partition(":")
        if not sep:
            raise ValueError(f"{path}:{number}: cannot parse {raw!r}")
        key, value = key.strip(), value.strip().strip('"')
        if key in KINDS and not raw[:1].isspace():
            blocks.append(Block(key, value))
        elif blocks:
            if key in blocks[-1].opts:
                raise ValueError(f"{path}:{number}: option {key} given twice in {blocks[-1].kind} {blocks[-1].name}")
            blocks[-1].opts[key] = value
        else:
            raise ValueError(f"{path}:{number}: option before the first block")
    return blocks


def render(blocks: list[Block], header: str = "") -> str:
    def val(k, v):
        return f'"{v}"' if k.endswith(QUOTED) else str(v)
    head = "".join(f"% {h}\n" for h in header.splitlines()) + ("\n" if header else "")
    return head + "\n\n".join(f'{b.kind}: "{b.name}"\n' + "\n".join(f"  {k}: {val(k, v)}" for k, v in b.opts.items())
                              for b in blocks) + "\n"


def first(blocks, kind) -> Block:
    found = [b for b in blocks if b.kind == kind]
    if len(found) != 1:
        raise ValueError(f"expected exactly one {kind} block, found {len(found)}")
    return found[0]


def listopt(value: str) -> list[str]:
    return [v.strip().strip('"') for v in str(value).split(",") if v.strip()]


def _dig(d: dict, dotted: str):
    for part in dotted.split("."):
        d = d[part]
    return float(d)


def acceptance_nps(spec: dict | None) -> list[dict]:
    """[{name, rel, title}] for a channel's acceptance uncertainties, read from its metadata file."""
    if not spec:
        return []
    meta = json.loads(repo_path(spec["source"]).read_text())
    out = []
    for np_ in spec["nps"]:
        rel = _dig(meta, np_["key"])
        if "divide_by" in np_:
            rel /= _dig(meta, np_["divide_by"])
        out.append({"name": np_["name"], "rel": rel, "title": np_["title"], "source": f'{spec["source"]}: {np_["key"]}'})
    return out


#: speed/size settings for the combination's own h/w/f runs; plots are made in python
QUIET = {"SystControlPlots": "FALSE", "DoSummaryPlot": "FALSE", "DoPieChartPlot": "FALSE", "DoSignalRegionsPlot": "FALSE",
         "DebugLevel": "1", "ImageFormat": "png"}


def _check_empty(histo_path: Path, job: Block, region: Block, blocks: list[Block], bins: list[int]) -> None:
    """Refuse to drop a bin unless data and every sample are exactly zero there (1-based indices)."""
    import uproot
    f = uproot.open(histo_path / f'{job.opts["HistoFile"]}.root')
    for b in blocks:
        if b.kind != "Sample":
            continue
        regions = listopt(b.opts.get("Regions", region.name))
        if region.name not in regions and "all" not in regions:
            continue
        name = region.opts.get("HistoName", region.name) + b.opts.get("HistoNameSuff", "")
        values = f[name].values()
        for i in bins:
            if values[i - 1] != 0:
                raise ValueError(f"{region.name} bin {i} is not empty in {name} ({values[i - 1]}); refusing to drop it")


def adapt_channel(key: str, spec: dict, poi: dict, poi_name: str, workdir: Path, histo_path: Path,
                  fit: dict | None = None) -> tuple[list[Block], dict]:
    """The channel config as it enters the combined likelihood.

    Changes, and only these:
      * Job renamed to the channel key, OutputDir -> workdir, HistoPath -> histo_path (absolute),
        plotting switched off (QUIET), the POI renamed to `poi_name`;
      * the channel POI NormFactor renamed to `poi_name`, range = poi["range"];
      * a constant NormFactor `xsref_<key>` = poi["reference_pb"] / sigma_reference(channel) on the same
        samples, so that poi_name x poi["reference_pb"] is the channel's sigma(60 < m < 120 GeV) in pb;
      * systematics renamed per spec["rename_systematics"] (NuisanceParameter only, histograms unchanged);
      * one OVERALL systematic per acceptance source (spec["acceptance"]) on the signal samples, when the
        channel quotes sigma(60-120) = sigma_fid / A outside its own fit;
      * MINOS on the POI(s) only, plus any extra Fit options given in `fit` (config/channels.json "fit");
      * DropBins for the bins in spec["drop_empty_bins"], after checking that data and every sample are
        exactly zero there;
      * spec["split_shape_norm"]: shape and normalisation of HISTO systematics as separate parameters
        (split_shape_norm below); spec["overall"] {name: size}: replace the normalisation part of such a
        split systematic by a symmetric OVERALL (used only by the diagnostic variation).
    """
    blocks = read(repo_path(spec["config"]))
    job, fit_block = first(blocks, "Job"), first(blocks, "Fit")
    old_poi = job.opts["POI"]
    poi_nf = [b for b in blocks if b.kind == "NormFactor" and b.name == old_poi]
    if len(poi_nf) != 1:
        raise ValueError(f"{key}: no NormFactor for the POI {old_poi}")
    signal = poi_nf[0].opts["Samples"]
    if set(listopt(signal)) != set(spec["signal"]):
        raise ValueError(f"{key}: POI acts on {signal}, the manifest says {spec['signal']}")

    job.name = key
    job.opts.update(POI=poi_name, OutputDir=str(workdir), HistoPath=str(histo_path), **QUIET)
    job.opts.pop("Lumi", None)
    fit_block.opts.update(UseMinos=poi_name, **(fit or {}))
    for b in blocks:
        if b.kind == "Sample" and "HistoPath" in b.opts:
            raise ValueError(f"{key}: per-sample HistoPath is not supported")
    lo, hi = poi["range"]
    poi_nf[0].name = poi_name
    poi_nf[0].opts.update(Title=poi.get("title", "#mu_{Z}"), Nominal="1", Min=str(lo), Max=str(hi))

    factor = poi["reference_pb"] / spec["sigma_reference_pb"]
    blocks.append(Block("NormFactor", f"xsref_{key}", {
        "Title": f"sigma_ref / sigma_ref({key})", "Nominal": repr(factor), "Min": repr(factor), "Max": repr(factor),
        "Constant": "TRUE", "Samples": signal}))

    renames = dict(spec.get("rename_systematics", {}))
    seen = set()
    for b in blocks:
        if b.kind == "Systematic":
            seen.add(b.name)
            if b.name in renames:
                b.opts["NuisanceParameter"] = renames[b.name]
    missing = set(renames) - seen
    if missing:
        raise ValueError(f"{key}: rename_systematics names systematics the config does not have: {sorted(missing)}")

    for region_name, bins in spec.get("drop_empty_bins", {}).items():
        region = [b for b in blocks if b.kind == "Region" and b.name == region_name]
        if len(region) != 1:
            raise ValueError(f"{key}: drop_empty_bins names unknown region {region_name}")
        if "DropBins" in region[0].opts:
            raise ValueError(f"{key}: {region_name} already has DropBins")
        _check_empty(histo_path, job, region[0], blocks, bins)
        region[0].opts["DropBins"] = ",".join(str(i) for i in bins)

    split = split_shape_norm(key, blocks, spec.get("split_shape_norm"))
    for name, size in (spec.get("overall") or {}).items():
        norm = [b for b in blocks if b.kind == "Systematic" and b.name == name and "DropShapeIn" in b.opts]
        if len(norm) != 1:
            raise ValueError(f"{key}: 'overall' needs the normalisation part of a split systematic, {name} is not one")
        keep = {k: norm[0].opts[k] for k in ("Samples", "Category", "SubCategory", "NuisanceParameter") if k in norm[0].opts}
        norm[0].opts = {"Title": f'{norm[0].opts.get("Title", name)} (norm. #pm{100 * size:g}%)', "Type": "OVERALL",
                        "OverallUp": repr(size), "OverallDown": repr(-size), **keep}

    acc = acceptance_nps(spec.get("acceptance"))
    for a in acc:
        blocks.append(Block("Systematic", a["name"], {
            "Title": a["title"], "Type": "OVERALL", "OverallUp": repr(a["rel"]), "OverallDown": repr(-a["rel"]),
            "Samples": signal, "Category": "Acceptance"}))
    info = {"poi": poi_name, "signal": listopt(signal), "xsref_factor": factor, "acceptance": acc,
            "renamed": renames, "config": spec["config"], "histo_path": str(histo_path),
            "dropped_empty_bins": spec.get("drop_empty_bins", {}), "split_shape_norm": split,
            "overall_override": spec.get("overall") or {}}
    return blocks, info


def split_shape_norm(key: str, blocks: list[Block], which) -> list[str]:
    """Give the shape part of HISTO systematics its own, channel-specific nuisance parameter.

    `which` is None/False (nothing), True/"all" (every HISTO systematic) or a list of names. Each selected
    systematic becomes two blocks built from the same templates:
      <name>          DropShapeIn: all regions -> normalisation only, NuisanceParameter unchanged
                      (so still correlated by name with the other channels)
      <name>_<key>Shape  DropNorm: all regions    -> shape only, NuisanceParameter <NP>_<key>Shape
                      (TRExFitter reserves "shape_", "alpha", "gamma" in NP names)
    """
    if not which:
        return []
    regions = [b.name for b in blocks if b.kind == "Region"]
    histo = [b for b in blocks if b.kind == "Systematic" and b.opts.get("Type", "").upper() == "HISTO"]
    if which not in (True, "all"):
        unknown = set(which) - {b.name for b in histo}
        if unknown:
            raise ValueError(f"{key}: split_shape_norm names non-HISTO or unknown systematics {sorted(unknown)}")
        histo = [b for b in histo if b.name in which]
    done = []
    for b in histo:
        if "DropShapeIn" in b.opts or "DropNorm" in b.opts:
            raise ValueError(f"{key}: {b.name} already has DropShapeIn/DropNorm")
        np_name = b.opts.get("NuisanceParameter", b.name)
        shape = Block("Systematic", f"{b.name}_{key}Shape", dict(b.opts))
        shape.opts.update(NuisanceParameter=f"{np_name}_{key}Shape", DropNorm=",".join(regions),
                          Title=f'{b.opts.get("Title", b.name)} (shape, {key})')
        b.opts["DropShapeIn"] = ",".join(regions)
        blocks.insert(blocks.index(b) + 1, shape)
        done.append(b.name)
    return done


def multifit(name: str, fits: list[dict], poi_names: list[str], workdir: Path, poi: dict, scan: bool = True) -> list[Block]:
    lo, hi = poi["scan_range"]
    mf = Block("MultiFit", name, {
        "Label": "Z #rightarrow ll", "OutputDir": str(workdir), "Combine": "TRUE", "Compare": "FALSE",
        "POIName": ",".join(poi_names) if len(poi_names) > 1 else poi_names[0],
        "POIRange": ",".join(f"{poi['range'][0]}:{poi['range'][1]}" for _ in poi_names),
        "DataName": "obsData", "CmeLabel": "13 TeV", "LumiLabel": "16.4 fb^{-1}", "NumCPU": "8",
        "PlotCombCorrMatrix": "FALSE", **{k: v for k, v in (poi.get("fit") or {}).items() if not k.startswith("_")}})
    if scan:
        mf.opts.update(doLHscan=poi_names[0], LHscanMin=str(lo), LHscanMax=str(hi), LHscanSteps="41")
    out = [mf]
    for f in fits:
        out.append(Block("Fit", f["job"], {"ConfigFile": f["config"], "Directory": f["directory"], "Label": f["label"]}))
    return out
