"""Skims for tau_h tau_h: NanoAOD parents -> small files with NanoAOD branch names (docs/03-skims.md).

One output file per parent file (resumable), written with uproot so the branches keep their
NanoAOD names (`nTau`, `Tau_pt`, ...) and any NanoAOD-aware code reads them.

Preselection
  data : certified lumisection, run 278820-284044, PV_npvsGood >= 1, di-tau trigger OR, and
         >= 2 tau candidates with pT > 35, |eta| < 2.3, decay mode 0/1/10/11 and DeepTau2017v2p1
         VSjet >= VVVLoose, VSe >= VVVLoose, VSmu >= VLoose.
  MC   : the same (the trigger is required as well; scale factors correct the simulated efficiency),
         plus generator information, and a `GenSums` tree with sums over *all* generated events of the
         file, before any selection (normalisation, pileup profile, acceptance and its uncertainties).
Only tau trigger objects (TrigObj id 15) are kept.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from . import config, gen, io

SKIM_VERSION = "tautau-v1"

EVENT = ["run", "luminosityBlock", "event", "PV_npvs", "PV_npvsGood", "fixedGridRhoFastjetAll",
         "MET_pt", "MET_phi", "MET_sumEt", "MET_significance", "MET_covXX", "MET_covXY", "MET_covYY",
         "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY", "PuppiMET_pt", "PuppiMET_phi"]
FLAGS = list(config.MET_FILTERS)
HLT = ["HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
       "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", "HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg",
       "HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg", "HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg"]
COLLECTIONS = {
    "Tau": ["pt", "eta", "phi", "mass", "charge", "decayMode", "dxy", "dz", "jetIdx",
            "idDecayModeOldDMs", "idDeepTau2017v2p1VSjet", "idDeepTau2017v2p1VSe", "idDeepTau2017v2p1VSmu",
            "rawDeepTau2017v2p1VSjet", "rawDeepTau2017v2p1VSe", "rawDeepTau2017v2p1VSmu",
            "chargedIso", "neutralIso", "photonsOutsideSignalCone", "leadTkPtOverTauPt"],
    "Muon": ["pt", "eta", "phi", "mass", "charge", "looseId", "mediumId", "tightId", "pfRelIso04_all", "dxy", "dz"],
    "Electron": ["pt", "eta", "phi", "mass", "charge", "cutBased", "mvaFall17V2noIso_WP90",
                 "mvaFall17V2Iso_WP90", "pfRelIso03_all", "dxy", "dz", "convVeto", "lostHits"],
    "Jet": ["pt", "eta", "phi", "mass", "jetId", "puId", "btagDeepFlavB"],
    "TrigObj": ["pt", "eta", "phi", "id", "filterBits"],
}
MC_COLLECTION_EXTRA = {"Tau": ["genPartFlav", "genPartIdx"], "Muon": ["genPartFlav"],
                       "Electron": ["genPartFlav"], "Jet": ["hadronFlavour"]}
MC_EVENT = ["genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up",
            "L1PreFiringWeight_Dn", "GenMET_pt", "GenMET_phi"]
MC_VECTORS = ["LHEScaleWeight", "PSWeight", "LHEPdfWeight"]
RUNS_BRANCHES = ["genEventCount", "genEventSumw", "genEventSumw2"]
PU_BINS = 100


def branch_list(is_mc: bool, keep_pdf: bool) -> list[str]:
    names = EVENT + FLAGS + HLT
    for coll, fields in COLLECTIONS.items():
        extra = MC_COLLECTION_EXTRA.get(coll, []) if is_mc else []
        names += [f"n{coll}"] + [f"{coll}_{f}" for f in fields + extra]
    if is_mc:
        names += MC_EVENT + gen.LHE_BRANCHES + ["nGenVisTau"] + gen.VISTAU_BRANCHES
        names += [v for v in MC_VECTORS if keep_pdf or v != "LHEPdfWeight"]
    return names


def n_tau_candidates(ev) -> np.ndarray:
    dm = ev.Tau_decayMode
    good = ((ev.Tau_pt > config.SKIM_TAU_PT) & (abs(ev.Tau_eta) < config.SKIM_TAU_ETA)
            & ((dm == 0) | (dm == 1) | (dm == 10) | (dm == 11))
            & ((ev.Tau_idDeepTau2017v2p1VSjet & config.SKIM_TAU_VSJET_BIT) > 0)
            & ((ev.Tau_idDeepTau2017v2p1VSe & config.SKIM_TAU_VSE_BIT) > 0)
            & ((ev.Tau_idDeepTau2017v2p1VSmu & config.SKIM_TAU_VSMU_BIT) > 0))
    return ak.to_numpy(ak.sum(good, axis=1))


# ------------------------------------------------------------------------------ generator sums
def blank_gensums():
    g = {"n_events": 0.0, "sumw": 0.0, "sumw2": 0.0, "pu_true": np.zeros(PU_BINS)}
    for fl in ("ee", "mumu", "tautau"):
        g[f"sumw_lhe_{fl}"] = 0.0
        g[f"sumw_lhe_{fl}_60_120"] = 0.0
    g["sumw_fid"] = 0.0
    g["sumw2_fid"] = 0.0
    return g


def _regular(arr):
    counts = ak.to_numpy(ak.num(arr, axis=1))
    if counts.size == 0 or counts.min() != counts.max() or counts.max() == 0:
        return None
    return ak.to_numpy(ak.to_regular(arr, axis=1)).astype(np.float64)


def accumulate_gensums(g, ev, flav, mll, fid):
    w = ak.to_numpy(ev.genWeight).astype(np.float64)
    g["n_events"] += len(ev)
    g["sumw"] += float(w.sum())
    g["sumw2"] += float((w ** 2).sum())
    g["pu_true"] += np.histogram(ak.to_numpy(ev.Pileup_nTrueInt), bins=PU_BINS, range=(0, PU_BINS), weights=w)[0]
    sels = {"all": np.ones(len(ev), dtype=bool), "fid": fid}
    for code, fl in ((11, "ee"), (13, "mumu"), (15, "tautau")):
        s = flav == code
        win = s & (mll > config.MASS_LO) & (mll < config.MASS_HI)
        g[f"sumw_lhe_{fl}"] += float(w[s].sum())
        g[f"sumw_lhe_{fl}_60_120"] += float(w[win].sum())
        if fl == "tautau":
            sels["lhe_tautau"] = s
            sels["lhe_tautau_60_120"] = win
    g["sumw_fid"] += float(w[fid].sum())
    g["sumw2_fid"] += float((w[fid] ** 2).sum())
    fields = set(ak.fields(ev))
    for branch, vname in (("LHEScaleWeight", "scale"), ("PSWeight", "ps"), ("LHEPdfWeight", "pdf")):
        if branch not in fields:
            continue
        v = _regular(ev[branch])
        if v is None:
            continue
        for sname, sel in sels.items():
            key = f"{vname}_{sname}"
            add = (w[sel, None] * v[sel]).sum(axis=0)
            g[key] = g[key] + add if key in g and np.shape(g[key]) == add.shape else add
    return g


# ------------------------------------------------------------------------------ one file
def build_output(ev, present, is_mc, sel, extra):
    e = ev[sel]
    out = {}
    for name in EVENT + FLAGS + HLT + (MC_EVENT if is_mc else []):
        if name in present:
            out[name] = ak.to_numpy(e[name])
    for coll, fields in COLLECTIONS.items():
        extra_f = MC_COLLECTION_EXTRA.get(coll, []) if is_mc else []
        cols = {f: e[f"{coll}_{f}"] for f in fields + extra_f if f"{coll}_{f}" in present}
        if not cols:
            continue
        rec = ak.zip(cols)
        if coll == "TrigObj":
            rec = rec[rec.id == config.TRIGOBJ_TAU_ID]
        out[coll] = rec
    if is_mc:
        for name in MC_VECTORS:
            if name in present:
                out[name] = e[name]
        for coll, branches in (("LHEPart", gen.LHE_BRANCHES), ("GenVisTau", gen.VISTAU_BRANCHES)):
            cols = {b.split("_", 1)[1]: e[b] for b in branches if b in present}
            if cols:
                out[coll] = ak.zip(cols)
    sel_np = ak.to_numpy(sel)
    for name, arr in extra.items():
        out[name] = np.asarray(arr)[sel_np]
    return out


def process_file(source: str, out_path: Path, sample: dict, chunk_size: str = config.CHUNK_SIZE) -> dict:
    """Skim one NanoAOD file into `out_path` (+ `<out_path>.json` provenance)."""
    is_mc = sample["is_mc"]
    t0 = time.time()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    grl = None if is_mc else io.load_grl()
    prov = {"source": str(source), "output": str(out_path), "sample": sample["key"], "is_mc": is_mc,
            "skim_version": SKIM_VERSION, "n_in": 0, "n_certified": 0, "n_trigger": 0, "n_out": 0,
            "missing_branches": [], "hostname": os.uname().nodename}
    gsums = blank_gensums() if is_mc else None

    fh = uproot.open(source, timeout=600)
    tree = fh["Events"]
    available = set(tree.keys())
    wanted = branch_list(is_mc, sample.get("keep_pdf", False))
    present = [b for b in wanted if b in available]
    prov["missing_branches"] = sorted(set(wanted) - available)
    runs = {}
    if is_mc and "Runs" in [k.split(";")[0] for k in fh.keys(recursive=False)]:
        rt = fh["Runs"]
        for b in RUNS_BRANCHES:
            if b in rt.keys():
                runs[b] = float(np.sum(rt[b].array(library="np")))
    fout = uproot.recreate(tmp_path, compression=uproot.ZSTD(5))
    have_tree = False
    for ev in tree.iterate(present, step_size=chunk_size, library="ak"):
        n = len(ev)
        prov["n_in"] += n
        for name in FLAGS + HLT:
            if name not in ak.fields(ev):
                ev = ak.with_field(ev, ak.Array(np.zeros(n, dtype=bool)), name)
        present_now = set(ak.fields(ev))
        extra = {}
        if is_mc:
            flav, mll = gen.lhe_flavour(ev)
            fid = gen.fiducial(ev, flav, mll)
            accumulate_gensums(gsums, ev, flav, mll, fid["fid"])
            keep = np.ones(n, dtype=bool)
            extra.update(gen_lhe_flavour=flav, gen_mll_lhe=mll.astype(np.float32), gen_fid=fid["fid"],
                         gen_n_vistau=fid["n_vis"], gen_vis1_pt=fid["vis1_pt"], gen_vis2_pt=fid["vis2_pt"])
        else:
            run = ak.to_numpy(ev.run)
            keep = (run >= config.RUN_MIN) & (run <= config.RUN_MAX) & io.lumi_mask(ev.run, ev.luminosityBlock, grl)
            prov["n_certified"] += int(keep.sum())
            extra["era"] = np.where(run >= config.ERA_H_FIRST_RUN, 1, 0).astype(np.uint8)
        keep &= io.trigger_or(ev, config.DITAU_TRIGGERS)
        prov["n_trigger"] += int(keep.sum())
        keep &= ak.to_numpy(ev.PV_npvsGood >= 1)
        keep &= n_tau_candidates(ev) >= 2
        n_sel = int(keep.sum())
        if n_sel == 0:
            continue
        out = build_output(ev, present_now, is_mc, ak.Array(keep), extra)
        if not have_tree:
            # mktree with arrays (not types) creates the tree *and* writes this first chunk
            fout.mktree("Events", out, counter_name=lambda c: "n" + c,
                        field_name=lambda o, i: i if o == "" else o + "_" + i)
            have_tree = True
        else:
            fout["Events"].extend(out)
        prov["n_out"] += n_sel
    if is_mc:
        # mktree, not fout["GenSums"] = {...}: uproot >= 5.6 writes a dict assignment as an RNTuple
        fout.mktree("GenSums", {k: (np.asarray([v], dtype=np.float64) if np.ndim(v) == 0
                                    else np.asarray(v, dtype=np.float64)[None, :]) for k, v in gsums.items()})
        prov["gensums"] = {k: float(v) for k, v in gsums.items() if np.ndim(v) == 0}
        prov["runs"] = runs
    fout.close()
    fh.close()
    tmp_path.rename(out_path)
    prov["has_events_tree"] = have_tree
    prov["seconds"] = round(time.time() - t0, 1)
    prov["output_bytes"] = out_path.stat().st_size
    with open(str(out_path) + ".json", "w") as f:
        json.dump(prov, f, indent=1)
    return prov


# ------------------------------------------------------------------------------ bookkeeping
def skim_files(key: str, base: Path | None = None) -> list[Path]:
    d = Path(base or config.SKIM_DIR) / key
    return sorted(p for p in d.glob("*.root") if Path(str(p) + ".json").exists())


def load_gensums(key: str, base: Path | None = None) -> dict:
    """Sum of the GenSums of all skim files of an MC sample (arrays for the vector sums)."""
    total = {}
    for p in skim_files(key, base):
        with uproot.open(p) as f:
            if "GenSums" not in [k.split(";")[0] for k in f.keys()]:
                continue
            for name, arr in f["GenSums"].arrays(library="np").items():
                val = np.sum(arr, axis=0)
                total[name] = total.get(name, 0.0) + val
    return total


def build_manifest(keys, base: Path | None = None) -> dict:
    from . import samples
    base = Path(base or config.SKIM_DIR)
    manifest = {"skim_version": SKIM_VERSION, "samples": {}}
    for key in keys:
        files = skim_files(key, base)
        s = samples.SAMPLES[key]
        entry = {"recid": s["recid"], "is_mc": s["is_mc"], "xsec_pb": s["xsec_pb"], "n_files": len(files),
                 "n_files_parent": len(samples.read_filelist(key)), "n_in": 0, "n_out": 0, "bytes": 0,
                 "sumw": 0.0, "files": []}
        for p in files:
            prov = json.load(open(str(p) + ".json"))
            entry["n_in"] += prov["n_in"]
            entry["n_out"] += prov["n_out"]
            entry["bytes"] += prov["output_bytes"]
            if prov["is_mc"]:
                entry["sumw"] += prov["gensums"]["sumw"]
            entry["files"].append(p.name)
        manifest["samples"][key] = entry
    with open(base / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    return manifest
