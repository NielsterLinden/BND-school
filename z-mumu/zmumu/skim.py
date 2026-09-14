"""v2 skims: NanoAOD parents -> small analysis-ready files with NanoAOD branch names.

One output file per input file (resumable), written with uproot (`mktree`/`extend`) so the
branches keep their NanoAOD names (`nMuon`, `Muon_pt`, ...) and any NanoAOD-aware code
reads them. See docs/10-skims.md.

Preselection (data): certified lumisections, run 278820-284044, PV_npvsGood >= 1,
HLT_IsoMu24 || HLT_IsoTkMu24, and at least one category:
    A (bit 1)  >= 2 loose muons (pT > 20) with >= 1 above 24 GeV    SR, SS, T&P, fake factor
    B (bit 2)  >= 1 loose muon (pT > 24) and >= 1 electron (pT > 20)  e-mu region
    C (bit 4)  exactly 1 loose muon (pT > 24), MET < 30, mT < 30, >= 1 jet   single-mu fake CR,
               prescaled 1:10 (event % 10 == 0, `skim_prescale` = 10) when it is the only category
MC: same categories, no trigger requirement (kept as branches), no lumi mask, plus generator
information and a per-file `GenSums` tree with everything the C factor and the acceptance
need (sums over *all* generated events, before any selection).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from . import config, gen, io, objects

SKIM_VERSION = "v2.0"
RUN_MIN, RUN_MAX = 278820, 284044          # Run2016G + Run2016H
ERA_H_FIRST_RUN = 281613                  # first run of Run2016H (280385 is the last of G)

# ---------------------------------------------------------------------------- branches
EVENT = ["run", "luminosityBlock", "event", "PV_npvs", "PV_npvsGood", "fixedGridRhoFastjetAll",
         "MET_pt", "MET_phi", "PuppiMET_pt", "PuppiMET_phi"]
FLAGS = ["Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
         "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
         "Flag_BadPFMuonFilter", "Flag_BadPFMuonDzFilter", "Flag_eeBadScFilter"]
HLT = ["HLT_IsoMu24", "HLT_IsoTkMu24", "HLT_Mu50", "HLT_TkMu50", "HLT_IsoMu27", "HLT_Mu27",
       "HLT_Mu45_eta2p1", "HLT_Mu17", "HLT_Mu8", "HLT_Mu17_TrkIsoVVL", "HLT_Mu8_TrkIsoVVL",
       "HLT_Mu3_PFJet40", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
       "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"]
COLLECTIONS = {
    "Muon": ["pt", "eta", "phi", "mass", "charge", "tightId", "mediumId", "looseId", "highPtId",
             "pfRelIso04_all", "pfRelIso03_all", "tkRelIso", "miniPFRelIso_all", "pfIsoId",
             "dxy", "dz", "dxyErr", "dzErr", "sip3d", "ip3d", "isGlobal", "isTracker", "isPFcand",
             "nTrackerLayers", "nStations", "jetIdx", "ptErr", "fsrPhotonIdx"],
    "Electron": ["pt", "eta", "phi", "mass", "charge", "cutBased", "pfRelIso03_all", "dxy", "dz",
                 "mvaFall17V2Iso_WP90", "convVeto", "lostHits", "deltaEtaSC"],
    "Jet": ["pt", "eta", "phi", "mass", "jetId", "puId", "btagDeepFlavB", "muonIdx1", "muonIdx2",
            "nMuons"],
    "FsrPhoton": ["pt", "eta", "phi", "relIso03", "dROverEt2", "muonIdx"],
    "TrigObj": ["pt", "eta", "phi", "id", "filterBits", "l1pt", "l2pt"],
}
MC_COLLECTION_EXTRA = {"Muon": ["genPartFlav", "genPartIdx"], "Electron": ["genPartFlav"],
                       "Jet": ["hadronFlavour"]}
MC_EVENT = ["genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up",
            "L1PreFiringWeight_Dn", "L1PreFiringWeight_Muon_Nom", "L1PreFiringWeight_Muon_StatUp",
            "L1PreFiringWeight_Muon_StatDn", "L1PreFiringWeight_Muon_SystUp",
            "L1PreFiringWeight_Muon_SystDn", "L1PreFiringWeight_ECAL_Nom",
            "L1PreFiringWeight_ECAL_Up", "L1PreFiringWeight_ECAL_Dn"]
MC_VECTORS = ["LHEScaleWeight", "PSWeight", "LHEPdfWeight"]     # jagged in NanoAOD
MC_GEN = gen.LHE_BRANCHES + gen.DRESSED_BRANCHES
RUNS_BRANCHES = ["genEventCount", "genEventSumw", "genEventSumw2", "LHEScaleSumw", "LHEPdfSumw"]

# skim-level object definitions (looser than the analysis)
MU_LOOSE_PT, MU_TAG_PT = 20.0, 24.0
EL_PT, EL_ETA, EL_ID = 20.0, 2.5, 1
JET_PT, JET_ETA, JET_ID, JET_DR = 30.0, 2.4, 2, 0.4
FF_MET_MAX, FF_MT_MAX = 30.0, 30.0
PRESCALE_C = 10
CAT_A, CAT_B, CAT_C = 1, 2, 4
PU_BINS = 100


def branch_list(is_mc: bool, keep_pdf: bool) -> list[str]:
    names = EVENT + FLAGS + HLT
    for coll, fields in COLLECTIONS.items():
        extra = MC_COLLECTION_EXTRA.get(coll, []) if is_mc else []
        names += [f"n{coll}"] + [f"{coll}_{f}" for f in fields + extra]
    if is_mc:
        names += MC_EVENT + MC_GEN + [v for v in MC_VECTORS if keep_pdf or v != "LHEPdfWeight"]
    return names


# ---------------------------------------------------------------------------- selection
def _mt(pt, phi, met, met_phi):
    return np.sqrt(np.maximum(2.0 * pt * met * (1.0 - np.cos(objects.delta_phi(phi, met_phi))), 0.0))


def categories(ev):
    """(category bitmask, n loose muons) per event; see the module docstring."""
    mu_loose = ((ev.Muon_isGlobal | ev.Muon_isTracker) & (ev.Muon_pt > MU_LOOSE_PT)
                & (abs(ev.Muon_eta) < config.MU_ETA_MAX))
    n_loose = ak.to_numpy(ak.sum(mu_loose, axis=1))
    n_tag = ak.to_numpy(ak.sum(mu_loose & (ev.Muon_pt > MU_TAG_PT), axis=1))
    el = (ev.Electron_pt > EL_PT) & (abs(ev.Electron_eta) < EL_ETA) & (ev.Electron_cutBased >= EL_ID)
    n_el = ak.to_numpy(ak.sum(el, axis=1))

    cat_a = (n_loose >= 2) & (n_tag >= 1)
    cat_b = (n_tag >= 1) & (n_el >= 1)

    # single-muon + jet region: leading loose muon, W suppressed by MET and mT
    mu = ak.zip({"pt": ev.Muon_pt, "eta": ev.Muon_eta, "phi": ev.Muon_phi})[mu_loose]
    mu = ak.pad_none(mu[ak.argsort(mu.pt, axis=1, ascending=False)], 1, axis=1)[:, 0]
    lead_pt = ak.to_numpy(ak.fill_none(mu.pt, 0.0))
    lead_eta = ak.fill_none(mu.eta, 0.0)
    lead_phi = ak.fill_none(mu.phi, 0.0)
    met = ak.to_numpy(ev.MET_pt)
    mt = _mt(lead_pt, ak.to_numpy(lead_phi), met, ak.to_numpy(ev.MET_phi))
    jet = ((ev.Jet_pt > JET_PT) & (abs(ev.Jet_eta) < JET_ETA) & (ev.Jet_jetId >= JET_ID)
           & (objects.delta_r(ev.Jet_eta, ev.Jet_phi, lead_eta, lead_phi) > JET_DR))
    n_jet = ak.to_numpy(ak.sum(jet, axis=1))
    cat_c = (n_loose == 1) & (n_tag == 1) & (met < FF_MET_MAX) & (mt < FF_MT_MAX) & (n_jet >= 1)

    cat = (cat_a.astype(np.uint8) * CAT_A + cat_b.astype(np.uint8) * CAT_B
           + cat_c.astype(np.uint8) * CAT_C)
    return cat, n_loose


# ---------------------------------------------------------------------------- gen sums
def blank_gensums():
    g = {"n_events": 0, "sumw": 0.0, "sumw2": 0.0,
         "pu_true": np.zeros(PU_BINS), "h_lhe_mll": np.zeros(200)}
    for fl in ("ee", "mumu", "tautau"):
        g[f"sumw_lhe_{fl}"] = 0.0
        g[f"sumw_lhe_{fl}_60_120"] = 0.0
    for name in ("fid_dressed", "fid_born"):
        g[f"sumw_{name}"] = 0.0
        g[f"sumw2_{name}"] = 0.0
    g["vector_sizes"] = {}
    return g


def _regular(arr):
    counts = ak.to_numpy(ak.num(arr, axis=1))
    if counts.size == 0 or counts.min() != counts.max() or counts.max() == 0:
        return None
    return ak.to_numpy(ak.to_regular(arr, axis=1)).astype(np.float64)


def accumulate_gensums(g, ev, flav, mll, flags):
    """Sums over all events of a chunk (called before any selection)."""
    w = ak.to_numpy(ev.genWeight).astype(np.float64)
    g["n_events"] += len(ev)
    g["sumw"] += float(w.sum())
    g["sumw2"] += float((w ** 2).sum())
    g["pu_true"] += np.histogram(ak.to_numpy(ev.Pileup_nTrueInt), bins=PU_BINS, range=(0, PU_BINS),
                                 weights=w)[0]
    g["h_lhe_mll"] += np.histogram(mll[flav == 13], bins=200, range=(0, 200), weights=w[flav == 13])[0]
    selections = {}
    for code, fl in ((11, "ee"), (13, "mumu"), (15, "tautau")):
        sel = flav == code
        g[f"sumw_lhe_{fl}"] += float(w[sel].sum())
        win = sel & (mll > config.MASS_LO) & (mll < config.MASS_HI)
        g[f"sumw_lhe_{fl}_60_120"] += float(w[win].sum())
        if fl == "mumu":
            selections["lhe_mumu"] = sel
            selections["lhe_mumu_60_120"] = win
    for name in ("fid_dressed", "fid_born"):
        sel = flags[name]
        g[f"sumw_{name}"] += float(w[sel].sum())
        g[f"sumw2_{name}"] += float((w[sel] ** 2).sum())
        selections[name] = sel
    selections["all"] = np.ones(len(ev), dtype=bool)
    fields = set(ak.fields(ev))
    for branch, vname in (("LHEPdfWeight", "pdf"), ("LHEScaleWeight", "scale"), ("PSWeight", "ps")):
        if branch not in fields:
            continue
        v = _regular(ev[branch])
        if v is None:
            continue
        g["vector_sizes"][vname] = v.shape[1]
        for sname, sel in selections.items():
            key = f"{vname}_{sname}"
            g[key] = g.get(key, 0.0) + (w[sel, None] * v[sel]).sum(axis=0)
    return g


# ---------------------------------------------------------------------------- one file
def build_output(ev, present, is_mc, sel, extra_scalars):
    """The dict of arrays written for the selected events of one chunk."""
    e = ev[sel]
    out = {}
    for name in EVENT + FLAGS + HLT + (MC_EVENT if is_mc else []):
        if name in present:
            out[name] = ak.to_numpy(e[name])
    for coll, fields in COLLECTIONS.items():
        extra = MC_COLLECTION_EXTRA.get(coll, []) if is_mc else []
        cols = {f: e[f"{coll}_{f}"] for f in fields + extra if f"{coll}_{f}" in present}
        if not cols:
            continue
        rec = ak.zip(cols)
        if coll == "TrigObj":
            rec = rec[rec.id == 13]          # only muon trigger objects are used
        out[coll] = rec
    if is_mc:
        for name in MC_VECTORS:
            if name in present:
                out[name] = e[name]
        for coll, fields in (("LHEPart", ["pt", "eta", "phi", "mass", "pdgId", "status"]),
                             ("GenDressedLepton", ["pt", "eta", "phi", "mass", "pdgId", "hasTauAnc"])):
            cols = {f: e[f"{coll}_{f}"] for f in fields if f"{coll}_{f}" in present}
            if cols:
                out[coll] = ak.zip(cols)
    sel_np = ak.to_numpy(sel)
    for name, arr in extra_scalars.items():
        out[name] = np.asarray(arr)[sel_np]
    return out


def process_file(source: str, out_path: Path, sample: dict, chunk_size: str = "200 MB") -> dict:
    """Skim one NanoAOD file into `out_path` (+ `<out_path>.json`); returns the provenance."""
    is_mc = sample["is_mc"]
    keep_pdf = sample.get("keep_pdf", False)
    t0 = time.time()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    grl = io.load_grl() if not is_mc else None

    prov = {"source": str(source), "output": str(out_path), "sample": sample["key"], "is_mc": is_mc,
            "skim_version": SKIM_VERSION, "n_in": 0, "n_certified": 0, "n_trigger": 0, "n_out": 0,
            "n_cat": {"A": 0, "B": 0, "C": 0, "C_only_prescaled_in": 0}, "missing_branches": [],
            "n_written_chunks": 0}
    gsums = blank_gensums() if is_mc else None
    runs = {}

    fh = uproot.open(source)
    tree = fh["Events"]
    available = set(tree.keys())
    wanted = branch_list(is_mc, keep_pdf)
    present = [b for b in wanted if b in available]
    prov["missing_branches"] = sorted(set(wanted) - available)
    if is_mc and "Runs" in fh:
        rt = fh["Runs"]
        for b in RUNS_BRANCHES:
            if b in rt.keys():
                arr = rt[b].array(library="ak")
                runs[b] = np.asarray(ak.sum(arr, axis=0))
    fout = uproot.recreate(tmp_path, compression=uproot.ZSTD(5))
    have_tree = False

    for ev in tree.iterate(present, step_size=chunk_size, library="ak"):
        n = len(ev)
        prov["n_in"] += n
        # missing scalar flags/triggers -> False so downstream code can rely on the names
        for name in FLAGS + HLT:
            if name not in ak.fields(ev):
                ev = ak.with_field(ev, ak.Array(np.zeros(n, dtype=bool)), name)
        present_now = set(ak.fields(ev))

        extra = {}
        if is_mc:
            flav, mll, _ = gen.lhe_flavour(ev)
            flags = gen.fiducial_flags(ev, flav, mll)
            accumulate_gensums(gsums, ev, flav, mll, flags)
            keep = np.ones(n, dtype=bool)
            extra.update(gen_lhe_flavour=flav.astype(np.uint8), gen_mll_lhe=mll.astype(np.float32),
                         gen_fid_dressed=flags["fid_dressed"], gen_fid_born=flags["fid_born"],
                         gen_m_dressed=flags["m_dressed"].astype(np.float32),
                         gen_pt1_dressed=flags["pt1_dressed"].astype(np.float32),
                         gen_pt2_dressed=flags["pt2_dressed"].astype(np.float32))
        else:
            run = ak.to_numpy(ev.run)
            keep = (run >= RUN_MIN) & (run <= RUN_MAX) & io.lumi_mask(ev.run, ev.luminosityBlock, grl)
            prov["n_certified"] += int(keep.sum())
            keep &= ak.to_numpy(io.trigger_or(ev, config.TRIGGERS))
            prov["n_trigger"] += int(keep.sum())
        keep &= ak.to_numpy(ev.PV_npvsGood >= 1)

        cat, _ = categories(ev)
        cat = np.where(keep, cat, 0).astype(np.uint8)
        only_c = cat == CAT_C
        prov["n_cat"]["C_only_prescaled_in"] += int(only_c.sum())
        prescale = np.ones(n, dtype=np.uint8)
        prescale[only_c] = PRESCALE_C
        pass_prescale = ~only_c | (ak.to_numpy(ev.event) % PRESCALE_C == 0)
        sel = (cat > 0) & pass_prescale
        for bit, name in ((CAT_A, "A"), (CAT_B, "B"), (CAT_C, "C")):
            prov["n_cat"][name] += int(((cat & bit) > 0)[sel].sum())
        run = ak.to_numpy(ev.run)
        era = (np.where(run >= ERA_H_FIRST_RUN, 1, 0).astype(np.uint8) if not is_mc
               else np.full(n, 255, np.uint8))
        extra.update(skim_cat=cat, skim_prescale=prescale, era=era)

        n_sel = int(sel.sum())
        if n_sel == 0:
            continue
        out = build_output(ev, present_now, is_mc, ak.Array(sel), extra)
        if not have_tree:
            fout.mktree("Events", out, counter_name=lambda c: "n" + c,
                        field_name=lambda o, i: i if o == "" else o + "_" + i)
            have_tree = True
        else:
            fout["Events"].extend(out)
        prov["n_out"] += n_sel
        prov["n_written_chunks"] += 1

    if is_mc:
        gs = {}
        for k, v in gsums.items():
            if k == "vector_sizes":
                continue
            gs[k] = np.asarray([v], dtype=np.float64) if np.ndim(v) == 0 else np.asarray(v, dtype=np.float64)[None, :]
        fout.mktree("GenSums", gs)            # explicit TTree (a dict assignment would write an RNTuple)
        if runs:
            fout.mktree("Runs", {k: (np.asarray([v], dtype=np.float64) if np.ndim(v) == 0
                                     else np.asarray(v, dtype=np.float64)[None, :]) for k, v in runs.items()})
        prov["gensums"] = {k: float(v) for k, v in gsums.items() if np.ndim(v) == 0 and k != "vector_sizes"}
        prov["gensums"]["vector_sizes"] = gsums["vector_sizes"]
        prov["runs"] = {k: float(v) for k, v in runs.items() if np.ndim(v) == 0}
    fout.close()
    fh.close()
    tmp_path.rename(out_path)
    prov["has_events_tree"] = have_tree
    prov["seconds"] = round(time.time() - t0, 1)
    prov["output_bytes"] = out_path.stat().st_size
    with open(str(out_path) + ".json", "w") as f:
        json.dump(prov, f, indent=1)
    return prov


# ---------------------------------------------------------------------------- bookkeeping
def skim_dir(base=None) -> Path:
    import os
    return Path(base or os.environ.get("BND_SKIM_DIR", "/data/atlas/users/nterlind/BND-school-cache/skims_v2"))


def skim_files(key: str, base=None) -> list[Path]:
    """Skim files of a sample that have a completed sidecar."""
    d = skim_dir(base) / key
    return sorted(p for p in d.glob("*.root") if Path(str(p) + ".json").exists())


def build_manifest(keys, base=None) -> dict:
    """Aggregate the sidecars of every skimmed file into <skim_dir>/manifest.json."""
    from . import samples
    base = skim_dir(base)
    manifest = {"skim_version": SKIM_VERSION, "samples": {}}
    for key in keys:
        files = skim_files(key, base)
        entry = {"recid": samples.SAMPLES[key]["recid"], "is_mc": samples.SAMPLES[key]["is_mc"],
                 "xsec_pb": samples.SAMPLES[key]["xsec_pb"], "n_files": len(files),
                 "n_in": 0, "n_out": 0, "bytes": 0, "sumw": 0.0, "genEventSumw": 0.0, "files": []}
        for p in files:
            prov = json.load(open(str(p) + ".json"))
            entry["n_in"] += prov["n_in"]
            entry["n_out"] += prov["n_out"]
            entry["bytes"] += prov["output_bytes"]
            if prov["is_mc"]:
                entry["sumw"] += prov["gensums"]["sumw"]
                entry["genEventSumw"] += prov.get("runs", {}).get("genEventSumw", 0.0)
            entry["files"].append(p.name)
        manifest["samples"][key] = entry
    with open(base / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    return manifest


def load_gensums(key: str, base=None) -> dict:
    """Sum of the GenSums trees of all skim files of an MC sample."""
    total = {}
    for p in skim_files(key, base):
        with uproot.open(p) as f:
            if "GenSums" not in f:
                continue
            arrays = f["GenSums"].arrays(library="ak")
            for name in arrays.fields:
                val = np.sum(np.asarray(arrays[name]), axis=0)
                total[name] = total.get(name, 0.0) + val
    return total
