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
SKIM_VERSION_V4 = "tautau-v4"

# v4 (docs/10-v4-plan.md): one skim for the mu tau_h, e tau_h and e mu channels and the trigger-efficiency
# side samples. Every event carries the bitmask `skim_cat` of the preselections it passes; a data stream only
# keeps its own categories (samples.SAMPLES[key]["stream"]), the simulation keeps every category.
CAT_TT, CAT_MT, CAT_ET, CAT_EM, CAT_EM_MU, CAT_EM_EL = 1, 2, 4, 8, 16, 32
STREAM_CATS = {"Tau": CAT_TT, "SingleMuon": CAT_MT | CAT_EM_MU, "SingleElectron": CAT_ET | CAT_EM_EL,
               "MuonEG": CAT_EM, None: CAT_TT | CAT_MT | CAT_ET | CAT_EM | CAT_EM_MU | CAT_EM_EL}

EVENT = ["run", "luminosityBlock", "event", "PV_npvs", "PV_npvsGood", "fixedGridRhoFastjetAll",
         "MET_pt", "MET_phi", "MET_sumEt", "MET_significance", "MET_covXX", "MET_covXY", "MET_covYY",
         "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY", "PuppiMET_pt", "PuppiMET_phi"]
FLAGS = list(config.MET_FILTERS)
HLT = ["HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
       "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", "HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg",
       "HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg", "HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg"]
HLT_V4 = HLT + config.MU_TRIGGERS + config.EL_TRIGGERS + config.EMU_TRIGGERS + [
    "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
    "HLT_Ele25_eta2p1_WPTight_Gsf", "HLT_Ele27_eta2p1_WPTight_Gsf", "HLT_Mu50", "HLT_IsoMu27"]
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
# v4 adds what the lepton scale factors, the electron energy scale and the b-tagging need
COLLECTIONS_V4_EXTRA = {"Muon": ["isGlobal", "isTracker", "isPFcand", "pfRelIso03_all"],
                        "Electron": ["deltaEtaSC", "mvaFall17V2noIso_WP80", "dEscaleUp", "dEscaleDown",
                                     "dEsigmaUp", "dEsigmaDown", "eCorr"],
                        "Jet": ["rawFactor", "area"]}
MC_COLLECTION_EXTRA = {"Tau": ["genPartFlav", "genPartIdx"], "Muon": ["genPartFlav"],
                       "Electron": ["genPartFlav"], "Jet": ["hadronFlavour"]}
GENPART_BRANCHES = ["GenPart_pt", "GenPart_pdgId", "GenPart_statusFlags"]     # top pT (ttbar samples only)
MC_EVENT = ["genWeight", "Pileup_nTrueInt", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up",
            "L1PreFiringWeight_Dn", "GenMET_pt", "GenMET_phi", "LHE_Njets", "LHE_NpNLO", "LHE_Vpt"]
# LHE_NpNLO = number of partons of the NLO matrix element (0, 1, 2): the jet-binned aMC@NLO samples are
# exclusive in it (LHE_Njets is not: the real-emission parton adds one). Generator-weight sums per bin
# for the jet-binned stitching (analysis.dy_norm).
NJET_STITCH_BINS = 3
MC_VECTORS = ["LHEScaleWeight", "PSWeight", "LHEPdfWeight"]
RUNS_BRANCHES = ["genEventCount", "genEventSumw", "genEventSumw2"]
PU_BINS = 100


def collection_fields(coll: str, is_mc: bool, v4: bool = False) -> list[str]:
    fields = list(COLLECTIONS[coll])
    if v4:
        fields += COLLECTIONS_V4_EXTRA.get(coll, [])
    if is_mc:
        fields += MC_COLLECTION_EXTRA.get(coll, [])
    return fields


def branch_list(is_mc: bool, keep_pdf: bool, v4: bool = False, top: bool = False) -> list[str]:
    names = EVENT + FLAGS + (HLT_V4 if v4 else HLT)
    for coll in COLLECTIONS:
        names += [f"n{coll}"] + [f"{coll}_{f}" for f in collection_fields(coll, is_mc, v4)]
    if is_mc:
        names += MC_EVENT + gen.LHE_BRANCHES + ["nGenVisTau"] + gen.VISTAU_BRANCHES
        names += [v for v in MC_VECTORS if keep_pdf or v != "LHEPdfWeight"]
        if top:
            names += ["nGenPart"] + GENPART_BRANCHES
    return names


def n_tau_candidates(ev) -> np.ndarray:
    dm = ev.Tau_decayMode
    good = ((ev.Tau_pt > config.SKIM_TAU_PT) & (abs(ev.Tau_eta) < config.SKIM_TAU_ETA)
            & ((dm == 0) | (dm == 1) | (dm == 10) | (dm == 11))
            & ((ev.Tau_idDeepTau2017v2p1VSjet & config.SKIM_TAU_VSJET_BIT) > 0)
            & ((ev.Tau_idDeepTau2017v2p1VSe & config.SKIM_TAU_VSE_BIT) > 0)
            & ((ev.Tau_idDeepTau2017v2p1VSmu & config.SKIM_TAU_VSMU_BIT) > 0))
    return ak.to_numpy(ak.sum(good, axis=1))


def _count(mask) -> np.ndarray:
    return ak.to_numpy(ak.sum(mask, axis=1))


def v4_categories(ev) -> np.ndarray:
    """Bitmask of the v4 preselections (CAT_*) per event (docs/10-v4-plan.md, "Skims").

    Objects are looser than the analysis ones so the isolation sidebands and energy-scale variations stay
    inside the skim: muons medium ID, pT > 22 (9 for the e mu trailing leg), I_rel(0.4) < 0.5; electrons MVA
    noIso 90%, pT > 25 (12), I_rel(0.3) < 0.5, conversion veto; tau_h pT > 18, |eta| < 2.3, decay mode
    0/1/10/11, VSjet >= VVVLoose, VSe >= VVLoose, VSmu >= VLoose.
    """
    dm = ev.Tau_decayMode
    tau = ((ev.Tau_pt > 18.0) & (abs(ev.Tau_eta) < 2.3) & ((dm == 0) | (dm == 1) | (dm == 10) | (dm == 11))
           & ((ev.Tau_idDeepTau2017v2p1VSjet & 1) > 0) & ((ev.Tau_idDeepTau2017v2p1VSe & 2) > 0)
           & ((ev.Tau_idDeepTau2017v2p1VSmu & 1) > 0))
    mu_base = (ev.Muon_mediumId & (abs(ev.Muon_eta) < 2.4) & (abs(ev.Muon_dxy) < 0.045) & (abs(ev.Muon_dz) < 0.2)
               & (ev.Muon_pfRelIso04_all < 0.5))
    el_base = (ev.Electron_mvaFall17V2noIso_WP90 & (abs(ev.Electron_eta) < 2.5) & (abs(ev.Electron_dxy) < 0.045)
               & (abs(ev.Electron_dz) < 0.2) & (ev.Electron_pfRelIso03_all < 0.5) & ev.Electron_convVeto
               & (ev.Electron_lostHits <= 1))
    n_tau = _count(tau)
    n_mu22, n_mu9 = _count(mu_base & (ev.Muon_pt > 22.0)), _count(mu_base & (ev.Muon_pt > 9.0))
    n_el25, n_el12 = _count(el_base & (ev.Electron_pt > 25.0)), _count(el_base & (ev.Electron_pt > 12.0))
    t_tt = io.trigger_or(ev, config.DITAU_TRIGGERS)
    t_mu = io.trigger_or(ev, config.MU_TRIGGERS)
    t_el = io.trigger_or(ev, config.EL_TRIGGERS)
    t_em = io.trigger_or(ev, config.EMU_TRIGGERS)
    cat = np.zeros(len(ev), dtype=np.uint8)
    cat |= np.where(t_tt & (n_tau_candidates(ev) >= 2), CAT_TT, 0).astype(np.uint8)
    cat |= np.where(t_mu & (n_mu22 >= 1) & (n_tau >= 1), CAT_MT, 0).astype(np.uint8)
    cat |= np.where(t_el & (n_el25 >= 1) & (n_tau >= 1), CAT_ET, 0).astype(np.uint8)
    cat |= np.where(t_em & (n_mu9 >= 1) & (n_el12 >= 1), CAT_EM, 0).astype(np.uint8)
    cat |= np.where(t_mu & (n_mu22 >= 1) & (n_el12 >= 1), CAT_EM_MU, 0).astype(np.uint8)
    cat |= np.where(t_el & (n_el25 >= 1) & (n_mu9 >= 1), CAT_EM_EL, 0).astype(np.uint8)
    return cat


def top_pt(ev):
    """(pT of the top, pT of the antitop) from the last-copy GenPart entries (0 where absent)."""
    if "GenPart_pdgId" not in ak.fields(ev):
        return np.zeros(len(ev), np.float32), np.zeros(len(ev), np.float32)
    last = (ev.GenPart_statusFlags & (1 << 13)) > 0
    out = []
    for pdg in (6, -6):
        sel = last & (ev.GenPart_pdgId == pdg)
        pt = ak.fill_none(ak.firsts(ev.GenPart_pt[sel]), 0.0)
        out.append(ak.to_numpy(pt).astype(np.float32))
    return out[0], out[1]


# ------------------------------------------------------------------------------ generator sums
def blank_gensums():
    g = {"n_events": 0.0, "sumw": 0.0, "sumw2": 0.0, "pu_true": np.zeros(PU_BINS)}
    for fl in ("ee", "mumu", "tautau"):
        g[f"sumw_lhe_{fl}"] = 0.0
        g[f"sumw_lhe_{fl}_60_120"] = 0.0
    g["sumw_fid"] = 0.0
    g["sumw2_fid"] = 0.0
    g["sumw_npnlo"] = np.zeros(NJET_STITCH_BINS)        # per LHE_NpNLO bin (0, 1, >= 2), all events
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
    if "LHE_NpNLO" in fields:
        nj = np.clip(ak.to_numpy(ev.LHE_NpNLO).astype(int), 0, NJET_STITCH_BINS - 1)
        g["sumw_npnlo"] += np.bincount(nj, weights=w, minlength=NJET_STITCH_BINS)
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
def build_output(ev, present, is_mc, sel, extra, v4: bool = False):
    e = ev[sel]
    out = {}
    for name in EVENT + FLAGS + (HLT_V4 if v4 else HLT) + (MC_EVENT if is_mc else []):
        if name in present:
            out[name] = ak.to_numpy(e[name])
    for coll in COLLECTIONS:
        cols = {f: e[f"{coll}_{f}"] for f in collection_fields(coll, is_mc, v4) if f"{coll}_{f}" in present}
        if not cols:
            continue
        rec = ak.zip(cols)
        if coll == "TrigObj":
            keep = rec.id == config.TRIGOBJ_TAU_ID
            if v4:
                keep = keep | (rec.id == config.TRIGOBJ_MU_ID) | (rec.id == config.TRIGOBJ_EL_ID)
            rec = rec[keep]
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


def process_file(source: str, out_path: Path, sample: dict, chunk_size: str = config.CHUNK_SIZE,
                 v4: bool = False) -> dict:
    """Skim one NanoAOD file into `out_path` (+ `<out_path>.json` provenance). `v4`: the four-channel
    preselection (`v4_categories`) instead of the tau_h tau_h one."""
    is_mc = sample["is_mc"]
    top = is_mc and sample.get("group") == "top"
    cats_keep = STREAM_CATS[None if is_mc else sample.get("stream")]
    t0 = time.time()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    grl = None if is_mc else io.load_grl()
    prov = {"source": str(source), "output": str(out_path), "sample": sample["key"], "is_mc": is_mc,
            "skim_version": SKIM_VERSION_V4 if v4 else SKIM_VERSION, "n_in": 0, "n_certified": 0, "n_trigger": 0,
            "n_out": 0, "missing_branches": [], "hostname": os.uname().nodename}
    if v4:
        prov["n_cat"] = {c: 0 for c in ("TT", "MT", "ET", "EM", "EM_MU", "EM_EL")}
    gsums = blank_gensums() if is_mc else None

    fh = uproot.open(source, timeout=600)
    tree = fh["Events"]
    available = set(tree.keys())
    wanted = branch_list(is_mc, sample.get("keep_pdf", False), v4=v4, top=top)
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
        for name in FLAGS + (HLT_V4 if v4 else HLT):
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
        keep &= ak.to_numpy(ev.PV_npvsGood >= 1)
        if v4:
            cat = v4_categories(ev) & cats_keep
            keep &= cat > 0
            prov["n_trigger"] += int(keep.sum())
            for c, bit in zip(("TT", "MT", "ET", "EM", "EM_MU", "EM_EL"), (CAT_TT, CAT_MT, CAT_ET, CAT_EM, CAT_EM_MU, CAT_EM_EL)):
                prov["n_cat"][c] += int(((cat & bit) > 0).sum())
            extra["skim_cat"] = cat
            if top:
                extra["gen_top_pt"], extra["gen_antitop_pt"] = top_pt(ev)
        else:
            keep &= io.trigger_or(ev, config.DITAU_TRIGGERS)
            prov["n_trigger"] += int(keep.sum())
            keep &= n_tau_candidates(ev) >= 2
        n_sel = int(keep.sum())
        if n_sel == 0:
            continue
        out = build_output(ev, present_now, is_mc, ak.Array(keep), extra, v4=v4)
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


def build_manifest(keys, base: Path | None = None, v4: bool = False) -> dict:
    from . import samples
    base = Path(base or (config.SKIM_DIR_V4 if v4 else config.SKIM_DIR))
    manifest = {"skim_version": SKIM_VERSION_V4 if v4 else SKIM_VERSION, "samples": {}}
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
            for c, n in prov.get("n_cat", {}).items():
                entry.setdefault("n_cat", {})[c] = entry.get("n_cat", {}).get(c, 0) + n
            entry["files"].append(p.name)
        manifest["samples"][key] = entry
    with open(base / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    return manifest
