#!/usr/bin/env python
"""Freeze the Z -> mumu signal-region m(mumu) stack in four correction stages (one new pass over the MC skims).

    source setup.sh
    python presentation/data/extract_zmumu_sr_stack.py --samples WW ST_tW_top --max-files 2 --check-fill-chunk \
           --parts-dir presentation/work/extract/smoke                       # smoke test (~20 s)
    python presentation/data/extract_zmumu_sr_stack.py --workers 8            # full run (~3-5 min)
    python presentation/data/extract_zmumu_sr_stack.py --check-only           # re-run the assertions on the JSON

Stages (prompt-prompt MC only in every stage, momentum calibration applied in every stage):
    raw        genWeight x sigma L / sum(genWeight)                          -- the uncorrected simulation
    pileup     raw x w_PU(nTrueInt)
    prefiring  pileup x L1PreFiringWeight_Nom
    nominal    prefiring x SF_ID(mu1) SF_ID(mu2) SF_iso(mu1) SF_iso(mu2) x SF_trigger(event)   == histograms.pkl
Block `recut` (17 Sep 2026): the same frozen weights in the order the talk tells them, tag-and-probe first, plus the
measured reconstruction SF that the frozen fit inputs carry (FREEZE.md; flat per event, applied at merge time):
    raw -> sf_id (x SF_ID SF_ID) -> sf_muon (x SF_iso SF_iso SF_trigger SF_reco) -> sf_pileup (x w_PU) -> final (x L1 prefiring)
`final` == nominal x SF_reco == the frozen fit input (mumu_SR_prefit.yaml, results_v2.json /yields/SR); asserted.
The worker is a copy of z-mumu/scripts/v2_4_histograms.py:process_file with BRANCHES and _FixedWeighter imported from
that script and the same iterate(filter_name, step_size="150 MB") chunking (the momentum smearing is seeded per chunk,
zmumu/momentum.py:122-124), so that stage `nominal` reproduces output/v2/histograms.pkl bin by bin.
Data and fakes are taken from histograms.pkl / fakes.json (no data pass). Reads z-mumu/ read-only; part files go to
presentation/work/extract/ (git-ignored); the only output is presentation/data/zmumu_sr_stack.json.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (FIT_RESULTS, V2, WORK, ZMUMU, Checker, add_zmumu_path, dump_json, load_json, provenance,  # noqa: E402
                             rebin, standard_args)

add_zmumu_path()
from zmumu import batch, histograms as H, momentum, pileup, regions, samples, skim, tnp  # noqa: E402
from v2_4_histograms import BRANCHES, _FixedWeighter  # noqa: E402  (z-mumu/scripts, read-only import)

STACK_KEYS = ["DY_NLO", "TTbar", "ST_tW_top", "ST_tW_antitop", "WW", "WZ_3LNu", "WZ_2Q2L", "ZZ_4L", "ZZ_2L2Nu", "ZZ_2Q2L"]
SAMPLES_OUT = ["DYmumu", "DYtautau", "TTbar", "SingleTop", "WW", "WZ", "ZZ"]
STACK_ORDER = ["Fakes", "WW", "WZ", "ZZ", "SingleTop", "TTbar", "DYtautau", "DYmumu"]     # zmumu/plotting.py:9, bottom-up
STAGES = ["raw", "pileup", "prefiring", "nominal"]
STAGE_LABELS = {
    "raw": "generator weight x sigma L / sum(genWeight) only: no pileup weight, no L1 prefiring weight, no scale factors "
           "(the 'uncorrected simulation'); prompt-prompt MC only -- the non-prompt layer is always the data-driven Fakes template",
    "pileup": "raw x pileup weight w_PU(nTrueInt) (N_PV-matched profile, scale 1.035, rel_smear 0.09)",
    "prefiring": "pileup x L1PreFiringWeight_Nom (muon system x ECAL)",
    "nominal": "prefiring x SF_ID(mu1) SF_ID(mu2) SF_iso(mu1) SF_iso(mu2) x SF_trigger(event) from tag-and-probe "
               "== output/v2/histograms.pkl <sample>|SR|mass_fit|nominal (the fit input before mu_Z)",
}
MOMENTUM_NOTE = ("the Z-peak momentum calibration (kappa per |eta| bin, extra smearing) is applied to the MC muon pT in ALL four "
                 "stages: it is a re-selection, so the stages share one event set and their ratios are clean mean weight factors")
SUM_NAMES = ["raw", "raw_pu", "pileup", "pileup_pref", "pileup_pref_mu", "pileup_pref_ecal", "prefiring", "prefiring_sf",
             "prefiring_id", "prefiring_iso", "prefiring_trig", "raw_id", "raw_iso", "raw_trig", "raw_sf", "raw_sf_pu", "n_entries"]
# ---- the order the talk tells it in (re-cut of 17 Sep 2026): tag-and-probe first, then pileup, then L1 prefiring. The weights
# are the same frozen per-event factors; only the order of multiplication differs, so `final` is the frozen fit input.
RECUT_STAGES = ["raw", "sf_id", "sf_muon", "sf_pileup", "final"]
RECUT_FILLED = {"sf_id": "recut_id", "sf_muon": "recut_sf", "sf_pileup": "recut_sf_pu"}      # worker histograms (before the reco SF)
RECUT_LABELS = {
    "raw": "identical to stage `raw`",
    "sf_id": "raw x SF_ID(mu1) SF_ID(mu2): the tag-and-probe tight-ID map alone",
    "sf_muon": "raw x SF_ID SF_ID x SF_iso SF_iso x SF_trigger(event) x SF_reco(event): every tag-and-probe scale factor; the "
               "reconstruction SF (fit meta reco_sf.sf_per_event, flat) is applied at merge time, which is exact for a flat factor",
    "sf_pileup": "sf_muon x pileup weight w_PU(nTrueInt)",
    "final": "sf_pileup x L1PreFiringWeight_Nom == stage `nominal` x SF_reco == the frozen fit input (FREEZE.md, mumu_SR_prefit.yaml)",
}
EDGES = H.edges("SR", "mass_fit")          # 60 x 1 GeV, 60-120
PARTS_DEFAULT = WORK / "sr_stack"


# ----------------------------------------------------------------------------- worker
def fill_stages(ev, out, key, weighter, sf, calib, split, fit_sample):
    """Mirror of zmumu/histograms.py:fill_chunk restricted to SR / mass_fit / nominal pT, with the weight built up in stages."""
    keep = regions.fired(ev) & regions.event_clean(ev)
    ev = ev[keep]
    if len(ev) == 0:
        return
    n = len(ev)
    pieces = weighter.pieces(ev)                                   # zmumu/weights.py:53-64
    gen, pu, pref = pieces["gen"], pieces["pu"], pieces["pref"]
    pref_mu = np.asarray(ev.L1PreFiringWeight_Muon_Nom, dtype=float)
    pref_ecal = np.asarray(ev.L1PreFiringWeight_ECAL_Nom, dtype=float)
    flav = np.asarray(ev.gen_lhe_flavour) if split else None
    e = regions.with_muon_pt(ev, calib.scaled_pt(ev, "nominal"))   # calibration on (seeded on the filtered chunk, as fill_chunk)
    masks = regions.muon_masks(e)
    matched = regions.trigger_match(e)
    ones = np.ones(n, dtype=bool)
    d = regions.dimuon_regions(e, masks, matched, ones)["SR"]
    if d is None or len(d["idx"]) == 0:
        return
    idx = d["idx"]
    pp = regions.is_prompt(d["flav1"]) & regions.is_prompt(d["flav2"])       # histograms.py:142
    pt1, eta1, pt2, eta2 = d["pt1"], d["eta1"], d["pt2"], d["eta2"]
    sfw = sf.event_weights(pt1, eta1, pt2, eta2)["nominal"]                  # tnp.py:330-349
    id_pair = sf.muon_sf("id", pt1, eta1) * sf.muon_sf("id", pt2, eta2)
    iso_pair = sf.muon_sf("iso", pt1, eta1) * sf.muon_sf("iso", pt2, eta2)
    d1, d2 = sf.trig_eff("data", pt1, eta1), sf.trig_eff("data", pt2, eta2)
    m1, m2 = sf.trig_eff("mc", pt1, eta1), sf.trig_eff("mc", pt2, eta2)
    num, den = 1.0 - (1.0 - d1) * (1.0 - d2), 1.0 - (1.0 - m1) * (1.0 - m2)
    trig = np.where(den > 0, num / np.where(den > 0, den, 1.0), 1.0)
    if not np.allclose(id_pair * iso_pair * trig, sfw, rtol=1e-12, atol=0):
        raise RuntimeError("SF decomposition does not reproduce ScaleFactors.event_weights")
    w_raw = np.where(pp, gen[idx], 0.0)
    w_pu = w_raw * pu[idx]
    w_pref = w_pu * pref[idx]
    w_nom = w_pref * sfw
    stage_w = {"raw": w_raw, "pileup": w_pu, "prefiring": w_pref, "nominal": w_nom,
               "recut_id": w_raw * id_pair, "recut_sf": w_raw * sfw, "recut_sf_pu": w_raw * sfw * pu[idx]}
    if flav is None:
        here = np.full(len(idx), fit_sample, dtype=object)
    else:
        here = np.array([H.FLAVOUR_SAMPLE.get(int(f), "DYother") for f in flav[idx]], dtype=object)
    mass = d["mass"]
    for smp in np.unique(here):
        m = here == smp
        for stage, w in stage_w.items():
            H._fill(out, f"{smp}|{stage}", mass[m], w[m], EDGES)
        sums = {"raw": w_raw[m].sum(), "raw_pu": (w_raw * pu[idx])[m].sum(), "pileup": w_pu[m].sum(),
                "pileup_pref": (w_pu * pref[idx])[m].sum(), "pileup_pref_mu": (w_pu * pref_mu[idx])[m].sum(),
                "pileup_pref_ecal": (w_pu * pref_ecal[idx])[m].sum(), "prefiring": w_pref[m].sum(),
                "prefiring_sf": w_nom[m].sum(), "prefiring_id": (w_pref * id_pair)[m].sum(),
                "prefiring_iso": (w_pref * iso_pair)[m].sum(), "prefiring_trig": (w_pref * trig)[m].sum(),
                "raw_id": (w_raw * id_pair)[m].sum(), "raw_iso": (w_raw * iso_pair)[m].sum(), "raw_trig": (w_raw * trig)[m].sum(),
                "raw_sf": (w_raw * sfw)[m].sum(), "raw_sf_pu": (w_raw * sfw * pu[idx])[m].sum(),
                "n_entries": float(np.count_nonzero(pp[m]))}
        for name, v in sums.items():
            k = f"{smp}|sum|{name}"
            out[k] = out.get(k, 0.0) + float(v)


def process_file(source, key, aux_dir, check_fill_chunk=False):
    """Copy of v2_4_histograms.process_file (:47-69) with fill_stages instead of H.fill_chunk."""
    import uproot
    s = samples.SAMPLES[key]
    assert s["is_mc"]
    aux = Path(aux_dir)
    pu = pileup.PileupWeights.from_json(aux / "tnp" / "pileup_weights.json")
    gens = json.load(open(aux / "gensums.json")).get(key)          # never skim.load_gensums in a worker
    gens = {k: (np.asarray(v) if isinstance(v, list) else v) for k, v in gens.items()}
    weighter = _FixedWeighter(key, gens, pu)
    sf = tnp.ScaleFactors(aux / "tnp" / "tnp_result.json")
    calib = momentum.MomentumCalibration(aux / "momentum.json")
    fit_sample = s["fit_sample"]
    split = bool(s.get("split_lhe"))
    out = {"n_files": 1}
    out2 = {} if check_fill_chunk else None
    n_chunks = 0
    with uproot.open(source) as f:
        if "Events" not in f:
            return out
        for ev in f["Events"].iterate(filter_name=BRANCHES, step_size="150 MB", library="ak"):
            n_chunks += 1
            fill_stages(ev, out, key, weighter, sf, calib, split, fit_sample)
            if out2 is not None:
                H.fill_chunk(ev, out2, key, True, weighter, sf, calib, split_flavour=split, theory=False,
                             fit_sample=fit_sample, variations=False)
    out["n_chunks"] = n_chunks
    if out2 is not None:
        smps = sorted({k.split("|")[0] for k in out2 if "|SR|mass_fit|nominal" in k} | {k.split("|")[0] for k in out if k.endswith("|nominal")})
        worst, mode = 0.0, "exact"
        for smp in smps:
            a = np.asarray(out.get(f"{smp}|nominal", np.zeros(60)))
            b = np.asarray(out2.get(f"{smp}|SR|mass_fit|nominal", np.zeros(60)))
            worst = max(worst, float(np.max(np.abs(a - b))) if a.size else 0.0)
            if not np.array_equal(a, b):
                mode = "allclose" if np.allclose(a, b, rtol=1e-9, atol=0) else "FAIL"
        out["fill_chunk_check"] = [{"file": Path(source).name, "samples": smps, "max_abs_diff": worst, "mode": mode}]
        if mode == "FAIL":
            raise RuntimeError(f"{source}: fill_stages nominal != H.fill_chunk nominal (max |d| = {worst})")
    return out


# ----------------------------------------------------------------------------- driver
def run(args):
    parts = Path(args.parts_dir)
    keys = args.samples or STACK_KEYS
    for key in keys:
        if key not in STACK_KEYS:
            sys.exit(f"{key} is not one of the SR stack samples {STACK_KEYS}")
    t0 = time.time()
    if not args.merge_only:
        for key in keys:
            files = skim.skim_files(key)[: args.max_files]
            if not files:
                print(f"[stack] {key}: no skim files, skipped"); continue
            tasks = [(str(f), f"{key}__{f.stem}") for f in files]
            print(f"[stack] {key}: {len(tasks)} files", flush=True)
            extra = ["--sample", key, "--aux", str(V2)] + (["--check-fill-chunk"] if args.check_fill_chunk else [])
            batch.run_files(tasks, Path(__file__), parts / key, args.workers, extra_args=extra, stall_s=2400, retries=1,
                            log=lambda s: print(s, flush=True))
    total, n_files = {}, {}
    for key in keys:
        d = parts / key
        if d.exists() and any(d.glob("*.pkl")):
            part = batch.load_parts(d, {})
            n_files[key] = int(part.pop("n_files", 0))
            part.pop("n_chunks", None)
            batch.merge(total, part)
    runtime = time.time() - t0
    print(f"[stack] pass over {sum(n_files.values())} files of {list(n_files)} in {runtime:.0f} s")
    return total, n_files, runtime


def build_json(total, n_files, runtime, args):
    with open(V2 / "histograms.pkl", "rb") as fh:
        hall = pickle.load(fh)
    fk = load_json(V2 / "fakes.json")
    data = np.asarray(hall["Data|SR|mass_fit|nominal"], dtype=float)
    assert np.all(data == np.round(data))
    ft = np.array(fk["templates"]["nominal"], dtype=float)
    fv = np.array(fk["templates"]["nominal_var"], dtype=float)
    fakes = {"counts": ft.reshape(60, 2).sum(1).tolist(), "var": (np.sqrt(fv).reshape(60, 2).sum(1) ** 2).tolist(), "total": float(ft.sum()),
             "source": "fakes.json templates.nominal (120 x 0.5 GeV) summed in pairs; var = (sum sqrt(var))^2 (fully correlated template error)"}
    mc, extra = {}, {}
    present = sorted({k.split("|")[0] for k in total if "|" in k})
    for smp in present:
        block = {}
        for st in STAGES:
            block[st] = np.asarray(total.get(f"{smp}|{st}", np.zeros(60)), dtype=float).tolist()
            block[st + "_w2"] = np.asarray(total.get(f"{smp}|{st}|w2", np.zeros(60)), dtype=float).tolist()
        block["totals"] = {st: float(np.sum(block[st])) for st in STAGES}
        block["sums"] = {nm: float(total.get(f"{smp}|sum|{nm}", 0.0)) for nm in SUM_NAMES}
        block["n_entries"] = int(round(block["sums"]["n_entries"]))
        (mc if smp in SAMPLES_OUT else extra)[smp] = block
    for smp in SAMPLES_OUT:
        if smp not in mc:
            print(f"[stack] WARNING: {smp} has no entries in this pass")

    def means(sums):
        f = lambda a, b: (sums[a] / sums[b]) if sums[b] else float("nan")
        return {"pileup": f("raw_pu", "raw"), "prefiring": f("pileup_pref", "pileup"), "prefiring_muon": f("pileup_pref_mu", "pileup"),
                "prefiring_ecal": f("pileup_pref_ecal", "pileup"), "sf_event": f("prefiring_sf", "prefiring"),
                "sf_id_pair": f("prefiring_id", "prefiring"), "sf_iso_pair": f("prefiring_iso", "prefiring"),
                "sf_trigger_event": f("prefiring_trig", "prefiring")}

    mean_factors = {smp: means(mc[smp]["sums"]) for smp in mc}
    all_sums = {nm: sum(mc[s]["sums"][nm] for s in mc) for nm in SUM_NAMES}
    mean_factors["all_mc"] = means(all_sums)
    totals, ratio = {}, {}
    for st in STAGES:
        mc_tot = np.sum([mc[s][st] for s in mc], axis=0) if mc else np.zeros(60)
        pred = mc_tot + np.array(fakes["counts"])
        totals[st] = {"mc": float(mc_tot.sum()), "mc_plus_fakes": float(pred.sum()), "data_over_pred": float(data.sum() / pred.sum()),
                      "bkg_plus_fakes": float(pred.sum() - (np.sum(mc["DYmumu"][st]) if "DYmumu" in mc else 0.0))}
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio[st] = np.where(pred > 0, data / pred, np.nan).tolist()
    # ---- re-cut order, with the frozen reconstruction SF on the prompt simulation (never on the data-driven fakes)
    if "DYmumu" in mc and "DYmumu|recut_sf" not in total:
        sys.exit("[stack] the part files carry no recut histograms (stale parts cache): delete the parts directory and rerun")
    reco = float(load_json(FIT_RESULTS / "zmumu_fit_result.json")["meta"]["reco_sf"]["sf_per_event"])

    def recut_hist(smp, st):
        if st == "raw":
            return np.asarray(mc[smp]["raw"], dtype=float)
        if st == "final":
            return np.asarray(mc[smp]["nominal"], dtype=float) * reco
        h = np.asarray(total.get(f"{smp}|{RECUT_FILLED[st]}", np.zeros(60)), dtype=float)
        return h if st == "sf_id" else h * reco

    recut_mc = {smp: {st: recut_hist(smp, st).tolist() for st in RECUT_STAGES} for smp in mc}
    recut_totals, recut_ratio = {}, {}
    for st in RECUT_STAGES:
        mc_tot = np.sum([recut_mc[s][st] for s in mc], axis=0) if mc else np.zeros(60)
        pred = mc_tot + np.array(fakes["counts"])
        recut_totals[st] = {"mc": float(mc_tot.sum()), "mc_plus_fakes": float(pred.sum()), "data_over_pred": float(data.sum() / pred.sum()),
                            "bkg_plus_fakes": float(pred.sum() - (np.sum(recut_mc["DYmumu"][st]) if "DYmumu" in mc else 0.0))}
        with np.errstate(divide="ignore", invalid="ignore"):
            recut_ratio[st] = np.where(pred > 0, data / pred, np.nan).tolist()
    n = {st: recut_totals[st]["mc"] for st in RECUT_STAGES}
    g = lambda a, b: (all_sums[a] / all_sums[b]) if all_sums[b] else float("nan")
    recut = {"stages": RECUT_STAGES, "stage_labels": RECUT_LABELS, "reco_sf_per_event": reco, "mc": recut_mc, "totals": recut_totals,
             "data_over_pred": recut_ratio,
             "step_factors": {"muon_id": n["sf_id"] / n["raw"], "muon_iso_trigger_reco": n["sf_muon"] / n["sf_id"],
                              "muon_efficiency": n["sf_muon"] / n["raw"], "pileup": n["sf_pileup"] / n["sf_muon"],
                              "prefiring": n["final"] / n["sf_pileup"], "all": n["final"] / n["raw"],
                              "note": "ratios of the total simulated SR yield (all prompt samples) between consecutive recut stages: the mean "
                                      "event weight of each step, in the order shown; they multiply to `all` exactly"},
             "raw_weighted_means": {"sf_id_pair": g("raw_id", "raw"), "sf_iso_pair": g("raw_iso", "raw"), "sf_trigger_event": g("raw_trig", "raw"),
                                    "sf_reco_event": reco, "sf_event": g("raw_sf", "raw"),
                                    "note": "mean pair / event scale factors over all prompt simulation, weighted with the raw event weight"}}
    out = {
        "provenance": provenance("extract_zmumu_sr_stack.py",
                                 [skim.skim_dir()] + [V2 / p for p in ("gensums.json", "momentum.json", "tnp/tnp_result.json", "tnp/pileup_weights.json", "histograms.pkl", "fakes.json")]
                                 + [FIT_RESULTS / "zmumu_fit_result.json"],
                                 skim_dir=str(skim.skim_dir()), skim_files=n_files, mc_keys=list(n_files),
                                 fit_sample_map={k: (samples.SAMPLES[k]["fit_sample"] or "DYmumu/DYee/DYtautau by gen_lhe_flavour") for k in n_files},
                                 selection=("HLT_IsoMu24 || HLT_IsoTkMu24, MET filters, PV_npvsGood >= 1; exactly two tight muons (tightId, pfRelIso04 < 0.15, "
                                            "pT > 20, |eta| < 2.4, |dxy| < 0.2 cm, |dz| < 0.5 cm), leading pT > 26, >= 1 trigger-matched, opposite sign, "
                                            "60 < m(mumu, FSR-recovered) < 120 GeV; MC: both muons prompt (genPartFlav 1 or 15)"),
                                 variable="m(mu mu) FSR-recovered [GeV], 60 x 1 GeV bins (histograms.py VARIABLES['mass_fit'])",
                                 momentum_calibration=MOMENTUM_NOTE, chunking="uproot iterate(filter_name=v2_4_histograms.BRANCHES, step_size='150 MB')",
                                 data_source="histograms.pkl Data|SR|mass_fit|nominal (Run2016G+H, prescale-weighted; SR prescale is 1)",
                                 runtime_s=round(runtime, 1), workers=args.workers, parts_dir=str(args.parts_dir),
                                 anchors=["z-mumu/handoff.md:35-36", "z-mumu/output/v2/RESULTS_v2.md yields", "z-mumu/docs/11-mc-weights.md:31",
                                          "z-mumu/output/v2/results_v2.json /yields/SR"]),
        "edges": EDGES.tolist(), "samples": SAMPLES_OUT, "stack_order": STACK_ORDER, "stages": STAGES, "stage_labels": STAGE_LABELS,
        "momentum_calibration": MOMENTUM_NOTE,
        "data": {"counts": [int(v) for v in data], "total": int(data.sum())},
        "mc": mc, "mc_extra": extra, "fakes": fakes, "totals": totals, "data_over_pred": ratio, "mean_factors": mean_factors,
        "recut": recut,
    }
    return out


def verify(d, ck: Checker):
    from extract_zmumu_fit import read_plot_yaml
    with open(V2 / "histograms.pkl", "rb") as fh:
        hall = pickle.load(fh)
    res = load_json(V2 / "results_v2.json")
    pre = read_plot_yaml(ZMUMU / "fit" / "results" / "zmumu" / "Plots" / "mumu_SR_prefit.yaml")
    fitres = load_json(FIT_RESULTS / "zmumu_fit_result.json")
    rc = d["recut"]
    ck.check("data total 10 378 567", d["data"]["total"], 10378567, 0, "z-mumu/handoff.md:56")
    ck.check("data counts sum == total", sum(d["data"]["counts"]), d["data"]["total"], 0)
    ck.check("data == histograms.pkl Data|SR|mass_fit|nominal", d["data"]["counts"], hall["Data|SR|mass_fit|nominal"], 0, "output/v2/histograms.pkl")
    ck.check("data 5 GeV rebin == prefit YAML Data", rebin(d["data"]["counts"], 5), pre["data"], 0, "fit/results/zmumu/Plots/mumu_SR_prefit.yaml")
    ck.check_true("edges 60 x 1 GeV, 60-120", len(d["edges"]) == 61 and d["edges"][0] == 60 and d["edges"][-1] == 120)
    ck.check_true("all 7 samples present", all(s in d["mc"] for s in SAMPLES_OUT), detail=str([s for s in SAMPLES_OUT if s not in d["mc"]]))
    ck.check_true("no unexpected samples with entries (DYee/DYother/WJets)", not any(abs(v["totals"]["nominal"]) > 0 for v in d.get("mc_extra", {}).values()),
                  detail=str({k: v["totals"]["nominal"] for k, v in d.get("mc_extra", {}).items()}))
    modes = {}
    for smp in SAMPLES_OUT:
        if smp not in d["mc"]:
            continue
        got, want = np.array(d["mc"][smp]["nominal"]), np.asarray(hall[f"{smp}|SR|mass_fit|nominal"], dtype=float)
        want_w2 = np.asarray(hall[f"{smp}|SR|mass_fit|nominal|w2"], dtype=float)
        if np.allclose(got, want, rtol=1e-6, atol=0):
            modes[smp] = "rtol 1e-6" + (" (bit-identical)" if np.array_equal(got, want) else "")
            ck.check(f"nominal {smp} == histograms.pkl (rtol 1e-6)", got, want, 1e-6, "output/v2/histograms.pkl", rel=True)
        else:
            sig = np.sqrt(np.maximum(want_w2, 1e-300))
            ok = bool(np.all(np.abs(got - want) <= sig)) and abs(got.sum() - want.sum()) <= 1e-4 * abs(want.sum())
            modes[smp] = "FALLBACK per-bin <= 1 sigma and total rel 1e-4" if ok else "FAIL"
            ck.check_true(f"nominal {smp} == histograms.pkl (fallback: per-bin <= 1 sigma, total rel 1e-4)", ok, "output/v2/histograms.pkl",
                          f"max |d|/sigma = {np.max(np.abs(got - want) / sig):.3g}, total rel = {abs(got.sum() - want.sum()) / abs(want.sum()):.3g}")
        ck.check(f"nominal w2 {smp} == histograms.pkl w2 (rtol 1e-6)", d["mc"][smp]["nominal_w2"], want_w2, 1e-6, rel=True)
    d.setdefault("checks_meta", {})["histograms_pkl_comparison_mode"] = modes
    print(f"[stack] histograms.pkl comparison modes: {modes}")
    for smp, want in (("DYmumu", 10373541.7), ("DYtautau", 11022.9), ("TTbar", 31607.2), ("SingleTop", 2935.7), ("WW", 3847.3), ("WZ", 9194.9), ("ZZ", 6308.7)):
        ck.check(f"nominal total {smp} {want:,.1f} (before the reconstruction SF)", d["mc"][smp]["totals"]["nominal"], want, 0.05, "output/v2/histograms.pkl")
        ck.check(f"recut final total {smp} == results_v2.json /yields/SR (frozen, incl. reconstruction SF)", sum(rc["mc"][smp]["final"]),
                 res["yields"]["SR"][smp], 1e-9, "results_v2.json /yields/SR", rel=True)
        ck.check(f"recut final {smp} 5 GeV rebin == frozen prefit YAML", rebin(rc["mc"][smp]["final"], 5), pre["samples"][smp], 1e-6,
                 "fit/results/zmumu/Plots/mumu_SR_prefit.yaml", rel=True)
    bkg = sum(d["mc"][s]["totals"]["nominal"] for s in SAMPLES_OUT if s != "DYmumu") + d["fakes"]["total"]
    ck.check("sum bkg + fakes 68 786.5 (nominal, before the reconstruction SF)", bkg, 68786.5, 0.05, "fit/results/zmumu_v2_15sep_fit_result.json counting n_bkg")
    ck.check("bkg_plus_fakes (nominal) == 15 Sep counting n_bkg 68786.525", d["totals"]["nominal"]["bkg_plus_fakes"], 68786.52541969114, 1e-6, "15 Sep fit result meta", rel=True)
    ck.check("recut final: bkg + fakes 68 799.0", rc["totals"]["final"]["bkg_plus_fakes"], 68799.0, 0.05, "z-mumu/handoff.md:57 / FREEZE.md")
    ck.check("recut final: bkg + fakes == frozen counting n_bkg", rc["totals"]["final"]["bkg_plus_fakes"], fitres["meta"]["counting"]["n_bkg"], 1e-6,
             "fit/results/zmumu_fit_result.json meta.counting", rel=True)
    ck.check("fakes total 3869.955", d["fakes"]["total"], 3869.955, 5e-4, "fakes.json")
    ck.check("fakes counts sum == total", sum(d["fakes"]["counts"]), d["fakes"]["total"], 1e-9, rel=True)
    ck.check_true("fakes no negative bins", min(d["fakes"]["counts"]) >= 0)
    ck.check("DYmumu nominal 5 GeV rebin x reco SF == prefit YAML DYmumu", np.asarray(rebin(d["mc"]["DYmumu"]["nominal"], 5)) * rc["reco_sf_per_event"],
             pre["samples"]["DYmumu"], 1e-6, "mumu_SR_prefit.yaml", rel=True)
    for smp in SAMPLES_OUT:
        for st, nm in (("raw", "raw"), ("pileup", "pileup"), ("prefiring", "prefiring"), ("nominal", "prefiring_sf")):
            ck.check(f"{smp} {st} histogram sum == weight sum", d["mc"][smp]["totals"][st], d["mc"][smp]["sums"][nm], 1e-9, rel=True)
    mf = d["mean_factors"]
    ck.check("<L1 prefiring> DYmumu SR 0.9803", mf["DYmumu"]["prefiring"], 0.9803, 5e-4, "z-mumu/docs/11-mc-weights.md:31")
    ck.check("<L1 prefiring muon> DYmumu 0.9818", mf["DYmumu"]["prefiring_muon"], 0.9818, 5e-4, "docs/11-mc-weights.md:31")
    ck.check("<L1 prefiring ECAL> DYmumu 0.9985", mf["DYmumu"]["prefiring_ecal"], 0.9985, 5e-4, "docs/11-mc-weights.md:31")
    ck.check("stage ratio nominal/prefiring DYmumu == <sf_event>", d["mc"]["DYmumu"]["totals"]["nominal"] / d["mc"]["DYmumu"]["totals"]["prefiring"], mf["DYmumu"]["sf_event"], 1e-9, rel=True)
    ck.check("stage ratio pileup/raw DYmumu == <pileup>", d["mc"]["DYmumu"]["totals"]["pileup"] / d["mc"]["DYmumu"]["totals"]["raw"], mf["DYmumu"]["pileup"], 1e-9, rel=True)
    r_raw, r_nom = d["totals"]["raw"]["data_over_pred"], d["totals"]["nominal"]["data_over_pred"]
    ck.check_true("raw stage: data / (raw MC + fakes) in [0.90, 0.98] (~0.94 expected)", 0.90 < r_raw < 0.98, "plan: raw simulation ~6% above data", f"{r_raw:.4f}")
    ck.check_true("nominal stage: data / (MC + fakes) rounds to 0.994", round(r_nom, 3) == 0.994, detail=f"{r_nom:.6f}")
    r_fin = rc["totals"]["final"]["data_over_pred"]
    ck.check("recut final: data / (MC + fakes) == frozen counting N_obs / (N_DYmumu + N_bkg)", r_fin,
             10378567 / (res["yields"]["SR"]["DYmumu"] + fitres["meta"]["counting"]["n_bkg"]), 1e-6, "fit result meta.counting", rel=True)
    ck.check("recut final 5 GeV data/pred == results_v2 lineshape (frozen)",
             rebin(d["data"]["counts"], 5) / (sum(np.asarray(rebin(rc["mc"][s]["final"], 5)) for s in SAMPLES_OUT) + rebin(d["fakes"]["counts"], 5)),
             res["lineshape"]["data_over_pred_prefit_5gev"], 1e-9, "results_v2.json /lineshape", rel=True)
    # ---- the re-cut order itself
    ck.check("reco SF per event == fit meta", rc["reco_sf_per_event"], fitres["meta"]["reco_sf"]["sf_per_event"], 0, "fit result meta.reco_sf")
    for smp in SAMPLES_OUT:
        ck.check(f"recut raw {smp} == stage raw", rc["mc"][smp]["raw"], d["mc"][smp]["raw"], 0)
        ck.check(f"recut sf_id {smp} histogram sum == weight sum", sum(rc["mc"][smp]["sf_id"]), d["mc"][smp]["sums"]["raw_id"], 1e-9, rel=True)
        ck.check(f"recut sf_muon {smp} histogram sum == weight sum x reco SF", sum(rc["mc"][smp]["sf_muon"]),
                 d["mc"][smp]["sums"]["raw_sf"] * rc["reco_sf_per_event"], 1e-9, rel=True)
        ck.check(f"recut sf_pileup {smp} histogram sum == weight sum x reco SF", sum(rc["mc"][smp]["sf_pileup"]),
                 d["mc"][smp]["sums"]["raw_sf_pu"] * rc["reco_sf_per_event"], 1e-9, rel=True)
    sfac = rc["step_factors"]
    ck.check("recut step factors multiply to the overall factor", sfac["muon_efficiency"] * sfac["pileup"] * sfac["prefiring"], sfac["all"], 1e-12, rel=True)
    ck.check("recut muon_efficiency == muon_id x (iso, trigger, reco)", sfac["muon_id"] * sfac["muon_iso_trigger_reco"], sfac["muon_efficiency"], 1e-12, rel=True)
    ladder = [round(rc["totals"][st]["data_over_pred"], 3) for st in RECUT_STAGES]
    ck.check("recut ladder data/pred: raw, sf_id, sf_muon, sf_pileup, final", ladder, [0.944, 0.971, 0.968, 0.974, 0.994], 1e-9,
             "this extractor (the numbers the clips 4-07, 4-14 ... 4-17 print)")
    ck.check_true("recut: every stage has a ratio inside the fixed panel range 0.88-1.12 in all 60 bins",
                  all(0.88 < v < 1.12 for st in RECUT_STAGES for v in rc["data_over_pred"][st]),
                  detail=str({st: (round(min(rc["data_over_pred"][st]), 4), round(max(rc["data_over_pred"][st]), 4)) for st in RECUT_STAGES}))
    ck.check("lumi_pb", d["provenance"]["lumi_pb"], 16393.381, 1e-6)
    return ck


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--samples", nargs="*", default=None, help=f"subset of {STACK_KEYS}")
    ap.add_argument("--check-fill-chunk", action="store_true", help="also run H.fill_chunk on every chunk and assert equality (smoke test)")
    ap.add_argument("--parts-dir", type=Path, default=PARTS_DEFAULT, help="per-file parts (git-ignored), resumable")
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument("--no-verify", action="store_true", help="write the JSON without the anchor assertions (smoke runs on a subset)")
    standard_args(ap, "zmumu_sr_stack.json")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)          # part file (zmumu/batch.py), never the JSON
    ap.add_argument("--sample", help=argparse.SUPPRESS)
    ap.add_argument("--aux", default=str(V2), help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:
        batch.write_part(process_file(args.one_file, args.sample, args.aux, args.check_fill_chunk), args.out)
        return
    if args.check_only:
        d = load_json(args.json)
        print(f"[zmumu_sr_stack] --check-only on {args.json}")
        verify(d, Checker("zmumu_sr_stack")).finish()
        return
    total, n_files, runtime = run(args)
    fcc = total.pop("fill_chunk_check", None)
    d = build_json(total, n_files, runtime, args)
    if fcc:
        d["provenance"]["fill_chunk_check"] = fcc
        print(f"[stack] --check-fill-chunk: {len(fcc)} files, modes {sorted({r['mode'] for r in fcc})}, max |d| {max(r['max_abs_diff'] for r in fcc):.3g}")
    subset = bool(args.samples) or args.max_files is not None
    if args.no_verify or subset:
        print(f"[stack] subset run ({list(n_files)}, max_files={args.max_files}): anchor assertions skipped, JSON -> {args.json}")
        for smp in d["mc"]:
            print(f"        {smp}: " + ", ".join(f"{st} {d['mc'][smp]['totals'][st]:,.2f}" for st in STAGES))
        print("        mean factors:", {s: {k: round(v, 5) for k, v in m.items()} for s, m in d["mean_factors"].items()})
        dump_json(args.json, d)
        return
    rows = verify(d, Checker("zmumu_sr_stack")).finish()
    d["checks"] = rows
    dump_json(args.json, d)


if __name__ == "__main__":
    main()
