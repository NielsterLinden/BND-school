#!/usr/bin/env python
"""Freeze real CMS 2016 Open Data events for the Z -> mumu event displays.

    source fitting/setup.sh && python presentation/data/extract_zmumu_events.py [--files 6] [--seed 20260915] [--json PATH] [--check-only]

Reads the first N data skim files of $BND_SKIM_DIR/data_2016G (SingleMuon Run2016G, NanoAODv9, record 30530) and
z-mumu/output/v2/fakes.json. Signal-region events use the analysis selection (zmumu/regions.py dimuon_regions, FSR-recovered
mass); the same-sign tag+probe pairs re-implement zmumu/fakes.py:85-101 at pair level (bare mass). Events are chosen for
DISPLAY with a seeded generator and documented rules -- this is not a physics result. Writes presentation/data/zmumu_events.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import V2, Checker, add_zmumu_path, finalize, load_json, provenance, standard_args  # noqa: E402

add_zmumu_path()
import awkward as ak  # noqa: E402
from zmumu import config, fakes, objects, regions, skim  # noqa: E402

BRANCHES = ["run", "luminosityBlock", "event", "PV_npvsGood", "era", "skim_cat", "skim_prescale", "HLT_IsoMu24", "HLT_IsoTkMu24",
            "Flag_*", "nMuon", "Muon_*", "TrigObj_*", "FsrPhoton_*", "Jet_*", "MET_pt", "MET_phi"]
RUN_G = (278820, 280385)                 # Run2016G (skim.py:32: 280385 is the last run of G)
N_MASSES = 500
N_SR_KEEP = 20
POOL_MAX = 120
JET_PT_MIN, JET_MAX = 20.0, 3
SR_RULES = {
    "textbook": "89 < m < 93 GeV, both |eta| < 0.9, back-to-back (|dphi| > 3.0), no FSR photon attached, no jet with pT > 30 GeV",
    "fsr": "FSR-recovered mass exceeds the bare mass by more than 1 GeV (>= 1 attached FSR photon)",
    "forward": "one muon with |eta| > 2.1",
    "boosted": "pT(mumu) > 60 GeV",
}
SS_RULE = ("9 anti-isolated probes (tightId, 0.20 < pfRelIso04 < 1.0, probe inside a jet preferred) + 1 isolated probe (tight): the display "
           "split reflects FF ~ 0.10 in the barrel cells (pT < 50, |eta| < 0.9: 0.093 / 0.103 / 0.087 / 0.118); unique events, seeded order "
           "with the isolated probe at a random position >= 1; the pool's raw n_tight/(n_tight+n_anti) is stored alongside")


def r5(x):
    return float(round(float(x), 5))


def r6(x):
    return float(round(float(x), 6))


def _muon(sub, i, matched, k):
    """Record of muon index i[k] of event k of `sub` (all values verbatim from the skim, rounded)."""
    def g(name):
        return objects.take(sub[name], i)[k]
    return {"idx": int(i[k]), "pt": r5(g("Muon_pt")), "eta": r6(g("Muon_eta")), "phi": r6(g("Muon_phi")), "mass": r6(g("Muon_mass")),
            "charge": int(g("Muon_charge")), "iso": r6(g("Muon_pfRelIso04_all")), "tightId": bool(g("Muon_tightId")),
            "dxy": r6(g("Muon_dxy")), "dz": r6(g("Muon_dz")), "matched": bool(objects.take(matched, i)[k]),
            "isGlobal": bool(g("Muon_isGlobal")), "isTracker": bool(g("Muon_isTracker")), "jetIdx": int(g("Muon_jetIdx"))}


def _jets(sub, k, mu_list):
    """Jets with pT > 20 of event k, dR to the listed muons and whether they contain them (Jet_muonIdx1/2)."""
    pt, eta, phi, mass = (np.asarray(sub["Jet_pt"][k]), np.asarray(sub["Jet_eta"][k]), np.asarray(sub["Jet_phi"][k]), np.asarray(sub["Jet_mass"][k]))
    jid, m1, m2 = np.asarray(sub["Jet_jetId"][k]), np.asarray(sub["Jet_muonIdx1"][k]), np.asarray(sub["Jet_muonIdx2"][k])
    order = np.argsort(-pt)
    out = []
    for j in order:
        if pt[j] < JET_PT_MIN or len(out) >= JET_MAX:
            continue
        rec = {"idx": int(j), "pt": r5(pt[j]), "eta": r6(eta[j]), "phi": r6(phi[j]), "mass": r5(mass[j]), "jetId": int(jid[j]),
               "muonIdx1": int(m1[j]), "muonIdx2": int(m2[j])}
        for name, mu in mu_list.items():
            rec[f"dr_{name}"] = r5(objects.delta_r(eta[j], phi[j], mu["eta"], mu["phi"]))
            rec[f"contains_{name}"] = bool(m1[j] == mu["idx"] or m2[j] == mu["idx"])
        out.append(rec)
    return out


def _event_head(sub, k):
    return {"run": int(sub["run"][k]), "lumi": int(sub["luminosityBlock"][k]), "event": int(sub["event"][k]), "era": "G" if int(sub["era"][k]) == 0 else "H",
            "npv": int(sub["PV_npvsGood"][k]), "skim_cat": int(sub["skim_cat"][k]), "skim_prescale": int(sub["skim_prescale"][k]),
            "hlt_isomu24": bool(sub["HLT_IsoMu24"][k]), "hlt_isotkmu24": bool(sub["HLT_IsoTkMu24"][k]),
            "n_muons": int(sub["nMuon"][k]), "met_pt": r5(sub["MET_pt"][k]), "met_phi": r6(sub["MET_phi"][k])}


class FileData:
    """One skim file: SR pair arrays and same-sign tag+probe pair arrays (all after fired & clean)."""

    def __init__(self, path, file_no, ffmap, fferr):
        import uproot
        self.path, self.file_no = Path(path), file_no
        with uproot.open(path) as f:
            ev = f["Events"].arrays(filter_name=BRANCHES, library="ak")
        keep = regions.fired(ev) & regions.event_clean(ev)
        e = ev[keep]
        self.e = e
        n = len(e)
        masks = fakes.muon_classes(e)                         # tight / anti / loose / anti_alt (regions.muon_masks + alt)
        matched = regions.trigger_match(e)
        self.matched = matched
        self.n_tight = ak.to_numpy(ak.sum(masks["tight"], axis=1))
        # ---- signal region (regions.dimuon_regions, FSR-recovered mass) -------------------------------------------
        sr = regions.dimuon_regions(e, masks, matched, np.ones(n, dtype=bool))["SR"]
        self.sr = sr
        if sr is not None and len(sr["idx"]):
            sub = e[sr["idx"]]
            i1, i2 = objects.leading_two(sub, masks["tight"][sr["idx"]])
            assert np.allclose(ak.to_numpy(objects.take(sub.Muon_pt, i1)), sr["pt1"]) and np.allclose(ak.to_numpy(objects.take(sub.Muon_pt, i2)), sr["pt2"])
            self.sr_sub, self.sr_i1, self.sr_i2 = sub, i1, i2
            fsr_att = objects.fsr_photon_mask(sub) & ((sub.FsrPhoton_muonIdx == i1) | (sub.FsrPhoton_muonIdx == i2))
            self.sr_nfsr = ak.to_numpy(ak.sum(fsr_att, axis=1))
            self.sr_fsr_att = fsr_att
            self.sr_njet30 = ak.to_numpy(ak.sum((sub.Jet_pt > 30) & (abs(sub.Jet_eta) < 2.4) & (sub.Jet_jetId >= 2), axis=1))
            self.sr_dphi = objects.delta_phi(sr["phi1"], sr["phi2"])
        else:
            self.sr_sub = None
        # ---- same-sign tag + probe pairs (zmumu/fakes.py:85-101 at pair level, bare mass) --------------------------
        tight = masks["tight"]
        is_tag = tight & (e.Muon_pt > config.MU_PT_LEAD) & matched
        idx = ak.local_index(e.Muon_pt, axis=1)
        p = ak.combinations(idx, 2, fields=["a", "b"])
        tag = ak.concatenate([p.a, p.b], axis=1)
        probe = ak.concatenate([p.b, p.a], axis=1)
        px1, py1, pz1, e1 = objects.p4(e.Muon_pt[tag], e.Muon_eta[tag], e.Muon_phi[tag], e.Muon_mass[tag])
        px2, py2, pz2, e2 = objects.p4(e.Muon_pt[probe], e.Muon_eta[probe], e.Muon_phi[probe], e.Muon_mass[probe])
        m = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
        same = (e.Muon_charge[tag] * e.Muon_charge[probe]) > 0
        base = is_tag[tag] & same & (m > fakes.MIN_PAIR_MASS)
        n_pairs = ak.to_numpy(ak.num(tag, axis=1))
        flat = lambda a: ak.to_numpy(ak.flatten(a))
        self.tp = {"ev": np.repeat(np.arange(n), n_pairs), "tag": flat(tag), "probe": flat(probe), "mass": flat(m),
                   "base": flat(base), "tight": flat(base & masks["tight"][probe]), "anti": flat(base & masks["anti"][probe]),
                   "ppt": flat(e.Muon_pt[probe]), "peta": flat(e.Muon_eta[probe]),
                   "pjet": flat(e.Muon_jetIdx[probe])}
        t = self.tp
        ipt = np.clip(np.digitize(np.clip(t["ppt"], fakes.PT_EDGES[0], fakes.PT_EDGES[-1] - 1e-3), fakes.PT_EDGES) - 1, 0, len(fakes.PT_EDGES) - 2)
        ieta = np.clip(np.digitize(np.abs(t["peta"]), fakes.ETA_EDGES) - 1, 0, len(fakes.ETA_EDGES) - 2)
        t["ipt"], t["ieta"], t["ff"], t["ff_err"] = ipt, ieta, ffmap[ipt, ieta], fferr[ipt, ieta]
        t["barrel"] = (np.abs(t["peta"]) < 0.9) & (t["ppt"] < 50.0)
        t["in_jet"] = t["pjet"] >= 0

    # ---- records -------------------------------------------------------------------------------------------------
    def sr_record(self, k, rule=None):
        sub, i1, i2, sr = self.sr_sub, self.sr_i1, self.sr_i2, self.sr
        mu1, mu2 = _muon(sub, i1, self.matched[sr["idx"]], k), _muon(sub, i2, self.matched[sr["idx"]], k)
        att = self.sr_fsr_att[k]
        fsr = [{"pt": r5(pt), "eta": r6(eta), "phi": r6(phi), "muonIdx": int(mi), "relIso03": r6(ri), "dROverEt2": r6(dr)}
               for pt, eta, phi, mi, ri, dr in zip(sub.FsrPhoton_pt[k][att], sub.FsrPhoton_eta[k][att], sub.FsrPhoton_phi[k][att],
                                                    sub.FsrPhoton_muonIdx[k][att], sub.FsrPhoton_relIso03[k][att], sub.FsrPhoton_dROverEt2[k][att])]
        rec = _event_head(sub, k)
        rec.update({"rule": rule, "n_tight": int(self.n_tight[sr["idx"][k]]), "muons": [mu1, mu2], "fsr_photons": fsr,
                    "n_fsr_all": int(ak.num(sub.FsrPhoton_pt, axis=1)[k]),
                    "mass": r5(sr["mass"][k]), "mass_bare": r5(sr["mass_bare"][k]), "zpt": r5(sr["zpt"][k]), "zy": r5(sr["zy"][k]),
                    "dphi": r6(self.sr_dphi[k]), "deta": r6(sr["eta1"][k] - sr["eta2"][k]),
                    "dr": r6(objects.delta_r(sr["eta1"][k], sr["phi1"][k], sr["eta2"][k], sr["phi2"][k])),
                    "n_jets30": int(self.sr_njet30[k]), "jets": _jets(sub, k, {"mu1": mu1, "mu2": mu2}),
                    "skim_file": self.path.name, "file_no": self.file_no})
        return rec

    def tp_record(self, j):
        t = self.tp
        k = int(t["ev"][j])
        e = self.e
        it, ip = np.array([int(t["tag"][j])]), np.array([int(t["probe"][j])])
        sub = e[k:k + 1]
        tag = _muon(sub, it, self.matched[k:k + 1], 0)
        probe = _muon(sub, ip, self.matched[k:k + 1], 0)
        cls = "tight" if t["tight"][j] else "anti"
        probe.update({"cls": cls, "in_jet": bool(t["in_jet"][j])})
        rec = _event_head(sub, 0)
        rec.update({"tag": tag, "probe": probe, "mass_bare": r5(t["mass"][j]), "passes": cls == "tight", "n_tight": int(self.n_tight[k]),
                    "ff_cell": {"ipt": int(t["ipt"][j]), "ieta": int(t["ieta"][j]),
                                "pt_range": [float(fakes.PT_EDGES[t["ipt"][j]]), float(fakes.PT_EDGES[t["ipt"][j] + 1])],
                                "eta_range": [float(fakes.ETA_EDGES[t["ieta"][j]]), float(fakes.ETA_EDGES[t["ieta"][j] + 1])],
                                "ff": float(t["ff"][j]), "ff_err": float(t["ff_err"][j])},
                    "dphi": r6(objects.delta_phi(tag["phi"], probe["phi"])), "jets": _jets(sub, 0, {"tag": tag, "probe": probe}),
                    "skim_file": self.path.name, "file_no": self.file_no})
        return rec


def build(args):
    fk = load_json(V2 / "fakes.json")
    ffmap, fferr = np.array(fk["maps"]["nominal"]), np.array(fk["maps"]["nominal_err"])
    files = skim.skim_files("data_2016G")[: args.files]
    if len(files) < args.files:
        sys.exit(f"only {len(files)} data_2016G skim files found")
    rng = np.random.default_rng(args.seed)
    print(f"[events] {len(files)} files of data_2016G, seed {args.seed}")
    fds = [FileData(f, i, ffmap, fferr) for i, f in enumerate(files)]
    # ---- SR: masses, rule candidates ---------------------------------------------------------------------------------
    masses, cands, n_sr = [], {r: [] for r in SR_RULES}, 0
    for fd in fds:
        if fd.sr_sub is None:
            continue
        sr = fd.sr
        n_sr += len(sr["idx"])
        masses.extend(sr["mass"].tolist())
        aeta1, aeta2 = np.abs(sr["eta1"]), np.abs(sr["eta2"])
        rules = {
            "textbook": (sr["mass"] > 89) & (sr["mass"] < 93) & (aeta1 < 0.9) & (aeta2 < 0.9) & (np.abs(fd.sr_dphi) > 3.0) & (fd.sr_nfsr == 0) & (fd.sr_njet30 == 0),
            "fsr": (sr["mass"] - sr["mass_bare"] > 1.0) & (fd.sr_nfsr >= 1),
            "forward": np.maximum(aeta1, aeta2) > 2.1,
            "boosted": sr["zpt"] > 60.0,
        }
        for r, mask in rules.items():
            cands[r].extend((fd.file_no, int(k)) for k in np.nonzero(mask)[0])
    print(f"[events] SR events in the pass: {n_sr:,}; rule candidates: " + ", ".join(f"{r} {len(v)}" for r, v in cands.items()))
    chosen, used = [], set()
    for r in SR_RULES:
        pool = [c for c in cands[r] if c not in used]
        if not pool:
            sys.exit(f"no candidate for rule {r}")
        fno, k = pool[int(rng.integers(len(pool)))]
        used.add((fno, k))
        chosen.append(fds[fno].sr_record(k, rule=r))
    extra = []
    all_sr = [(fd.file_no, k) for fd in fds if fd.sr_sub is not None for k in range(len(fd.sr["idx"]))]
    for j in rng.permutation(len(all_sr)):
        c = all_sr[j]
        if c in used:
            continue
        used.add(c)
        extra.append(fds[c[0]].sr_record(c[1], rule="random"))
        if len(extra) >= N_SR_KEEP - len(chosen):
            break
    # ---- SS tag + probe pool (barrel probes) --------------------------------------------------------------------------
    tight_all, anti_all = [], []
    for fd in fds:
        t = fd.tp
        tight_all.extend((fd.file_no, int(j)) for j in np.nonzero(t["tight"] & t["barrel"])[0])
        anti_all.extend((fd.file_no, int(j)) for j in np.nonzero(t["anti"] & t["barrel"])[0])
    n_t, n_a = len(tight_all), len(anti_all)
    print(f"[events] barrel same-sign tag+probe pairs: {n_t} tight, {n_a} anti-isolated (raw fraction {n_t / max(n_t + n_a, 1):.4f})")
    if n_t < 3 or n_a < 27:
        sys.exit(f"barrel same-sign pool too small ({n_t} tight / {n_a} anti): raise --files")
    keep_t = [tight_all[j] for j in rng.permutation(n_t)[: min(n_t, 40)]]
    keep_a = [anti_all[j] for j in rng.permutation(n_a)[: max(0, POOL_MAX - len(keep_t))]]
    pool = [fds[f].tp_record(j) for f, j in keep_t] + [fds[f].tp_record(j) for f, j in keep_a]
    # chosen 10: 9 anti (in-jet preferred) + 1 tight, unique events, seeded order, tight at a random position >= 1
    def uniq(recs, n, prefer=None):
        out, seen = [], set()
        order = rng.permutation(len(recs))
        for pref in ([True, False] if prefer else [None]):
            for j in order:
                r = recs[j]
                if pref is not None and bool(r["probe"]["in_jet"]) != pref:
                    continue
                key = (r["run"], r["lumi"], r["event"])
                if key in seen:
                    continue
                seen.add(key); out.append(r)
                if len(out) >= n:
                    return out
        return out
    anti_recs = [r for r in pool if not r["passes"]]
    tight_recs = [r for r in pool if r["passes"]]
    nine = uniq(anti_recs, 9, prefer=True)
    one = uniq([r for r in tight_recs if (r["run"], r["lumi"], r["event"]) not in {(x["run"], x["lumi"], x["event"]) for x in nine}], 1)
    if len(nine) < 9 or len(one) < 1:
        sys.exit("not enough unique same-sign pairs for the ten")
    pos = int(rng.integers(1, 10))
    ten = nine[:pos] + one + nine[pos:]
    for j, r in enumerate(ten):
        r["display_index"] = j
    sidecars = {f.name: json.load(open(str(f) + ".json")) for f in files}
    out = {
        "provenance": provenance("extract_zmumu_events.py", files + [V2 / "fakes.json"],
                                 statement=("real CMS 2016 Open Data events (SingleMuon Run2016G, NanoAODv9, record 30530), run/lumi/event listed; "
                                            f"selected for display with seed {args.seed} and the rules stated below; not a physics result"),
                                 seed=args.seed, n_files=len(files), skim_dir=str(skim.skim_dir()), skim_sidecars=sidecars,
                                 selection_sr=("HLT_IsoMu24 || HLT_IsoTkMu24, MET filters, PV_npvsGood >= 1; exactly two tight muons, leading pT > 26, "
                                               ">= 1 trigger-matched, opposite sign, 60 < m(FSR-recovered) < 120 GeV (zmumu/regions.py:dimuon_regions)"),
                                 selection_ss_tp=("same-sign tag (tight, pT > 26, trigger-matched) + probe (tight or anti-isolated: tightId, 0.20 < iso < 1.0), "
                                                  "bare pair mass > 12 GeV, both orientations of every muon pair (zmumu/fakes.py:85-101); "
                                                  "pool restricted to barrel probes |eta| < 0.9, pT < 50 GeV"),
                                 mass_definitions={"sr.mass": "FSR-recovered (attached FsrPhoton candidates added)", "sr.mass_bare": "two muons only",
                                                   "ss_tp.mass_bare": "two muons only (fake-factor definition)"},
                                 units="pT, mass, MET in GeV; dxy, dz in cm; angles in rad; floats rounded to 1e-5 (pT/mass) and 1e-6 (eta/phi/iso)"),
        "sr": {"rules": SR_RULES, "chosen": chosen, "events": chosen + extra, "n_events_in_pass": n_sr,
               "n_candidates": {r: len(v) for r, v in cands.items()},
               "masses_first500": [r5(m) for m in masses[:N_MASSES]], "masses_note": "the first 500 SR masses in file order (for the rain)"},
        "ss_tp": {"chosen": ten, "chosen_rule": SS_RULE, "pool": pool, "pool_n_tight": n_t, "pool_n_anti": n_a,
                  "pool_fraction_tight": n_t / (n_t + n_a), "barrel_ff_cells": [float(ffmap[i, 0]) for i in range(4)],
                  "pool_note": f"pool = a seeded sample of barrel tight probes (<= 40) + barrel anti-isolated probes (<= {POOL_MAX} records in total)"},
    }
    return out


def _mass_of(parts):
    px = py = pz = en = 0.0
    for pt, eta, phi, m in parts:
        a, b, c, d = objects.p4(pt, eta, phi, m)
        px += a; py += b; pz += c; en += d
    return float(objects.invariant_mass(px, py, pz, en))


def _tight_ok(mu):
    return (mu["tightId"] and mu["iso"] < regions.ISO_TIGHT and mu["pt"] > regions.MU_PT_MIN and abs(mu["eta"]) < regions.MU_ETA
            and abs(mu["dxy"]) < config.MU_DXY_MAX and abs(mu["dz"]) < config.MU_DZ_MAX)


def _anti_ok(mu):
    return (mu["tightId"] and regions.ANTI_ISO_LO < mu["iso"] < regions.ANTI_ISO_HI and mu["pt"] > regions.MU_PT_MIN
            and abs(mu["eta"]) < regions.MU_ETA and abs(mu["dxy"]) < config.MU_DXY_MAX and abs(mu["dz"]) < config.MU_DZ_MAX)


def verify(d, ck: Checker):
    fk = load_json(V2 / "fakes.json")
    ffmap = np.array(fk["maps"]["nominal"])
    sr_events, ten, pool = d["sr"]["events"], d["ss_tp"]["chosen"], d["ss_tp"]["pool"]
    allrecs = sr_events + ten + pool
    ck.check_true("run in Run2016G [278820, 280385]", all(RUN_G[0] <= r["run"] <= RUN_G[1] for r in allrecs))
    ck.check_true("skim_prescale == 1", all(r["skim_prescale"] == 1 for r in allrecs))
    ck.check_true("skim category A (bit 1)", all(r["skim_cat"] & 1 for r in allrecs))
    ck.check_true("trigger fired (IsoMu24 || IsoTkMu24)", all(r["hlt_isomu24"] or r["hlt_isotkmu24"] for r in allrecs))
    ck.check_true("era G", all(r["era"] == "G" for r in allrecs))
    ck.check_true("npv >= 1 (clean)", all(r["npv"] >= 1 for r in allrecs))
    # SR re-selection on the stored numbers
    ok_sel, ok_mass, ok_bare, worst = True, True, True, 0.0
    for r in sr_events:
        m1, m2 = r["muons"]
        sel = (_tight_ok(m1) and _tight_ok(m2) and r["n_tight"] == 2 and m1["pt"] > config.MU_PT_LEAD and (m1["matched"] or m2["matched"])
               and m1["charge"] * m2["charge"] < 0 and config.MASS_LO < r["mass"] < config.MASS_HI and m1["pt"] >= m2["pt"])
        ok_sel &= sel
        parts = [(m["pt"], m["eta"], m["phi"], m["mass"]) for m in r["muons"]] + [(g["pt"], g["eta"], g["phi"], 0.0) for g in r["fsr_photons"]]
        dm = abs(_mass_of(parts) - r["mass"]); worst = max(worst, dm)
        ok_mass &= dm < 1e-3
        ok_bare &= abs(_mass_of(parts[:2]) - r["mass_bare"]) < 1e-3
    ck.check_true("SR events re-select under the stored selection", ok_sel, "zmumu/regions.py")
    ck.check_true("SR mass == muons + attached FSR photons (1e-3 GeV)", ok_mass, detail=f"max |dm| = {worst:.2e}")
    ck.check_true("SR mass_bare == two muons (1e-3 GeV)", ok_bare)
    ck.check_true("20 SR events, 4 chosen with unique run/lumi/event", len(sr_events) == 20 and len(d["sr"]["chosen"]) == 4
                  and len({(r["run"], r["lumi"], r["event"]) for r in sr_events}) == 20)
    ch = {r["rule"]: r for r in d["sr"]["chosen"]}
    ck.check_true("chosen rules textbook/fsr/forward/boosted", list(ch) == list(SR_RULES))
    t = ch["textbook"]
    ck.check_true("textbook: 89 < m < 93, |eta| < 0.9, |dphi| > 3.0, no FSR, no jet > 30", 89 < t["mass"] < 93 and all(abs(m["eta"]) < 0.9 for m in t["muons"])
                  and abs(t["dphi"]) > 3.0 and not t["fsr_photons"] and t["n_jets30"] == 0, detail=f"m={t['mass']:.2f}")
    ck.check_true("fsr: m - m_bare > 1 GeV with >= 1 photon", ch["fsr"]["mass"] - ch["fsr"]["mass_bare"] > 1.0 and len(ch["fsr"]["fsr_photons"]) >= 1)
    ck.check_true("forward: one |eta| > 2.1", max(abs(m["eta"]) for m in ch["forward"]["muons"]) > 2.1)
    ck.check_true("boosted: pT(mumu) > 60", ch["boosted"]["zpt"] > 60)
    ck.check_true("500 SR masses in 60-120", len(d["sr"]["masses_first500"]) == 500 and all(60 <= m <= 120 for m in d["sr"]["masses_first500"]))
    # same-sign tag + probe
    ck.check_true("ten chosen: exactly 1 passes, 9 anti", len(ten) == 10 and sum(r["passes"] for r in ten) == 1)
    ck.check_true("ten chosen: unique events", len({(r["run"], r["lumi"], r["event"]) for r in ten}) == 10)
    ck.check_true("ten chosen: first is anti-isolated with the probe inside a jet", (not ten[0]["passes"]) and ten[0]["probe"]["in_jet"])
    ok = True
    for r in ten + pool:
        tg, pr = r["tag"], r["probe"]
        ok &= (tg["charge"] == pr["charge"] and _tight_ok(tg) and tg["pt"] > config.MU_PT_LEAD and tg["matched"] and r["mass_bare"] > fakes.MIN_PAIR_MASS
               and (_tight_ok(pr) if r["passes"] else _anti_ok(pr)) and pr["cls"] == ("tight" if r["passes"] else "anti")
               and abs(pr["eta"]) < 0.9 and pr["pt"] < 50)
        ok &= abs(_mass_of([(tg["pt"], tg["eta"], tg["phi"], tg["mass"]), (pr["pt"], pr["eta"], pr["phi"], pr["mass"])]) - r["mass_bare"]) < 1e-3
        c = r["ff_cell"]
        ok &= abs(c["ff"] - ffmap[c["ipt"], c["ieta"]]) < 1e-12 and c["ieta"] == 0 and c["ipt"] <= 3
        ok &= c["pt_range"][0] <= pr["pt"] < c["pt_range"][1]
    ck.check_true("SS pairs: same sign, tag tight pT > 26 matched, m_bare > 12, probe class consistent, barrel FF cell correct, mass recomputed", ok)
    ck.check_true(f"pool >= 3 tight and >= 27 anti (got {d['ss_tp']['pool_n_tight']} / {d['ss_tp']['pool_n_anti']})",
                  d["ss_tp"]["pool_n_tight"] >= 3 and d["ss_tp"]["pool_n_anti"] >= 27)
    ck.check("barrel FF cells", d["ss_tp"]["barrel_ff_cells"], ffmap[:4, 0], 1e-12, "fakes.json maps.nominal")
    ck.check_true(f"pool size <= {POOL_MAX}", len(pool) <= POOL_MAX)
    return ck


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--files", type=int, default=6, help="number of data_2016G skim files to read (default 6)")
    ap.add_argument("--seed", type=int, default=20260915)
    standard_args(ap, "zmumu_events.json")
    args = ap.parse_args()
    finalize(args, lambda: build(args), verify, "zmumu_events", indent=None)   # compact: ~20 events + 120 pairs


if __name__ == "__main__":
    main()
