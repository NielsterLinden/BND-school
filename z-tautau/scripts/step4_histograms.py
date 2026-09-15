#!/usr/bin/env python
"""Step 4 -- histograms: TRExFitter inputs with all variations, control plots, yields (docs/07, docs/08, docs/09).

    python scripts/step4_histograms.py [--ff-variant mcsub|nosub] [--no-plots]

Fit inputs (fitting/CONVENTIONS.md): fit/fitinputs/ztautau[_nosub].root with, for every BDT category
region tautau_SR<k> (config.REGIONS),
    tautau_SR<k>__Data, tautau_SR<k>__<sample>, tautau_SR<k>__<sample>__<syst>Up/Down
of the fit variable (config.FIT_VARIABLE). Samples: DYtautau (fiducial signal), DYtautau_nonfid, DYee,
DYmumu, DYlowmass, WJets, TTbar, SingleTop, WW, WZ, ZZ (simulation, leading tau not a jet) and Fakes
(data x FF, simulation with a genuine leading tau subtracted). The fake-factor statistics enter the Fakes
template variance per bin; the residual same-sign non-closure per category and mass region is one
nuisance parameter each (FakeClosure_tautau_c<k>_lo/hi).
Also writes output/data/yields[_nosub].json and output/plots/step4_*.png (data vs prediction in the SR,
per category, and in the fake-factor regions for many variables).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import hist
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root: fitting/
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fitting import trexhist  # noqa: E402
from ztautau import analysis, config, fakes, plotting, samples  # noqa: E402

REGIONS = config.REGIONS
NCAT = len(REGIONS)
CONTROL_VARS = {
    "m_tt": (config.FIT_BINS, r"$m_{\tau\tau}$ (MET likelihood) [GeV]", True),
    "m_vis": (np.arange(0, 305, 10), r"$m_{vis}$ [GeV]", False),
    "m_col": (np.arange(0, 405, 20), r"$m_{col}$ (collinear, valid events) [GeV]", False),
    "t1_pt": (np.array([40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 120, 150, 200, 300]), r"$p_T(\tau_1)$ [GeV]", True),
    "t2_pt": (np.array([40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 120, 150]), r"$p_T(\tau_2)$ [GeV]", True),
    "t1_eta": (np.linspace(-2.1, 2.1, 22), r"$\eta(\tau_1)$", False),
    "t2_eta": (np.linspace(-2.1, 2.1, 22), r"$\eta(\tau_2)$", False),
    "t1_dm": (np.array([-0.5, 0.5, 1.5, 9.5, 10.5, 11.5]), r"decay mode $\tau_1$ (0, 1, 10, 11)", False),
    "t2_dm": (np.array([-0.5, 0.5, 1.5, 9.5, 10.5, 11.5]), r"decay mode $\tau_2$ (0, 1, 10, 11)", False),
    "met": (np.arange(0, 155, 10), r"$p_T^{miss}$ [GeV]", False),
    "pt_tt": (np.arange(0, 205, 10), r"$p_T^{\tau\tau}$ (incl. MET) [GeV]", False),
    "dr_tt": (np.linspace(0.5, 5.0, 19), r"$\Delta R(\tau_1, \tau_2)$", False),
    "njets": (np.arange(-0.5, 6.5, 1), r"$N_{jets}$ ($p_T$ > 30 GeV)", False),
    "nbjets": (np.arange(-0.5, 3.5, 1), r"$N_{b-jets}$", False),
    "npv": (np.arange(0, 60, 3), r"$N_{PV}$", False),
    "mt_tot": (np.arange(0, 305, 15), r"$m_T^{tot}$ [GeV]", False),
    "bdt": (np.linspace(0, 1, 21), "BDT score", True),
}
CATEGORY_VARS = ["m_tt", "m_vis", "t1_pt", "t2_pt", "dr_tt", "met", "t1_eta", "njets"]
STACK_GROUPS = [("DYll", ["DYee", "DYmumu"]), ("DYlowmass", ["DYlowmass"]), ("Diboson", ["WW", "WZ", "ZZ"]),
                ("Top", ["TTbar", "SingleTop"]), ("WJets", ["WJets"]), ("Fakes", ["Fakes"]),
                ("DYtautau_nonfid", ["DYtautau_nonfid"]), ("DYtautau", ["DYtautau"])]
FIT_SAMPLES = [samples.SIGNAL] + samples.FIT_BACKGROUNDS + ["Fakes"]
DY_COMPS = (samples.SIGNAL, samples.SIGNAL_NONFID)


def h1(x, w, edges, w2=None):
    """2 x nbins: sum of weights and sum of squares (of `w2` if given, else of `w`)."""
    edges = np.asarray(edges, dtype=float)
    xc = np.clip(x, edges[0], edges[-1] - 1e-6)
    v, _ = np.histogram(xc, bins=edges, weights=w)
    v2, _ = np.histogram(xc, bins=edges, weights=(w if w2 is None else w2) ** 2)
    return np.stack([v, v2])


def to_hist(vv, edges):
    h = hist.Hist(hist.axis.Variable(np.asarray(edges, dtype=float)), storage=hist.storage.Weight())
    h.view().value[...] = vv[0]
    h.view().variance[...] = vv[1]
    return h


class Book:
    """name -> 2 x nbins array (sum w, sum w^2), added on repeated fills."""

    def __init__(self):
        self.h = {}

    def add(self, name, vv):
        self.h[name] = self.h[name] + vv if name in self.h else vv.copy()

    def get(self, name, n):
        return self.h.get(name, np.zeros((2, n)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    nominal_variant = "mcsub" if config.FF_SUBTRACT_MC else "nosub"
    ap.add_argument("--ff-variant", choices=["mcsub", "nosub"], default=nominal_variant)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    variant = args.ff_variant
    subtract = variant == "mcsub"
    suffix = "" if variant == nominal_variant else f"_{variant}"
    ffres = json.loads((config.DATA_DIR / "fakefactors.json").read_text())
    table = fakes.from_json(ffres[variant]["ff"])
    edges = np.asarray(config.FIT_BINS)
    nb = len(edges) - 1
    var = config.FIT_VARIABLE
    fit = Book()
    ctrl = {rg: {v: Book() for v in CONTROL_VARS} for rg in ("SR", "AR", "SS_T", "OSAI_T")}
    ctrl.update({rg: {v: Book() for v in CATEGORY_VARS} for rg in REGIONS})

    def val(d, v, kin=None, sc=None):
        if v == "bdt":
            return sc
        return (kin or d)[v] if v in analysis.KIN_KEYS else d[v]

    # ------------------------------------------------------------------ data, simulation, C per category
    data = analysis.load_data()
    reg = analysis.regions(data)
    sc_d = analysis.scores(data, "data")
    cat_d = analysis.categories(data, "data")
    data["_cat"] = cat_d
    ones = np.ones(len(data["run"]))
    loaded = []
    for key in analysis.available_mc():
        d, meta = analysis.load(key)
        if not len(d):
            print(f"  {key}: empty ntuple")
            continue
        d = dict(d)
        d["_cat"] = analysis.categories(d, key)
        loaded.append((key, d, analysis.regions(d, is_mc=True), analysis.weights(d, key)))
    # the OS/SS correction per (era, N_jets, BDT category): the charge correlation of the two jets depends on
    # the topology the BDT selects (inclusive C under-predicts the signal-like categories by 3-6%, docs/05)
    sub_c = [(d, r, analysis.subtraction_weights(k, w)) for k, d, r, w in loaded if k not in analysis.SUBTRACT_EXCLUDE] if subtract else []
    osss = fakes.osss_correction(data, reg, sub_c, n_cat=NCAT)
    c_osss = np.asarray(osss["C"])                                      # (era, N_jets, category)
    c_rel = np.hypot(np.asarray(osss["stat"]) / c_osss, config.FF_OSSS_SYST)
    if osss.get("replaced_by_njet_inclusive"):
        print(f"  C_OS/SS bins (era, N_jets, category) replaced by the N_jets-inclusive value (stat > 10%): {osss['replaced_by_njet_inclusive']}")
    for k in range(NCAT):
        print(f"  C_OS/SS category {k}: " + "  ".join(f"era {'GH'[e]}: " + ", ".join(f"{c_osss[e, j, k]:.3f}+-{np.asarray(osss['stat'])[e, j, k]:.3f}" for j in range(c_osss.shape[1])) for e in range(c_osss.shape[0])))
    wf, wf2 = fakes.fake_weights(data, reg["AR"], table, c_osss, with_err=True)
    ar = reg["AR"]
    ss_book = {k: Book() for k in range(NCAT)}            # same-sign closure per category, MC subtracted
    wss = fakes.fake_weights(data, reg["SS_L"], table, 1.0)
    for k, region in enumerate(REGIONS):
        sr_k, ar_k = reg["SR"] & (cat_d == k), ar & (cat_d == k)
        fit.add(f"{region}__Data", h1(data[var][sr_k], ones[sr_k], edges))
        fit.add(f"{region}__Fakes", h1(data[var][ar_k], wf[ar_k], edges, wf2[ar_k]))
        for d_, sgn in (("Up", 1), ("Down", -1)):
            rel = fakes.per_event(c_rel, data)
            fit.add(f"{region}__Fakes__FakeOSSS_tautau_c{k}{d_}", h1(data[var][ar_k], (wf * (1 + sgn * rel))[ar_k], edges))
        ss_book[k].add("obs", h1(data[var][reg["SS_T"] & (cat_d == k)], ones[reg["SS_T"] & (cat_d == k)], edges))
        ss_book[k].add("pred", h1(data[var][reg["SS_L"] & (cat_d == k)], wss[reg["SS_L"] & (cat_d == k)], edges))
        for v in CATEGORY_VARS:
            ed = CONTROL_VARS[v][0]
            ctrl[region][v].add("Data", h1(val(data, v)[sr_k], ones[sr_k], ed))
            ctrl[region][v].add("Fakes", h1(val(data, v)[ar_k], wf[ar_k], ed, wf2[ar_k]))
    for rg, book in ((r, ctrl[r]) for r in ("SR", "AR", "SS_T", "OSAI_T")):
        for v, (ed, _, _) in CONTROL_VARS.items():
            x = val(data, v, sc=sc_d)
            book[v].add("Data", h1(x[reg[rg]], ones[reg[rg]], ed))
            if rg == "SR":
                book[v].add("Fakes", h1(x[ar], wf[ar], ed, wf2[ar]))
            if rg == "SS_T":
                book[v].add("Fakes", h1(x[reg["SS_L"]], wss[reg["SS_L"]], ed))
            if rg == "OSAI_T":
                tai = osss["table_ai"]
                wai = fakes.fake_weights(data, reg["OSAI_L"], tai, c_osss)
                book[v].add("Fakes", h1(x[reg["OSAI_L"]], wai[reg["OSAI_L"]], ed))

    # ------------------------------------------------------------------ simulation
    mc_syst_names = analysis.WEIGHT_SYSTS + analysis.KINEMATIC_SYSTS
    relaxed = {}
    dy_keys = analysis.dy_stitch_keys()
    for key, d, r, w in loaded:
        comps = analysis.mc_components(key)
        wsub = analysis.subtraction_weights(key, w)       # W+jets: uniform weights in the FF regions
        sc_m = analysis.scores(d, key)
        cat_m = d["_cat"]
        for name, cm in comps:
            for k, region in enumerate(REGIONS):
                sr = r["SR"] & cm & (cat_m == k)
                fit.add(f"{region}__{name}", h1(d[var][sr], w[sr], edges))
                if name in analysis.SMOOTHED_SAMPLES:
                    m = r["SR_relaxed"] & cm & (cat_m == k)
                    relaxed[(name, region)] = relaxed.get((name, region), 0) + h1(d[var][m], w[m], edges)
                for v in CATEGORY_VARS:
                    ctrl[region][v].add(name, h1(val(d, v)[sr], w[sr], CONTROL_VARS[v][0]))
                if subtract and key not in analysis.SUBTRACT_EXCLUDE:   # genuine taus out of the FF regions
                    m = r["AR"] & cm & (cat_m == k)
                    wm = fakes.fake_weights(d, m, table, c_osss) * wsub
                    neg = -h1(d[var][m], wm[m], edges) * np.array([[1], [0]])
                    fit.add(f"{region}__Fakes", neg)
                    for d_ in ("Up", "Down"):
                        fit.add(f"{region}__Fakes__FakeOSSS_tautau_c{k}{d_}", neg)
                    for v in CATEGORY_VARS:
                        ctrl[region][v].add("Fakes", -h1(val(d, v)[m], wm[m], CONTROL_VARS[v][0]) * np.array([[1], [0]]))
                    mss = r["SS_T"] & cm & (cat_m == k)
                    ss_book[k].add("obs_mc", h1(d[var][mss], wsub[mss], edges))
                    msl = r["SS_L"] & cm & (cat_m == k)
                    wm = fakes.fake_weights(d, msl, table, 1.0) * wsub
                    ss_book[k].add("pred_mc", h1(d[var][msl], wm[msl], edges))
            for rg, book in ((rr, ctrl[rr]) for rr in ("SR", "AR", "SS_T", "OSAI_T")):
                for v, (ed, _, _) in CONTROL_VARS.items():
                    x = val(d, v, sc=sc_m)
                    book[v].add(name, h1(x[r[rg] & cm], w[r[rg] & cm], ed))
                    if rg == "SR" and subtract and key not in analysis.SUBTRACT_EXCLUDE:
                        m = r["AR"] & cm
                        wm = fakes.fake_weights(d, m, table, c_osss) * wsub
                        book[v].add("Fakes", -h1(x[m], wm[m], ed) * np.array([[1], [0]]))
        for syst in analysis.WEIGHT_SYSTS:
            for d_ in ("Up", "Down"):
                ws = analysis.weights(d, key, syst, d_)
                for name, cm in comps:
                    for k, region in enumerate(REGIONS):
                        sr = r["SR"] & cm & (cat_m == k)
                        fit.add(f"{region}__{name}__{syst}{d_}", h1(d[var][sr], ws[sr], edges))
        for syst in analysis.KINEMATIC_SYSTS:
            for d_ in ("Up", "Down"):
                kin = analysis.kinematics(d, key, syst, d_)
                rs = analysis.regions(d, kin, is_mc=True)
                ws = analysis.weights(d, key, kin=kin)
                cs = analysis.categories(d, key, kin, syst, d_)
                for name, cm in comps:
                    for k, region in enumerate(REGIONS):
                        sr = rs["SR"] & cm & (cs == k)
                        fit.add(f"{region}__{name}__{syst}{d_}", h1(kin[var][sr], ws[sr], edges))
        if key in dy_keys:
            for name, cm in comps:
                if name not in DY_COMPS:
                    continue
                for syst in analysis.THEORY_SYSTS:
                    fac = analysis.theory_weights(d, key, syst)
                    if fac is None:
                        continue
                    for k, region in enumerate(REGIONS):
                        sr = r["SR"] & cm & (cat_m == k)
                        members = np.stack([h1(d[var][sr], (w * fac[:, i])[sr], edges)[0] for i in range(fac.shape[1])])
                        # accumulate the member histograms over the stitched samples, combine at the end
                        fit.add(f"__members__{region}__{name}__{syst}", np.concatenate([members, np.zeros((1, nb))]))
        print(f"  {key}: " + ", ".join(f"{n} {sum(fit.get(f'{rg}__{n}', nb)[0].sum() for rg in REGIONS):.1f}" for n, _ in comps), flush=True)
    # theory variations: envelope / quadrature on the summed member histograms
    for name in list(fit.h):
        if not name.startswith("__members__"):
            continue
        _, _, region, sample, syst = name.split("__")
        members = fit.h.pop(name)[:-1]
        nominal = fit.get(f"{region}__{sample}", nb)
        up, dn = analysis.combine_theory(syst, nominal[0], members)
        fit.add(f"{region}__{sample}__{syst}Up", np.stack([up, nominal[1]]))
        fit.add(f"{region}__{sample}__{syst}Down", np.stack([dn, nominal[1]]))

    # smoothed templates: shape from the relaxed selection, normalisation of each variation from the SR;
    # the statistical uncertainty of the SR normalisation becomes one OVERALL nuisance parameter
    smooth_info = {}
    for sname in analysis.SMOOTHED_SAMPLES:
        tot_nom = sum(fit.get(f"{rg}__{sname}", nb) for rg in REGIONS)
        if tot_nom[0].sum() <= 0:
            continue
        smooth_info[sname] = {"sr_yield": float(tot_nom[0].sum()), "sr_stat_rel": float(np.sqrt(tot_nom[1].sum()) / tot_nom[0].sum())}
        for region in REGIONS:
            shape = relaxed.get((sname, region))
            nom = fit.h.get(f"{region}__{sname}")
            if shape is None or nom is None or shape[0].sum() <= 0 or nom[0].sum() <= 0:
                continue
            norm_shape = np.maximum(shape[0], 0) / np.maximum(shape[0], 0).sum()
            for name in [n for n in fit.h if n == f"{region}__{sname}" or n.startswith(f"{region}__{sname}__")]:
                tot = fit.h[name][0].sum()
                fit.h[name] = np.stack([norm_shape * tot, (norm_shape * tot) ** 2 * (shape[1].sum() / shape[0].sum() ** 2)])
        print(f"  {sname}: shape from SR_relaxed per category, SR yield {tot_nom[0].sum():.0f} "
              f"+- {100 * smooth_info[sname]['sr_stat_rel']:.0f}% (MC stat)")

    # alternative generator (LO madgraph) for the *fiducial* signal C factor, normalised to the NLO fiducial
    # cross section: a one-sided, normalisation-only variation SigModel_tautau (this channel only)
    sigmodel = {}
    try:
        dlo, _ = analysis.load("DY_LO")
        dlo = dict(dlo)
        nlo, lo = analysis.signal_prediction(samples.DY_INCLUSIVE), analysis.signal_prediction("DY_LO")
        rlo = analysis.regions(dlo, is_mc=True)
        clo = analysis.categories(dlo, "DY_LO")
        wlo = analysis.weights(dlo, "DY_LO") * nlo["sigma_fid_pb"] / lo["sigma_fid_pb"]
        for k, region in enumerate(REGIONS):
            m = rlo["SR"] & (dlo["gen_lhe_flavour"] == 15) & dlo["gen_fid"].astype(bool) & (clo == k)
            fit.add(f"{region}__{samples.SIGNAL}__SigModel_tautauUp", h1(dlo[var][m], wlo[m], edges))
        tot_lo = sum(fit.get(f"{rg}__{samples.SIGNAL}__SigModel_tautauUp", nb)[0].sum() for rg in REGIONS)
        tot_nlo = sum(fit.get(f"{rg}__{samples.SIGNAL}", nb)[0].sum() for rg in REGIONS)
        sigmodel = {"C_LO_over_NLO_fiducial": tot_lo / tot_nlo}
        print(f"  SigModel_tautau: fiducial C(LO) / C(NLO) = {tot_lo / tot_nlo:.4f}")
    except FileNotFoundError:
        print("  DY_LO ntuple missing: no SigModel_tautau variation")

    # residual same-sign non-closure per category and mass region -> one nuisance parameter each
    closure_nps = {}
    lo_bins = edges[:-1] < config.FF_CLOSURE_MASS_SPLIT
    for k, region in enumerate(REGIONS):
        obs = ss_book[k].get("obs", nb)[0] - ss_book[k].get("obs_mc", nb)[0]
        pred = ss_book[k].get("pred", nb)[0] - ss_book[k].get("pred_mc", nb)[0]
        pred_var = ss_book[k].get("pred", nb)[1]
        nom = fit.get(f"{region}__Fakes", nb)
        for tag, sel in (("lo", lo_bins), ("hi", ~lo_bins)):
            o, p, pv = obs[sel].sum(), pred[sel].sum(), pred_var[sel].sum()
            ratio = o / p if p > 0 else 1.0
            stat = ratio * np.sqrt(1 / max(o, 1) + pv / max(p, 1e-9) ** 2)
            delta = float(np.hypot(ratio - 1, stat))
            closure_nps[f"FakeClosure_tautau_c{k}_{tag}"] = {"region": region, "ratio": float(ratio), "stat": float(stat),
                                                             "delta": delta, "obs": float(o), "pred": float(p)}
            for d_, sgn in (("Up", 1), ("Down", -1)):
                f = np.where(sel, 1 + sgn * delta, 1.0)
                fit.add(f"{region}__Fakes__FakeClosure_tautau_c{k}_{tag}{d_}", np.stack([nom[0] * f, nom[1]]))
            print(f"  {region} closure {tag}: obs/pred = {ratio:.3f} +- {stat:.3f} -> NP {100 * delta:.1f}%")

    # ------------------------------------------------------------------ write fit inputs
    names = {}
    for name, vv in fit.h.items():
        names[name] = to_hist(vv, edges)
    for region in REGIONS:
        for s in FIT_SAMPLES:
            nm = f"{region}__{s}"
            if nm not in names:      # samples without any event still need a (tiny) template for TRExFitter
                names[nm] = to_hist(np.stack([np.full(nb, 1e-6), np.zeros(nb)]), edges)
    out = config.FIT_DIR / "fitinputs" / f"{config.JOB}{suffix}.root"
    sig = analysis.signal_prediction(samples.DY_INCLUSIVE)
    meta = {"variant": variant, "fit_variable": var, "bins": edges.tolist(), "regions": REGIONS,
            "region_labels": config.REGION_LABELS, "category_edges": config.BDT_CATEGORY_EDGES, "lumi_pb": config.LUMI_PB,
            "C_osss": c_osss.tolist(), "C_osss_stat": osss["stat"], "C_osss_rel_unc": c_rel.tolist(), "signal_prediction": sig,
            "has_sigmodel": bool(sigmodel) and config.SIGMODEL_IN_FIT, "sigmodel": sigmodel, "mc_systs": mc_syst_names,
            "theory_systs": analysis.THEORY_SYSTS, "fake_systs": [f"FakeOSSS_tautau_c{k}" for k in range(NCAT)] + list(closure_nps),
            "osss_nps": {f"FakeOSSS_tautau_c{k}": {"region": REGIONS[k], "rel_unc": c_rel[:, :, k].tolist()} for k in range(NCAT)},
            "C_osss_replaced": osss.get("replaced_by_njet_inclusive", []), "C_osss_njet_inclusive": osss.get("C_njet_inclusive"),
            "closure_nps": closure_nps, "smoothed": smooth_info, "dy_samples": dy_keys}
    rep = trexhist.write_fitinputs(out, names, meta=meta)
    print(f"fit inputs -> {out}: {rep['n_hists']} histograms, clipped {len(rep['clipped'])}")

    yields = {"regions": {}, "yields": {}}
    for region in REGIONS:
        yields["regions"][region] = {s: {"value": float(fit.get(f"{region}__{s}", nb)[0].sum()),
                                         "stat": float(np.sqrt(fit.get(f"{region}__{s}", nb)[1].sum()))} for s in FIT_SAMPLES}
        yields["regions"][region]["Data"] = {"value": float(fit.get(f"{region}__Data", nb)[0].sum())}
    for s in FIT_SAMPLES + ["Data"]:
        yields["yields"][s] = {"value": sum(yields["regions"][r][s]["value"] for r in REGIONS)}
        if s != "Data":
            yields["yields"][s]["stat"] = float(np.sqrt(sum(yields["regions"][r][s]["stat"] ** 2 for r in REGIONS)))
    yields["yields"]["Total"] = {"value": sum(yields["yields"][s]["value"] for s in FIT_SAMPLES)}
    impacts = {}
    for name, vv in fit.h.items():
        parts = name.split("__")
        if len(parts) == 3:
            nom = fit.h.get(f"{parts[0]}__{parts[1]}")
            if nom is not None and nom[0].sum() > 0:
                impacts.setdefault(parts[1], {}).setdefault(parts[2], {})[parts[0]] = float(vv[0].sum() / nom[0].sum() - 1)
    # stat-only sensitivity of the m_tt shape, backgrounds fixed: inclusive vs categories
    sens = {}
    tot_s = sum(fit.get(f"{rg}__{samples.SIGNAL}", nb)[0] for rg in REGIONS)
    tot_b = sum(fit.get(f"{rg}__{s}", nb)[0] for rg in REGIONS for s in FIT_SAMPLES if s != samples.SIGNAL)

    def _sens(s_, b_):
        ok = (s_ + b_) > 0
        return float(np.sum(s_[ok] ** 2 / (s_[ok] + b_[ok])))
    sens["inclusive"] = 1 / np.sqrt(_sens(tot_s, tot_b))
    sens["categories"] = 1 / np.sqrt(sum(_sens(fit.get(f"{rg}__{samples.SIGNAL}", nb)[0],
                                              sum(fit.get(f"{rg}__{s}", nb)[0] for s in FIT_SAMPLES if s != samples.SIGNAL)) for rg in REGIONS))
    (config.DATA_DIR / f"yields{suffix}.json").write_text(json.dumps({**yields, "prefit_norm_effects": impacts, "meta": meta,
                                                                        "stat_only_sensitivity": sens}, indent=1, default=float))
    print("SR yields: " + ", ".join(f"{s} {v['value']:.0f}" for s, v in yields["yields"].items()))
    for region in REGIONS:
        y = yields["regions"][region]
        print(f"  {region}: data {y['Data']['value']:.0f}, fakes {y['Fakes']['value']:.0f}, {samples.SIGNAL} {y[samples.SIGNAL]['value']:.0f}, "
              f"{samples.SIGNAL_NONFID} {y[samples.SIGNAL_NONFID]['value']:.0f}")
    print(f"stat-only relative sigma(mu), backgrounds fixed: inclusive {sens['inclusive']:.4f}, categories {sens['categories']:.4f}")

    if args.no_plots:
        return
    titles = {"SR": "signal region (OS, both Medium)", "AR": r"application region (OS, $\tau_1$ fails Medium)",
              "SS_T": "same sign, both Medium (FF closure)", "OSAI_T": r"OS, $\tau_2$ anti-isolated (C$_{OS/SS}$ check)"}
    titles.update({rg: f"signal region, {lab}" for rg, lab in zip(REGIONS, config.REGION_LABELS)})
    for rg, books in ctrl.items():
        for v, book_ in books.items():
            ed, xl, logy = CONTROL_VARS[v]
            book = book_.h
            if "Data" not in book:
                continue
            stack = []
            for group, members in STACK_GROUPS:
                parts = [book[m] for m in members if m in book]
                if parts:
                    tot = sum(parts)
                    stack.append((group, np.maximum(tot[0], 0), tot[1]))
            dv = book["Data"][0]
            variable_bins = len(set(np.round(np.diff(ed), 6))) > 1 and v not in ("t1_dm", "t2_dm")
            plotting.stack_plot(config.PLOT_DIR / f"step4_{rg}_{v}{suffix}.png", ed, (dv, np.sqrt(dv)), stack, xl,
                                title=titles[rg], logy=False, density=variable_bins)
            if logy and rg in ("SR",) + tuple(REGIONS):
                plotting.stack_plot(config.PLOT_DIR / f"step4_{rg}_{v}_log{suffix}.png", ed, (dv, np.sqrt(dv)), stack,
                                    xl, title=titles[rg], logy=True)


if __name__ == "__main__":
    main()
