#!/usr/bin/env python
"""Step 4 -- histograms: TRExFitter inputs with all variations, control plots, yields (docs/07, docs/08).

    python scripts/step4_histograms.py [--ff-variant nominal|mcsub]

Fit inputs (fitting/CONVENTIONS.md): fit/fitinputs/ztautau[_mcsub].root with
    tautau_SR__Data, tautau_SR__<sample>, tautau_SR__<sample>__<syst>Up/Down
for the fit variable (config.FIT_VARIABLE) in the signal region. Samples: DYtautau (signal), DYee,
DYmumu, WJets, TTbar, SingleTop, WW, WZ, ZZ (simulation, leading tau not a jet) and Fakes (data x FF).
Also writes output/data/yields[_mcsub].json and output/plots/step4_*.png (data vs prediction in the SR
and in the fake-factor regions for many variables) and a stat-only sensitivity comparison of the
di-tau mass estimators.
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

REGION = config.REGION
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
}
STACK_GROUPS = [("DYll", ["DYee", "DYmumu"]), ("DYlowmass", ["DYlowmass"]), ("Diboson", ["WW", "WZ", "ZZ"]), ("Top", ["TTbar", "SingleTop"]),
                ("WJets", ["WJets"]), ("Fakes", ["Fakes"]), ("DYtautau", ["DYtautau"])]
FIT_SAMPLES = [samples.SIGNAL] + samples.FIT_BACKGROUNDS + ["Fakes"]


def h1(x, w, edges):
    edges = np.asarray(edges, dtype=float)
    xc = np.clip(x, edges[0], edges[-1] - 1e-6)
    v, _ = np.histogram(xc, bins=edges, weights=w)
    v2, _ = np.histogram(xc, bins=edges, weights=w ** 2)
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ff-variant", choices=["nominal", "mcsub"], default="mcsub" if config.FF_SUBTRACT_MC else "nominal")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    variant = args.ff_variant
    suffix = "" if variant == ("mcsub" if config.FF_SUBTRACT_MC else "nominal") else f"_{variant}"
    ffres = json.loads((config.DATA_DIR / "fakefactors.json").read_text())
    table = fakes.from_json(ffres[variant]["ff"])
    c_osss = np.asarray(ffres[variant]["osss"]["C"])                     # per (era, jet-multiplicity category)
    c_rel = np.hypot(np.asarray(ffres[variant]["osss"]["stat"]) / c_osss, config.FF_OSSS_SYST)
    nonclosure = np.clip(np.asarray(ffres[variant]["nonclosure_m_tt"]["ratio"]), 0.5, 1.5)
    edges = config.FIT_BINS
    var = config.FIT_VARIABLE
    fit = Book()
    ctrl = {rg: {v: Book() for v in CONTROL_VARS} for rg in ("SR", "AR", "SS_T", "OSAI_T")}

    # ------------------------------------------------------------------ data and fakes
    data = analysis.load_data()
    reg = analysis.regions(data)
    ones = np.ones(len(data["run"]))
    fit.add(f"{REGION}__Data", h1(data[var][reg["SR"]], ones[reg["SR"]], edges))
    wf = fakes.fake_weights(data, reg["AR"], table, c_osss)
    ar = reg["AR"]
    fit.add(f"{REGION}__Fakes", h1(data[var][ar], wf[ar], edges))
    for dm in fakes.DMS:
        for d_, sgn in (("Up", 1), ("Down", -1)):
            w = fakes.fake_weights(data, ar, table, c_osss, dm, sgn)
            fit.add(f"{REGION}__Fakes__FakeStat_tautau_DM{dm}{d_}", h1(data[var][ar], w[ar], edges))
    for d_, sgn in (("Up", 1), ("Down", -1)):
        rel = fakes.per_event(c_rel, data)
        fit.add(f"{REGION}__Fakes__FakeOSSS_tautau{d_}", h1(data[var][ar], (wf * (1 + sgn * rel))[ar], edges))
    idx = np.clip(np.searchsorted(np.asarray(edges), data[var], side="right") - 1, 0, len(edges) - 2)
    for d_, f in (("Up", nonclosure), ("Down", 1.0 / nonclosure)):
        fit.add(f"{REGION}__Fakes__FakeClosure_tautau{d_}", h1(data[var][ar], (wf * f[idx])[ar], edges))
    for rg, book in ctrl.items():
        for v, (ed, _, _) in CONTROL_VARS.items():
            book[v].add("Data", h1(data[v][reg[rg]], ones[reg[rg]], ed))
            if rg == "SR":
                book[v].add("Fakes", h1(data[v][ar], wf[ar], ed))
            if rg == "SS_T":
                wss = fakes.fake_weights(data, reg["SS_L"], table, 1.0)
                book[v].add("Fakes", h1(data[v][reg["SS_L"]], wss[reg["SS_L"]], ed))
            if rg == "OSAI_T":
                tai = fakes.from_json(ffres[variant]["ff_ai"])
                wai = fakes.fake_weights(data, reg["OSAI_L"], tai, c_osss)
                book[v].add("Fakes", h1(data[v][reg["OSAI_L"]], wai[reg["OSAI_L"]], ed))
    fake_syst_names = [f"FakeStat_tautau_DM{dm}" for dm in fakes.DMS] + ["FakeOSSS_tautau", "FakeClosure_tautau"]

    # ------------------------------------------------------------------ simulation
    mc_syst_names = analysis.WEIGHT_SYSTS + analysis.KINEMATIC_SYSTS
    relaxed = {}
    for key in analysis.available_mc():
        d, meta = analysis.load(key)
        if not len(d):
            print(f"  {key}: empty ntuple")
            continue
        comps = analysis.mc_components(key)
        r = analysis.regions(d, is_mc=True)
        w = analysis.weights(d, key)
        for name, cm in comps:
            sr = r["SR"] & cm
            fit.add(f"{REGION}__{name}", h1(d[var][sr], w[sr], edges))
            if name in analysis.SMOOTHED_SAMPLES:
                m = r["SR_relaxed"] & cm
                relaxed[name] = relaxed.get(name, 0) + h1(d[var][m], w[m], edges)
            for rg, book in ctrl.items():
                for v, (ed, _, _) in CONTROL_VARS.items():
                    book[v].add(name, h1(d[v][r[rg] & cm], w[r[rg] & cm], ed))
            if variant == "mcsub":        # remove genuine taus from the FF application regions
                wm = fakes.fake_weights(d, r["AR"] & cm, table, c_osss) * w
                m = r["AR"] & cm
                fit.add(f"{REGION}__Fakes", -h1(d[var][m], wm[m], edges) * np.array([[1], [0]]))
                for fs in fake_syst_names:
                    for d_ in ("Up", "Down"):
                        fit.add(f"{REGION}__Fakes__{fs}{d_}", -h1(d[var][m], wm[m], edges) * np.array([[1], [0]]))
                book = ctrl["SR"]
                for v, (ed, _, _) in CONTROL_VARS.items():
                    book[v].add("Fakes", -h1(d[v][m], wm[m], ed) * np.array([[1], [0]]))
        for syst in analysis.WEIGHT_SYSTS:
            for d_ in ("Up", "Down"):
                ws = analysis.weights(d, key, syst, d_)
                for name, cm in comps:
                    sr = r["SR"] & cm
                    fit.add(f"{REGION}__{name}__{syst}{d_}", h1(d[var][sr], ws[sr], edges))
        for syst in analysis.KINEMATIC_SYSTS:
            for d_ in ("Up", "Down"):
                kin = analysis.kinematics(d, key, syst, d_)
                rs = analysis.regions(d, kin, is_mc=True)
                ws = analysis.weights(d, key, kin=kin)
                for name, cm in comps:
                    sr = rs["SR"] & cm
                    fit.add(f"{REGION}__{name}__{syst}{d_}", h1(kin[var][sr], ws[sr], edges))
        if key == "DY_NLO":
            sr = r["SR"] & (d["gen_lhe_flavour"] == 15)
            nominal = h1(d[var][sr], w[sr], edges)
            for syst in analysis.THEORY_SYSTS:
                fac = analysis.theory_weights(d, key, syst)
                if fac is None:
                    continue
                members = [h1(d[var][sr], (w * fac[:, i])[sr], edges)[0] for i in range(fac.shape[1])]
                up, dn = analysis.combine_theory(syst, nominal[0], members)
                fit.add(f"{REGION}__{samples.SIGNAL}__{syst}Up", np.stack([up, nominal[1]]))
                fit.add(f"{REGION}__{samples.SIGNAL}__{syst}Down", np.stack([dn, nominal[1]]))
        print(f"  {key}: " + ", ".join(f"{n} {fit.h[f'{REGION}__{n}'][0].sum():.1f}" for n, _ in comps
                                        if f"{REGION}__{n}" in fit.h), flush=True)

    # smoothed templates: shape from the relaxed selection, normalisation of each variation from the SR;
    # the statistical uncertainty of the SR normalisation becomes one OVERALL nuisance parameter
    smooth_info = {}
    for sname in analysis.SMOOTHED_SAMPLES:
        shape = relaxed.get(sname)
        nom = fit.h.get(f"{REGION}__{sname}")
        if shape is None or nom is None or shape[0].sum() <= 0 or nom[0].sum() <= 0:
            continue
        norm_shape = np.maximum(shape[0], 0) / np.maximum(shape[0], 0).sum()
        smooth_info[sname] = {"sr_yield": float(nom[0].sum()), "sr_stat_rel": float(np.sqrt(nom[1].sum()) / nom[0].sum()),
                              "relaxed_yield": float(shape[0].sum())}
        for name in [n for n in fit.h if n == f"{REGION}__{sname}" or n.startswith(f"{REGION}__{sname}__")]:
            tot = fit.h[name][0].sum()
            fit.h[name] = np.stack([norm_shape * tot, (norm_shape * tot) ** 2 * (shape[1].sum() / shape[0].sum() ** 2)])
        print(f"  {sname}: SR shape from SR_relaxed ({shape[0].sum():.0f} events), SR yield {nom[0].sum():.0f} "
              f"+- {100 * smooth_info[sname]['sr_stat_rel']:.0f}% (MC stat)")

    # alternative generator (LO madgraph) for the signal C factor, normalised to the NLO fiducial yield
    try:
        dlo, _ = analysis.load("DY_LO")
        nlo, lo = analysis.signal_prediction("DY_NLO"), analysis.signal_prediction("DY_LO")
        rlo = analysis.regions(dlo, is_mc=True)
        wlo = analysis.weights(dlo, "DY_LO") * nlo["sigma_fid_pb"] / lo["sigma_fid_pb"]
        m = rlo["SR"] & (dlo["gen_lhe_flavour"] == 15)
        fit.add(f"{REGION}__{samples.SIGNAL}__SigModelUp", h1(dlo[var][m], wlo[m], edges))
        has_lo = True
    except FileNotFoundError:
        print("  DY_LO ntuple missing: no SigModel variation")
        has_lo = False

    # ------------------------------------------------------------------ write fit inputs
    names = {}
    for name, vv in fit.h.items():
        names[name] = to_hist(vv, edges)
    # samples without any event still need a (tiny) template for TRExFitter
    for s in FIT_SAMPLES:
        nm = f"{REGION}__{s}"
        if nm not in names:
            names[nm] = to_hist(np.stack([np.full(len(edges) - 1, 1e-6), np.zeros(len(edges) - 1)]), edges)
    out = config.FIT_DIR / "fitinputs" / f"{config.JOB}{suffix}.root"
    sig = analysis.signal_prediction("DY_NLO")
    meta = {"variant": variant, "fit_variable": var, "bins": edges, "lumi_pb": config.LUMI_PB, "C_osss": c_osss.tolist(),
            "C_osss_rel_unc": c_rel.tolist(), "signal_prediction": sig, "has_sigmodel": has_lo,
            "mc_systs": mc_syst_names, "theory_systs": analysis.THEORY_SYSTS, "fake_systs": fake_syst_names,
            "smoothed": smooth_info}
    rep = trexhist.write_fitinputs(out, names, meta=meta)
    print(f"fit inputs -> {out}: {rep['n_hists']} histograms, clipped {len(rep['clipped'])}")

    yields = {s: {"value": float(fit.h[f'{REGION}__{s}'][0].sum()) if f'{REGION}__{s}' in fit.h else 0.0,
                  "stat": float(np.sqrt(fit.h[f'{REGION}__{s}'][1].sum())) if f'{REGION}__{s}' in fit.h else 0.0}
              for s in FIT_SAMPLES}
    yields["Data"] = {"value": float(fit.h[f"{REGION}__Data"][0].sum())}
    pred = sum(yields[s]["value"] for s in FIT_SAMPLES)
    yields["Total"] = {"value": pred}
    # relative size of every variation on every sample (prefit)
    impacts = {}
    for name, vv in fit.h.items():
        parts = name.split("__")
        if len(parts) == 3:
            nom = fit.h.get(f"{parts[0]}__{parts[1]}")
            if nom is not None and nom[0].sum() > 0:
                impacts.setdefault(parts[1], {})[parts[2]] = float(vv[0].sum() / nom[0].sum() - 1)
    (config.DATA_DIR / f"yields{suffix}.json").write_text(json.dumps({"yields": yields, "prefit_norm_effects": impacts,
                                                                        "meta": meta}, indent=1, default=float))
    print("SR yields: " + ", ".join(f"{s} {v['value']:.0f}" for s, v in yields.items()))

    # ------------------------------------------------------------------ sensitivity of the mass variables
    sens = {}
    for v in ("m_tt", "m_vis"):
        book = ctrl["SR"][v]
        s_ = book.h.get("DYtautau", np.zeros((2, 1)))[0]
        b_ = sum(book.h[k][0] for k in book.h if k not in ("Data", "DYtautau"))
        ok = (s_ + b_) > 0
        sens[v] = float(1.0 / np.sqrt(np.sum(s_[ok] ** 2 / (s_[ok] + b_[ok]))))
    print("stat-only relative sigma(mu) from the SR shape (backgrounds fixed): " + ", ".join(f"{k} {v:.4f}" for k, v in sens.items()))
    js = json.loads((config.DATA_DIR / f"yields{suffix}.json").read_text())
    js["stat_only_sensitivity"] = sens
    (config.DATA_DIR / f"yields{suffix}.json").write_text(json.dumps(js, indent=1))

    if args.no_plots:
        return
    titles = {"SR": "signal region (OS, both Medium)", "AR": r"application region (OS, $\tau_1$ fails Medium)",
              "SS_T": "same sign, both Medium (FF closure)", "OSAI_T": r"OS, $\tau_2$ anti-isolated (C$_{OS/SS}$ check)"}
    for rg, books in ctrl.items():
        for v, (ed, xl, logy) in CONTROL_VARS.items():
            book = books[v].h
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
            if logy and rg == "SR":
                plotting.stack_plot(config.PLOT_DIR / f"step4_{rg}_{v}_log{suffix}.png", ed, (dv, np.sqrt(dv)), stack,
                                    xl, title=titles[rg], logy=True)


if __name__ == "__main__":
    main()
