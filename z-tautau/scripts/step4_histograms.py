#!/usr/bin/env python
"""Step 4 (v4) -- fit inputs of the four-channel fit (docs/10-v4-plan.md sections 8-9).

    python scripts/step4_histograms.py [--channels tautau mutau etau emu] [--no-plots]

Regions: tautau_SR0/1/2 (the v3 categories: Data and Fakes templates copied from fit/fitinputs/ztautau.root, the
simulation rebuilt with the v4 signal definition, no VSjet scale factor, the 3% tau energy-scale prior and the
decay-mode split), mutau_SR, etau_SR, emu_SR, emu_CRtt.
Templates: <region>__Data, <region>__Fakes (+ <region>__Fakes__<np>Up/Down), <region>__<sample>_tDM<key>
(+ systematics), key = 'none' or the decay mode(s) of the genuine tau_h legs (docs/10 section 8).
Writes fit_v4/fitinputs/ztautau_v4.root (+ .meta.json with the systematic registry that step 5 turns into the
TRExFitter config), output_v4/data/yields_v4.json and control plots output_v4/plots/step4_<region>_<var>.png.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import hist
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fitting import trexhist  # noqa: E402
from ztautau import analysis, analysis_v4 as an, config, corrections, fakes, fakes_v4 as fk, plotting, samples  # noqa: E402

BINS = {r: config.FIT_BINS for r in config.REGIONS}
BINS.update({r: config.FIT_BINS_LTAU for r in config.LTAU_REGIONS + config.LTAU_PTSPLIT_REGIONS})
BINS.update({"emu_SR": config.FIT_BINS_EMU, "emu_CRtt": config.FIT_BINS_CRTT})
XSEC_UNC = {**samples.XSEC_UNC, config.V4_SIGNAL_OUT: 0.05}
TITLES = {"DYtautau": "Z/#gamma*#rightarrow#tau#tau (60-120)", "DYtautau_out": "Z/#gamma*#rightarrow#tau#tau (outside)",
          "DYee": "Z#rightarrowee", "DYmumu": "Z#rightarrow#mu#mu", "DYlowmass": "Z/#gamma*#rightarrowll (m<50)",
          "WJets": "W+jets", "TTbar": "t#bar{t}", "SingleTop": "single t", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ", "Fakes": "jet#rightarrow#tau_{h} / multijet"}
CONTROL = {
    "mutau": {"m_tt": (config.FIT_BINS_LTAU, r"$m_{\tau\tau}$ [GeV]"), "m_vis": (np.arange(0, 205, 10), r"$m_{vis}$ [GeV]"),
              "l_pt": (np.array([26, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150]), r"$p_T(\mu)$ [GeV]"),
              "t_pt": (np.array([30, 35, 40, 45, 50, 60, 70, 80, 100, 150]), r"$p_T(\tau_h)$ [GeV]"),
              "mt_1": (np.arange(0, 42, 4), r"$m_T(\mu, MET)$ [GeV]"), "met": (np.arange(0, 105, 7), r"$p_T^{miss}$ [GeV]"),
              "njets": (np.arange(-0.5, 5.5, 1), r"$N_{jets}$"), "t_dm": (np.array([-0.5, 0.5, 1.5, 9.5, 10.5, 11.5]), "decay mode")},
    "emu": {"m_tt": (config.FIT_BINS_EMU, r"$m_{\tau\tau}$ [GeV]"), "m_vis": (np.arange(0, 205, 10), r"$m_{vis}$ [GeV]"),
            "mu_pt": (np.array([10, 15, 20, 25, 30, 40, 50, 70, 100]), r"$p_T(\mu)$ [GeV]"), "el_pt": (np.array([13, 20, 25, 30, 40, 50, 70, 100]), r"$p_T(e)$ [GeV]"),
            "dzeta": (np.arange(-100, 105, 10), r"$D_\zeta$ [GeV]"), "met": (np.arange(0, 105, 7), r"$p_T^{miss}$ [GeV]"),
            "njets": (np.arange(-0.5, 5.5, 1), r"$N_{jets}$"), "dr": (np.linspace(0.3, 5, 20), r"$\Delta R(e,\mu)$")},
}
CONTROL["etau"] = {**CONTROL["mutau"], "l_pt": (np.array([29, 32, 35, 40, 45, 50, 60, 70, 80, 100, 150]), r"$p_T(e)$ [GeV]"), "mt_1": (np.arange(0, 42, 4), r"$m_T(e, MET)$ [GeV]")}
STACK_GROUPS = [("DYll", ["DYee", "DYmumu"]), ("DYlowmass", ["DYlowmass"]), ("Diboson", ["WW", "WZ", "ZZ"]), ("Top", ["TTbar", "SingleTop"]),
                ("WJets", ["WJets"]), ("Fakes", ["Fakes"]), ("DYtautau_out", ["DYtautau_out"]), ("DYtautau", ["DYtautau"])]


def h1(x, w, edges, w2=None):
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
    def __init__(self):
        self.h = {}

    def add(self, name, vv):
        self.h[name] = self.h[name] + vv if name in self.h else vv.copy()

    def get(self, name, n):
        return self.h.get(name, np.zeros((2, n)))


class Registry:
    """Systematic-parameter registry -> step 5 config (name -> dict(type, samples, regions, category, ...))."""

    def __init__(self):
        self.systs = {}
        self.samples = {}         # template name -> dict(base, key, regions, is_signal)

    def syst(self, name, samples_, regions_, category, kind="HISTO", subcategory=None, smoothing=None, up=None, down=None, title=None):
        s = self.systs.setdefault(name, {"name": name, "type": kind, "samples": [], "regions": [], "category": category,
                                         "subcategory": subcategory, "smoothing": smoothing, "up": up, "down": down, "title": title or name})
        for x in samples_:
            if x not in s["samples"]:
                s["samples"].append(x)
        for r in regions_:
            if r not in s["regions"]:
                s["regions"].append(r)

    def sample(self, tname, base, key, region):
        e = self.samples.setdefault(tname, {"base": base, "key": key, "regions": [], "is_signal": base == config.V4_SIGNAL})
        if region not in e["regions"]:
            e["regions"].append(region)


def tname(base, key):
    return f"{base}_tDM{key}"


# ------------------------------------------------------------------------------ simulation, generic
def fill_mc(book, reg_names, channel, key, d, w, comps, keys, region_masks, kin_fn, regs_fn, var, edges_of, registry, ctrl=None, ctrl_vars=None):
    """Nominal + every variation of one MC ntuple into the regions `reg_names` (region -> mask under nominal
    kinematics via `region_masks`; `regs_fn(kin, syst, dir)` returns the masks under a kinematic variation)."""
    nb = {r: len(edges_of[r]) - 1 for r in reg_names}
    for name, cm in comps:
        for r in reg_names:
            m = region_masks[r] & cm
            for k in np.unique(keys[m]):
                mk = m & (keys == k)
                t = tname(name, k)
                book.add(f"{r}__{t}", h1(d[var][mk], w[mk], edges_of[r]))
                registry.sample(t, name, k, r)
                if ctrl is not None and r in ctrl:
                    for v, (ed, _) in ctrl_vars.items():
                        ctrl[r][v].add(name, h1(d[v][mk], w[mk], ed))
    mc_all = sorted({tname(n, k) for n, _ in comps for k in np.unique(keys)})
    for syst in an.WEIGHT_SYSTS[channel]:
        if syst == "TopPt" and samples.SAMPLES[key]["group"] != "top":
            continue
        for d_ in ("Up", "Down"):
            ws = an.weights(d, key, channel, syst, d_)
            for name, cm in comps:
                for r in reg_names:
                    m = region_masks[r] & cm
                    for k in np.unique(keys[m]):
                        mk = m & (keys == k)
                        book.add(f"{r}__{tname(name, k)}__{syst}{d_}", h1(d[var][mk], ws[mk], edges_of[r]))
        registry.syst(syst, mc_all, reg_names, *SYST_INFO[syst])
    for syst in an.KINEMATIC_SYSTS[channel]:
        for d_ in ("Up", "Down"):
            kin, masks = kin_fn(syst, d_)
            ws = an.weights(d, key, channel, kin=kin)
            for name, cm in comps:
                for r in reg_names:
                    m = masks[r] & cm
                    for k in np.unique(keys[m]):
                        mk = m & (keys == k)
                        book.add(f"{r}__{tname(name, k)}__{syst}{d_}", h1(kin[var][mk], ws[mk], edges_of[r]))
        registry.syst(syst, mc_all, reg_names, *SYST_INFO[syst])
    if key in an.dy_stitch_keys(channel):
        for name, cm in comps:
            if name not in (config.V4_SIGNAL, config.V4_SIGNAL_OUT):
                continue
            for syst in an.THEORY_SYSTS:
                fac = an.theory_weights(d, key, channel, syst)
                if fac is None:
                    continue
                for r in reg_names:
                    m = region_masks[r] & cm
                    for k in np.unique(keys[m]):
                        mk = m & (keys == k)
                        members = np.stack([h1(d[var][mk], (w * fac[:, i])[mk], edges_of[r])[0] for i in range(fac.shape[1])])
                        book.add(f"__members__{r}__{tname(name, k)}__{syst}", np.concatenate([members, np.zeros((1, nb[r]))]))
                registry.syst(syst, [tname(name, k) for k in np.unique(keys[cm])], reg_names, "Signal modelling", title=syst)


SYST_INFO = {   # name -> (category, kind, subcategory, smoothing, up, down, title)
    "Pileup": ("Pileup", "HISTO", None, None, None, None, "Pileup"), "L1Prefiring": ("L1 prefiring", "HISTO", None, None, None, None, "L1 prefiring"),
    "TauFakeEle": ("Tau", "HISTO", "Tau ID", None, None, None, "e#rightarrow#tau_{h} ID"), "TauFakeMu": ("Tau", "HISTO", "Tau ID", None, None, None, "#mu#rightarrow#tau_{h} ID"),
    "MuonID": ("Muon efficiency", "HISTO", None, None, None, None, "muon ID"), "MuonIso": ("Muon efficiency", "HISTO", None, None, None, None, "muon isolation"),
    "MuonTrigger": ("Muon efficiency", "HISTO", None, None, None, None, "IsoMu24 trigger"), "ElectronReco": ("Electron efficiency", "HISTO", None, None, None, None, "electron reconstruction"),
    "ElectronID": ("Electron efficiency", "HISTO", None, None, None, None, "electron ID"), "ElectronTrigger": ("Electron trigger", "HISTO", None, None, None, None, "Ele27 trigger (in situ, plateau)"),
    "ElectronTrigger_lowpt": ("Electron trigger", "HISTO", None, None, None, None, "Ele27 trigger (in situ, 29-35 GeV turn-on)"),
    "EmuTrigger": ("Emu trigger", "HISTO", None, None, None, None, "e#mu cross trigger (in situ)"), "BTag": ("b tagging", "HISTO", None, None, None, None, "b-tag efficiency"),
    "TopPt": ("Background modelling", "HISTO", None, None, None, None, "top p_{T} reweighting"), "MET_Unclustered": ("MET", "HISTO", None, 40, None, None, "MET unclustered energy"),
    "JES": ("Jets", "HISTO", None, 40, None, None, "jet energy scale (total)"), "MuonScale": ("Muon momentum", "HISTO", None, 40, None, None, "muon momentum scale"),
    "ElectronScale": ("Electron energy", "HISTO", None, 40, None, None, "electron energy scale"),
}
for _d in config.TAU_DMS:
    SYST_INFO[f"TauES_DM{_d}"] = ("Tau", "HISTO", "Tau energy scale", 40, None, None, f"#tau_{{h}} energy scale DM{_d} (3% prior)")
    SYST_INFO[f"TauTrigger_DM{_d}"] = ("Tau", "HISTO", "Tau trigger", None, None, None, f"di-#tau trigger DM{_d}")


def combine_members(book, edges_of):
    for name in list(book.h):
        if not name.startswith("__members__"):
            continue
        _, _, region, sample, syst = name.split("__")
        members = book.h.pop(name)[:-1]
        nominal = book.get(f"{region}__{sample}", len(edges_of[region]) - 1)
        up, dn = an.combine_theory(syst, nominal[0], members)
        book.add(f"{region}__{sample}__{syst}Up", np.stack([up, nominal[1]]))
        book.add(f"{region}__{sample}__{syst}Down", np.stack([dn, nominal[1]]))


# ------------------------------------------------------------------------------ tau_h tau_h
def tautau(book, registry, yields):
    regs = config.REGIONS
    REGION_SET.update({r: "nominal" for r in regs})
    var = config.FIT_VARIABLE
    edges = np.asarray(config.FIT_BINS)
    nb = len(edges) - 1
    src = config.FIT_DIR / "fitinputs" / f"{config.JOB}.root"
    meta3 = json.loads((str(src) + ".meta.json").read_text() if False else Path(str(src) + ".meta.json").read_text())
    with uproot.open(src) as f:
        for k in f.keys():
            k = k.split(";")[0]
            parts = k.split("__")
            if len(parts) < 2 or parts[1] not in ("Data", "Fakes"):
                continue
            h = f[k]
            book.add(k, np.stack([h.values(), h.variances()]))
    for name in meta3["fake_systs"]:
        info = meta3["osss_nps"].get(name) or meta3["closure_nps"].get(name)
        r = info["region"]
        if name.endswith("_c0_lo"):
            continue
        registry.syst(name, ["Fakes"], [r], "Fakes", title=name.replace("_tautau", ""))
    # the paper's 30% on the W+jets / ttbar part of the tau_h tau_h fakes: jet-faked leading taus of those samples in the AR
    data = analysis.load_data()
    reg_d = analysis.regions(data)
    cat_d = analysis.categories(data, "data")
    frac = np.zeros(len(regs))
    for key in ("WJets", "TTTo2L2Nu", "TTToSemiLeptonic", "ST_tW_top", "ST_tW_antitop"):
        if not analysis.has_ntuple(key):
            continue
        d, _ = analysis.load(key)
        if not len(d):
            continue
        d = dict(d)
        r = analysis.regions(d, is_mc=False)
        w = analysis.subtraction_weights(key, analysis.weights(d, key))
        cat = analysis.categories(d, key)
        for k in range(len(regs)):
            m = r["AR"] & (d["t1_genflav"] == 0) & (cat == k)
            frac[k] += float(w[m].sum())
    for k, region in enumerate(regs):
        n_ar = float((reg_d["AR"] & (cat_d == k)).sum())
        f = frac[k] / max(n_ar, 1)
        nom = book.get(f"{region}__Fakes", nb)
        for d_, sgn in (("Up", 1), ("Down", -1)):
            book.add(f"{region}__Fakes__FakeWTT_tautau{d_}", np.stack([nom[0] * (1 + sgn * config.LTAU_FF_WTT_SYST * f), nom[1]]))
        yields.setdefault("tautau_wtt_fraction", {})[region] = f
    registry.syst("FakeWTT_tautau", ["Fakes"], regs, "Fakes", title="W/t#bar{t} part of the #tau_{h}#tau_{h} fakes (30%)")
    # simulation
    relaxed = {}
    for key in analysis.available_mc():
        d, _ = analysis.load(key)
        if not len(d):
            continue
        d = dict(d)
        d["_cat"] = analysis.categories(d, key)
        r = analysis.regions(d, is_mc=True)
        w = an.weights(d, key, "tautau")
        comps = an.mc_components(key, "tautau")
        keys = an.dm_key(d, "tautau")
        masks = {region: r["SR"] & (d["_cat"] == k) for k, region in enumerate(regs)}

        def kin_fn(syst, d_, d=d, key=key):
            kin = an.kinematics(d, key, "tautau", syst, d_)
            rs = analysis.regions(d, kin, is_mc=True)
            cs = analysis.categories(d, key, kin, syst, d_)
            return kin, {region: rs["SR"] & (cs == k) for k, region in enumerate(regs)}

        fill_mc(book, regs, "tautau", key, d, w, comps, keys, masks, kin_fn, None, var, BINS, registry)
        for name, cm in comps:
            if name in analysis.SMOOTHED_SAMPLES:
                for k, region in enumerate(regs):
                    for kk in np.unique(keys):
                        m = r["SR_relaxed"] & cm & (d["_cat"] == k) & (keys == kk)
                        relaxed[(tname(name, kk), region)] = relaxed.get((tname(name, kk), region), 0) + h1(d[var][m], w[m], edges)
        print(f"  tautau {key}: " + ", ".join(f"{n} {sum(book.get(f'{rg}__{tname(n, kk)}', nb)[0].sum() for rg in regs for kk in np.unique(keys)):.1f}" for n, _ in comps), flush=True)
    combine_members(book, BINS)
    # W+jets smoothing (docs CLAUDE.md pitfall 2), per DM template
    for (t, region), shape in relaxed.items():
        nom = book.h.get(f"{region}__{t}")
        if nom is None or shape[0].sum() <= 0 or nom[0].sum() <= 0:
            continue
        norm_shape = np.maximum(shape[0], 0) / np.maximum(shape[0], 0).sum()
        for name in [n for n in book.h if n == f"{region}__{t}" or n.startswith(f"{region}__{t}__")]:
            tot = book.h[name][0].sum()
            book.h[name] = np.stack([norm_shape * tot, (norm_shape * tot) ** 2 * (shape[1].sum() / shape[0].sum() ** 2)])
    wj = [t for t in registry.samples if t.startswith("WJets_")]
    tot = sum(book.get(f"{rg}__{t}", nb) for rg in regs for t in wj)
    if tot[0].sum() > 0:
        u = float(np.sqrt(tot[1].sum()) / tot[0].sum())
        registry.syst("MCStatNorm_WJets_tautau", wj, regs, "Background normalisation", kind="OVERALL", up=u, down=-u, title="W+jets MC stat. (norm.)")
        yields["wjets_smoothed"] = {"sr_yield": float(tot[0].sum()), "sr_stat_rel": u, "samples": wj}


# ------------------------------------------------------------------------------ lepton channels
REGION_SET = {}          # region -> 'nominal' (the regions of the measurement) | 'ptsplit' (cross-check only)


def dm_regions(channel):
    """[(region, decay mode or None, pT(tau_h) window, region set)] of a lepton channel.

    'nominal': one region per tau_h decay mode (docs/10 section 8) -- these are the regions of the fit.
    'ptsplit': the same regions split at config.LTAU_PT_SPLIT. They are always filled but never fitted
    together with the nominal ones (they hold the same events); the cross-check job `ztautau_ptsplit`
    (step 5 --region-set ptsplit) fits them with a separate tau_h ID scale factor below 40 GeV and so
    tests the one assumption the tau_h tau_h / (l tau_h)^2 lever rests on (review/REVIEW_v4.md finding 4).
    """
    if not config.LTAU_DM_REGIONS:
        return [(f"{channel}_SR", None, (0.0, np.inf), "nominal")]
    out = [(f"{channel}_SR_dm{dm}", dm, (0.0, np.inf), "nominal") for dm in config.TAU_DMS]
    for dm in config.TAU_DMS:
        out.append((f"{channel}_SRlo_dm{dm}", dm, (0.0, config.LTAU_PT_SPLIT), "ptsplit"))
        out.append((f"{channel}_SRhi_dm{dm}", dm, (config.LTAU_PT_SPLIT, np.inf), "ptsplit"))
    return out


def ltau_mask(d, kin, dm, win):
    """Events of one (decay mode, pT(tau_h) window) region; `kin` carries the varied tau_h pT."""
    m = np.ones(len(d["run"]), bool) if dm is None else (d["t_dm"] == dm)
    if win[0] > 0.0 or np.isfinite(win[1]):
        pt = (kin if kin is not None else d)["t_pt"]
        m = m & (pt >= win[0]) & (pt < win[1])
    return m


def ltau(book, registry, yields, channel, ctrl, plots):
    regs = dm_regions(channel)
    reg_names = [r for r, *_ in regs]
    REGION_SET.update({r: rset for r, _, _, rset in regs})
    edges = np.asarray(BINS[reg_names[0]]); nb = len(edges) - 1
    var = "m_tt"
    ff = json.loads((config.DATA_DIR_V4 / f"fakes_{channel}.json").read_text())
    tables = fk.from_json(ff["tables"])
    fr, osss, rw = ff["fractions"], ff["osss"], ff["w_mt"]
    data = an.load_data(channel)
    reg = an.regions(data, channel)
    ones = np.ones(len(data["run"]))
    ctrl_vars = CONTROL[channel]
    ctrl_region = f"{channel}_SR"
    ctrl[ctrl_region] = {v: Book() for v in ctrl_vars}

    wf, wf2 = fk.fake_weights(data, reg["AR"], tables, fr, osss, rw, with_err=True)
    fake_vars = {"osss": f"FakeOSSS_{channel}", "wmt": f"FakeWmT_{channel}", "frac": f"FakeFrac_{channel}", "tt": f"FakeTT_{channel}"}
    wvar = {(v, sgn): fk.fake_weights(data, reg["AR"], tables, fr, osss, rw, variation=v, direction=sgn) for v in fake_vars for sgn in (1, -1)}
    for region, dm, win, _ in regs:
        sel = ltau_mask(data, None, dm, win)
        sr, ar = reg["SR"] & sel, reg["AR"] & sel
        book.add(f"{region}__Data", h1(data[var][sr], ones[sr], edges))
        book.add(f"{region}__Fakes", h1(data[var][ar], wf[ar], edges, wf2[ar]))
        for v, name in fake_vars.items():
            for d_, sgn in (("Up", 1), ("Down", -1)):
                book.add(f"{region}__Fakes__{name}{d_}", h1(data[var][ar], wvar[(v, sgn)][ar], edges))
    for v, name in fake_vars.items():
        registry.syst(name, ["Fakes"], reg_names, "Fakes", title={"osss": "multijet FF OS/SS", "wmt": "W FF m_{T} extrapolation", "frac": f"AR fractions (W #pm {100 * config.LTAU_FF_FRAC_SYST:.0f}%)", "tt": "t#bar{t} FF (30%)"}[v] + f" {channel}")
    for v, (ed, _) in ctrl_vars.items():
        ctrl[ctrl_region][v].add("Data", h1(data[v][reg["SR"]], ones[reg["SR"]], ed))
        ctrl[ctrl_region][v].add("Fakes", h1(data[v][reg["AR"]], wf[reg["AR"]], ed, wf2[reg["AR"]]))
    subtracted = {r: np.zeros((2, nb)) for r in reg_names}
    for key in an.available_mc(channel):
        d, _ = an.load(key, channel)
        d = dict(d)
        r = an.regions(d, channel, is_mc=True)
        w = an.weights(d, key, channel)
        comps = an.mc_components(key, channel)
        keys = an.dm_key(d, channel)
        cand = an.fit_candidates(d, channel, is_mc=True)

        masks = {region: r["SR"] & ltau_mask(d, None, dm, win) for region, dm, win, _ in regs}

        def kin_fn(syst, d_, d=d, key=key, cand=cand):
            kin = an.kinematics(d, key, channel, syst, d_, fit_only=cand)
            sr = an.regions(d, channel, kin, is_mc=True)["SR"]
            return kin, {region: sr & ltau_mask(d, kin, dm, win) for region, dm, win, _ in regs}

        ctrl_book = {ctrl_region: ctrl[ctrl_region]}
        fill_mc(book, reg_names, channel, key, d, w, comps, keys, masks, kin_fn, None, var, BINS, registry, None, None)
        # control plots (channel-inclusive) and the genuine / lepton-faked tau_h of the AR, subtracted from the fake estimate
        for name, cm in comps:
            m_sr = r["SR"] & cm
            for v, (ed, _) in ctrl_vars.items():
                ctrl[ctrl_region][v].add(name, h1(d[v][m_sr], w[m_sr], ed))
            m_all = r["AR"] & cm
            wm = fk.fake_weights(d, m_all, tables, fr, osss, rw) * w
            wv_all = {(v, sgn): fk.fake_weights(d, m_all, tables, fr, osss, rw, variation=v, direction=sgn) * w for v in fake_vars for sgn in (1, -1)}
            for region, dm, win, _ in regs:
                m = m_all & ltau_mask(d, None, dm, win)
                neg = -h1(d[var][m], wm[m], edges) * np.array([[1], [0]])
                book.add(f"{region}__Fakes", neg)
                subtracted[region] += -neg
                for v, name_ in fake_vars.items():
                    for d_, sgn in (("Up", 1), ("Down", -1)):
                        book.add(f"{region}__Fakes__{name_}{d_}", -h1(d[var][m], wv_all[(v, sgn)][m], edges) * np.array([[1], [0]]))
            for v, (ed, _) in ctrl_vars.items():
                ctrl[ctrl_region][v].add("Fakes", -h1(d[v][m_all], wm[m_all], ed) * np.array([[1], [0]]))
        print(f"  {channel} {key}: " + ", ".join(f"{n} {sum(book.get(f'{rg}__{tname(n, kk)}', nb)[0].sum() for rg in reg_names for kk in np.unique(keys)):.1f}" for n, _ in comps), flush=True)
    combine_members(book, BINS)
    # subtraction uncertainty (10% of the subtracted simulation) and the same-sign non-closure, per region
    delta = ff["closure_ss"]["delta"]
    for region in reg_names:
        nom = book.get(f"{region}__Fakes", nb)
        for d_, sgn in (("Up", 1), ("Down", -1)):
            book.add(f"{region}__Fakes__FakeMCSub_{channel}{d_}", np.stack([nom[0] - sgn * 0.1 * subtracted[region][0], nom[1]]))
            book.add(f"{region}__Fakes__FakeClosure_{channel}{d_}", np.stack([nom[0] * (1 + sgn * delta), nom[1]]))
    registry.syst(f"FakeMCSub_{channel}", ["Fakes"], reg_names, "Fakes", title=f"genuine-#tau_{{h}} subtraction {channel}")
    registry.syst(f"FakeClosure_{channel}", ["Fakes"], reg_names, "Fakes", title=f"FF same-sign non-closure {channel}")
    yields[f"{channel}_closure"] = ff["closure_ss"]
    if plots:
        control_plots(ctrl_region, ctrl[ctrl_region], ctrl_vars, f"{channel} signal region (all decay modes)")


def emu(book, registry, yields, ctrl, plots):
    channel = "emu"
    ff = json.loads((config.DATA_DIR_V4 / "fakes_emu.json").read_text())
    dr_edges = np.asarray(ff["osss"]["dr_edges"])
    ratio, rel = np.asarray(ff["osss"]["ratio"]), np.asarray(ff["osss"]["rel_unc"])
    data = an.load_data(channel)
    reg = an.regions(data, channel)
    ones = np.ones(len(data["run"]))
    ctrl_vars = CONTROL["emu"]
    regs = {"emu_SR": ("SR", "SS"), "emu_CRtt": ("CRtt", "CRtt_SS")}
    REGION_SET.update({r: "nominal" for r in regs})
    var = "m_tt"

    def c_of(d, m, shift=0):
        ib = np.clip(np.searchsorted(dr_edges, d["dr"][m], side="right") - 1, 0, len(dr_edges) - 2)
        return ratio[ib] * (1 + shift * rel[ib])

    for region, (sr, ss) in regs.items():
        edges = np.asarray(BINS[region])
        ctrl[region] = {v: Book() for v in ctrl_vars}
        book.add(f"{region}__Data", h1(data[var][reg[sr]], ones[reg[sr]], edges))
        c = c_of(data, reg[ss])
        book.add(f"{region}__Fakes", h1(data[var][reg[ss]], c, edges))
        for d_, sgn in (("Up", 1), ("Down", -1)):
            book.add(f"{region}__Fakes__QCDOSSS_emu{d_}", h1(data[var][reg[ss]], c_of(data, reg[ss], sgn), edges))
        for v, (ed, _) in ctrl_vars.items():
            ctrl[region][v].add("Data", h1(data[v][reg[sr]], ones[reg[sr]], ed))
            ctrl[region][v].add("Fakes", h1(data[v][reg[ss]], c, ed))
    registry.syst("QCDOSSS_emu", ["Fakes"], list(regs), "Fakes", title="multijet OS/SS (e#mu)")
    for key in an.available_mc(channel):
        d, _ = an.load(key, channel)
        d = dict(d)
        r = an.regions(d, channel, is_mc=True)
        w = an.weights(d, key, channel)
        comps = an.mc_components(key, channel)
        keys = an.dm_key(d, channel)
        cand = an.fit_candidates(d, channel, is_mc=True)

        def kin_fn(syst, d_, d=d, key=key, cand=cand):
            kin = an.kinematics(d, key, channel, syst, d_, fit_only=cand)
            rs = an.regions(d, channel, kin, is_mc=True)
            return kin, {region: rs[sr] for region, (sr, _) in regs.items()}

        fill_mc(book, list(regs), channel, key, d, w, comps, keys, {region: r[sr] for region, (sr, _) in regs.items()}, kin_fn, None, var, BINS, registry, ctrl, ctrl_vars)
        for region, (sr, ss) in regs.items():
            edges = np.asarray(BINS[region])
            for name, cm in comps:
                m = r[ss] & cm
                wm = c_of(d, m) * w[m]
                neg = -h1(d[var][m], wm, edges) * np.array([[1], [0]])
                book.add(f"{region}__Fakes", neg)
                for d_, sgn in (("Up", 1), ("Down", -1)):
                    book.add(f"{region}__Fakes__QCDOSSS_emu{d_}", -h1(d[var][m], c_of(d, m, sgn) * w[m], edges) * np.array([[1], [0]]))
                for v, (ed, _) in ctrl_vars.items():
                    ctrl[region][v].add("Fakes", -h1(d[v][m], wm, ed) * np.array([[1], [0]]))
        print(f"  emu {key}: " + ", ".join(f"{n} {sum(book.get(f'{rg}__{tname(n, kk)}', len(BINS[rg]) - 1)[0].sum() for rg in regs for kk in np.unique(keys)):.1f}" for n, _ in comps), flush=True)
    combine_members(book, BINS)
    if plots:
        for region in regs:
            control_plots(region, ctrl[region], ctrl_vars, {"emu_SR": "e#mu signal region", "emu_CRtt": "e#mu t#bar{t} control region"}[region].replace("#", "\\"))


def control_plots(region, books, ctrl_vars, title):
    for v, book_ in books.items():
        ed, xl = ctrl_vars[v]
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
        variable = len(set(np.round(np.diff(ed), 6))) > 1 and v not in ("t_dm",)
        plotting.stack_plot(config.PLOT_DIR_V4 / f"step4_{region}_{v}.png", ed, (dv, np.sqrt(dv)), stack, xl, title=title, density=variable)
        if v == "m_tt":
            plotting.stack_plot(config.PLOT_DIR_V4 / f"step4_{region}_{v}_log.png", ed, (dv, np.sqrt(dv)), stack, xl, title=title, logy=True)


# ------------------------------------------------------------------------------ completeness
def complete_variations(names, registry, edges_of):
    """Give every (region, template, HISTO systematic) the registry promises both of its variations.

    A kinematic variation (JES, MET_Unclustered, TauES_DM*) is filled from the events that pass the region
    *under that variation*. For a template with a handful of simulated events one direction can move all of
    them out, and then no histogram was booked at all -- the file promised a systematic it did not contain.
    TRExFitter reads that as "no variation" (with HistoChecks NOCRASH), but anything stricter, including the
    combination, fails on it.

    The missing side is filled by mirroring the side that exists about the nominal (the usual one-sided
    symmetrisation, clipped at zero), which is both complete and stable: taking the variation literally
    would turn the migration of one or two simulated events into a 100% uncertainty on that template.
    If neither side exists the nominal is copied, i.e. the variation genuinely does nothing here.
    """
    added = []
    for name, s in registry.systs.items():
        if s["type"] != "HISTO":
            continue
        for r in s["regions"]:
            for smp in s["samples"]:
                nom = names.get(f"{r}__{smp}")
                if nom is None:
                    continue
                up, dn = names.get(f"{r}__{smp}__{name}Up"), names.get(f"{r}__{smp}__{name}Down")
                if up is not None and dn is not None:
                    continue
                nv = nom.view().value
                for missing, present in (("Up", dn), ("Down", up)):
                    if names.get(f"{r}__{smp}__{name}{missing}") is not None:
                        continue
                    h = hist.Hist(hist.axis.Variable(np.asarray(edges_of[r], dtype=float)), storage=hist.storage.Weight())
                    h.view().value[...] = nv if present is None else np.maximum(2 * nv - present.view().value, 0.0)
                    h.view().variance[...] = nom.view().variance
                    names[f"{r}__{smp}__{name}{missing}"] = h
                    added.append(f"{r}__{smp}__{name}{missing}")
    return added


# ------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--channels", nargs="*", default=config.CHANNELS)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    config.DATA_DIR_V4.mkdir(parents=True, exist_ok=True)
    config.PLOT_DIR_V4.mkdir(parents=True, exist_ok=True)
    book, registry, yields, ctrl = Book(), Registry(), {}, {}
    if "tautau" in args.channels:
        tautau(book, registry, yields)
    for ch in ("mutau", "etau"):
        if ch in args.channels:
            ltau(book, registry, yields, ch, ctrl, not args.no_plots)
    if "emu" in args.channels:
        emu(book, registry, yields, ctrl, not args.no_plots)
    # region-wide bookkeeping: every template that appears in a region gets a (tiny) nominal there if a
    # systematic variation exists for it (TRExFitter needs the nominal)
    regions = [r for r in config.REGIONS_V4 if any(n.startswith(r + "__") for n in book.h)]
    names = {}
    for name, vv in book.h.items():
        names[name] = to_hist(vv, BINS[name.split("__")[0]])
    for t, info in registry.samples.items():
        for r in info["regions"]:
            if f"{r}__{t}" not in names:
                nb = len(BINS[r]) - 1
                names[f"{r}__{t}"] = to_hist(np.stack([np.full(nb, 1e-6), np.zeros(nb)]), BINS[r])
    for r in regions:
        if f"{r}__Fakes" not in names:
            nb = len(BINS[r]) - 1
            names[f"{r}__Fakes"] = to_hist(np.stack([np.full(nb, 1e-6), np.zeros(nb)]), BINS[r])
    added = complete_variations(names, registry, BINS)
    if added:
        print(f"completed {len(added)} systematic variation(s) that the region emptied: " + ", ".join(added[:6])
              + (" ..." if len(added) > 6 else ""))
    # normalisation NPs (registry only; templates are the nominal ones)
    mc_templates = sorted(registry.samples)
    registry.syst("Lumi", mc_templates, regions, "Luminosity", kind="OVERALL", up=config.LUMI_REL_UNC, down=-config.LUMI_REL_UNC, title="Luminosity")
    registry.syst("XS_DYll", [t for t in mc_templates if t.startswith(("DYee_", "DYmumu_"))], regions, "Background normalisation", kind="OVERALL", up=0.05, down=-0.05, title="#sigma(Z#rightarrowee/#mu#mu)")
    for base in (config.V4_SIGNAL_OUT, "DYlowmass", "WJets", "SingleTop", "WW", "WZ", "ZZ"):
        ts = [t for t in mc_templates if t.startswith(base + "_tDM")]
        if ts:
            u = XSEC_UNC[base]
            registry.syst(f"XS_{base}", ts, regions, "Background normalisation", kind="OVERALL", up=u, down=-u, title=f"#sigma({TITLES[base]})")
    pog_sf = {str(dm): corrections.pog()["id_vsjet_dm"][config.TAU_WP][str(dm)] for dm in config.TAU_DMS}
    out = config.FIT_DIR_V4 / "fitinputs" / f"{config.JOB_V4}.root"
    meta = {"regions": regions, "region_sets": {r: REGION_SET.get(r, "nominal") for r in regions},
            "bins": {r: list(BINS[r]) for r in regions}, "lumi_pb": config.LUMI_PB, "channels": args.channels,
            "samples": registry.samples, "systs": registry.systs, "tau_id_sf_pog": pog_sf, "tes_prior": config.TES_PRIOR_V4,
            "signal": config.V4_SIGNAL, "signal_out": config.V4_SIGNAL_OUT, "sideband_mtt_min": config.SIDEBAND_REGION_MTT_MIN,
            # for a combination (docs/11-combination-inputs.md): the signal is split by the decay mode of the
            # genuine tau_h legs, so `mu_Z` (and any per-channel reference factor) goes on *these* templates,
            # not on a sample called "DYtautau"; and there is no acceptance uncertainty to add outside the fit.
            "signal_samples": sorted(t for t, i in registry.samples.items() if i["is_signal"]),
            "acceptance_in_fit": True,
            "acceptance_note": ("A x epsilon is profiled inside the fit: every theory variation member is "
                                "renormalised to the same sigma(60 < m_LHE < 120) (docs/CONVENTIONS.md section 3), "
                                "so PDF / QCDScale / PS_ISR / PS_FSR vary only the acceptance, per region. There is "
                                "no A_unc block here on purpose -- adding Acc_* parameters on top would double count."),
            "signal_prediction": an.signal_prediction("tautau" if "tautau" in args.channels else "mutau"),
            "yields_extra": yields, "titles": TITLES}
    rep = trexhist.write_fitinputs(out, names, meta=meta)
    print(f"fit inputs -> {out}: {rep['n_hists']} histograms, clipped {len(rep['clipped'])}")
    # yields per region and base sample
    y = {}
    for r in regions:
        nb = len(BINS[r]) - 1
        y[r] = {"Data": float(book.get(f"{r}__Data", nb)[0].sum()), "Fakes": float(book.get(f"{r}__Fakes", nb)[0].sum())}
        for t, info in registry.samples.items():
            if r in info["regions"]:
                y[r][info["base"]] = y[r].get(info["base"], 0.0) + float(book.get(f"{r}__{t}", nb)[0].sum())
        print(f"  {r}: " + ", ".join(f"{k} {v:.0f}" for k, v in y[r].items()))
    (config.DATA_DIR_V4 / "yields_v4.json").write_text(json.dumps({"regions": y, "extra": yields, "n_templates": len(registry.samples),
                                                                     "n_systs": len(registry.systs)}, indent=1, default=float))


if __name__ == "__main__":
    main()
