"""Histograms of the v2 analysis for every sample, region, variable and systematic variation.

One pass per skim file (scripts/v2_4_histograms.py) fills numpy histograms keyed by
"<fit sample>|<region>|<variable>|<variation>" (and "...|w2" for the sum of squared weights),
merged over files with zmumu.batch. The fit templates are the `mass_fit` histograms of the
`SR` and `CRemu` regions; everything else is for the data/MC comparison plots.

MC enters the SR/SS only with *prompt* muons (genPartFlav 1 or 15); non-prompt muons are
estimated from data with the fake factor (zmumu.fakes). In the e-mu regions (CRemu, SSemu)
*all* MC events are kept: nothing data-driven replaces the non-prompt electrons there
(W+jets jet -> e, Z -> mumu + conversion, Z -> tautau tau_h -> e, top b -> e), which are 9% of
the OS e-mu region -- see review/README.md F1. DY_NLO is split by LHE flavour into DYmumu / DYee /
DYtautau. DY_powheg (raw generator weights) provides the SigModel template; for that
comparison the DYmumu nominal is also filled with the powheg generator window
50 < m_LHE < 120 GeV ("lhe50120") together with the fiducial sums needed to normalise it.
Jets are lepton-cleaned (dR > 0.4 to the selected leptons).

Variations filled for MC: nominal, PileupUp/Down, L1PrefiringUp/Down, MuonID/MuonIso/
MuonTriggerUp/Down (scale factors), MuonScaleUp/Down and MuonResUp/Down (re-selection with
varied muon pT) and, for the signal, the PDF / QCD-scale / parton-shower weight members
("theory_pdf", "theory_scale", "theory_ps": arrays of shape (n_members, n_bins) for mass_fit).
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config, regions

VARIABLES = {
    "mass_fit": (60, 60.0, 120.0), "mass_fine": (120, 60.0, 120.0),
    "pt1": (36, 20.0, 200.0), "pt2": (36, 20.0, 200.0), "eta1": (24, -2.4, 2.4), "eta2": (24, -2.4, 2.4),
    "phi1": (16, -np.pi, np.pi), "zpt": (40, 0.0, 200.0), "zy": (24, -2.4, 2.4), "npv": (50, 0.0, 50.0),
    "met": (30, 0.0, 150.0), "puppimet": (30, 0.0, 150.0), "njet": (8, 0.0, 8.0), "iso1": (30, 0.0, 0.15), "iso2": (30, 0.0, 0.15),
    "nfsr": (4, 0.0, 4.0),
}
CREMU_MASS_FIT = (12, 60.0, 120.0)
EMU_VARS = ["mass_fit", "mass_fine", "pt1", "pt_el", "met", "puppimet", "njet", "npv"]
REGION_VARS = {"SR": list(VARIABLES), "SS": ["mass_fit", "mass_fine", "pt1", "pt2"], "CRemu": EMU_VARS, "SSemu": EMU_VARS}
EMU_REGIONS = ("CRemu", "SSemu")
LHE_WINDOW = (50.0, 120.0)        # generator cut of the powheg ZToMuMu_M-50To120 sample
WEIGHT_VARIATIONS = ["PileupUp", "PileupDown", "L1PrefiringUp", "L1PrefiringDown",
                     "MuonIDUp", "MuonIDDown", "MuonIsoUp", "MuonIsoDown", "MuonTriggerUp", "MuonTriggerDown"]
PT_VARIATIONS = ["MuonScaleUp", "MuonScaleDown", "MuonResUp", "MuonResDown"]
THEORY = {"pdf": 103, "scale": 9, "ps": 4}
FLAVOUR_SAMPLE = {13: "DYmumu", 11: "DYee", 15: "DYtautau"}


def edges(region, var):
    if var == "mass_fit" and region in EMU_REGIONS:
        n, lo, hi = CREMU_MASS_FIT
    elif var == "pt_el":
        n, lo, hi = 36, 20.0, 200.0
    else:
        n, lo, hi = VARIABLES[var]
    return np.linspace(lo, hi, n + 1)


def _fill(out, key, x, w, e):
    h = np.histogram(x, bins=e, weights=w)[0]
    out[key] = out.get(key, 0) + h
    out[key + "|w2"] = out.get(key + "|w2", 0) + np.histogram(x, bins=e, weights=w * w)[0]


def _region_values(region, d, ev_idx_arrays):
    """Variable -> per-event array for a region selection dict `d`."""
    npv, met, puppimet, njet, nfsr = ev_idx_arrays
    vals = {"mass_fit": d["mass"], "mass_fine": d["mass"], "pt1": d.get("pt1"), "npv": npv[d["idx"]],
            "met": met[d["idx"]], "puppimet": puppimet[d["idx"]], "njet": njet[d["idx"]]}
    if region in ("SR", "SS"):
        vals.update(pt2=d["pt2"], eta1=d["eta1"], eta2=d["eta2"], phi1=d["phi1"], zpt=d["zpt"], zy=d["zy"],
                    iso1=d["iso1"], iso2=d["iso2"], nfsr=nfsr[d["idx"]])
    if region in EMU_REGIONS:
        vals["pt_el"] = d["pt_el"]
    return vals


def fill_chunk(ev, out, key, is_mc, weighter=None, sf=None, calib=None, split_flavour=False,
               theory=False, fit_sample=None, variations=True):
    """Fill all histograms of one chunk of skim events."""
    fired = regions.fired(ev)
    clean = regions.event_clean(ev)
    keep = fired & clean
    ev = ev[keep]
    if len(ev) == 0:
        return
    n = len(ev)
    npv = np.asarray(ev.PV_npvsGood, dtype=float)
    met = np.asarray(ev.MET_pt, dtype=float)
    puppimet = np.asarray(ev.PuppiMET_pt, dtype=float)
    njet = regions.clean_jet_count(ev, regions.muon_masks(ev))
    nfsr = np.asarray(ak.num(ev.FsrPhoton_pt, axis=1), dtype=float)
    extras = (npv, met, puppimet, njet, nfsr)

    if is_mc:
        pieces = weighter.pieces(ev)
        base = pieces["gen"] * pieces["pu"] * pieces["pref"]
        wvar = {
            "PileupUp": pieces["gen"] * pieces["pu_up"] * pieces["pref"],
            "PileupDown": pieces["gen"] * pieces["pu_down"] * pieces["pref"],
            "L1PrefiringUp": pieces["gen"] * pieces["pu"] * pieces["pref_up"],
            "L1PrefiringDown": pieces["gen"] * pieces["pu"] * pieces["pref_down"],
        }
        flav = np.asarray(ev.gen_lhe_flavour) if split_flavour else None
        in_window = fid = None
        if "gen_mll_lhe" in ak.fields(ev):
            m_lhe = np.asarray(ev.gen_mll_lhe, dtype=float)
            in_window = (m_lhe > LHE_WINDOW[0]) & (m_lhe < LHE_WINDOW[1])
            fid = np.asarray(ev.gen_fid_dressed, dtype=bool)
    else:
        base = np.asarray(ev.skim_prescale, dtype=float)
        wvar, flav = {}, None

    def sample_of(idx):
        if flav is None:
            return np.full(len(idx), fit_sample, dtype=object)
        return np.array([FLAVOUR_SAMPLE.get(int(f), "DYother") for f in flav[idx]], dtype=object)

    pt_variants = [("nominal", None)]
    if is_mc and calib is not None:
        pt_variants = [("nominal", calib.scaled_pt(ev, "nominal"))]
        if variations:
            pt_variants += [(v, calib.scaled_pt(ev, v)) for v in PT_VARIATIONS]

    for pt_name, pt_arr in pt_variants:
        e = ev if pt_arr is None else regions.with_muon_pt(ev, pt_arr)
        masks = regions.muon_masks(e)
        matched = regions.trigger_match(e)
        ones = np.ones(n, dtype=bool)
        sels = regions.dimuon_regions(e, masks, matched, ones)
        sels["CRemu"] = regions.emu_region(e, masks, matched, ones)
        sels["SSemu"] = regions.emu_region(e, masks, matched, ones, same_sign=True)
        for region, d in sels.items():
            if d is None or len(d["idx"]) == 0:
                continue
            idx = d["idx"]
            w0 = base[idx]
            pp = np.ones(len(idx), dtype=bool)
            if is_mc:
                if region not in EMU_REGIONS:       # e-mu regions keep the non-prompt MC (no replacement exists)
                    pp = regions.is_prompt(d["flav1"]) & regions.is_prompt(d["flav2"])
                w0 = np.where(pp, w0, 0.0)
                if sf is not None:
                    sfw = (sf.single_muon_weights(d["pt1"], d["eta1"]) if region in EMU_REGIONS
                           else sf.event_weights(d["pt1"], d["eta1"], d["pt2"], d["eta2"]))
                else:
                    sfw = {"nominal": np.ones(len(idx))}
            else:
                sfw = {"nominal": np.ones(len(idx))}
            samples_here = sample_of(idx)
            values = _region_values(region, d, extras)
            for smp in np.unique(samples_here):
                m = samples_here == smp
                if pt_name != "nominal":       # only the fit templates for the pT variations
                    for var in ("mass_fit", "mass_fine"):
                        _fill(out, f"{smp}|{region}|{var}|{pt_name}", values[var][m], (w0 * sfw["nominal"])[m], edges(region, var))
                    continue
                for var in REGION_VARS[region]:
                    _fill(out, f"{smp}|{region}|{var}|nominal", values[var][m], (w0 * sfw["nominal"])[m], edges(region, var))
                if is_mc and smp == "DYmumu" and region == "SR" and in_window is not None:
                    # nominal restricted to the powheg generator window, and the fiducial sums (all / in window)
                    wn = (w0 * sfw["nominal"])[m]
                    for var in ("mass_fit", "mass_fine"):
                        _fill(out, f"{smp}|{region}|{var}|lhe50120", values[var][m], wn * in_window[idx][m], edges(region, var))
                    f_ = fid[idx][m]
                    _fill(out, f"{smp}|{region}|fidsum|all", np.full(int(f_.sum()), 0.5), wn[f_], np.array([0.0, 1.0]))
                    _fill(out, f"{smp}|{region}|fidsum|lhe50120", np.full(int(f_.sum()), 0.5), (wn * in_window[idx][m])[f_], np.array([0.0, 1.0]))
                if not (is_mc and variations):
                    continue
                for var in ("mass_fit", "mass_fine"):
                    e_ = edges(region, var)
                    for vname, warr in wvar.items():
                        _fill(out, f"{smp}|{region}|{var}|{vname}", values[var][m],
                              (np.where(pp, warr[idx], 0.0) * sfw["nominal"])[m], e_)
                    for vname, s_ in sfw.items():
                        if vname == "nominal":
                            continue
                        _fill(out, f"{smp}|{region}|{var}|{vname}", values[var][m], (w0 * s_)[m], e_)
                if theory and region in ("SR", "CRemu") and smp == "DYmumu":
                    e_ = edges(region, "mass_fit")
                    bins = np.clip(np.digitize(values["mass_fit"][m], e_) - 1, 0, len(e_) - 2)
                    wn = (w0 * sfw["nominal"])[m]
                    for kind in THEORY:
                        v = weighter.theory_weights(ev[idx][m], kind, renorm=True)
                        if v is None:
                            continue
                        h = np.zeros((v.shape[1], len(e_) - 1))
                        for k in range(v.shape[1]):
                            h[k] = np.bincount(bins, weights=wn * v[:, k], minlength=len(e_) - 1)
                        key_ = f"{smp}|{region}|mass_fit|theory_{kind}"
                        out[key_] = out.get(key_, 0) + h
