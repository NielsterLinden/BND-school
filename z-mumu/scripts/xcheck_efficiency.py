#!/usr/bin/env python
"""Cross-check -- trigger efficiency with trigger-object matching, and the
background under the ID / isolation tag-and-probe.

Two concerns from reviewing step 2, both testable on the *parent* NanoAOD
(`config.PARENT_DIR`), which unlike the skim still has the `TrigObj_*` branches.

1. **Trigger.** Step 2 measures IsoMu24 || IsoTkMu24 relative to HLT_IsoMu27,
   HLT_Mu27, HLT_Mu50 and HLT_Mu45_eta2p1. Those share the L1 seed and the L3
   muon reconstruction with the measured paths, so the ratio is close to one by
   construction and mostly measures the online isolation. Here the per-muon
   efficiency is measured with tag-and-probe instead:

       tag    selected muon, pT > 26 GeV, matched to an IsoMu24/IsoTkMu24 object
       probe  selected muon, pT > 20 GeV, opposite sign, 70 < m < 110 GeV
       pass   probe matched as well

   and folded into the event efficiency 1 - (1 - eps1)(1 - eps2) over the
   signal muons, weighting each event by 1/eps (the observed sample is the
   triggered one).

2. **ID / isolation.** Step 2 counts pass and fail probes in 70-110 GeV with no
   background subtraction. Non-Z pairs end up mostly in the *fail* sample, so
   both efficiencies come out low. The pass/fail mass spectra are kept here,
   for opposite- and same-sign pairs and with and without a trigger-matched
   tag, so the background can be fitted.

As a closure test, step 2's own `_trigger_efficiency` and `_tag_and_probe` are
run on the same events; they must reproduce the step 2 numbers, since the
parent contains every event of the skim.

Nothing in steps 1-6 is changed. The output is a comparison.

Outputs
    output/data/xcheck_efficiency.pkl
    output/plots/xcheck_*.png

**Input.** The parent NanoAOD is read from CERN Open Data EOS over xrootd by
default (`--source eos`, file lists from `cernopendata-client` in
`config.FILELIST_DIR`), or from the dCache copy (`--source dcache`). The dCache
files match the catalogue checksums, but its NFS reads stalled repeatedly under
parallel access (14 Sep 2026), so EOS is the robust choice. Each file runs in
its own subprocess via `zmumu.batch`, with its result in
output/data/xcheck_parts/, so an interrupted run resumes.

Run:  .venv/bin/python scripts/xcheck_efficiency.py            # full parent from EOS, ~10 min
      .venv/bin/python scripts/xcheck_efficiency.py --summarise-only
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import awkward as ak
import numpy as np

import step2_efficiency as step2
from zmumu import batch, config, hists, io, objects, stats

BRANCHES = sorted(set(step2.BRANCHES + [
    "FsrPhoton_pt", "FsrPhoton_eta", "FsrPhoton_phi",
    "FsrPhoton_relIso03", "FsrPhoton_dROverEt2", "FsrPhoton_muonIdx",
    "TrigObj_id", "TrigObj_pt", "TrigObj_eta", "TrigObj_phi",
    "TrigObj_l1pt", "TrigObj_filterBits",
]))

N_ETA = len(config.EFF_ETA_BINS) - 1
N_TRIG_CELLS = (len(config.TRIG_PT_BINS) - 1) * N_ETA
N_EFF_CELLS = (len(config.EFF_PT_BINS) - 1) * N_ETA
TP_MASS_AXIS = (120, 60.0, 120.0)


# --------------------------------------------------------------------------
# Booking and merging
# --------------------------------------------------------------------------
def blank():
    out = {"n_files": 0, "n_signal": 0}
    for variant in ("nominal", "l1veto"):
        for era in config.ERAS:
            for kind in ("pass", "tot"):
                out[f"trig_{variant}_{era}_{kind}"] = hists.hist2d(
                    config.TRIG_PT_BINS, config.EFF_ETA_BINS,
                    r"probe $p_T$ [GeV]", r"probe $\eta$")
    for tag in ("any", "matched"):
        for charge in ("os", "ss"):
            for cat in ("id_pass", "id_fail", "iso_pass", "iso_fail"):
                h = hists.hist2d(config.EFF_PT_BINS, np.linspace(*TP_MASS_AXIS[1:], TP_MASS_AXIS[0] + 1),
                                 r"probe $p_T$ [GeV]", r"$m_{\mu\mu}$ [GeV]")
                out[f"tp_{tag}_{charge}_{cat}"] = h
    # (cell of muon 1, cell of muon 2) counts of the selected signal events.
    out["sig_trig_cells"] = np.zeros(N_TRIG_CELLS * N_TRIG_CELLS)
    # Same, split by how many of the two muons carry a trigger-object match.
    out["sig_trig_cells_both"] = np.zeros(N_TRIG_CELLS * N_TRIG_CELLS)
    out["sig_trig_cells_none"] = np.zeros(N_TRIG_CELLS * N_TRIG_CELLS)
    out["sig_eff_cells"] = np.zeros(N_EFF_CELLS * N_EFF_CELLS)
    # Step 2's own functions, run on the parent (closure test).
    out["step2"] = step2.blank()
    return out


# --------------------------------------------------------------------------
# Per-file processing
# --------------------------------------------------------------------------
def trigger_match(ev, veto_l1_band: bool = False) -> ak.Array:
    """(events, muons) mask: muon within dR of an analysis-path HLT muon object."""
    obj = ((ev.TrigObj_id == 13)
           & ((ev.TrigObj_filterBits & config.TRIG_OBJ_BITS) > 0)
           & (ev.TrigObj_pt > config.TRIG_OBJ_PT_MIN))
    if veto_l1_band:
        lo, hi = config.TRIG_L1_VETO_BAND
        obj = obj & ~((ev.TrigObj_l1pt >= lo) & (ev.TrigObj_l1pt < hi))
    mu = ak.zip({"eta": ev.Muon_eta, "phi": ev.Muon_phi})
    to = ak.zip({"eta": ev.TrigObj_eta[obj], "phi": ev.TrigObj_phi[obj]})
    m, t = ak.unzip(ak.cartesian([mu, to], axis=1, nested=True))
    dr = objects.delta_r(m.eta, m.phi, t.eta, t.phi)
    return ak.any(dr < config.TRIG_MATCH_DR, axis=2)


def _pairs(ev):
    """Both orientations of every muon pair: (tag idx, probe idx, mass, opposite sign)."""
    idx = ak.local_index(ev.Muon_pt, axis=1)
    p = ak.combinations(idx, 2, fields=["a", "b"])
    tag = ak.concatenate([p.a, p.b], axis=1)
    probe = ak.concatenate([p.b, p.a], axis=1)
    pt, eta, phi, mass = ev.Muon_pt, ev.Muon_eta, ev.Muon_phi, ev.Muon_mass
    px1, py1, pz1, e1 = objects.p4(pt[tag], eta[tag], phi[tag], mass[tag])
    px2, py2, pz2, e2 = objects.p4(pt[probe], eta[probe], phi[probe], mass[probe])
    pair_mass = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
    opposite = (ev.Muon_charge[tag] * ev.Muon_charge[probe]) < 0
    return tag, probe, pair_mass, opposite


def _flat(x):
    return ak.to_numpy(ak.flatten(x))


def _clip(pt, edges):
    return np.clip(pt, edges[0], edges[-1] - 1e-3)


def _cell(pt, eta, pt_edges):
    ipt = np.digitize(_clip(pt, pt_edges), pt_edges) - 1
    ieta = np.clip(np.digitize(eta, config.EFF_ETA_BINS) - 1, 0, N_ETA - 1)
    return ipt * N_ETA + ieta


def _trigger_tnp(ev, good, era, out):
    tag, probe, mass, opposite = _pairs(ev)
    window = opposite & (mass > config.TP_MASS_LO) & (mass < config.TP_MASS_HI)
    for variant, veto in (("nominal", False), ("l1veto", True)):
        matched = trigger_match(ev, veto_l1_band=veto)
        is_tag = good & (ev.Muon_pt > config.TAG_PT_MIN) & matched
        valid = is_tag[tag] & good[probe] & window
        ppt = _clip(_flat(ev.Muon_pt[probe][valid]), config.TRIG_PT_BINS)
        peta = _flat(ev.Muon_eta[probe][valid])
        passed = _flat(matched[probe][valid])
        out[f"trig_{variant}_{era}_tot"].fill(ppt, peta)
        out[f"trig_{variant}_{era}_pass"].fill(ppt[passed], peta[passed])
    return tag, probe, mass, opposite


def _idiso_tnp(ev, good, pairs, out):
    tag, probe, mass, opposite = pairs
    loose = objects.loose_muon_mask(ev)
    matched = trigger_match(ev)
    pass_id = (ev[config.MU_ID_BRANCH] & (abs(ev.Muon_dxy) < config.MU_DXY_MAX)
               & (abs(ev.Muon_dz) < config.MU_DZ_MAX))
    pass_iso = ev.Muon_pfRelIso04_all < config.MU_ISO_MAX
    in_axis = (mass > TP_MASS_AXIS[1]) & (mass < TP_MASS_AXIS[2])
    for tagdef in ("any", "matched"):
        is_tag = good & (ev.Muon_pt > config.TAG_PT_MIN)
        if tagdef == "matched":
            is_tag = is_tag & matched
        base = is_tag[tag] & loose[probe] & in_axis
        for charge, sign in (("os", opposite), ("ss", ~opposite)):
            valid = base & sign
            ppt = _clip(_flat(ev.Muon_pt[probe][valid]), config.EFF_PT_BINS)
            m = _flat(mass[valid])
            pid = _flat(pass_id[probe][valid])
            piso = _flat(pass_iso[probe][valid])
            key = f"tp_{tagdef}_{charge}"
            out[f"{key}_id_pass"].fill(ppt[pid], m[pid])
            out[f"{key}_id_fail"].fill(ppt[~pid], m[~pid])
            out[f"{key}_iso_pass"].fill(ppt[pid & piso], m[pid & piso])
            out[f"{key}_iso_fail"].fill(ppt[pid & ~piso], m[pid & ~piso])


def _signal_cells(ev, good, matched, out):
    """Step 1's signal selection (after trigger), recorded as muon-cell pairs."""
    keep = ak.sum(good, axis=1) == config.N_MUONS_REQUIRED
    ev, good, matched = ev[keep], good[keep], matched[keep]
    if len(ev) == 0:
        return
    i1, i2 = objects.leading_two(ev, good)
    pt1 = ak.to_numpy(objects.take(ev.Muon_pt, i1))
    q = ak.to_numpy(objects.take(ev.Muon_charge, i1) * objects.take(ev.Muon_charge, i2))
    sel = (pt1 > config.MU_PT_LEAD) & (q < 0)
    ev, matched, i1, i2 = ev[sel], matched[sel], i1[sel], i2[sel]
    if len(ev) == 0:
        return
    m = ak.to_numpy(objects.dimuon_mass(ev, i1, i2, with_fsr=config.FSR_ENABLED))
    win = (m > config.MASS_LO) & (m < config.MASS_HI)
    ev, matched, i1, i2 = ev[win], matched[win], i1[win], i2[win]
    if len(ev) == 0:
        return
    m1 = ak.to_numpy(objects.take(matched, i1))
    m2 = ak.to_numpy(objects.take(matched, i2))
    pt1, eta1 = ak.to_numpy(objects.take(ev.Muon_pt, i1)), ak.to_numpy(objects.take(ev.Muon_eta, i1))
    pt2, eta2 = ak.to_numpy(objects.take(ev.Muon_pt, i2)), ak.to_numpy(objects.take(ev.Muon_eta, i2))
    out["n_signal"] += len(ev)
    for key, edges, ncell in (("sig_trig_cells", config.TRIG_PT_BINS, N_TRIG_CELLS),
                              ("sig_eff_cells", config.EFF_PT_BINS, N_EFF_CELLS)):
        pair = _cell(pt1, eta1, edges) * ncell + _cell(pt2, eta2, edges)
        out[key] += np.bincount(pair, minlength=ncell * ncell)
        if key == "sig_trig_cells":
            out[key + "_both"] += np.bincount(pair[m1 & m2], minlength=ncell * ncell)
            out[key + "_none"] += np.bincount(pair[~m1 & ~m2], minlength=ncell * ncell)


def process_file(source, era: str) -> dict:
    grl = io.load_grl()
    out = blank()
    out["n_files"] = 1
    n_all = n_skim = 0
    for events, _ in io.read_chunks(source, BRANCHES, chunk_size="200 MB"):
        n_all += len(events)
        events = events[ak.num(events.Muon_pt, axis=1) >= 2]          # the skim
        n_skim += len(events)
        events = events[io.lumi_mask(events.run, events.luminosityBlock, grl)]
        events = events[io.pass_met_filters(events) & (events.PV_npvsGood >= 1)]
        if len(events) == 0:
            continue
        step2._trigger_efficiency(events, out["step2"])
        fired = events[io.trigger_or(events, config.TRIGGERS)]
        if len(fired) == 0:
            continue
        step2._tag_and_probe(fired, out["step2"])
        good = objects.good_muon_mask(fired)
        pairs = _trigger_tnp(fired, good, era, out)
        _idiso_tnp(fired, good, pairs, out)
        _signal_cells(fired, good, trigger_match(fired), out)
    out["file_counts"] = [(era, Path(source).name, n_all, n_skim)]
    return out


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
def _eff(pass_h, tot_h):
    e, lo, hi = stats.clopper_pearson(pass_h.values(), tot_h.values())
    return e, 0.5 * (lo + hi)


def fold_trigger(eff_map, cells):
    """Event trigger efficiency from a per-muon map, over the signal cell pairs.

    Returns (mean over the observed events, efficiency-corrected mean). The
    second is N_obs / sum(1/eps_event), which is what the cross section needs.
    """
    e = eff_map.ravel()
    ev_eff = 1.0 - np.outer(1.0 - e, 1.0 - e).ravel()
    n = cells
    ok = n > 0
    observed = float((n[ok] * ev_eff[ok]).sum() / n[ok].sum())
    corrected = float(n[ok].sum() / (n[ok] / ev_eff[ok]).sum())
    return observed, corrected


def fit_fail_signal(fail_m, pass_m, centres, window):
    """Signal fail probes in `window`: fail = s * pass-shape + exponential.

    Template fit in the full 60-120 GeV spectrum. Returns (signal count in the
    window, background count in the window).
    """
    from scipy.optimize import curve_fit

    tmpl = pass_m / pass_m.sum() if pass_m.sum() > 0 else np.zeros_like(pass_m)

    def model(x, s, b, c):
        return s * tmpl + b * np.exp(-c * (x - 60.0))

    err = np.sqrt(np.maximum(fail_m, 1.0))
    p0 = [max(fail_m.sum() * 0.6, 1.0), max(fail_m[:10].mean(), 1.0), 0.02]
    popt, _ = curve_fit(model, centres, fail_m, p0=p0, sigma=err,
                        bounds=([0, 0, -0.2], [np.inf, np.inf, 0.5]), maxfev=20000)
    s, b, c = popt
    sig = float((s * tmpl)[window].sum())
    bkg = float((b * np.exp(-c * (centres - 60.0)))[window].sum())
    return sig, bkg


def summarise(total, karel_step1, karel_step2, karel_step5):
    res = {}
    print(f"\n[xcheck] files {total['n_files']}, signal events (step 1 selection) "
          f"{total['n_signal']:,}")

    # ---- closure: step 2's own methods on the parent ---------------------
    s2 = total["step2"]
    trig_ref = s2["trig_pass"] / s2["trig_total"]
    id_ref = s2["id_pass_2d"].values().sum() / s2["id_tot_2d"].values().sum()
    iso_ref = s2["iso_pass_2d"].values().sum() / s2["iso_tot_2d"].values().sum()
    res.update(closure_trig=trig_ref, closure_id=id_ref, closure_iso=iso_ref,
               closure_pairs=s2["n_pairs"])
    print(f"  closure  step-2 trigger  {trig_ref:.4f}   (step 2: {karel_step2['trig_event']:.4f})")
    print(f"  closure  step-2 ID       {id_ref:.4f}   (step 2: {karel_step2['id_inclusive']:.4f})")
    print(f"  closure  step-2 ISO      {iso_ref:.4f}   (step 2: {karel_step2['iso_inclusive']:.4f})")
    print(f"  closure  T&P pairs       {s2['n_pairs']:,}   (step 2: {karel_step2['n_pairs']:,})")

    # ---- trigger: per-muon maps and folded event efficiency --------------
    print("\n  Trigger (IsoMu24 || IsoTkMu24), tag-and-probe with trigger-object matching")
    shape = (len(config.TRIG_PT_BINS) - 1, N_ETA)
    for variant in ("nominal", "l1veto"):
        p = sum((total[f"trig_{variant}_{e}_pass"] for e in config.ERAS[1:]),
                total[f"trig_{variant}_{config.ERAS[0]}_pass"])
        t = sum((total[f"trig_{variant}_{e}_tot"] for e in config.ERAS[1:]),
                total[f"trig_{variant}_{config.ERAS[0]}_tot"])
        e_map, e_err = _eff(p, t)
        obs, corr = fold_trigger(e_map, total["sig_trig_cells"])
        plateau = t.values()[4:].sum()
        plateau_eff = p.values()[4:].sum() / plateau
        res[f"trig_{variant}_map"] = e_map
        res[f"trig_{variant}_map_err"] = e_err
        res[f"trig_{variant}_event_observed"] = obs
        res[f"trig_{variant}_event"] = corr
        res[f"trig_{variant}_per_muon_plateau"] = plateau_eff
        print(f"    {variant:8s} per-muon eff, pT > 26 GeV : {plateau_eff:.4f}"
              f"   event eff: {corr:.4f}  (mean over observed events {obs:.4f})")
        for era in config.ERAS:
            pe, te = total[f"trig_{variant}_{era}_pass"], total[f"trig_{variant}_{era}_tot"]
            if te.values().sum() == 0:
                continue
            em, _ = _eff(pe, te)
            print(f"             {era}: per-muon plateau {pe.values()[4:].sum()/te.values()[4:].sum():.4f}"
                  f"   event {fold_trigger(em, total['sig_trig_cells'])[1]:.4f}")
    # Model check: with independent per-muon efficiencies, the fraction of
    # triggered signal events in which *both* muons are matched is
    # eps1*eps2 / (1 - (1-eps1)(1-eps2)). Compare with what is observed.
    e = res["trig_nominal_map"].ravel()
    ev_eff = 1.0 - np.outer(1.0 - e, 1.0 - e).ravel()
    both_pred = np.divide(np.outer(e, e).ravel(), ev_eff, out=np.zeros_like(ev_eff), where=ev_eff > 0)
    n = total["sig_trig_cells"]
    res["trig_both_matched_observed"] = float(total["sig_trig_cells_both"].sum() / n.sum())
    res["trig_both_matched_predicted"] = float((n * both_pred).sum() / n.sum())
    res["trig_none_matched_observed"] = float(total["sig_trig_cells_none"].sum() / n.sum())
    print(f"    model check: both muons matched, observed {res['trig_both_matched_observed']:.4f}"
          f" vs predicted {res['trig_both_matched_predicted']:.4f};"
          f" no muon matched in {100*res['trig_none_matched_observed']:.3f}% of triggered signal events")
    trig_new = res["trig_nominal_event"]
    trig_syst = abs(res["trig_nominal_event"] - res["trig_l1veto_event"])
    res["trig_event_syst"] = trig_syst

    # ---- ID / isolation: background in the probe samples -----------------
    print("\n  ID / isolation tag-and-probe, background in 70-110 GeV")
    centres = np.linspace(60.25, 119.75, 120)
    window = (centres > config.TP_MASS_LO) & (centres < config.TP_MASS_HI)
    npt = len(config.EFF_PT_BINS) - 1
    for tagdef in ("any", "matched"):
        for eff in ("id", "iso"):
            P = total[f"tp_{tagdef}_os_{eff}_pass"].values()
            F = total[f"tp_{tagdef}_os_{eff}_fail"].values()
            Pss = total[f"tp_{tagdef}_ss_{eff}_pass"].values()
            Fss = total[f"tp_{tagdef}_ss_{eff}_fail"].values()
            raw, ss_sub, fit = np.zeros(npt), np.zeros(npt), np.zeros(npt)
            fit_bkg_frac = np.zeros(npt)
            for i in range(npt):
                p, f = P[i][window].sum(), F[i][window].sum()
                pss, fss = Pss[i][window].sum(), Fss[i][window].sum()
                raw[i] = p / (p + f)
                ss_sub[i] = (p - pss) / ((p - pss) + (f - fss))
                f_sig, f_bkg = fit_fail_signal(F[i], P[i], centres, window)
                fit[i] = (p - pss) / ((p - pss) + f_sig)
                fit_bkg_frac[i] = f_bkg / f if f else 0.0
            p_all, f_all = P[:, window].sum(), F[:, window].sum()
            key = f"{eff}_{tagdef}"
            res[f"{key}_raw_pt"], res[f"{key}_ss_pt"], res[f"{key}_fit_pt"] = raw, ss_sub, fit
            res[f"{key}_fail_bkg_frac_pt"] = fit_bkg_frac
            # Inclusive numbers, weighting the per-bin results by the pass+fail counts.
            w = (P[:, window].sum(axis=1) + F[:, window].sum(axis=1))
            res[f"{key}_raw"] = p_all / (p_all + f_all)
            res[f"{key}_ss"] = float((ss_sub * w).sum() / w.sum())
            res[f"{key}_fit"] = float((fit * w).sum() / w.sum())
            print(f"    {eff.upper():3s} tag={tagdef:7s} raw {res[f'{key}_raw']:.4f}"
                  f"   SS-subtracted {res[f'{key}_ss']:.4f}   template fit {res[f'{key}_fit']:.4f}")

    # ---- impact on the cross section ------------------------------------
    # Step 2's (pT, eta) maps are rescaled by the per-pT-bin ratio of the
    # background-corrected to the raw efficiency, then folded over the signal
    # muons with a per-event 1/eps weight. The two background estimates bracket
    # the truth (same-sign pairs undercount the opposite-sign background; the
    # template fit also absorbs any extra width of the failing-probe peak), so
    # their midpoint is used, with half the difference as a systematic.
    def scaled_map(eff, tagdef, method):
        ratio = res[f"{eff}_{tagdef}_{method}_pt"] / res[f"{eff}_any_raw_pt"]
        return np.clip(karel_step2[f"{eff}_map"] * ratio[:, None], 0.0, 1.0)

    def step5_fold(m, kin):
        w = kin.values()
        use = m > 0
        return float((np.where(use, w, 0) * m).sum() / np.where(use, w, 0).sum())

    def event_fold(id_map, iso_map):
        per_mu = (id_map * iso_map).ravel()
        ev = np.outer(per_mu, per_mu).ravel()
        n = total["sig_eff_cells"]
        ok = (n > 0) & (ev > 0)
        return float(n[ok].sum() / (n[ok] / ev[ok]).sum())

    reco2 = config.MU_RECO_EFF ** 2
    k1, k2 = karel_step1["h_mu1_kin"], karel_step1["h_mu2_kin"]
    idm, isom = karel_step2["id_map"], karel_step2["iso_map"]
    idiso_step5 = step5_fold(idm, k1) * step5_fold(isom, k1) * step5_fold(idm, k2) * step5_fold(isom, k2)
    idiso_event = event_fold(idm, isom)
    idiso_ss = event_fold(scaled_map("id", "matched", "ss"), scaled_map("iso", "matched", "ss"))
    idiso_fit = event_fold(scaled_map("id", "matched", "fit"), scaled_map("iso", "matched", "fit"))
    idiso_mid = 0.5 * (idiso_ss + idiso_fit)

    n_signal = karel_step5["n_signal"]
    eff0 = karel_step5["eff_total"]
    rows = [
        ("step 5 as committed", eff0, config.LUMI_PB),
        ("+ per-event 1/eps folding", reco2 * idiso_event * karel_step2["trig_event"], config.LUMI_PB),
        ("+ trigger: TrigObj tag-and-probe", reco2 * idiso_event * trig_new, config.LUMI_PB),
        ("  (ID/iso bkg: same-sign subtracted)", reco2 * idiso_ss * trig_new, config.LUMI_PB),
        ("  (ID/iso bkg: template fit)", reco2 * idiso_fit * trig_new, config.LUMI_PB),
        ("+ ID/iso bkg: midpoint, matched tag", reco2 * idiso_mid * trig_new, config.LUMI_PB),
        ("+ luminosity with normtag", reco2 * idiso_mid * trig_new, config.LUMI_PB_NORMTAG),
    ]
    print("\n  Effect on sigma_fid (cumulative; rows in brackets are the two bracketing variants)")
    res["sigma_steps"] = []
    sigma0 = n_signal / (eff0 * config.LUMI_PB)
    for label, eff_tot, lumi in rows:
        sigma = n_signal / (eff_tot * lumi)
        res["sigma_steps"].append((label, eff_tot, sigma, lumi))
        print(f"    {label:38s} eff_total {eff_tot:.4f}  L {lumi:9.1f}  sigma_fid {sigma:7.1f} pb"
              f"   ({100*(sigma/sigma0-1):+.2f}%)")

    # Uncertainties that change with the revision (relative, on sigma).
    res["rel_unc_revised"] = {
        "eff: trigger (L1-band veto variation)": trig_syst / trig_new,
        "eff: ID/iso background method": 0.5 * abs(idiso_fit - idiso_ss) / idiso_mid,
    }
    print("    new systematics: " + ", ".join(f"{k} {100*v:.2f}%" for k, v in res["rel_unc_revised"].items()))
    res.update(idiso_step5=idiso_step5, idiso_event=idiso_event, idiso_ss=idiso_ss,
               idiso_fit=idiso_fit, idiso_mid=idiso_mid)
    return res


def plot(total, res, karel_step2):
    import matplotlib.pyplot as plt
    import mplhep as hep

    # Trigger efficiency vs pT, eta-inclusive, both eras and the step-2 value.
    fig, ax = plt.subplots(figsize=(9, 7))
    edges = np.asarray(config.TRIG_PT_BINS)
    centres, half = 0.5 * (edges[1:] + edges[:-1]), 0.5 * np.diff(edges)
    for (era, variant), colour in zip([(e, "nominal") for e in config.ERAS], ["#1f77b4", "#d62728"]):
        p = total[f"trig_{variant}_{era}_pass"].values().sum(axis=1)
        t = total[f"trig_{variant}_{era}_tot"].values().sum(axis=1)
        e, lo, hi = stats.clopper_pearson(p, t)
        ax.errorbar(centres, e, yerr=[lo, hi], xerr=half, fmt="o", color=colour,
                    markersize=4, label=f"T&P, TrigObj matched, {era.split('__')[0]}")
    ax.axhline(karel_step2["trig_event"], color="black", linestyle="--",
               label=f"step 2 event eff. (reference method) {karel_step2['trig_event']:.4f}")
    ax.axhline(res["trig_nominal_event"], color="#2ca02c", linestyle="-.",
               label=f"T&P folded event eff. {res['trig_nominal_event']:.4f}")
    ax.set_xscale("log")
    ax.set_xlabel(r"probe muon $p_T$ [GeV]")
    ax.set_ylabel("IsoMu24 || IsoTkMu24 efficiency")
    ax.set_ylim(0.0, 1.08)
    ax.legend(fontsize=11, loc="lower right")
    hists._decorate(ax)
    hists._title(ax, "Per-muon trigger efficiency (points) vs event efficiency (lines)")
    hists.save_fig(fig, "xcheck_trigger_tnp_vs_pt.png")

    # Per-muon efficiency vs eta on the plateau.
    fig, ax = plt.subplots(figsize=(9, 7))
    eedges = np.asarray(config.EFF_ETA_BINS)
    ec, eh = 0.5 * (eedges[1:] + eedges[:-1]), 0.5 * np.diff(eedges)
    for era, colour in zip(config.ERAS, ["#1f77b4", "#d62728"]):
        p = total[f"trig_nominal_{era}_pass"].values()[4:].sum(axis=0)
        t = total[f"trig_nominal_{era}_tot"].values()[4:].sum(axis=0)
        e, lo, hi = stats.clopper_pearson(p, t)
        ax.errorbar(ec, e, yerr=[lo, hi], xerr=eh, fmt="o", color=colour, markersize=4,
                    label=era.split("__")[0])
    ax.set_xlabel(r"probe muon $\eta$")
    ax.set_ylabel("per-muon efficiency, $p_T$ > 26 GeV")
    ax.set_ylim(0.70, 1.02)
    ax.legend(fontsize=12)
    hists._decorate(ax)
    hists._title(ax, "IsoMu24 || IsoTkMu24, tag-and-probe with TrigObj matching")
    hists.save_fig(fig, "xcheck_trigger_tnp_vs_eta.png")

    # Fail-probe mass spectra, lowest pT bin, with same-sign overlay.
    for eff in ("id", "iso"):
        fig, ax = plt.subplots(figsize=(9, 7))
        mc = np.linspace(60.25, 119.75, 120)
        for tagdef, charge, style in (("any", "os", "-"), ("any", "ss", ":"), ("matched", "os", "--")):
            v = total[f"tp_{tagdef}_{charge}_{eff}_fail"].values()[0]
            ax.step(mc, v, where="mid", linestyle=style,
                    label=f"fail probes, {charge.upper()}, tag {tagdef}")
        ax.set_yscale("log")
        ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
        ax.set_ylabel("Pairs / 0.5 GeV")
        ax.legend(fontsize=12)
        hists._decorate(ax)
        hists._title(ax, f"{eff.upper()} fail probes, 20 < pT < 25 GeV")
        hists.save_fig(fig, f"xcheck_{eff}_fail_mass_lowpt.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--source", choices=["eos", "dcache"], default="eos")
    ap.add_argument("--filelists", type=Path, default=config.FILELIST_DIR,
                    help="directory with the cernopendata-client *.sizes.tsv lists")
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--karel-data", type=Path, default=config.DATA_DIR,
                    help="directory with the step1/step2 pickles to compare against")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:
        era = args.key.split("__", 2)[0] + "__" + args.key.split("__", 2)[1]
        batch.write_part(process_file(args.one_file, era), args.out)
        return

    parts = config.DATA_DIR / "xcheck_parts"
    if not args.summarise_only:
        tasks = []
        for era in config.ERAS:
            if args.source == "dcache":
                paths = sorted((config.PARENT_DIR / era).glob("*.root"))
            else:
                recid = era.split("__")[1]
                tsv = next(args.filelists.glob(f"SingleMuon_*_{recid}.sizes.tsv"))
                paths = sorted(batch.read_catalogue(tsv).values(), key=lambda v: Path(v[0]).name)
                paths = [v[0] for v in paths]
            tasks += [(src, f"{era}__{Path(src).name}") for src in paths]
        tasks = tasks[: args.max_files]
        print(f"[xcheck] parent NanoAOD from {args.source}: {len(tasks)} files", flush=True)
        notes = batch.run_files(tasks, Path(__file__), parts, args.workers, retries=2,
                                log=lambda line: print(line, flush=True))
        total = batch.load_parts(parts, blank())
        total["notes"] = notes
        hists.save(total, "xcheck_efficiency.pkl")
    total = hists.load("xcheck_efficiency.pkl")

    import pickle
    with open(args.karel_data / "step1_selection.pkl", "rb") as fh:
        k1 = pickle.load(fh)
    with open(args.karel_data / "step2_efficiency.pkl", "rb") as fh:
        k2 = pickle.load(fh)["result"]
    with open(args.karel_data / "step5_crosssection.pkl", "rb") as fh:
        k5 = pickle.load(fh)
    res = summarise(total, k1, k2, k5)
    plot(total, res, k2)
    hists.save(res, "xcheck_efficiency_result.pkl")


if __name__ == "__main__":
    main()
