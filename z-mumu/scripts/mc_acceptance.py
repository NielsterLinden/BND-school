#!/usr/bin/env python
"""Acceptance of the fiducial volume, and closure of the data-driven method, in DY simulation.

The measurement in steps 1-6 is fiducial and uses no simulation. Two things
need simulation, and the DY samples the z-ee subgroup uses (recids 35669 NLO
and 35671 LO) can be streamed from CERN Open Data EOS over xrootd:

**1. The acceptance, for the combination.** The combination notebook wants a
total cross section per flavour, sigma(Z/gamma* -> mu mu, m > 50 GeV) =
6077.22/3 pb, so this channel needs

    A = sigma_fid(predicted) / (6077.22 / 3)  =  3 * sumw(fiducial) / sumw(all events)

with the fiducial volume at generator level: exactly two GenDressedLepton muons
(no tau ancestor) with pT > 20 GeV and |eta| < 2.4, leading pT > 26 GeV,
opposite charge, 60 < m < 120 GeV, in an event whose LHE record is Z -> mu mu.
The alternative denominators (all LHE mu mu events; LHE mu mu with
60 < m < 120 GeV) are recorded as well, because the convention has to be agreed
with the combination group. PDF (NNPDF3.1 NNLO Hessian, LHA 325300: square
root of the summed eigenvector shifts), alpha_s (members 101/102, 0.116/0.120,
scaled to +/-0.0015) and scale (7-point envelope) uncertainties are propagated
through the ratio.

**2. Closure of the data-driven method.** Running the analysis pieces on
simulation, where the truth is known, tests each one without any background:

    trigger     truth: fraction of offline-selected events that fired
                vs step 2's reference-trigger method
                vs tag-and-probe with trigger-object matching (xcheck_efficiency.py)
    ID, iso     truth: gen-matched muons  vs  step 2's tag-and-probe
    "A = 1"     reco events with two *prompt* muons passing only the kinematic
                cuts, divided by the generator-level fiducial events. Anything
                other than 1 (after the reconstruction efficiency) is migration
                across the volume boundary that the fiducial measurement ignores.

Outputs
    output/data/mc_acceptance_<sample>.pkl
    output/plots/mc_*.png

Run:  .venv/bin/python scripts/mc_acceptance.py --sample nlo [--max-files N]
      .venv/bin/python scripts/mc_acceptance.py --summarise-only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import awkward as ak
import numpy as np

import step2_efficiency as step2
import xcheck_efficiency as xcheck
from zmumu import batch, config, hists, io, objects

GEN_BRANCHES = [
    "genWeight",
    "LHEPart_pt", "LHEPart_eta", "LHEPart_phi", "LHEPart_mass",
    "LHEPart_pdgId", "LHEPart_status",
    "GenDressedLepton_pt", "GenDressedLepton_eta", "GenDressedLepton_phi",
    "GenDressedLepton_mass", "GenDressedLepton_pdgId", "GenDressedLepton_hasTauAnc",
    "LHEPdfWeight", "LHEScaleWeight", "Muon_genPartFlav",
]
BRANCHES = sorted(set(xcheck.BRANCHES + GEN_BRANCHES))
CATEGORIES = ("all", "lhe_mumu", "lhe_mumu_60_120", "fid")


# --------------------------------------------------------------------------
# Per-file processing
# --------------------------------------------------------------------------
def blank():
    out = {"n_files": 0, "n_events": 0, "sumw_runs": 0.0, "titles": []}
    for c in CATEGORIES:
        out[f"sumw_{c}"] = 0.0
        out[f"sumw2_{c}"] = 0.0
    for key in ("offline", "offline_fired", "offline_fired_fid", "kin", "kin_fid",
                "gen_mu", "gen_mu_reco", "gen_mu_acc", "gen_mu_id", "gen_mu_id_iso"):
        out[f"w_{key}"] = 0.0
    # Step 2's methods, and the trigger-object tag-and-probe, on simulation.
    out["step2"] = step2.blank()
    for variant in ("nominal", "l1veto"):
        for kind in ("pass", "tot"):
            out[f"trig_{variant}_mc_{kind}"] = hists.hist2d(
                config.TRIG_PT_BINS, config.EFF_ETA_BINS, r"probe $p_T$ [GeV]", r"probe $\eta$")
    out["sig_trig_cells"] = np.zeros(xcheck.N_TRIG_CELLS ** 2)
    out["sig_trig_cells_both"] = np.zeros(xcheck.N_TRIG_CELLS ** 2)
    out["sig_trig_cells_none"] = np.zeros(xcheck.N_TRIG_CELLS ** 2)
    out["sig_eff_cells"] = np.zeros(xcheck.N_EFF_CELLS ** 2)
    out["n_signal"] = 0
    out["h_m_dressed"] = hists.mass_hist()
    out["h_m_lhe"] = hists.mass_hist()
    return out


def _leading_pair(pt, eta, phi, mass, pdg, mask):
    """Two leading objects passing `mask`, padded with None where absent."""
    rec = ak.zip({"pt": pt, "eta": eta, "phi": phi, "mass": mass, "pdg": pdg})[mask]
    rec = rec[ak.argsort(rec.pt, axis=1, ascending=False)]
    n = ak.to_numpy(ak.num(rec, axis=1))
    rec = ak.pad_none(rec, 2, axis=1)
    a, b = rec[:, 0], rec[:, 1]
    f = lambda x: ak.to_numpy(ak.fill_none(x, 0.0))
    px1, py1, pz1, e1 = objects.p4(f(a.pt), f(a.eta), f(a.phi), f(a.mass))
    px2, py2, pz2, e2 = objects.p4(f(b.pt), f(b.eta), f(b.phi), f(b.mass))
    m = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)
    return n, f(a.pt), m, f(a.pdg) * f(b.pdg)


def _dimuon_selection(ev, mask) -> np.ndarray:
    """Step 1's pairing on reco muons passing `mask`: exactly two, leading pT,
    opposite sign, FSR-recovered mass in the window. Per-event boolean."""
    n = ak.to_numpy(ak.sum(mask, axis=1))
    keep = np.zeros(len(ev), dtype=bool)
    idx = np.nonzero(n == config.N_MUONS_REQUIRED)[0]
    if idx.size == 0:
        return keep
    sub, m = ev[idx], mask[idx]
    i1, i2 = objects.leading_two(sub, m)
    pt1 = ak.to_numpy(objects.take(sub.Muon_pt, i1))
    q = ak.to_numpy(objects.take(sub.Muon_charge, i1) * objects.take(sub.Muon_charge, i2))
    mass = ak.to_numpy(objects.dimuon_mass(sub, i1, i2, with_fsr=config.FSR_ENABLED))
    ok = (pt1 > config.MU_PT_LEAD) & (q < 0) & (mass > config.MASS_LO) & (mass < config.MASS_HI)
    keep[idx[ok]] = True
    return keep


def _weight_sums(out, name, w, sel, variations):
    out[f"sumw_{name}"] += float(w[sel].sum())
    out[f"sumw2_{name}"] += float((w[sel] ** 2).sum())
    for vname, v in variations.items():
        key = f"{vname}_{name}"
        vals = (w[sel, None] * v[sel]).sum(axis=0)
        out[key] = out.get(key, 0.0) + vals


def _regular(arr):
    counts = ak.to_numpy(ak.num(arr, axis=1))
    if counts.size == 0 or counts.min() != counts.max() or counts.max() == 0:
        return None
    return ak.to_numpy(ak.to_regular(arr, axis=1)).astype(float)


def process_file(source) -> dict:
    import uproot

    out = blank()
    out["n_files"] = 1
    with uproot.open(str(source)) as fh:
        out["sumw_runs"] = float(fh["Runs"]["genEventSumw"].array(library="np").sum())
        tree = fh["Events"]
        out["titles"] = [{b: tree[b].title for b in ("LHEPdfWeight", "LHEScaleWeight",
                                                     "GenDressedLepton_pt") if b in tree.keys()}]

    for ev, _ in io.read_chunks(source, BRANCHES, chunk_size="200 MB"):
        out["n_events"] += len(ev)
        w = ak.to_numpy(ev.genWeight).astype(float)
        variations = {}
        for branch, vname in (("LHEPdfWeight", "pdf"), ("LHEScaleWeight", "scale")):
            v = _regular(ev[branch])
            if v is not None:
                variations[vname] = v

        # ---- generator level -------------------------------------------
        lhe_mu = (ev.LHEPart_status == 1) & (abs(ev.LHEPart_pdgId) == 13)
        n_lhe, _, m_lhe, _ = _leading_pair(ev.LHEPart_pt, ev.LHEPart_eta, ev.LHEPart_phi,
                                           ev.LHEPart_mass, ev.LHEPart_pdgId, lhe_mu)
        lhe_mumu = n_lhe == 2
        lhe_6012 = lhe_mumu & (m_lhe > config.MASS_LO) & (m_lhe < config.MASS_HI)

        dressed = ((abs(ev.GenDressedLepton_pdgId) == 13) & ~ev.GenDressedLepton_hasTauAnc
                   & (ev.GenDressedLepton_pt > config.MU_PT_SUBLEAD)
                   & (abs(ev.GenDressedLepton_eta) < config.MU_ETA_MAX))
        n_d, pt_lead, m_d, charge_prod = _leading_pair(
            ev.GenDressedLepton_pt, ev.GenDressedLepton_eta, ev.GenDressedLepton_phi,
            ev.GenDressedLepton_mass, ev.GenDressedLepton_pdgId, dressed)
        fid = (lhe_mumu & (n_d == config.N_MUONS_REQUIRED) & (pt_lead > config.MU_PT_LEAD)
               & (charge_prod < 0) & (m_d > config.MASS_LO) & (m_d < config.MASS_HI))

        everything = np.ones(len(ev), dtype=bool)
        for name, sel in (("all", everything), ("lhe_mumu", lhe_mumu),
                          ("lhe_mumu_60_120", lhe_6012), ("fid", fid)):
            _weight_sums(out, name, w, sel, variations)
        out["h_m_dressed"].fill(m_d[fid], weight=w[fid])
        out["h_m_lhe"].fill(m_lhe[lhe_6012], weight=w[lhe_6012])

        # ---- truth efficiencies of the fiducial muons -------------------
        _truth_muon_efficiencies(ev[fid], w[fid], dressed[fid], out)

        # ---- reco level, as in data (no lumi mask) ----------------------
        clean = ak.to_numpy(io.pass_met_filters(ev) & (ev.PV_npvsGood >= 1))
        ev_c, w_c, fid_c = ev[clean], w[clean], fid[clean]
        good = objects.good_muon_mask(ev_c)
        fired = ak.to_numpy(io.trigger_or(ev_c, config.TRIGGERS))
        offline = _dimuon_selection(ev_c, good)
        out["w_offline"] += float(w_c[offline].sum())
        out["w_offline_fired"] += float(w_c[offline & fired].sum())
        out["w_offline_fired_fid"] += float(w_c[offline & fired & fid_c].sum())

        prompt = ((ev_c.Muon_genPartFlav == 1) & (ev_c.Muon_pt > config.MU_PT_SUBLEAD)
                  & (abs(ev_c.Muon_eta) < config.MU_ETA_MAX))
        kin = _dimuon_selection(ev_c, prompt)
        out["w_kin"] += float(w_c[kin].sum())
        out["w_kin_fid"] += float(w_c[kin & fid_c].sum())

        # ---- the data-driven methods, unweighted as in data -------------
        step2._trigger_efficiency(ev_c, out["step2"])
        ev_f = ev_c[fired]
        if len(ev_f):
            step2._tag_and_probe(ev_f, out["step2"])
            good_f = objects.good_muon_mask(ev_f)
            xcheck._trigger_tnp(ev_f, good_f, "mc", out)
            xcheck._signal_cells(ev_f, good_f, xcheck.trigger_match(ev_f), out)
    return out


def _truth_muon_efficiencies(ev, w, dressed, out):
    """Per-muon reco / ID / iso efficiency for generator-level fiducial muons."""
    if len(ev) == 0:
        return
    gen = ak.zip({"eta": ev.GenDressedLepton_eta[dressed], "phi": ev.GenDressedLepton_phi[dressed]})
    reco = ak.zip({"eta": ev.Muon_eta, "phi": ev.Muon_phi})
    g, r = ak.unzip(ak.cartesian([gen, reco], axis=1, nested=True))
    dr = objects.delta_r(g.eta, g.phi, r.eta, r.phi)
    best = ak.argmin(dr, axis=2, keepdims=True)
    dr_best = ak.fill_none(ak.flatten(dr[best], axis=2), 99.0)
    idx = ak.flatten(best, axis=2)

    def at(branch):
        return ak.fill_none(ev[branch][idx], 0)

    loose = ak.fill_none((ev.Muon_isTracker | ev.Muon_isGlobal)[idx], False)
    matched = (dr_best < config.GEN_MATCH_DR) & loose
    acc = matched & (at("Muon_pt") > config.MU_PT_SUBLEAD) & (abs(at("Muon_eta")) < config.MU_ETA_MAX)
    pass_id = (ak.fill_none(ev[config.MU_ID_BRANCH][idx], False)
               & (abs(at("Muon_dxy")) < config.MU_DXY_MAX) & (abs(at("Muon_dz")) < config.MU_DZ_MAX))
    pass_iso = at("Muon_pfRelIso04_all") < config.MU_ISO_MAX

    def wsum(mask):
        return float((ak.to_numpy(ak.sum(mask, axis=1)) * w).sum())

    ones = ak.ones_like(dr_best, dtype=bool)
    out["w_gen_mu"] += wsum(ones)
    out["w_gen_mu_reco"] += wsum(matched)
    out["w_gen_mu_acc"] += wsum(acc)
    out["w_gen_mu_id"] += wsum(acc & pass_id)
    out["w_gen_mu_id_iso"] += wsum(acc & pass_id & pass_iso)


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
def _ratio_unc(sumw_num, sumw2_num, sumw_den, sumw2_den):
    """Uncertainty of sumw_num/sumw_den when the numerator is a subset of the denominator."""
    a = sumw_num / sumw_den
    var = (sumw2_num * (1 - a) ** 2 + (sumw2_den - sumw2_num) * a ** 2) / sumw_den ** 2
    return a, float(np.sqrt(max(var, 0.0)))


def summarise(t, sample):
    r = {"sample": sample, "label": config.DY_SAMPLES[sample]["label"],
         "n_files": t["n_files"], "n_events": t["n_events"]}
    print(f"\n[mc] {r['label']}: {t['n_files']} files, {t['n_events']:,} events")
    print(f"  sum genWeight {t['sumw_all']:.6g}  vs Runs genEventSumw {t['sumw_runs']:.6g}"
          f"  (ratio {t['sumw_all']/t['sumw_runs']:.6f})")
    r["f_mumu"] = t["sumw_lhe_mumu"] / t["sumw_all"]

    frac, frac_err = _ratio_unc(t["sumw_fid"], t["sumw2_fid"], t["sumw_all"], t["sumw2_all"])
    r["sigma_fid_pred"] = config.DY_XSEC_PB * frac
    r["A"] = 3.0 * frac
    r["A_stat"] = 3.0 * frac_err
    r["A_lhe_mumu"] = t["sumw_fid"] / t["sumw_lhe_mumu"]
    r["A_lhe_mumu_60_120"] = t["sumw_fid"] / t["sumw_lhe_mumu_60_120"]

    def variation_ratios(vname):
        if f"{vname}_fid" not in t:
            return None
        return (t[f"{vname}_fid"] / t[f"{vname}_all"]) / frac

    pdf = variation_ratios("pdf")
    if pdf is not None and len(pdf) >= 101:
        # Hessian eigenvectors: sqrt of the summed squared shifts (not an RMS).
        r["A_pdf_rel"] = float(np.sqrt(np.sum((pdf[1:101] - pdf[0]) ** 2)))
        # alpha_s members are +/-0.002; PDF4LHC recommends +/-0.0015.
        r["A_alphas_rel"] = float(0.75 * 0.5 * abs(pdf[102] - pdf[101])) if len(pdf) >= 103 else 0.0
    else:
        r["A_pdf_rel"], r["A_alphas_rel"] = float("nan"), float("nan")
    scale = variation_ratios("scale")
    if scale is not None and len(scale) == 9:
        # NanoAOD order: [muR,muF] = 0.5,0.5 / 0.5,1 / 0.5,2 / 1,0.5 / 1,1 / 1,2 / 2,0.5 / 2,1 / 2,2
        seven = scale[[0, 1, 3, 4, 5, 7, 8]] / scale[4]
        r["A_scale_rel"] = float(np.max(np.abs(seven - 1.0)))
    else:
        r["A_scale_rel"] = float("nan")

    print(f"  LHE mu mu fraction of all events       : {r['f_mumu']:.5f}   (1/3 = 0.33333)")
    print(f"  predicted sigma_fid                    : {r['sigma_fid_pred']:.1f} pb"
          f"   (sigma_DY = {config.DY_XSEC_PB} pb x fiducial fraction)")
    print(f"  A = sigma_fid / (sigma_DY/3)            : {r['A']:.5f} +/- {r['A_stat']:.5f} (MC stat)"
          f"  +/- {100*r['A_pdf_rel']:.2f}% (PDF)  +/- {100*r['A_alphas_rel']:.2f}% (alpha_s)"
          f"  +/- {100*r['A_scale_rel']:.2f}% (scale)")
    print(f"  A w.r.t. LHE mu mu (m > 50)            : {r['A_lhe_mumu']:.5f}")
    print(f"  A w.r.t. LHE mu mu, 60 < m < 120 GeV   : {r['A_lhe_mumu_60_120']:.5f}")

    # ---- closure ---------------------------------------------------------
    s2 = t["step2"]
    r["trig_truth"] = t["w_offline_fired"] / t["w_offline"]
    r["trig_reference"] = s2["trig_pass"] / s2["trig_total"] if s2["trig_total"] else float("nan")
    tp = {}
    for variant in ("nominal", "l1veto"):
        e_map, _ = xcheck._eff(t[f"trig_{variant}_mc_pass"], t[f"trig_{variant}_mc_tot"])
        tp[variant] = xcheck.fold_trigger(e_map, t["sig_trig_cells"])[1]
    r["trig_tnp"], r["trig_tnp_l1veto"] = tp["nominal"], tp["l1veto"]
    r["reco_truth"] = t["w_gen_mu_reco"] / t["w_gen_mu"]
    r["id_truth"] = t["w_gen_mu_id"] / t["w_gen_mu_acc"]
    r["iso_truth"] = t["w_gen_mu_id_iso"] / t["w_gen_mu_id"]
    r["id_tnp"] = s2["id_pass_2d"].values().sum() / s2["id_tot_2d"].values().sum()
    r["iso_tnp"] = s2["iso_pass_2d"].values().sum() / s2["iso_tot_2d"].values().sum()
    r["C"] = t["w_offline_fired"] / t["sumw_fid"]
    r["purity"] = t["w_offline_fired_fid"] / t["w_offline_fired"]
    r["M_kin"] = t["w_kin"] / t["sumw_fid"]
    r["M_kin_purity"] = t["w_kin_fid"] / t["w_kin"]
    r["migration"] = r["M_kin"] / r["reco_truth"] ** 2

    print("\n  Closure in simulation (truth vs the data-driven method)")
    print(f"    trigger, event level : truth {r['trig_truth']:.4f}   step-2 reference method "
          f"{r['trig_reference']:.4f}   TrigObj tag-and-probe {r['trig_tnp']:.4f}"
          f" (L1-band veto {r['trig_tnp_l1veto']:.4f})")
    print(f"    ID  per muon         : truth {r['id_truth']:.4f}   step-2 tag-and-probe {r['id_tnp']:.4f}")
    print(f"    iso per muon (|ID)   : truth {r['iso_truth']:.4f}   step-2 tag-and-probe {r['iso_tnp']:.4f}")
    print(f"    reco per muon        : truth {r['reco_truth']:.4f}   (data analysis uses {config.MU_RECO_EFF})")
    print(f"    C = N_reco / N_fid   : {r['C']:.4f}   purity (reco events inside the gen volume) {r['purity']:.4f}")
    print(f"    kinematic migration  : N(2 prompt reco muons, kinematic cuts only) / N_fid = {r['M_kin']:.4f};"
          f" / reco_eff^2 = {r['migration']:.4f}")

    # The whole data-driven efficiency chain applied to simulation, with
    # per-muon inclusive numbers: how far is N_reco / eps_method from N_fid?
    reco2 = r["reco_truth"] ** 2
    idiso_tnp = (r["id_tnp"] * r["iso_tnp"]) ** 2
    r["closure_step2"] = r["C"] / (reco2 * idiso_tnp * r["trig_reference"])
    r["closure_trigobj"] = r["C"] / (reco2 * idiso_tnp * r["trig_tnp"])
    r["closure_trigobj_migration"] = r["closure_trigobj"] / r["migration"]
    print(f"    method closure N_reco/(eps_method N_fid): step-2 chain {r['closure_step2']:.4f};"
          f" TrigObj trigger {r['closure_trigobj']:.4f}; and migration-corrected"
          f" {r['closure_trigobj_migration']:.4f}")
    return r


def revised_result(results, karel_step5, xres):
    """Revised sigma_fid, its uncertainty budget, and the combination inputs."""
    nlo = next(r for r in results if r["sample"] == "nlo")
    lo = next((r for r in results if r["sample"] == "lo"), None)
    label, eff, sigma_lumi_eff, lumi = xres["sigma_steps"][-1]
    migration = nlo["migration"]
    sigma = sigma_lumi_eff / migration

    rel = dict(karel_step5["relative_uncertainties"])
    rel.pop("eff: trigger (method bias)")
    rel.update(xres["rel_unc_revised"])
    rel["eff: trigger T&P closure (MC)"] = abs(nlo["trig_tnp"] / nlo["trig_truth"] - 1)
    ref = lo or nlo
    rel["eff: ID/iso T&P closure (MC)"] = 2 * np.hypot(ref["id_tnp"] / ref["id_truth"] - 1,
                                                     ref["iso_tnp"] / ref["iso_truth"] - 1)
    rel["migration (NLO vs LO)"] = abs(nlo["migration"] / lo["migration"] - 1) if lo else 0.0
    stat = rel.pop("statistical (data)")
    syst = float(np.sqrt(sum(v ** 2 for v in rel.values())))

    A = nlo["A"]
    A_rel = float(np.sqrt((nlo["A_stat"] / A) ** 2 + nlo["A_pdf_rel"] ** 2
                          + nlo["A_alphas_rel"] ** 2 + nlo["A_scale_rel"] ** 2))
    out = {
        "sigma_fid": sigma, "stat_rel": stat, "syst_rel": syst, "rel": rel,
        "eff_total": eff, "lumi_pb": lumi, "migration": migration, "n_signal": karel_step5["n_signal"],
        "A": A, "A_rel": A_rel, "A_lo": lo["A"] if lo else None,
        "sigma_tot": sigma / A, "acc_eff": A * eff * migration,
        "sigma_fid_pred_nlo": nlo["sigma_fid_pred"],
    }
    print("\n[mc] Revised fiducial cross section (dressed muons)")
    print(f"  N_signal {out['n_signal']:,.0f}   eff_total {eff:.4f}   L {lumi:.1f} pb^-1"
          f"   migration {migration:.4f}")
    print(f"  sigma_fid = {sigma:.1f} +/- {sigma*stat:.1f} (stat) +/- {sigma*syst:.1f} (syst) pb"
          f"   [committed: {karel_step5['sigma_fid_pb']:.1f} +/- {karel_step5['sigma_syst_pb']:.1f} (syst)]")
    for k, v in sorted(rel.items(), key=lambda kv: -kv[1]):
        print(f"    {k:40s} {100*v:6.3f}%")
    print(f"  NLO prediction {nlo['sigma_fid_pred']:.1f} pb -> measured/predicted {sigma/nlo['sigma_fid_pred']:.3f}")
    print(f"  total: sigma(Z/gamma* -> mu mu, m > 50) = sigma_fid / A = {out['sigma_tot']:.1f} pb"
          f"  (+/- {100*A_rel:.2f}% from A)   target {config.DY_XSEC_PB/3:.1f} pb")
    print(f"  combination inputs: n_obs {karel_step5['n_observed']:,.0f}, n_bkg {karel_step5['n_background']:,.1f},"
          f" acc_eff = A * eff * migration = {out['acc_eff']:.4f}, lumi_pb {lumi:.3f}")
    return out


def plot(results, rev, karel_step5):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5.5))
    rows = []
    for r in results:
        rel = np.sqrt((r["A_stat"] / r["A"]) ** 2 + r["A_pdf_rel"] ** 2
                      + r["A_alphas_rel"] ** 2 + r["A_scale_rel"] ** 2)
        rows.append((f"DY prediction, {r['label']}", r["sigma_fid_pred"],
                     r["sigma_fid_pred"] * rel, "#1f77b4"))
    rows.append(("measured, step 5 as committed", karel_step5["sigma_fid_pb"],
                 karel_step5["sigma_total_pb"], "grey"))
    rows.append(("measured, revised", rev["sigma_fid"],
                 rev["sigma_fid"] * np.hypot(rev["stat_rel"], rev["syst_rel"]), "black"))
    for y, (label, value, err, colour) in enumerate(rows[::-1]):
        ax.errorbar(value, y, xerr=err, fmt="o", color=colour, capsize=4, markersize=7)
        ax.text(value, y + 0.2, f"{value:.1f} $\\pm$ {err:.1f} pb", ha="center", fontsize=11)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([row[0] for row in rows[::-1]], fontsize=12)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_xlabel(r"$\sigma_{fid}(pp \to Z \to \mu\mu)$ [pb]")
    hists._decorate(ax, lumi_fb=rev["lumi_pb"] / 1000)
    ax.text(0.99, 0.02, "predictions: 6077.22 pb x MC fiducial fraction (stat, PDF, $\\alpha_s$, scale);\n"
            "the uncertainty of the 6077.22 pb normalisation is not included",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10, style="italic")
    hists.save_fig(fig, "mc_fiducial_prediction.png")


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", choices=list(config.DY_SAMPLES), action="append")
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--filelists", type=Path, default=config.FILELIST_DIR)
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--karel-data", type=Path, default=config.DATA_DIR,
                    help="directory with the step 5 pickle")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.one_file:
        batch.write_part(process_file(args.one_file), args.out)
        return

    samples = args.sample or list(config.DY_SAMPLES)
    if not args.summarise_only:
        for sample in samples:
            cat = batch.read_catalogue(args.filelists / config.DY_SAMPLES[sample]["filelist"])
            urls = sorted(v[0] for v in cat.values())[: args.max_files]
            tasks = [(u, f"{sample}__{Path(u).stem}") for u in urls]
            print(f"[mc] {sample}: streaming {len(tasks)} files from EOS", flush=True)
            batch.run_files(tasks, Path(__file__), config.DATA_DIR / "mc_parts" / sample,
                            args.workers, retries=2, log=lambda s: print(s, flush=True))
            total = batch.load_parts(config.DATA_DIR / "mc_parts" / sample, blank())
            hists.save(total, f"mc_acceptance_{sample}.pkl")

    results = []
    for sample in samples:
        try:
            total = hists.load(f"mc_acceptance_{sample}.pkl")
        except FileNotFoundError:
            continue
        results.append(summarise(total, sample))
        hists.save(results[-1], f"mc_acceptance_{sample}_result.pkl")

    import pickle
    with open(args.karel_data / "step5_crosssection.pkl", "rb") as fh:
        karel_step5 = pickle.load(fh)
    try:
        xres = hists.load("xcheck_efficiency_result.pkl")
    except FileNotFoundError:
        print("[mc] run scripts/xcheck_efficiency.py first for the revised result")
        return
    if any(r["sample"] == "nlo" for r in results):
        rev = revised_result(results, karel_step5, xres)
        hists.save(rev, "revised_result.pkl")
        plot(results, rev, karel_step5)


if __name__ == "__main__":
    main()
