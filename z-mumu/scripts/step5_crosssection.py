#!/usr/bin/env python
"""Step 5 -- the fiducial cross section.

                     N_obs - N_bkg
    sigma_fid  =  ---------------------
                     eff_total * L

**Why fiducial and not total.** The total (inclusive) Z production cross
section requires extrapolating from the measured phase space to all pT and
rapidity, and that extrapolation factor -- the acceptance A -- can only come
from a theory calculation or simulation. No MC is available here, so instead
of inventing an acceptance we quote the cross section *in the phase space we
actually measured*, where A = 1 by construction:

    exactly two opposite-sign muons
    leading pT > 26 GeV, subleading pT > 20 GeV
    |eta| < 2.4
    60 < m(mumu) < 120 GeV

Anyone with a DY sample can multiply this number by 1/A to get the inclusive
cross section; nothing in the measurement has to be redone. See
docs/06-cross-section.md.

**The efficiency.** eff_total is built per muon and then per event:

    eff_total = eff_reco^2 * <eff_ID>_1 <eff_ISO>_1 * <eff_ID>_2 <eff_ISO>_2
                * eff_trigger(event)

The ID and isolation terms are the tag-and-probe maps from step 2, averaged
over the *observed* (pT, eta) distribution of the signal muons from step 1, so
a muon in a poorly-performing corner of the detector is weighted by how often
it actually occurs. Correlations between the two muons' efficiencies are
neglected; they are small because the two muons are well separated.

Outputs
    output/data/step5_crosssection.pkl   the result and full systematics table
    output/plots/step5_systematics.png   uncertainty breakdown

Run:  .venv/bin/python scripts/step5_crosssection.py   (needs steps 1-4)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from zmumu import config, hists, stats


def weighted_map_average(eff_map, eff_err, kin_hist):
    """Average an efficiency map over an observed (pT, eta) distribution.

    Bins where the tag-and-probe sample is empty would otherwise contribute a
    spurious zero, so they are dropped from both numerator and denominator and
    the fraction of signal muons they hold is reported.
    """
    weights = kin_hist.values()
    usable = eff_map > 0
    w = np.where(usable, weights, 0.0)
    total = w.sum()
    if total <= 0:
        raise ValueError("no overlap between efficiency map and signal kinematics")
    mean = float((w * eff_map).sum() / total)
    # Propagate the per-bin errors as a weighted quadrature sum.
    err = float(np.sqrt((w**2 * eff_err**2).sum()) / total)
    uncovered = float(weights[~usable].sum() / weights.sum()) if weights.sum() else 0.0
    return mean, err, uncovered


def main():
    step1 = hists.load("step1_selection.pkl")
    step2 = hists.load("step2_efficiency.pkl")["result"]
    step3 = hists.load("step3_backgrounds.pkl")
    step4 = hists.load("step4_signal.pkl")

    # ---- efficiency ------------------------------------------------------
    id1, id1_err, unc1 = weighted_map_average(step2["id_map"], step2["id_map_err"],
                                              step1["h_mu1_kin"])
    id2, id2_err, unc2 = weighted_map_average(step2["id_map"], step2["id_map_err"],
                                              step1["h_mu2_kin"])
    iso1, iso1_err, _ = weighted_map_average(step2["iso_map"], step2["iso_map_err"],
                                             step1["h_mu1_kin"])
    iso2, iso2_err, _ = weighted_map_average(step2["iso_map"], step2["iso_map_err"],
                                             step1["h_mu2_kin"])
    trig = step2["trig_event"]
    reco = step2["reco_eff"]

    eff_total = (reco**2) * id1 * id2 * iso1 * iso2 * trig

    # ---- yields ----------------------------------------------------------
    n_sig = step4["n_signal_counting"]
    n_obs = step4["n_observed"]
    lumi = config.LUMI_PB

    sigma = n_sig / (eff_total * lumi)

    # ---- uncertainties, as relative contributions ------------------------
    rel = {}
    rel["statistical (data)"] = np.sqrt(n_obs) / n_sig
    rel["background estimate"] = step3["n_bkg_err"] / n_sig
    rel["eff: ID (tag-and-probe stat)"] = np.hypot(id1_err / id1, id2_err / id2)
    rel["eff: isolation (T&P stat)"] = np.hypot(iso1_err / iso1, iso2_err / iso2)
    rel["eff: trigger (stat)"] = step2["trig_event_err"] / trig
    rel["eff: trigger (method bias)"] = config.TRIG_METHOD_REL_UNC
    rel["eff: reconstruction (external)"] = 2 * step2["reco_eff_err"] / reco

    # FSR: how much the yield in the window moves if FSR photons are not added
    # back. This is the genuine signal-extraction ambiguity (see step 4).
    n_bare = float(step1["h_os"].values().sum())
    n_fsr = float(step1["h_os_fsr"].values().sum())
    rel["FSR recovery"] = abs(n_fsr - n_bare) / n_fsr if n_fsr else 0.0

    rel["muon momentum scale"] = config.MUON_SCALE_REL_UNC
    rel["luminosity"] = config.LUMI_REL_UNC

    stat_rel = rel["statistical (data)"]
    syst_rel = stats.combine_in_quadrature(
        *[v for k, v in rel.items() if k != "statistical (data)"])
    total_rel = stats.combine_in_quadrature(stat_rel, syst_rel)

    result = {
        "sigma_fid_pb": sigma,
        "sigma_stat_pb": sigma * stat_rel,
        "sigma_syst_pb": sigma * syst_rel,
        "sigma_total_pb": sigma * total_rel,
        "eff_total": eff_total,
        "eff_id_mu1": id1, "eff_id_mu2": id2,
        "eff_iso_mu1": iso1, "eff_iso_mu2": iso2,
        "eff_trigger": trig, "eff_reco": reco,
        "n_signal": n_sig, "n_observed": n_obs,
        "n_background": step3["n_bkg"],
        "lumi_pb": lumi,
        "relative_uncertainties": rel,
        "stat_rel": stat_rel, "syst_rel": syst_rel, "total_rel": total_rel,
        "uncovered_kin_fraction": max(unc1, unc2),
        "fiducial_definition": {
            "n_muons": config.N_MUONS_REQUIRED,
            "charge": "opposite sign",
            "pt_lead_gev": config.MU_PT_LEAD,
            "pt_sublead_gev": config.MU_PT_SUBLEAD,
            "abs_eta_max": config.MU_ETA_MAX,
            "mass_window_gev": [config.MASS_LO, config.MASS_HI],
            "muon_id": "Muon_mediumId",
            "muon_iso": f"pfRelIso04_all < {config.MU_ISO_MAX}",
        },
    }
    _report(result)
    _plot(result)
    path = hists.save(result, "step5_crosssection.pkl")
    print(f"[step5] wrote {path}")


def _report(r):
    print("\n[step5] Fiducial cross section\n")
    print("  Efficiency breakdown")
    print(f"    eff(reco)  per muon, external : {r['eff_reco']:.4f}")
    print(f"    eff(ID)    leading / subleading: {r['eff_id_mu1']:.4f} / {r['eff_id_mu2']:.4f}")
    print(f"    eff(ISO)   leading / subleading: {r['eff_iso_mu1']:.4f} / {r['eff_iso_mu2']:.4f}")
    print(f"    eff(trigger) event level       : {r['eff_trigger']:.4f}")
    print(f"    eff(total) event               : {r['eff_total']:.4f}")
    if r["uncovered_kin_fraction"] > 0.001:
        print(f"    NOTE: {100*r['uncovered_kin_fraction']:.2f}% of signal muons fall in "
              "(pT,eta) bins with no tag-and-probe statistics")
    print()
    print("  Inputs")
    print(f"    N(observed)                    : {r['n_observed']:>14,.0f}")
    print(f"    N(background)                  : {r['n_background']:>14,.1f}")
    print(f"    N(signal)                      : {r['n_signal']:>14,.0f}")
    print(f"    integrated luminosity          : {r['lumi_pb']:>14,.1f} pb^-1")
    print()
    print("  Relative uncertainties")
    for k, v in sorted(r["relative_uncertainties"].items(), key=lambda x: -x[1]):
        print(f"    {k:<34s} {100*v:>8.3f}%")
    print(f"    {'-'*34} {'-'*9}")
    print(f"    {'statistical':<34s} {100*r['stat_rel']:>8.3f}%")
    print(f"    {'systematic':<34s} {100*r['syst_rel']:>8.3f}%")
    print(f"    {'TOTAL':<34s} {100*r['total_rel']:>8.3f}%")
    print()
    print("  " + "=" * 66)
    print(f"   sigma_fid(pp -> Z -> mu mu) = {r['sigma_fid_pb']:.1f}"
          f"  +/- {r['sigma_stat_pb']:.1f} (stat)"
          f"  +/- {r['sigma_syst_pb']:.1f} (syst) pb")
    print("  " + "=" * 66)


def _plot(r):
    import matplotlib.pyplot as plt
    import mplhep as hep

    items = sorted(r["relative_uncertainties"].items(), key=lambda x: x[1])
    labels = [k for k, _ in items]
    values = [100 * v for _, v in items]

    fig, ax = plt.subplots(figsize=(10, 7))
    ypos = np.arange(len(labels))
    ax.barh(ypos, values, color="#3f7fb5", height=0.65)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel("Relative uncertainty on $\\sigma_{fid}$ [%]")
    ax.set_xscale("log")
    for y, v in zip(ypos, values):
        ax.text(v * 1.12, y, f"{v:.3f}%", va="center", fontsize=10)
    ax.axvline(100 * r["total_rel"], color="crimson", linestyle="--",
               label=f"total = {100*r['total_rel']:.2f}%")
    ax.legend(fontsize=12)
    hep.cms.label("Open Data", data=True, lumi=round(config.LUMI_PB/1000, 1), year=2016, ax=ax)
    hists._title(ax, "Uncertainty breakdown")
    hists.save_fig(fig, "step5_systematics.png")


if __name__ == "__main__":
    main()
