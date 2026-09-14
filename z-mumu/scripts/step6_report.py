#!/usr/bin/env python
"""Step 6 -- collect everything into a report.

Reads the pickles written by steps 1-5 and produces:

    output/results.json      every number, machine readable
    output/RESULTS.md        the same numbers as a readable summary
    output/plots/step6_summary.png   the headline plot

Nothing is recomputed here -- if a number looks wrong, it is wrong in the step
that produced it.

Run:  .venv/bin/python scripts/step6_report.py   (needs steps 1-5)
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from zmumu import config, hists


def main():
    s1 = hists.load("step1_selection.pkl")
    s2 = hists.load("step2_efficiency.pkl")["result"]
    s3 = hists.load("step3_backgrounds.pkl")
    s4 = hists.load("step4_signal.pkl")
    s5 = hists.load("step5_crosssection.pkl")

    results = {
        "generated": str(date.today()),
        "dataset": {
            "path": str(config.SKIM_DIR),
            "eras": config.ERAS,
            "luminosity_pb": config.LUMI_PB,
            "luminosity_rel_unc": config.LUMI_REL_UNC,
        },
        "fiducial_definition": s5["fiducial_definition"],
        "cutflow": {k: int(v) for k, v in s1["cutflow"].items()},
        "efficiencies": {
            "reco_per_muon_external": s2["reco_eff"],
            "id_inclusive": s2["id_inclusive"],
            "iso_inclusive": s2["iso_inclusive"],
            "trigger_event": s2["trig_event"],
            "id_weighted_mu1": s5["eff_id_mu1"], "id_weighted_mu2": s5["eff_id_mu2"],
            "iso_weighted_mu1": s5["eff_iso_mu1"], "iso_weighted_mu2": s5["eff_iso_mu2"],
            "total_event": s5["eff_total"],
            "tag_probe_pairs": int(s2["n_pairs"]),
        },
        "yields": {
            "observed_os": s4["n_observed"],
            "same_sign": s3["n_ss"],
            "emu": s3["n_emu"],
            "background_fake": s3["n_fake"],
            "background_flavour_symmetric": s3["n_fs"],
            "background_total": s3["n_bkg"],
            "background_fraction": s3["bkg_fraction"],
            "signal_counting": s4["n_signal_counting"],
            "signal_fit_crosscheck": s4["n_signal_fit"],
        },
        "fit_validation": {
            "peak_gev": s4["fit_peak"], "peak_err_gev": s4["fit_peak_err"],
            "resolution_gev": s4["fit_sigma"], "resolution_err_gev": s4["fit_sigma_err"],
            "chi2": s4["chi2"], "ndf": s4["ndf"],
        },
        "cross_section": {
            "sigma_fid_pb": s5["sigma_fid_pb"],
            "stat_pb": s5["sigma_stat_pb"],
            "syst_pb": s5["sigma_syst_pb"],
            "total_pb": s5["sigma_total_pb"],
            "relative_uncertainties": s5["relative_uncertainties"],
        },
    }

    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = config.OUTPUT_DIR / "results.json"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)

    md_path = config.OUTPUT_DIR / "RESULTS.md"
    md_path.write_text(_markdown(results, s5))

    _summary_plot(s1, s3, s5)
    print(f"[step6] wrote {json_path}")
    print(f"[step6] wrote {md_path}")
    print("\n" + "=" * 70)
    print(f"  sigma_fid(pp -> Z -> mu mu) = {s5['sigma_fid_pb']:.1f} "
          f"+/- {s5['sigma_stat_pb']:.1f} (stat) "
          f"+/- {s5['sigma_syst_pb']:.1f} (syst) pb")
    print("=" * 70)


def _markdown(r, s5):
    fid = r["fiducial_definition"]
    cs = r["cross_section"]
    lines = [
        "# Z -> mu+ mu- fiducial cross section -- results",
        "",
        f"Generated {r['generated']} from CMS 2016 Open Data (SingleMuon, Run2016G+H).",
        "",
        "## Result",
        "",
        f"**sigma_fid = {cs['sigma_fid_pb']:.1f} +/- {cs['stat_pb']:.1f} (stat) "
        f"+/- {cs['syst_pb']:.1f} (syst) pb**",
        "",
        f"Total uncertainty {100*s5['total_rel']:.2f}% "
        f"({100*s5['stat_rel']:.2f}% stat, {100*s5['syst_rel']:.2f}% syst).",
        "",
        "## Fiducial volume",
        "",
        "| Requirement | Value |",
        "|---|---|",
        f"| muons | exactly {fid['n_muons']}, {fid['charge']} |",
        f"| leading muon pT | > {fid['pt_lead_gev']} GeV |",
        f"| subleading muon pT | > {fid['pt_sublead_gev']} GeV |",
        f"| muon abs(eta) | < {fid['abs_eta_max']} |",
        f"| muon ID | {fid['muon_id']} |",
        f"| muon isolation | {fid['muon_iso']} |",
        f"| dimuon mass | {fid['mass_window_gev'][0]}-{fid['mass_window_gev'][1]} GeV |",
        "",
        "Acceptance is 1 by construction: this is the phase space that was measured.",
        "",
        "## Cutflow",
        "",
        "| Cut | Events | Fraction |",
        "|---|---:|---:|",
    ]
    first = list(r["cutflow"].values())[0]
    for cut, n in r["cutflow"].items():
        lines.append(f"| {cut} | {n:,} | {100*n/first:.3f}% |")

    eff = r["efficiencies"]
    lines += [
        "",
        "## Efficiencies",
        "",
        "| Term | Value | Source |",
        "|---|---:|---|",
        f"| reconstruction (per muon) | {eff['reco_per_muon_external']:.4f} | "
        "external (CMS Muon POG) -- not measurable in NanoAOD |",
        f"| medium ID (inclusive) | {eff['id_inclusive']:.4f} | tag-and-probe |",
        f"| isolation (inclusive) | {eff['iso_inclusive']:.4f} | tag-and-probe |",
        f"| trigger (per event) | {eff['trigger_event']:.4f} | reference-trigger method |",
        f"| **total (per event)** | **{eff['total_event']:.4f}** | product, kinematics-weighted |",
        "",
        f"Tag-and-probe pairs used: {eff['tag_probe_pairs']:,}.",
        "",
        "## Yields",
        "",
        "| Quantity | Events |",
        "|---|---:|",
    ]
    y = r["yields"]
    lines += [
        f"| observed (opposite sign) | {y['observed_os']:,.0f} |",
        f"| same-sign control | {y['same_sign']:,.0f} |",
        f"| e-mu control | {y['emu']:,.0f} |",
        f"| background: non-prompt | {y['background_fake']:,.1f} |",
        f"| background: flavour-symmetric | {y['background_flavour_symmetric']:,.1f} |",
        f"| **background: total** | **{y['background_total']:,.1f}** "
        f"({100*y['background_fraction']:.3f}%) |",
        f"| **signal (counting, nominal)** | **{y['signal_counting']:,.0f}** |",
        f"| signal (fit, cross-check only) | {y['signal_fit_crosscheck']:,.0f} |",
        "",
        "## Uncertainty breakdown",
        "",
        "| Source | Relative |",
        "|---|---:|",
    ]
    for k, v in sorted(cs["relative_uncertainties"].items(), key=lambda x: -x[1]):
        lines.append(f"| {k} | {100*v:.3f}% |")
    lines += [
        f"| **total** | **{100*s5['total_rel']:.3f}%** |",
        "",
        "## Fit validation",
        "",
        f"- Fitted peak: {r['fit_validation']['peak_gev']:.3f} +/- "
        f"{r['fit_validation']['peak_err_gev']:.3f} GeV (PDG 91.1876)",
        f"- Fitted resolution: {r['fit_validation']['resolution_gev']:.3f} GeV",
        f"- chi2/ndf = {r['fit_validation']['chi2']:.0f}/{r['fit_validation']['ndf']} "
        f"= {r['fit_validation']['chi2']/r['fit_validation']['ndf']:.2f} "
        "(poor by construction: a single Voigt profile has no FSR tail)",
        "",
        "See `docs/` for the method behind each number.",
    ]
    return "\n".join(lines) + "\n"


def _summary_plot(s1, s3, s5):
    import matplotlib.pyplot as plt
    import mplhep as hep

    fig, ax = plt.subplots(figsize=(10, 7.5))
    hep.histplot(s1["h_os_fsr"], ax=ax, histtype="errorbar", color="black",
                 markersize=3, label="Data (opposite sign, FSR recovered)")
    hep.histplot(s3["h_bkg"], ax=ax, histtype="fill", color="#d98a3a", alpha=0.85,
                 label="Data-driven background")
    ax.set_yscale("log")
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    ax.set_ylabel("Events / 0.5 GeV")
    ax.legend(fontsize=12, loc="upper right")
    ax.set_ylim(top=ax.get_ylim()[1] * 60)
    ax.text(0.04, 0.94,
            r"$\sigma_{\mathrm{fid}} = %.1f \pm %.1f\,\mathrm{(stat)} \pm %.1f\,\mathrm{(syst)}$ pb"
            % (s5["sigma_fid_pb"], s5["sigma_stat_pb"], s5["sigma_syst_pb"]),
            transform=ax.transAxes, fontsize=15, va="top")
    ax.text(0.04, 0.86,
            "\n".join([
                r"exactly 2 OS $\mu$, $p_T > %g/%g$ GeV" % (config.MU_PT_LEAD, config.MU_PT_SUBLEAD),
                r"$|\eta| < %g$,  $%g < m_{\mu\mu} < %g$ GeV" % (
                    config.MU_ETA_MAX, config.MASS_LO, config.MASS_HI),
                r"$N_{\mathrm{sig}} = %s$,  $\epsilon = %.3f$" % (
                    f"{s5['n_signal']:,.0f}", s5["eff_total"]),
            ]),
            transform=ax.transAxes, fontsize=12, va="top")
    hep.cms.label("Open Data", data=True, lumi=round(config.LUMI_PB/1000, 1), year=2016, ax=ax)
    hists.save_fig(fig, "step6_summary.png")


if __name__ == "__main__":
    main()
