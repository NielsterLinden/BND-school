#!/usr/bin/env python
"""v2 step 9 -- uncertainty parity with CMS-SMP-20-004 (arXiv:2408.03744) and the comparison (docs/16).

    python scripts/v2_9_cms_parity.py                   # after steps 7 (reco T&P), 8 (generator study) and
                                                        # the fit with the measured reco SF (--tag zmumu_recosf)

Maps every row of the CMS systematic table (Table 7: total inclusive Z -> ll at 13 TeV) onto this
analysis, adds what the comparison showed to be missing -- measured where possible, estimated with three
options where not -- and writes

    output/v2/cms_parity/parity_zmumu.json      the parity file (fitting/uncertainty_parity.py format)
    output/v2/cms_parity/parity_zmumu.md        the table for docs/16
    output/v2/cms_parity/parity_zmumu_summary.json
    output/v2/cms_parity/inputs.json            every intermediate number, with its source
    output/v2/plots/cms_parity_option_{a,b,c}.png, cms_parity_options.png, cms_parity_estimate_inputs.png

Nothing here changes the channel's baseline result (fit/results/zmumu_fit_result.json).
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ROOT))

import numpy as np

from fitting import uncertainty_parity as up
from zmumu import config, hists, weights

OUT = config.OUTPUT_DIR / "v2" / "cms_parity"
PLOTS = config.OUTPUT_DIR / "v2" / "plots"
FIT = config.REPO_DIR / "fit" / "results"

# CMS-SMP-20-004, JHEP 04 (2025) 162, arXiv:2408.03744v2. Z -> l+l- (e and mu fitted together), 60-120 GeV, 13 TeV.
CMS = {
    "label": "CMS-SMP-20-004",
    "cite": "arXiv:2408.03744, Tables 4, 7, 10, 13 (Z -> l+l-, 13 TeV, 206 pb^-1, 2017 low-pileup runs)",
    "total_xsec": {"value": 1952.0, "stat": 4.0, "syst": 18.0, "lumi": 45.0},          # Table 13
    "fiducial_xsec": {"value": 754.0, "stat": 2.0, "syst": 3.0, "lumi": 17.0},        # Table 10
    # Table 7 (total inclusive), Z -> l+l- column, percent; Table 4 (fiducial) in brackets in the comments
    "table7": {"Total": 0.90, "muR and muF scales": 0.66, "PDF + alphaS": 0.43, "Resum. + FSR": 0.12,
               "Trigger prefire correction": 0.34, "Efficiency (stat)": 0.28, "QCD multijet (syst)": 0.14,
               "MC sim. stat": 0.16, "Efficiency (syst)": 0.16, "EW + tt cross section": 0.04,
               "QCD multijet (stat)": 0.07, "Hadronic recoil calibration": 0.05},
    "table4": {"Total": 0.40, "Efficiency (stat)": 0.26, "Trigger prefire correction": 0.26, "QCD multijet (syst)": 0.14,
               "EW + tt cross section": 0.01, "MC sim. stat": 0.12, "Efficiency (syst)": 0.15, "PDF + alphaS": 0.03,
               "muR and muF scales": 0.07, "QCD multijet (stat)": 0.02, "Hadronic recoil calibration": 0.02},
    "lumi_pct": 2.3,
}


# ----------------------------------------------------------------------------- inputs
def np_impact_pct(fit, name):
    """Post-fit impact of one nuisance parameter on mu, in percent (mean of the up/down shifts)."""
    for r in fit.get("ranking", []):
        if r["name"] == name:
            return 100 * 0.5 * (abs(r["dpoi_up_post"]) + abs(r["dpoi_down_post"])) / fit["mu"]
    return 0.0


def group_pct(fit, group):
    return 100 * fit["grouped_impacts_mu"].get(group, 0.0) / fit["mu"]


def acceptance_theory(gens):
    """Theory uncertainties of A(60-120) = fid(dressed) / LHE mumu (60 < m < 120), from the full-sample sums.

    The baseline fit result carries the same quantities for A(m > 50) (scripts/mc_acceptance.py, the v1
    review); the 60-120 GeV cross section needs the 60-120 GeV denominator."""
    g, pw = gens["DY_NLO"], gens["DY_powheg"]
    num, den = "fid_dressed", "lhe_mumu_60_120"
    A = g[f"sumw_{num}"] / g[f"sumw_{den}"]
    pdf = (np.array(g[f"pdf_{num}"]) / np.array(g[f"pdf_{den}"]))
    pdf = pdf / pdf[0]
    scale = np.array(g[f"scale_{num}"]) / np.array(g[f"scale_{den}"])
    scale = scale / scale[weights.SCALE_NOMINAL]
    ps = (np.array(g[f"ps_{num}"]) / np.array(g[f"ps_{den}"])) / A
    A_pw = pw[f"sumw_{num}"] / pw[f"sumw_{den}"]
    # binomial-like MC statistical uncertainty of a subset ratio
    n2, d2 = g[f"sumw2_{num}"], g["sumw2"]
    A_stat = np.sqrt(n2) / g[f"sumw_{num}"] * A           # dominated by the numerator
    out = {"A_60_120": A,
           "pdf_pct": 100 * float(np.sqrt(np.sum((pdf[weights.PDF_MEMBERS] - 1) ** 2))),
           "alphas_pct": 100 * float(0.75 * 0.5 * abs(pdf[102] - pdf[101])),
           "scale_pct": 100 * float(np.max(np.abs(scale[weights.SCALE_7POINT] - 1))),
           "mc_stat_pct": 100 * float(A_stat / A),
           "ps_isr_pct": 100 * float(0.5 * (abs(ps[weights.PS_ISR_UP] - 1) + abs(ps[weights.PS_ISR_DOWN] - 1))),
           "ps_fsr_pct": 100 * float(0.5 * (abs(ps[weights.PS_FSR_UP] - 1) + abs(ps[weights.PS_FSR_DOWN] - 1))),
           "generator_pct": 100 * float(abs(A_pw / A - 1)), "A_powheg": A_pw,
           "source": "output/v2/gensums.json (all 71.8 M aMC@NLO events, 2.96 M powheg events)"}
    return out


def ptz_reweighting(hists_all, theory):
    """Data/MC ratio of the reconstructed pT(mumu) in the SR (normalised), folded with A(pT(Z))."""
    mc_samples = ["DYmumu", "DYtautau", "DYee", "TTbar", "SingleTop", "WW", "WZ", "ZZ"]
    data = np.asarray(hists_all["Data|SR|zpt|nominal"], dtype=float)
    mc = sum(np.asarray(hists_all[f"{s}|SR|zpt|nominal"], dtype=float) for s in mc_samples
             if f"{s}|SR|zpt|nominal" in hists_all)
    pw = np.asarray(hists_all["DYmumu_powheg|SR|zpt|nominal"], dtype=float)
    sig = np.asarray(hists_all["DYmumu|SR|zpt|nominal"], dtype=float)
    reco_edges = np.linspace(0, 200, len(data) + 1)
    ratio = (data / data.sum()) / (mc / mc.sum())
    ratio_pw = (pw / pw.sum()) / (sig / sig.sum())
    e = np.array(theory["ptz_edges"])
    c = 0.5 * (e[1:] + e[:-1])
    idx = np.clip(np.digitize(c, reco_edges) - 1, 0, len(ratio) - 1)
    w = np.where(c < 200, ratio[idx], 1.0)
    w_pw = np.where(c < 200, ratio_pw[idx], 1.0)
    return {"reco_edges": reco_edges.tolist(), "ratio_data_mc": ratio.tolist(), "ratio_powheg_amcnlo": ratio_pw.tolist(),
            "gen_weights_data": w.tolist(), "gen_weights_powheg": w_pw.tolist()}


def fsr_photons(hists_all):
    mc_samples = ["DYmumu", "DYtautau", "DYee", "TTbar", "SingleTop", "WW", "WZ", "ZZ"]
    d = np.asarray(hists_all["Data|SR|nfsr|nominal"], dtype=float)
    m = sum(np.asarray(hists_all[f"{s}|SR|nfsr|nominal"], dtype=float) for s in mc_samples if f"{s}|SR|nfsr|nominal" in hists_all)
    fd, fm = d[1:].sum() / d.sum(), m[1:].sum() / m.sum()
    return {"frac_with_fsr_photon_data": fd, "frac_with_fsr_photon_mc": fm, "ratio": fd / fm,
            "n_data": d.tolist(), "n_mc": m.tolist()}


# ----------------------------------------------------------------------------- parity file
def build(base, fit, reco, acc, theory, ptz, fsr, dA_ptz):
    mu = fit["mu"]
    sigma = fit["sigma_60_120_pb"]
    stat_pct = 100 * fit["mu_stat_only_fit"] / mu
    t7 = CMS["table7"]
    reco_meta = reco["meta"]
    old_reco = np_impact_pct(base, "MuonReco")
    rows = []

    def row(**kw):
        kw.setdefault("rho_reference", 0.0)
        kw.setdefault("estimate", None)
        rows.append(kw)

    row(id="lumi", source="Luminosity", group="luminosity", reference_pct=CMS["lumi_pct"], ours_pct=100 * config.LUMI_REL_UNC,
        status="better", treatment="external",
        ours_how="1.2 % for the 2016 legacy calibration (CMS-LUM-17-003), normtag luminosity 16393 pb^-1",
        reference_how="2.3 % for the 2017 low-pileup runs (LUM-17-004)")
    scale_c = np_impact_pct(fit, "QCDScale")
    row(id="scale", source="μR and μF scales", group="theory", reference_pct=t7["muR and muF scales"],
        ours_pct=float(np.hypot(acc["scale_pct"], scale_c)), status="similar", treatment="measured", rho_reference=1.0,
        components={"acceptance (7-point, 60-120 GeV denominator)": acc["scale_pct"], "C factor (fit, QCDScale)": scale_c},
        ours_how=f"7-point envelope on A(60-120) {acc['scale_pct']:.2f} % ⊕ in the fit on C {scale_c:.2f} %; "
                 f"a CMS-like 25/25 GeV volume gives {theory['volumes']['dressed2525']['scale_7pt_rel']*100:.2f} %",
        reference_how="largest acceptance shift of the 6 non-extreme variations")
    pdf_c = float(np.hypot(np_impact_pct(fit, "PDF"), np_impact_pct(fit, "AlphaS")))
    row(id="pdf", source="PDF + αs", group="theory", reference_pct=t7["PDF + alphaS"],
        ours_pct=float(np.sqrt(acc["pdf_pct"] ** 2 + acc["alphas_pct"] ** 2 + pdf_c ** 2)), status="similar",
        treatment="measured", rho_reference=1.0,
        components={"acceptance PDF (NNPDF3.1 Hessian)": acc["pdf_pct"], "acceptance αs": acc["alphas_pct"],
                    "C factor (fit, PDF ⊕ AlphaS)": pdf_c},
        ours_how="NNPDF3.1 Hessian eigenvectors and αs ± 0.0015 on A(60-120) and in the fit on C",
        reference_how="NNPDF3.1 prescription on the acceptance")
    pt_part = max(abs(dA_ptz["data"]) * 100, acc["ps_isr_pct"])
    shower_c = float(np.sqrt(np_impact_pct(fit, "SigModel") ** 2 + np_impact_pct(fit, "PS_ISR") ** 2 + np_impact_pct(fit, "PS_FSR") ** 2))
    row(id="resum", source="Boson pT (resummation) and generator on A", group="theory",
        reference_pct=t7["Resum. + FSR"], status="measured-now", treatment="measured",
        ours_pct=float(np.sqrt(pt_part ** 2 + acc["generator_pct"] ** 2 + acc["ps_fsr_pct"] ** 2)),
        components={"pT(Z) spectrum reweighted to our data (max with PS ISR)": pt_part,
                    "powheg vs aMC@NLO": acc["generator_pct"], "PS FSR (QCD)": acc["ps_fsr_pct"]},
        ours_how=f"new: A with the pT(Z) spectrum reweighted to our measured pT(μμ) ({100*dA_ptz['data']:+.2f} %; "
                 f"PS ISR {acc['ps_isr_pct']:.2f} %), powheg vs aMC@NLO ({acc['generator_pct']:.2f} %), PS FSR "
                 f"({acc['ps_fsr_pct']:.2f} %)",
        reference_how="DYTURBO NNLO+NNLL vs aMC@NLO on A, together with PYTHIA vs PHOTOS")
    row(id="genC", source="Generator and parton shower on C", group="theory", reference_pct=None,
        ours_pct=shower_c, status="ours-only", treatment="measured",
        components={n: np_impact_pct(fit, n) for n in ("SigModel", "PS_ISR", "PS_FSR")},
        ours_how="in the fit: powheg vs aMC@NLO lineshape (SigModel), PS ISR/FSR weights, all on C",
        reference_how="no generator or shower term on the efficiency in the table")
    fsr_size = abs(theory["fsr"]["A_bare_over_dressed"] - 1) * 100
    r_fsr = abs(fsr["ratio"] - 1)
    row(id="qedfsr", source="QED FSR model (PHOTOS vs PYTHIA)", group="theory", reference_pct=None, status="partial",
        treatment="estimate", ours_pct=None,
        ours_how="not measurable here: estimated (three options)", reference_how="inside Resum. + FSR (0.12 %)",
        estimate={"why_not_measured": "every 2016 UL Drell-Yan sample in the Open Data (aMC@NLO FxFx, madgraph MLM, powheg) "
                                      "showers QED radiation with PYTHIA 8; a PHOTOS sample needs generation, GEANT4 simulation "
                                      "and reconstruction in CMSSW, which this cluster cannot run, and a generator-level PHOTOS run "
                                      "would only change A, while the sensitivity sits in C (reconstructed bare-muon pT cuts). "
                                      "The data rate of recovered FSR photons mixes the FSR model with the photon reconstruction "
                                      "efficiency, so it can size the effect but not measure it.",
                  "why_not_cited": "CMS-SMP-20-004 quotes PYTHIA vs PHOTOS only merged with the resummation uncertainty "
                                   "(0.12 %), for electrons and muons fitted together and for a Born-level fiducial volume with "
                                   "symmetric 25 GeV cuts; there is no muon-only or FSR-only number, and a Born-level volume is far "
                                   "more FSR-sensitive than our dressed one, so the number neither separates nor transfers. "
                                   "The generator-group PHOTOS/PYTHIA comparisons are internal.",
                  "options": {"a": t7["Resum. + FSR"], "b": 0.0, "c": float(r_fsr * fsr_size)},
                  "derivation": f"(size of the FSR effect) x (plausible model difference) = |A_bare/A_dressed - 1| "
                                f"= {fsr_size:.2f} % (generator level, bounds what the reconstructed bare-muon pT cuts can see) "
                                f"x |data/MC - 1| of the fraction of SR events with a recovered FSR photon = {100*r_fsr:.1f} % "
                                f"(attributing all of it to the FSR model, none to the photon efficiency) = {r_fsr*fsr_size:.3f} %"})
    row(id="prefire", source="L1 trigger prefiring", group="experimental", reference_pct=t7["Trigger prefire correction"],
        ours_pct=np_impact_pct(fit, "L1Prefiring"), status="similar", treatment="measured",
        ours_how="NanoAOD L1PreFiringWeight Up/Dn (muon + ECAL maps, stat ⊕ 20 %)",
        reference_how="prefiring maps, stat ⊕ ±20 % of the probability")
    eff = float(np.sqrt(sum(np_impact_pct(fit, n) ** 2 for n in ("MuonID", "MuonIso", "MuonTrigger"))))
    row(id="eff", source="Lepton ID, isolation, trigger efficiency", group="experimental",
        reference_pct=float(np.hypot(t7["Efficiency (stat)"], t7["Efficiency (syst)"])), ours_pct=eff, status="similar",
        treatment="measured",
        components={n: np_impact_pct(fit, n) for n in ("MuonID", "MuonIso", "MuonTrigger")},
        ours_how="own T&P fits, stat ⊕ syst (bkg shape, fit range, tag, powheg template, closure) shifted coherently "
                 "over cells (conservative for the stat part); no charge binning (second order for Z: SF(μ+)SF(μ-))",
        reference_how="T&P in charge × pT × η bins, stat uncorrelated per bin; syst: fit model, PHOTOS FSR in templates, tag")
    row(id="reco", source="Muon reconstruction (tracking, stand-alone)", group="experimental", reference_pct=None,
        ours_pct=np_impact_pct(fit, "MuonReco"), status="measured-now", treatment="measured",
        components={"assigned before (0.4 %/muon)": old_reco},
        ours_how=f"new: T&P on the unskimmed NanoAOD (stand-alone and isolated-track probes), SF per muon "
                 f"{reco_meta['sf_reco_per_muon']:.4f} ± {reco_meta['sf_reco_per_muon_err']:.4f}; was 1 ± 0.004 assigned "
                 f"({old_reco:.2f} % on μ)",
        reference_how="T&P, inside Efficiency (stat)/(syst)")
    row(id="fakes", source="QCD multijet / non-prompt", group="experimental",
        reference_pct=float(np.hypot(t7["QCD multijet (syst)"], t7["QCD multijet (stat)"])), ours_pct=group_pct(fit, "Fakes"),
        status="similar", treatment="measured", ours_how="fake factor (anti-isolated muons), stat + 30 % method",
        reference_how="misidentification factor (iso/anti-iso), stat + ±10/20 % prompt contamination")
    row(id="mcstat", source="MC statistics", group="experimental", reference_pct=t7["MC sim. stat"],
        ours_pct=float(np.hypot(group_pct(fit, "Gammas"), acc["mc_stat_pct"])), status="similar", treatment="measured",
        ours_how="Barlow-Beeston lite gammas per bin ⊕ MC stat on A", reference_how="Barlow-Beeston")
    row(id="bkgxs", source="EW + tt̄ cross sections", group="experimental", reference_pct=t7["EW + tt cross section"],
        ours_pct=group_pct(fit, "Background normalisation"), status="similar", treatment="measured",
        ours_how="6 % tt̄, 10 % tW/WW/WZ/ZZ, 5 % Z→ττ", reference_how="5 % top and diboson; τ backgrounds via theory ratios")
    row(id="recoil", source="Hadronic recoil calibration", group="experimental", reference_pct=t7["Hadronic recoil calibration"],
        ours_pct=0.0, status="not-applicable", treatment="none",
        ours_how="no p_T^miss in the Z selection", reference_how="W recoil calibration, enters Z via the joint W+Z fit")
    row(id="momentum", source="Muon momentum scale and resolution", group="experimental", reference_pct=None,
        ours_pct=group_pct(fit, "Muon momentum"), status="similar", treatment="measured",
        ours_how="Z-peak calibration per |η| (Rochester corrections not obtainable), shape NPs",
        reference_how="Rochester corrections; shape NP, 'small', not tabulated")
    row(id="pileup", source="Pileup", group="experimental", reference_pct=None, ours_pct=group_pct(fit, "Pileup"),
        status="ours-only", treatment="measured", ours_how="N_PV-matched profile, σ_minbias ± 4.6 %",
        reference_how="low-pileup runs (<μ> ≈ 3): not needed")
    # the proper measured systematic: the fit's correlated total without luminosity (MINOS, grouped) and the
    # acceptance rows outside the fit in quadrature -- not the sum of individual impacts
    fit_syst = 100 * fit["mu_syst_excl_lumi"] / mu
    acc_syst = float(np.sqrt(acc["scale_pct"] ** 2 + acc["pdf_pct"] ** 2 + acc["alphas_pct"] ** 2 + acc["mc_stat_pct"] ** 2
                             + pt_part ** 2 + acc["generator_pct"] ** 2 + acc["ps_fsr_pct"] ** 2))
    syst_measured = float(np.hypot(fit_syst, acc_syst))
    return {
        "schema": "uncertainty-parity/1", "channel": "zmumu", "channel_label": "Z→μμ (this work)", "colour": "#2a78d6",
        "observable": "σ(Z/γ*→ℓℓ, 60 < m < 120 GeV)", "unit": "pb",
        "reference": {"label": CMS["label"], "cite": CMS["cite"], "value": CMS["total_xsec"]["value"],
                      "stat_pct": 100 * CMS["total_xsec"]["stat"] / CMS["total_xsec"]["value"],
                      "lumi_pct": CMS["lumi_pct"],
                      "total_pct": 100 * float(np.sqrt(sum(CMS["total_xsec"][k] ** 2 for k in ("stat", "syst", "lumi"))))
                      / CMS["total_xsec"]["value"],
                      "note": "e and μ fitted together (lepton universality); 2017 low-pileup runs, 13 TeV"},
        "measurement": {"value": sigma, "stat_pct": stat_pct, "fit": fit.get("_path"),
                        "syst_pct_measured": syst_measured, "syst_pct_fit": fit_syst, "syst_pct_acceptance": acc_syst,
                        "label": "σ_fid / A(60-120), measured reconstruction SF applied"},
        "rows": rows,
    }


@up._default_style
def plot_estimate_inputs(ptz, theory, fsr, dA, path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    e = np.array(ptz["reco_edges"])
    c = 0.5 * (e[1:] + e[:-1])
    ax = axes[0]
    ax.step(e[:-1], ptz["ratio_data_mc"], where="post", color=up.INK, label="data / aMC@NLO (+ bkg), normalised")
    ax.step(e[:-1], ptz["ratio_powheg_amcnlo"], where="post", color="#eb6834", label="powheg / aMC@NLO (signal)")
    ax.axhline(1, color="#9a9a9a", lw=0.8)
    ax.set_xlim(0, 100)
    ax.set_ylim(0.75, 1.25)
    ax.set_xlabel(r"reconstructed $p_T(\mu\mu)$ [GeV]")
    ax.set_ylabel("shape ratio")
    ax.set_title(f"pT(Z) shape: reweighting A to data gives {100*dA['data']:+.2f} %", loc="left", fontsize=12)
    ax.legend(fontsize=10, frameon=False)
    ax = axes[1]
    ge = np.array(theory["ptz_edges"])
    ge[-1] = 200
    gc = 0.5 * (ge[1:] + ge[:-1])
    for v, col, lab in (("dressed", "#2a78d6", "dressed, 26/20 GeV (ours)"), ("bare", "#9a9a9a", "bare muons, 26/20 GeV"),
                        ("dressed2525", "#1baf7a", "dressed, 25/25 GeV (CMS-like)")):
        ax.plot(gc, theory["volumes"][v]["A_ptz"], "o-", color=col, label=lab, markersize=4, lw=1.5)
    ax.set_xscale("log")
    ax.set_xlim(1, 200)
    ax.set_xlabel(r"generator $p_T(Z)$ [GeV]")
    ax.set_ylabel(r"A($p_T$(Z)) = fiducial / LHE $\mu\mu$ 60-120")
    ax.set_title(f"acceptance vs pT(Z); bare/dressed = {theory['fsr']['A_bare_over_dressed']:.4f}", loc="left", fontsize=12)
    ax.legend(fontsize=10, frameon=False)
    ax = axes[2]
    nd, nm = np.array(fsr["n_data"]), np.array(fsr["n_mc"])
    x = np.arange(len(nd))
    ax.bar(x - 0.2, nd / nd.sum(), width=0.38, color=up.INK, label="data")
    ax.bar(x + 0.2, nm / nm.sum(), width=0.38, color="#2a78d6", label="simulation")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xlabel("recovered FSR photons per SR event")
    ax.set_ylabel("fraction of events")
    ax.set_title(f"≥ 1 FSR photon: data/MC = {fsr['ratio']:.3f}", loc="left", fontsize=12)
    ax.legend(fontsize=10, frameon=False)
    for ax in axes:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle("Inputs of the resummation measurement and of the QED FSR estimate (docs/16)", fontsize=14, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight", facecolor=up.BG)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", type=Path, default=FIT / "zmumu_recosf_fit_result.json")
    ap.add_argument("--baseline", type=Path, default=FIT / "zmumu_fit_result.json")
    args = ap.parse_args()
    base = json.load(open(args.baseline))
    fit = json.load(open(args.fit))
    fit["_path"] = str(args.fit.resolve().relative_to(config.REPO_DIR.parent.resolve()))
    reco = json.load(open(config.OUTPUT_DIR / "v2" / "tnp" / "reco_result.json"))
    gens = json.load(open(config.OUTPUT_DIR / "v2" / "gensums.json"))
    theory = json.load(open(OUT / "theory_acceptance.json"))
    hists_all = pickle.load(open(config.OUTPUT_DIR / "v2" / "histograms.pkl", "rb"))

    acc = acceptance_theory(gens)
    ptz = ptz_reweighting(hists_all, theory)
    fsr = fsr_photons(hists_all)
    # fold the reweighting into A(pT(Z)); the denominator spectrum is the generator one of step 8
    import importlib
    step8 = importlib.import_module("v2_8_theory_acceptance")
    den = np.array(theory["den_ptz"])
    dA = {}
    for name, key, vol in (("data", "gen_weights_data", "dressed"), ("powheg", "gen_weights_powheg", "dressed"),
                           ("data_2525", "gen_weights_data", "dressed2525")):
        A_bin = np.array(theory["volumes"][vol]["A_ptz"])
        dA[name] = step8.reweight_acceptance(A_bin, den, np.array(ptz[key])) / theory["volumes"][vol]["A"] - 1
    parity = build(base, fit, reco, acc, theory, ptz, fsr, dA)
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = {"acceptance_60_120": acc, "ptz_reweighting": {**ptz, "dA_rel": dA}, "fsr_photons": fsr,
              "reco_sf": reco["meta"], "fit": str(args.fit), "baseline_fit": str(args.baseline),
              "sigma_60_120_baseline_pb": base["sigma_60_120_pb"], "sigma_60_120_pb": fit["sigma_60_120_pb"],
              "cms": CMS}
    (OUT / "inputs.json").write_text(json.dumps(inputs, indent=1, default=float))
    path = OUT / "parity_zmumu.json"
    path.write_text(json.dumps(parity, indent=1, default=float))
    p = up.load(path)
    for o in up.OPTIONS:
        up.plot_option(p, o, PLOTS / f"cms_parity_option_{o}.png")
    up.plot_overview(p, PLOTS / f"cms_parity_options.png")
    plot_estimate_inputs(ptz, theory, fsr, dA, PLOTS / "cms_parity_estimate_inputs.png")
    md = up.markdown(p)
    path.with_suffix(".md").write_text(md + "\n")
    (OUT / "parity_zmumu_summary.json").write_text(json.dumps(up.summary(p), indent=1, default=float))
    print(md)
    print(f"\n[parity] acceptance (60-120): {json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in acc.items()})}")
    print(f"[parity] pT(Z) reweighting: dA/A = {dA}")
    print(f"[parity] FSR photons: {fsr['frac_with_fsr_photon_data']:.4f} data vs {fsr['frac_with_fsr_photon_mc']:.4f} MC")


if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    main()
