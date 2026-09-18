#!/usr/bin/env python
"""Freeze every number printed in section 1 (theory): the prediction the final plot carries, the published predictions and
measurement it is set against, the generator-level m_ll spectrum of our Drell-Yan sample, and the literature constants of the
lepton-universality and background clips.

    source fitting/setup.sh && python presentation/data/extract_theory_reference.py [--json PATH] [--check-only]

Reads (read-only):
  combination/combLieke/output/result.json    prediction: aMC@NLO sigma(Z/gamma* -> ll, 60 < m_LHE < 120) and its scale / PDF /
                                               alpha_s components (combLieke/mf/prediction.py) = the band of the final plot
  combination/combLieke/config/references.json the published CMS 13 TeV measurement (CMS-SMP-20-004)
  z-mumu/output/v2/gensums.json                DY_NLO: h_lhe_mll (1 GeV, 0-200 GeV, one lepton flavour), sumw, per-flavour 60-120 sums
Literature values are typed in below with their reference and checked against each other where possible (PDF-set predictions
re-read from arXiv:2408.03744 Table 5 on 17 Sep 2026). Writes presentation/data/theory_reference.json.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import REPO, Checker, finalize, load_json, provenance, standard_args  # noqa: E402

RESULT_JSON = REPO / "combination" / "combLieke" / "output" / "result.json"
REFERENCES_JSON = REPO / "combination" / "combLieke" / "config" / "references.json"
GENSUMS = REPO / "z-mumu" / "output" / "v2" / "gensums.json"
DY_XSEC_PB = 6077.22                       # sigma(Z/gamma* -> ll, m > 50), NNLO, all three flavours (repository CLAUDE.md)

# NNLO+NNLL (DYturbo v1.3.2) sigma(Z -> ll, 60 < m < 120), 13 TeV, per PDF set: arXiv:2408.03744 (CMS-SMP-20-004), Table 5
PDF_SETS = [
    {"label": "CT18", "value": 1921.0, "up": 30.0, "down": 33.0},
    {"label": "MSHT20", "value": 1935.0, "up": 23.0, "down": 27.0},
    {"label": "NNPDF3.1", "value": 1940.0, "up": 15.0, "down": 21.0},
    {"label": "NNPDF4.0", "value": 1970.0, "up": 11.0, "down": 14.0},
]
PDF_SETS_REF = "CMS-SMP-20-004, JHEP 04 (2025) 162, arXiv:2408.03744, Table 5 (NNLO+NNLL, DYturbo; stat, PDF, alpha_s, muR/muF)"

# the perturbative orders of the inclusive Drell-Yan cross section
ORDERS = [
    {"order": "LO", "year": 1970, "reference": "S. D. Drell, T.-M. Yan, PRL 25 (1970) 316"},
    {"order": "NLO", "year": 1979, "reference": "G. Altarelli, R. K. Ellis, G. Martinelli, NPB 157 (1979) 461"},
    {"order": "NNLO", "year": 1991, "reference": "R. Hamberg, W. L. van Neerven, T. Matsuura, NPB 359 (1991) 343"},
    {"order": "N3LO", "year": 2022, "reference": "C. Duhr, B. Mistlberger, JHEP 03 (2022) 116 (neutral current, gamma* + Z)"},
]

LEPTON_MASS_MEV = {"e": 0.51099895, "mu": 105.6583755, "tau": 1776.93}          # PDG 2024
BR_LL = {"value": 3.3658, "err": 0.0023, "unit": "%", "reference": "PDG 2024, Z boson, Gamma(l+ l-)/Gamma_total (lepton universality)"}
LEP_RATIOS = {"mumu_over_ee": {"value": 1.0009, "err": 0.0028}, "tautau_over_ee": {"value": 1.0019, "err": 0.0032}, "correlation": 0.63,
              "reference": "ALEPH, DELPHI, L3, OPAL, SLD, Phys. Rept. 427 (2006) 257 (hep-ex/0509008); PDG Z boson review"}
SIN2_EFF = {"value": 0.23153, "err": 0.00016, "reference": "Phys. Rept. 427 (2006) 257, LEP + SLD average"}
# (T3, Q) of the left-handed member; g_V = T3 - 2 Q sin^2, g_A = T3
FERMIONS = {"nu": (0.5, 0.0), "l": (-0.5, -1.0), "u": (0.5, 2.0 / 3.0), "d": (-0.5, -1.0 / 3.0)}

M_Z_GEV, SQRT_S_GEV, Y_MAX = 91.1876, 13000.0, 2.4                               # PDG; the LHC at 13 TeV; central boson rapidity
# LHC Higgs WG YR4 (CERN-2017-002-M), m_H = 125.09 GeV, 13 TeV: ggF N3LO, VBF, WH, ZH, ttH, bbH, tHq, tHW [pb]; BR(H -> tautau)
HIGGS_XS_PB = {"ggF": 48.58, "VBF": 3.782, "WH": 1.373, "ZH": 0.8839, "ttH": 0.5071, "bbH": 0.4880, "tHq": 0.07425, "tHW": 0.01517}
HIGGS_BR_TAUTAU = 0.06272
HIGGS_REF = "LHC Higgs Cross Section WG, Yellow Report 4 (CERN-2017-002-M), m_H = 125.09 GeV, 13 TeV"

SPEC_LO, SPEC_HI, SPEC_REBIN = 50, 200, 2         # the m_ll spectrum of clip mass_window: 2 GeV bins, 50-200 GeV


def build():
    res = load_json(RESULT_JSON)
    pred = res["prediction"]
    comp = pred["components_rel"]
    s0 = float(pred["sigma_pb"])
    steps = []                                     # the band grown component by component, in quadrature
    up2 = down2 = 0.0
    for name, u, d in (("scale", comp["scale_up"], -comp["scale_down"]), ("alphas", comp["alphas"], comp["alphas"]), ("pdf", comp["pdf"], comp["pdf"])):
        up2 += u * u
        down2 += d * d
        steps.append({"added": name, "rel_up": math.sqrt(up2), "rel_down": math.sqrt(down2),
                      "lo_pb": s0 * (1 - math.sqrt(down2)), "hi_pb": s0 * (1 + math.sqrt(up2))})

    refs = load_json(REFERENCES_JSON)
    cms = next(r for r in refs["combined"] if r["label"] == "CMS")
    cms_out = {"value": cms["value"], **cms["errors"], "total": math.sqrt(sum(v * v for v in cms["errors"].values())),
               "window": cms["window"], "detail": cms["detail"], "citation": cms["citation"]}

    g = load_json(GENSUMS)["DY_NLO"]
    h = np.asarray(g["h_lhe_mll"], dtype=float)    # 1 GeV bins from 0 GeV
    to_pb = DY_XSEC_PB / float(g["sumw"])
    edges = list(range(SPEC_LO, SPEC_HI + 1, SPEC_REBIN))
    dens = [float(h[a:a + SPEC_REBIN].sum() * to_pb / SPEC_REBIN) for a in edges[:-1]]
    sig_60_120 = float(h[60:120].sum() * to_pb)
    sig_m50 = float(h.sum() * to_pb)

    c = M_Z_GEV / SQRT_S_GEV
    higgs_xs = sum(HIGGS_XS_PB.values())

    gv = {f: t3 - 2 * q * SIN2_EFF["value"] for f, (t3, q) in FERMIONS.items()}
    return {
        "provenance": provenance("extract_theory_reference.py", [RESULT_JSON, REFERENCES_JSON, GENSUMS],
                                 dataset="aMC@NLO DYJetsToLL_M-50 generator sums (z-mumu v2 skim); literature values with references",
                                 note="section 1 (theory) of the BND-school talk; every number printed there comes from this file"),
        "unit": "pb",
        "prediction": {"label": pred["label"], "value": s0, "err_up": float(pred["err_up_pb"]), "err_down": float(pred["err_down_pb"]),
                       "components_rel": comp, "band_steps": steps, "source": "combination/combLieke/output/result.json prediction (" + pred["source"] + ")"},
        "pdf_sets": {"sets": PDF_SETS, "reference": PDF_SETS_REF, "spread_pb": max(p["value"] for p in PDF_SETS) - min(p["value"] for p in PDF_SETS)},
        "cms_2024": cms_out,
        "orders": ORDERS,
        "x_range": {"m_Z_GeV": M_Z_GEV, "sqrt_s_GeV": SQRT_S_GEV, "y_max": Y_MAX, "x_central": c,
                    "x_min": c * math.exp(-Y_MAX), "x_max": c * math.exp(Y_MAX), "tau": c * c,
                    "note": "leading order: x_{1,2} = (m_Z / sqrt(s)) exp(+-y), x1 x2 = m_Z^2 / s"},
        "leptons": {"mass_MeV": LEPTON_MASS_MEV, "mass_reference": "PDG 2024",
                    "ratio_mu_e": LEPTON_MASS_MEV["mu"] / LEPTON_MASS_MEV["e"], "ratio_tau_mu": LEPTON_MASS_MEV["tau"] / LEPTON_MASS_MEV["mu"],
                    "ratio_tau_e": LEPTON_MASS_MEV["tau"] / LEPTON_MASS_MEV["e"],
                    "br_ll": BR_LL, "lep_ratios": LEP_RATIOS, "sin2_eff": SIN2_EFF,
                    "couplings": {f: {"T3": t3, "Q": q, "g_V": gv[f], "g_A": t3} for f, (t3, q) in FERMIONS.items()},
                    "tau_phase_space_axial": (1 - 4 * (LEPTON_MASS_MEV["tau"] / 1e3) ** 2 / M_Z_GEV ** 2) ** 1.5},
        "spectrum": {"edges_GeV": edges, "dsigma_dm_pb_per_GeV": dens, "flavour": "one lepton flavour (h_lhe_mll)",
                     "sigma_60_120_pb": sig_60_120, "sigma_m50_pb": sig_m50, "fraction_60_120": sig_60_120 / sig_m50,
                     "source": "z-mumu/output/v2/gensums.json DY_NLO.h_lhe_mll x 6077.22 pb / sumw"},
        "higgs": {"xs_pb": HIGGS_XS_PB, "xs_total_pb": higgs_xs, "br_tautau": HIGGS_BR_TAUTAU, "xs_br_pb": higgs_xs * HIGGS_BR_TAUTAU,
                  "z_over_h_tautau": s0 / (higgs_xs * HIGGS_BR_TAUTAU), "reference": HIGGS_REF,
                  "note": "same final state: sigma(Z/gamma* -> tautau, 60-120) = the prediction value (lepton universality) vs sigma(pp -> H) x BR(H -> tautau)"},
    }


def verify(d, ck: Checker):
    R = "combination/combLieke/output/result.json"
    ck.check("prediction 1953.93 +56.06 -81.66 pb", [d["prediction"]["value"], d["prediction"]["err_up"], d["prediction"]["err_down"]], [1953.93, 56.06, 81.66], 0.01, R)
    last = d["prediction"]["band_steps"][-1]
    ck.check("band after the last step == the quoted band", [last["hi_pb"] - d["prediction"]["value"], d["prediction"]["value"] - last["lo_pb"]],
             [d["prediction"]["err_up"], d["prediction"]["err_down"]], 1e-6, "consistency")
    ck.check("scale +2.46 / -3.91 %, alpha_s 1.28 %, PDF 0.74 %", [d["prediction"]["components_rel"][k] for k in ("scale_up", "scale_down", "alphas", "pdf")],
             [0.0246, -0.0391, 0.0128, 0.0074], 0.0001, R)
    ck.check("PDF sets 1921 / 1935 / 1940 / 1970 pb", [p["value"] for p in d["pdf_sets"]["sets"]], [1921, 1935, 1940, 1970], 0, "arXiv:2408.03744 Table 5")
    ck.check("PDF-set spread 49 pb", d["pdf_sets"]["spread_pb"], 49, 0, "consistency")
    ck.check("CMS 1952 +- 4 +- 18 +- 45 pb", [d["cms_2024"][k] for k in ("value", "stat", "syst", "lumi")], [1952, 4, 18, 45], 0, "combination/combLieke/config/references.json")
    ck.check("m_ll spectrum: sigma(60-120) == the prediction", d["spectrum"]["sigma_60_120_pb"], d["prediction"]["value"], 1e-9, "same generator sums", rel=True)
    ck.check("sigma(m > 50) per flavour ~ 6077.22 / 3 pb", d["spectrum"]["sigma_m50_pb"], DY_XSEC_PB / 3, 0.002, "repository CLAUDE.md", rel=True)
    ck.check("fraction 60-120 = 0.965", d["spectrum"]["fraction_60_120"], 0.965, 0.0005, "gensums.json")
    spec = d["spectrum"]
    ck.check("spectrum bins x width over 60-120 == sigma(60-120)",
             sum(v * 2 for e, v in zip(spec["edges_GeV"][:-1], spec["dsigma_dm_pb_per_GeV"]) if 60 <= e < 120), spec["sigma_60_120_pb"], 1e-6, "consistency", rel=True)
    ck.check("x range 6.4e-4 ... 0.077 at 13 TeV, |y| < 2.4", [d["x_range"]["x_min"], d["x_range"]["x_max"]], [6.36e-4, 0.0773], 0.002, "computed", rel=True)
    ck.check("m_mu/m_e = 206.77, m_tau/m_e = 3477", [d["leptons"]["ratio_mu_e"], d["leptons"]["ratio_tau_e"]], [206.77, 3477.4], 0.001, "PDG", rel=True)
    ck.check("g_V(l) = -0.037", d["leptons"]["couplings"]["l"]["g_V"], -0.03694, 1e-5, "computed")
    ck.check("tau phase space (axial) -0.23 %", d["leptons"]["tau_phase_space_axial"], 0.99772, 1e-4, "computed")
    ck.check("H -> tautau 3.49 pb, Z/H ~ 560", [d["higgs"]["xs_br_pb"], d["higgs"]["z_over_h_tautau"]], [3.494, 559.2], 0.01, "YR4", rel=True)
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "theory_reference.json")
    finalize(ap.parse_args(), build, verify, "theory_reference")


if __name__ == "__main__":
    main()
