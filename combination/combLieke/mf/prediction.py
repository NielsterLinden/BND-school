"""The aMC@NLO prediction of sigma(pp -> Z/gamma* -> ll, 60 < m_LHE < 120 GeV) and its uncertainty.

Central value: 6077.22 pb (the NNLO-normalised cross section every channel uses for DYJetsToLL_M-50
amcatnloFXFX) x sum w(LHE flavour, 60-120) / sum w, from the generator sums of the z-mumu v2 skim
(z-mumu/output/v2/gensums.json, DY_NLO). This is exactly poi.reference_pb.

Uncertainty, from the generator's own weights evaluated on sigma(60 < m_LHE < 120), with the
prescription z-mumu uses for its templates (z-mumu/zmumu/weights.py):
    muR/muF   7-point envelope (LHEScaleWeight 0,1,3,4,5,7,8), asymmetric
    PDF       NNPDF3.1 Hessian, sqrt(sum_i (sigma_i - sigma_0)^2) over LHEPdfWeight 1..100
    alpha_s   0.75 x half the difference of LHEPdfWeight 101/102 (alpha_s 0.116/0.120 -> +-0.0015)
added in quadrature. The uncertainty of the NNLO normalisation 6077.22 pb itself is not public
(the CMS cross-section database is internal) and is not included.
"""

from __future__ import annotations

import json

import numpy as np

from .paths import repo_path

GENSUMS = "z-mumu/output/v2/gensums.json"
DY_XSEC_PB = 6077.22
SCALE_7POINT = [0, 1, 3, 4, 5, 7, 8]


def amcatnlo(flavour: str = "mumu") -> dict:
    g = json.loads(repo_path(GENSUMS).read_text())["DY_NLO"]
    nom = g[f"sumw_lhe_{flavour}_60_120"]
    scale = np.asarray(g[f"scale_lhe_{flavour}_60_120"])[SCALE_7POINT] / nom - 1.0
    pdf = np.asarray(g[f"pdf_lhe_{flavour}_60_120"])
    if not np.isclose(pdf[0], nom, rtol=1e-6):
        raise ValueError("PDF member 0 is not the nominal weight sum")
    pdf_rel = float(np.sqrt(np.sum((pdf[1:101] - pdf[0]) ** 2)) / pdf[0])
    alphas_rel = float(0.75 * 0.5 * abs(pdf[102] - pdf[101]) / pdf[0])
    sigma = DY_XSEC_PB * nom / g["sumw"]
    up = float(np.sqrt(max(scale.max(), 0) ** 2 + pdf_rel ** 2 + alphas_rel ** 2))
    down = float(np.sqrt(min(scale.min(), 0) ** 2 + pdf_rel ** 2 + alphas_rel ** 2))
    return {"label": "aMC@NLO", "flavour": flavour, "sigma_pb": sigma,
            "err_up_pb": up * sigma, "err_down_pb": down * sigma,
            "components_rel": {"scale_up": float(scale.max()), "scale_down": float(scale.min()),
                               "pdf": pdf_rel, "alphas": alphas_rel},
            "source": f"{GENSUMS} (DY_NLO); normalisation {DY_XSEC_PB} pb (NNLO), weights: see mf/prediction.py"}


def window_ratio(lo=(66, 116), hi=(60, 120)) -> float:
    """sigma(60-120) / sigma(66-116) at LHE level (1 GeV histogram h_lhe_mll of the same sample)."""
    h = np.asarray(json.loads(repo_path(GENSUMS).read_text())["DY_NLO"]["h_lhe_mll"])
    return float(h[hi[0]:hi[1]].sum() / h[lo[0]:lo[1]].sum())
