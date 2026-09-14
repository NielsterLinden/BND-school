"""Muon momentum scale and resolution from the Z peak (docs/14-fit-and-systematics.md).

Rochester corrections for UL2016 are not obtainable from this cluster, so the residual
data/MC difference is calibrated on the Z peak itself: in categories where both muons fall in
the same |eta| bin, the dimuon mass peak is fitted in data and in simulation with a
Breit-Wigner (PDG Z mass and width) convolved with a Gaussian (shift, width) on a linear
background. The MC muon pT is then scaled by kappa_j = 1 + (mu_data - mu_MC)/m_Z and smeared
by an extra relative resolution s_j = sqrt(2 (sigma_data^2 - sigma_MC^2)) / m_Z.

Nuisance parameters: `MuonScale` = kappa +/- its uncertainty (fit error, at least 0.05%),
`MuonRes` = smearing removed (Down) / doubled (Up).
"""

from __future__ import annotations

import json

import awkward as ak
import numpy as np

from . import objects

M_Z, G_Z = 91.1876, 2.4952
ETA_EDGES = np.array([0.0, 0.9, 1.2, 2.1, 2.4])
FIT_LO, FIT_HI, NBINS = 80.0, 102.0, 88
EDGES = np.linspace(FIT_LO, FIT_HI, NBINS + 1)
CENTRES = 0.5 * (EDGES[1:] + EDGES[:-1])
KAPPA_ERR_FLOOR = 0.0005


def blank():
    n = len(ETA_EDGES) - 1
    return {f"m_eta{j}": np.zeros(NBINS) for j in range(n)} | {f"m_eta{j}_w2": np.zeros(NBINS) for j in range(n)}


def fill(out, sr, weight):
    """`sr` is a regions.dimuon_regions()['SR'] dict; both muons in the same |eta| bin."""
    if sr is None or len(sr["idx"]) == 0:
        return
    j1 = np.digitize(np.abs(sr["eta1"]), ETA_EDGES) - 1
    j2 = np.digitize(np.abs(sr["eta2"]), ETA_EDGES) - 1
    w = weight
    for j in range(len(ETA_EDGES) - 1):
        sel = (j1 == j) & (j2 == j)
        out[f"m_eta{j}"] += np.histogram(sr["mass"][sel], bins=EDGES, weights=w[sel])[0]
        out[f"m_eta{j}_w2"] += np.histogram(sr["mass"][sel], bins=EDGES, weights=w[sel] ** 2)[0]


def _voigt_shifted(x, mu, sigma):
    from scipy.special import voigt_profile
    return voigt_profile(x - mu, sigma, G_Z / 2.0)


def fit_peak(h, var=None):
    """Fit (mu, sigma) of BW(m_Z) x Gauss(mu - m_Z, sigma) + linear background. Returns dict."""
    from iminuit import Minuit
    n = h
    v = var if var is not None else np.maximum(h, 1.0)
    k = v.sum() / max(n.sum(), 1e-9)
    x = CENTRES

    def model(N, mu, sigma, b0, b1):
        s = _voigt_shifted(x, mu, sigma)
        s = s / s.sum()
        b = np.clip(b0 + b1 * (x - 91.0), 0, None)
        return N * s + b

    def nll(N, mu, sigma, b0, b1):
        m = np.clip(model(N, mu, sigma, b0, b1), 1e-9, None) / k
        a = n / k
        return float(np.sum(m - a * np.log(m)))

    mn = Minuit(nll, N=n.sum(), mu=M_Z - 0.3, sigma=1.5, b0=max(n[:5].mean(), 0.1), b1=0.0)
    mn.errordef = Minuit.LIKELIHOOD
    mn.limits["N"] = (0, None); mn.limits["mu"] = (85, 95); mn.limits["sigma"] = (0.3, 5.0)
    mn.limits["b0"] = (0, None)
    mn.migrad(); mn.hesse()
    return {"mu": mn.values["mu"], "mu_err": mn.errors["mu"], "sigma": mn.values["sigma"],
            "sigma_err": mn.errors["sigma"], "valid": bool(mn.valid), "model": model(*mn.values), "n": n.sum()}


def calibrate(data, mc) -> dict:
    """Per |eta| bin: kappa (MC pT scale factor) and extra relative smearing."""
    res = {"eta_edges": ETA_EDGES.tolist(), "kappa": [], "kappa_err": [], "smear": [], "fits": {}}
    for j in range(len(ETA_EDGES) - 1):
        fd = fit_peak(data[f"m_eta{j}"], data[f"m_eta{j}_w2"])
        fm = fit_peak(mc[f"m_eta{j}"], mc[f"m_eta{j}_w2"])
        kappa = 1.0 + (fd["mu"] - fm["mu"]) / M_Z
        kerr = max(np.hypot(fd["mu_err"], fm["mu_err"]) / M_Z, KAPPA_ERR_FLOOR)
        smear = np.sqrt(max(2.0 * (fd["sigma"] ** 2 - fm["sigma"] ** 2), 0.0)) / M_Z
        res["kappa"].append(float(kappa)); res["kappa_err"].append(float(kerr)); res["smear"].append(float(smear))
        res["fits"][j] = {"data": {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in fd.items()},
                          "mc": {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in fm.items()}}
        print(f"  |eta| bin {j}: data peak {fd['mu']:.3f}+-{fd['mu_err']:.3f} width {fd['sigma']:.3f} | "
              f"MC {fm['mu']:.3f}+-{fm['mu_err']:.3f} width {fm['sigma']:.3f} -> kappa {kappa:.5f}+-{kerr:.5f}, "
              f"smear {100*smear:.2f}%")
    return res


class MomentumCalibration:
    def __init__(self, path):
        with open(path) as fh:
            d = json.load(fh)
        self.eta_edges = np.array(d["eta_edges"])
        self.kappa = np.array(d["kappa"]); self.kappa_err = np.array(d["kappa_err"]); self.smear = np.array(d["smear"])

    def scaled_pt(self, ev, variation="nominal"):
        """Muon_pt of a MC chunk with the calibration applied (deterministic smearing)."""
        eta = ak.to_numpy(ak.flatten(ev.Muon_eta))
        pt = ak.to_numpy(ak.flatten(ev.Muon_pt)).astype(np.float64)
        j = np.clip(np.digitize(np.abs(eta), self.eta_edges) - 1, 0, len(self.eta_edges) - 2)
        k = self.kappa[j]
        if variation == "MuonScaleUp":
            k = k + self.kappa_err[j]
        elif variation == "MuonScaleDown":
            k = k - self.kappa_err[j]
        s = self.smear[j]
        if variation == "MuonResUp":
            s = 2.0 * s
        elif variation == "MuonResDown":
            s = 0.0 * s
        seed = int(ak.to_numpy(ev.event)[0]) if len(ev) else 0
        rng = np.random.default_rng(seed)
        new = pt * k * (1.0 + s * rng.standard_normal(len(pt)))
        counts = ak.num(ev.Muon_pt, axis=1)
        return ak.values_astype(ak.unflatten(new, counts), np.float32)
