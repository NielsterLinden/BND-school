"""Pileup reweighting from the CMS Open Data luminosity record (docs/11-mc-weights.md).

No official `puWeights.json` is reachable from this cluster, so the data profile is built
from the per-lumisection table of record 1059 (`pp_2016lumibyls.csv`, brilcalc with the
PHYSICS normtag):

  * only the certified lumisections of runs 278820-284044 (Run2016G+H) enter, weighted by
    the recorded luminosity; the sum of that luminosity reproduces 16393.381 pb^-1 exactly;
  * brilcalc's `avgpu` column is computed with sigma_minbias = 80 mb, while the CMS
    Ultra-Legacy recommendation is 69.2 mb: the profile is `avgpu x 69.2/80`;
  * the CSV gives one average per lumisection and nothing about the bunch-to-bunch spread
    that the official `pileupCalc.py` ("true" mode) adds from the per-bunch luminosity.
    Two parameters restore what is missing: a global `scale` of the per-lumisection mean and a
    Gaussian `rel_smear` (sigma / mean) applied to every lumisection. They are fixed by
    requiring the reweighted simulation to reproduce the N_PV distribution of the Z -> mu mu
    signal region (`match_npv`, run by scripts/v2_2_pileup.py); the result (scale 1.035,
    rel_smear 0.09 with the v2 skims) is the "N_PV-matched" profile;
  * the +/-4.6% variation of the minimum-bias cross section (72.4 / 66.0 mb, applied on top
    of the matched scale) gives the `Pileup` Up/Down weights.

The MC profile is the genWeight-weighted `Pileup_nTrueInt` histogram of *all* generated
events, stored in the `GenSums` tree of every skim file.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from . import io, skim

CSV_PATH = Path("/data/atlas/users/nterlind/BND-school-cache/external/pp_2016lumibyls.csv")
CSV_URL = "root://eospublic.cern.ch//eos/opendata/cms/luminosity/2016/pp_2016lumibyls.csv"
SIGMA_MB_BRIL = 80.0          # what brilcalc used for `avgpu`
SIGMA_MB_NOMINAL = 69.2       # CMS UL recommendation, +/- 4.6%
SIGMA_MB_UNC = 0.046
NBINS = skim.PU_BINS
MIN_SMEAR = 0.3               # interactions; keeps the per-lumisection Gaussian well defined


def _lumisections(csv_path=CSV_PATH, run_min=skim.RUN_MIN, run_max=skim.RUN_MAX):
    """(recorded lumi [pb^-1], avgpu at 80 mb) arrays of the certified lumisections."""
    grl = io.load_grl()
    rec, mu = [], []
    with open(csv_path) as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith("#"):
                continue
            run = int(row[0].split(":")[0])
            if run < run_min or run > run_max:
                continue
            ls = int(row[1].split(":")[0])
            blocks = grl.get(run)
            if blocks is None or not any(lo <= ls <= hi for lo, hi in blocks):
                continue
            rec.append(float(row[6]) * 1e3)          # /fb -> /pb
            mu.append(float(row[7]))
    # aggregate lumisections with the same avgpu (rounded to 0.01 interactions at 80 mb): ~150 k
    # lumisections -> ~1000 weighted entries, which makes the smeared profiles ~100x faster
    mu = np.round(np.asarray(mu), 2)
    uniq, inv = np.unique(mu, return_inverse=True)
    rec_agg = np.zeros(len(uniq))
    np.add.at(rec_agg, inv, np.asarray(rec))
    return rec_agg, uniq


def profile_from_lumisections(rec, avgpu, sigma_mb=SIGMA_MB_NOMINAL, scale=1.0, rel_smear=0.0):
    """Lumi-weighted true-pileup histogram (bins [i, i+1)) of per-lumisection means."""
    mu = avgpu * sigma_mb / SIGMA_MB_BRIL * scale
    hist = np.zeros(NBINS)
    if rel_smear <= 0:
        b = np.clip(mu.astype(int), 0, NBINS - 1)
        np.add.at(hist, b, rec)
    else:
        from scipy.stats import norm
        sig = np.maximum(rel_smear * mu, MIN_SMEAR)
        edges = np.arange(NBINS + 1, dtype=float)
        cdf_prev = norm.cdf((edges[0] - mu) / sig)
        for i in range(NBINS):
            cdf = norm.cdf((edges[i + 1] - mu) / sig)
            hist[i] = np.sum(rec * (cdf - cdf_prev))
            cdf_prev = cdf
    return hist


def data_profile(csv_path=CSV_PATH, sigma_mb=SIGMA_MB_NOMINAL, run_min=skim.RUN_MIN,
                 run_max=skim.RUN_MAX, scale=1.0, rel_smear=0.0, lumisections=None) -> dict:
    """{'hist': 100-bin lumi-weighted true-pileup histogram, 'lumi_pb', 'n_ls', 'mean', 'rms',
    'sigma_mb', 'scale', 'rel_smear'} for the certified lumisections in [run_min, run_max]."""
    rec, avgpu = lumisections if lumisections is not None else _lumisections(csv_path, run_min, run_max)
    hist = profile_from_lumisections(rec, avgpu, sigma_mb, scale, rel_smear)
    x = np.arange(NBINS) + 0.5
    mean = float(np.sum(hist * x) / hist.sum()) if hist.sum() else 0.0
    rms = float(np.sqrt(np.sum(hist * (x - mean) ** 2) / hist.sum())) if hist.sum() else 0.0
    return {"hist": hist, "lumi_pb": float(rec.sum()), "n_ls": int(len(rec)), "mean": mean, "rms": rms,   # n_ls: distinct avgpu values
            "sigma_mb": sigma_mb, "scale": scale, "rel_smear": rel_smear}


def mc_profile(keys=("DY_NLO",), base=None) -> np.ndarray:
    """Summed GenSums.pu_true of the given samples (all UL16 samples share the PU scenario)."""
    total = np.zeros(NBINS)
    for key in keys:
        g = skim.load_gensums(key, base)
        if "pu_true" in g:
            total += np.asarray(g["pu_true"], dtype=float)
    return total


def _weights(data_hist, mc_norm):
    d = data_hist / data_hist.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(mc_norm > 0, d / mc_norm, 0.0)
    w = w / np.sum(w * mc_norm)          # unit mean over the MC profile (keeps the MC yield fixed)
    return np.clip(w, 0.0, 10.0)


def match_npv(mc_ntrue_npv, data_npv, mc_hist, csv_path=CSV_PATH, scales=None, smears=None, min_count=50):
    """Choose (scale, rel_smear) so that the reweighted MC N_PV distribution matches the data.

    mc_ntrue_npv: 2D histogram (NBINS true-pileup bins x N_PV bins) of the selected MC events
    weighted with everything except the pileup weight; data_npv: 1D data histogram with the
    same N_PV bins. Returns the best point and the chi2 grid."""
    scales = np.arange(0.96, 1.16, 0.005) if scales is None else np.asarray(scales)
    smears = np.arange(0.0, 0.31, 0.01) if smears is None else np.asarray(smears)
    rec, avgpu = _lumisections(csv_path)
    mc_norm = np.asarray(mc_hist, dtype=float) / np.sum(mc_hist)
    data_npv = np.asarray(data_npv, dtype=float)
    sel = data_npv > min_count
    grid = np.zeros((len(scales), len(smears)))
    for i, s in enumerate(scales):
        for j, r in enumerate(smears):
            w = _weights(profile_from_lumisections(rec, avgpu, SIGMA_MB_NOMINAL, s, r), mc_norm)
            pred = (w[:, None] * mc_ntrue_npv).sum(axis=0)
            pred = pred * data_npv.sum() / pred.sum()
            grid[i, j] = np.sum((data_npv[sel] - pred[sel]) ** 2 / data_npv[sel])
    i, j = np.unravel_index(np.argmin(grid), grid.shape)
    return {"scale": float(scales[i]), "rel_smear": float(smears[j]), "chi2": float(grid[i, j]),
            "chi2_unmatched": float(grid[np.argmin(abs(scales - 1.0)), np.argmin(abs(smears))]),
            "ndf": int(sel.sum()), "scales": scales.tolist(), "smears": smears.tolist(), "grid": grid.tolist()}


class PileupWeights:
    """w(nTrueInt) for the nominal and the two sigma_mb variations.

    `scale` and `rel_smear` are the N_PV-matching parameters of the data profile (1, 0 = the
    raw per-lumisection CSV profile)."""

    def __init__(self, mc_hist, csv_path=CSV_PATH, scale=1.0, rel_smear=0.0):
        self.mc = np.asarray(mc_hist, dtype=float)
        self.mc_norm = self.mc / self.mc.sum()
        self.scale, self.rel_smear = float(scale), float(rel_smear)
        self.data = {}
        self.weights = {}
        ls = _lumisections(csv_path)
        for name, sigma in (("nominal", SIGMA_MB_NOMINAL), ("up", SIGMA_MB_NOMINAL * (1 + SIGMA_MB_UNC)),
                            ("down", SIGMA_MB_NOMINAL * (1 - SIGMA_MB_UNC))):
            prof = data_profile(csv_path, sigma_mb=sigma, scale=scale, rel_smear=rel_smear, lumisections=ls)
            self.data[name] = prof
            self.weights[name] = _weights(prof["hist"], self.mc_norm)
        # the raw CSV profile, kept for reference and plots
        self.data["csv_raw"] = data_profile(csv_path, lumisections=ls)

    def __call__(self, n_true, variation="nominal"):
        b = np.clip(np.asarray(n_true, dtype=float).astype(int), 0, NBINS - 1)
        return self.weights[variation][b]

    @classmethod
    def from_json(cls, path):
        obj = cls.__new__(cls)
        with open(path) as fh:
            d = json.load(fh)
        obj.mc = np.asarray(d["mc_profile"], dtype=float)
        obj.mc_norm = obj.mc / obj.mc.sum()
        obj.scale, obj.rel_smear = float(d.get("scale", 1.0)), float(d.get("rel_smear", 0.0))
        obj.data = {k: {kk: (np.asarray(vv) if kk == "hist" else vv) for kk, vv in v.items()} for k, v in d["data"].items()}
        obj.weights = {k: np.asarray(v, dtype=float) for k, v in d["weights"].items()}
        return obj

    def to_json(self, path, extra=None):
        with open(path, "w") as fh:
            json.dump({"mc_profile": self.mc.tolist(), "scale": self.scale, "rel_smear": self.rel_smear,
                       "data": {k: {kk: (vv.tolist() if hasattr(vv, "tolist") else vv) for kk, vv in v.items()}
                                for k, v in self.data.items()},
                       "weights": {k: v.tolist() for k, v in self.weights.items()}, **(extra or {})}, fh)
