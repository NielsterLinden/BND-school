"""Pileup reweighting from the CMS Open Data luminosity record (docs/11-mc-weights.md).

No official `puWeights.json` is reachable from this cluster, so the data profile is built
from the per-lumisection table of record 1059 (`pp_2016lumibyls.csv`, brilcalc with the
PHYSICS normtag):

  * only the certified lumisections of runs 278820-284044 (Run2016G+H) enter, weighted by
    the recorded luminosity; the sum of that luminosity reproduces 16393.381 pb^-1 exactly;
  * brilcalc's `avgpu` column is computed with sigma_minbias = 80 mb, while the CMS
    Ultra-Legacy recommendation is 69.2 mb: the profile is `avgpu x 69.2/80` (the ratio of the
    full-2016 mean to the official UL2016 histogram is 0.87, consistent with 69.2/80 = 0.865);
  * the +/-4.6% variation of the minimum-bias cross section (72.4 / 66.0 mb) gives the
    `Pileup` Up/Down weights.

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


def data_profile(csv_path=CSV_PATH, sigma_mb=SIGMA_MB_NOMINAL, run_min=skim.RUN_MIN,
                 run_max=skim.RUN_MAX) -> dict:
    """{'hist': 100-bin lumi-weighted true-pileup histogram, 'lumi_pb', 'n_ls', 'mean'} for the
    certified lumisections in [run_min, run_max]."""
    grl = io.load_grl()
    hist = np.zeros(NBINS)
    lumi = 0.0
    n_ls = 0
    mean_num = 0.0
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
            rec = float(row[6]) * 1e3          # /fb -> /pb
            mu = float(row[7]) * sigma_mb / SIGMA_MB_BRIL
            b = min(int(mu), NBINS - 1)
            hist[b] += rec
            lumi += rec
            mean_num += rec * mu
            n_ls += 1
    return {"hist": hist, "lumi_pb": lumi, "n_ls": n_ls, "mean": mean_num / lumi if lumi else 0.0,
            "sigma_mb": sigma_mb}


def mc_profile(keys=("DY_NLO",), base=None) -> np.ndarray:
    """Summed GenSums.pu_true of the given samples (all UL16 samples share the PU scenario)."""
    total = np.zeros(NBINS)
    for key in keys:
        g = skim.load_gensums(key, base)
        if "pu_true" in g:
            total += np.asarray(g["pu_true"], dtype=float)
    return total


class PileupWeights:
    """w(nTrueInt) for the nominal and the two sigma_mb variations."""

    def __init__(self, mc_hist, csv_path=CSV_PATH):
        self.mc = np.asarray(mc_hist, dtype=float)
        self.mc_norm = self.mc / self.mc.sum()
        self.data = {}
        self.weights = {}
        for name, sigma in (("nominal", SIGMA_MB_NOMINAL), ("up", SIGMA_MB_NOMINAL * (1 + SIGMA_MB_UNC)),
                            ("down", SIGMA_MB_NOMINAL * (1 - SIGMA_MB_UNC))):
            prof = data_profile(csv_path, sigma_mb=sigma)
            d = prof["hist"] / prof["hist"].sum()
            with np.errstate(divide="ignore", invalid="ignore"):
                w = np.where(self.mc_norm > 0, d / self.mc_norm, 0.0)
            # renormalise to unit mean over the MC profile (keeps the MC yield fixed)
            w = w / np.sum(w * self.mc_norm)
            self.data[name] = prof
            self.weights[name] = np.clip(w, 0.0, 10.0)

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
        obj.data = {k: {"hist": np.asarray(v["hist"]), "lumi_pb": v["lumi_pb"], "mean": v["mean"],
                        "sigma_mb": v["sigma_mb"]} for k, v in d["data"].items()}
        obj.weights = {k: np.asarray(v, dtype=float) for k, v in d["weights"].items()}
        return obj

    def to_json(self, path):
        with open(path, "w") as fh:
            json.dump({"mc_profile": self.mc.tolist(),
                       "data": {k: {"hist": v["hist"].tolist(), "lumi_pb": v["lumi_pb"],
                                    "mean": v["mean"], "sigma_mb": v["sigma_mb"]}
                                for k, v in self.data.items()},
                       "weights": {k: v.tolist() for k, v in self.weights.items()}}, fh)
