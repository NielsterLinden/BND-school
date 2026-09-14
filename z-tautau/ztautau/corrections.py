"""TauPOG corrections for the simulation and the pileup reweighting (docs/07-corrections-and-systematics.md).

All numbers come from external/tau_pog_UL2016postVFP.json (written by scripts/step0_external.py from
the public TauPOG repositories) and from the CMS Open Data luminosity table; nothing here needs ROOT.

Per tau, by generator match (NanoAOD `Tau_genPartFlav`):
    5   genuine tau_h            DeepTau VSjet Medium ID SF (per decay mode), energy scale (per DM)
    1,3 prompt e / tau -> e      VSe VVLoose SF (per |eta|)
    2,4 prompt mu / tau -> mu    VSmu VLoose SF (per |eta|)
    0   jet (only in the tau_2 leg of the simulation, see fakes.py)   no ID SF
Every leg (whatever its origin) gets the di-tau trigger leg SF(pT, DM) = eff_data / eff_MC of the
TauPOG fitted efficiencies. Scale-factor uncertainties are uncorrelated between decay modes and
fully correlated between the two legs of an event (one nuisance parameter per DM).
"""

from __future__ import annotations

import csv
import json
from functools import lru_cache

import numpy as np

from . import config, io

POG_JSON = config.EXTERNAL_DIR / "tau_pog_UL2016postVFP.json"


@lru_cache(maxsize=1)
def pog() -> dict:
    return json.loads(POG_JSON.read_text())


def tes_nominal() -> dict:
    return {int(dm): v[0] for dm, v in pog()["tes_dm"].items()}


def tes_uncertainty() -> dict:
    return {int(dm): v[1] for dm, v in pog()["tes_dm"].items()}


def _binned(table, x):
    edges = np.asarray(table["edges"])
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(table["values"]) - 1)
    return np.asarray(table["values"])[idx], np.asarray(table["errors"])[idx]


def id_sf(dm, genflav, abseta, shift_dm: int | None = None, direction: int = 0, kind: str = "id"):
    """Per-tau ID scale factor (numpy arrays in, array out).

    `kind` selects which uncertainty `shift_dm/direction` moves: 'id' (VSjet, one NP per DM),
    'vse', 'vsmu' (shift_dm ignored). direction = +1 / -1 / 0.
    """
    p = pog()
    sf = np.ones(len(dm))
    genuine = genflav == 5
    for d in config.TAU_DMS:
        val, err = p["id_vsjet_dm"]["Medium"][str(d)]
        sel = genuine & (dm == d)
        shift = direction * err if (kind == "id" and shift_dm == d) else 0.0
        sf[sel] = val + shift
    ele = (genflav == 1) | (genflav == 3)
    v, e = _binned(p["vse"]["VVLoose"], abseta)
    sf = np.where(ele, v + (direction * e if kind == "vse" else 0.0), sf)
    muo = (genflav == 2) | (genflav == 4)
    v, e = _binned(p["vsmu"]["VLoose"], abseta)
    sf = np.where(muo, v + (direction * e if kind == "vsmu" else 0.0), sf)
    return sf


def trigger_sf(dm, pt, shift_dm: int | None = None, direction: int = 0):
    """Di-tau trigger leg scale factor SF(pT, DM) with an optional +-1 sigma shift of one DM."""
    t = pog()["trigger_ditau_Medium"]
    grid = np.asarray(t["pt"])
    idx = np.clip(np.searchsorted(grid, pt, side="right") - 1, 0, len(grid) - 1)
    sf = np.ones(len(pt))
    for d in config.TAU_DMS:
        sel = dm == d
        vals = np.asarray(t[str(d)]["sf"])[idx[sel]]
        if shift_dm == d and direction:
            vals = vals + direction * np.asarray(t[str(d)]["sf_err"])[idx[sel]]
        sf[sel] = vals
    return sf


# ------------------------------------------------------------------------------ pileup
# Same construction as z-mumu/zmumu/pileup.py (v2), so the `Pileup` nuisance parameter means the same
# thing in both channels: data profile from the certified lumisections of runs 278820-284044 in the
# Open Data lumi-by-LS table (record 1059), avgpu rescaled from 80 mb to 69.2 mb (UL recommendation),
# +-4.6% on the minimum-bias cross section for the variations; MC profile = genWeight-weighted
# Pileup_nTrueInt of all generated events (GenSums).
SIGMA_MB_BRIL, SIGMA_MB_NOMINAL, SIGMA_MB_UNC = 80.0, 69.2, 0.046
PU_BINS = 100


def data_pileup_profile(sigma_mb=SIGMA_MB_NOMINAL):
    grl = io.load_grl()
    hist = np.zeros(PU_BINS)
    lumi = 0.0
    with open(config.LUMIBYLS_CSV) as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith("#"):
                continue
            run = int(row[0].split(":")[0])
            if run < config.RUN_MIN or run > config.RUN_MAX:
                continue
            ls = int(row[1].split(":")[0])
            blocks = grl.get(run)
            if blocks is None or not any(lo <= ls <= hi for lo, hi in blocks):
                continue
            rec = float(row[6]) * 1e3
            mu = float(row[7]) * sigma_mb / SIGMA_MB_BRIL
            hist[min(int(mu), PU_BINS - 1)] += rec
            lumi += rec
    return hist, lumi


class PileupWeights:
    def __init__(self, mc_profile):
        mc = np.asarray(mc_profile, dtype=float)
        self.mc = mc / mc.sum()
        self.weights, self.data, self.lumi = {}, {}, {}
        for name, sigma in (("nominal", SIGMA_MB_NOMINAL), ("up", SIGMA_MB_NOMINAL * (1 + SIGMA_MB_UNC)),
                            ("down", SIGMA_MB_NOMINAL * (1 - SIGMA_MB_UNC))):
            hist, lumi = data_pileup_profile(sigma)
            d = hist / hist.sum()
            with np.errstate(divide="ignore", invalid="ignore"):
                w = np.where(self.mc > 0, d / self.mc, 0.0)
            w /= np.sum(w * self.mc)
            self.weights[name] = np.clip(w, 0.0, 10.0)
            self.data[name] = d
            self.lumi[name] = lumi

    def __call__(self, n_true, variation="nominal"):
        b = np.clip(np.asarray(n_true).astype(int), 0, PU_BINS - 1)
        return self.weights[variation][b]

    def to_json(self):
        return {"mc": self.mc.tolist(), "data": {k: v.tolist() for k, v in self.data.items()},
                "weights": {k: v.tolist() for k, v in self.weights.items()}, "lumi_pb": self.lumi}
