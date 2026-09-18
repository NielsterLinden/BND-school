#!/usr/bin/env python
"""Freeze real CMS 2016 Open Data events for the Z -> tautau event displays: one signal-region event with its MET-likelihood
posterior on the (x1, x2) grid, and ten same-sign pairs for the fake-factor tally.

    source fitting/setup.sh && python presentation/data/extract_ztautau_events.py [--seed 20260916] [--json PATH] [--check-only]

Reads (read-only, uproot): the data ntuples data_2016G.root and data_2016H.root of $BND_TAUTAU_CACHE/ntuples_v1 (flat TTree
'ntuple', candidates from 38 GeV: the analysis cuts at 40). The signal region is re-applied from the DeepTau bitmasks
(config.VSJET_BITS: Tight = 32, VVVLoose = 1): SR = OS, both pT > 40, both Tight -> must count 21 160. The posterior of the
chosen event replicates ztautau/mass.py likelihood_mass on its 60 x 60 grid and is checked against the ntuple m_tt.
Events are chosen for DISPLAY with a seeded generator and documented rules -- this is not a physics result.
Writes presentation/data/ztautau_events.json. Anchors: docs/02 cutflow, output/RESULTS.md, docs/04, docs/05.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import (DATASET_TAUTAU, TAUTAU_NTUPLES, VERSION_TAUTAU, ZTAUTAU, Checker, add_ztautau_path,  # noqa: E402
                             finalize, load_json, provenance, standard_args)

M_Z = 91.1876
FF_JSON = ZTAUTAU / "output" / "data" / "fakefactors.json"
BRANCHES = ["run", "lumi", "event", "npv", "era", "os", "njets", "jet1_pt", "dr_tt", "met", "met_x", "met_y", "met_covxx", "met_covxy", "met_covyy",
            "m_vis", "m_col", "m_tt", "pt_tt"] + [f"t{i}_{b}" for i in (1, 2) for b in ("pt", "eta", "phi", "mass", "dm", "charge", "vsjet", "rawvsjet", "genflav")]
SR_RULE = "OS, both pT > 40 GeV, both DeepTau VSjet Tight; display candidates: |m_tt - m_Z| < 5 GeV, one DM in (0, 1) and the other in (10, 11), 10 < pT_miss < 60 GeV"
SS_RULE = ("same-sign pairs (os == 0), both pT > 40, tau2 Tight; tau1 Tight ('passes') or VVVLoose-and-not-Tight ('fails'); pool restricted to N_jets == 0 and "
           "Run2016G; 2 passing + 8 failing drawn with the seeded generator, chosen[0] a failing pair (the display pair); the pool's raw n_T / (n_T + n_L) is stored")
N_BLOCK = 12


def r4(x):
    return float(round(float(x), 4))


def tau(a, i, k):
    dm = int(a[f"t{i}_dm"][k])
    return {"pt": r4(a[f"t{i}_pt"][k]), "eta": r4(a[f"t{i}_eta"][k]), "phi": r4(a[f"t{i}_phi"][k]), "mass": r4(a[f"t{i}_mass"][k]), "dm": dm,
            "charge": int(a[f"t{i}_charge"][k]), "prongs": 1 if dm in (0, 1) else 3, "pi0": dm in (1, 11),
            "vsjet_bits": int(a[f"t{i}_vsjet"][k]), "tight": bool(a[f"t{i}_vsjet"][k] & 32), "vvvloose": bool(a[f"t{i}_vsjet"][k] & 1),
            "raw_vsjet": r4(a[f"t{i}_rawvsjet"][k])}


def head(a, k):
    return {"run": int(a["run"][k]), "lumi": int(a["lumi"][k]), "event": int(a["event"][k]), "era": "G" if int(a["era"][k]) == 0 else "H", "npv": int(a["npv"][k]),
            "njets": int(a["njets"][k]), "jet1_pt": r4(a["jet1_pt"][k]), "dr_tt": r4(a["dr_tt"][k])}


def posterior(ev, n=60, prior_pow=2.0, m_tau=1.77686, cov_min=1.0):
    """likelihood_mass (ztautau/mass.py) for one event, returning the grid, the normalised weights and the posterior median."""
    pt1, phi1, m1 = ev["t1"]["pt"], ev["t1"]["phi"], ev["t1"]["mass"]
    pt2, phi2, m2 = ev["t2"]["pt"], ev["t2"]["phi"], ev["t2"]["mass"]
    u = (np.arange(n) + 0.5) / n
    xmin1 = float(np.clip((m1 / m_tau) ** 2, 1e-4, 0.99))
    xmin2 = float(np.clip((m2 / m_tau) ** 2, 1e-4, 0.99))
    x1 = xmin1 + (1 - xmin1) * u
    x2 = xmin2 + (1 - xmin2) * u
    a1, a2 = 1.0 / x1 - 1.0, 1.0 / x2 - 1.0
    nux = a1[:, None] * pt1 * np.cos(phi1) + a2[None, :] * pt2 * np.cos(phi2)
    nuy = a1[:, None] * pt1 * np.sin(phi1) + a2[None, :] * pt2 * np.sin(phi2)
    rx, ry = ev["met"]["x"] - nux, ev["met"]["y"] - nuy
    cxx, cyy, cxy = max(ev["met"]["covxx"], cov_min), max(ev["met"]["covyy"], cov_min), ev["met"]["covxy"]
    det = cxx * cyy - cxy ** 2
    if det <= 0:
        cxy, det = 0.0, cxx * cyy
    ixx, iyy, ixy = cyy / det, cxx / det, -cxy / det
    chi2 = ixx * rx ** 2 + iyy * ry ** 2 + 2 * ixy * rx * ry
    chi2 -= chi2.min()
    mgrid = ev["m_vis"] / np.sqrt(x1[:, None] * x2[None, :])
    like = np.exp(-0.5 * chi2) * mgrid ** (-prior_pow)
    lf, mf = like.ravel(), mgrid.ravel()
    o = np.argsort(mf)
    cdf = np.cumsum(lf[o]); cdf /= cdf[-1]
    med = float(mf[o][min(int((cdf < 0.5).sum()), len(o) - 1)])
    w = like / like.sum()
    b = n // N_BLOCK
    blocks = w.reshape(N_BLOCK, b, N_BLOCK, b).sum(axis=(1, 3))
    e1 = (xmin1 + (1 - xmin1) * np.arange(N_BLOCK + 1) / N_BLOCK).tolist()
    e2 = (xmin2 + (1 - xmin2) * np.arange(N_BLOCK + 1) / N_BLOCK).tolist()
    imax = np.unravel_index(int(np.argmax(w)), w.shape)
    chi2min_idx = np.unravel_index(int(np.argmin(chi2)), chi2.shape)
    return {"x1_edges": e1, "x2_edges": e2, "weights": blocks.tolist(), "x_lo": [xmin1, xmin2], "m_median": med, "grid_n": n,
            "block": b, "mode_x": [float(x1[imax[0]]), float(x2[imax[1]])], "mode_m": float(mgrid[imax]),
            "chi2min_x": [float(x1[chi2min_idx[0]]), float(x2[chi2min_idx[1]])], "m_at_chi2min": float(mgrid[chi2min_idx]),
            "m_block_median": (ev["m_vis"] / np.sqrt(np.outer(0.5 * (np.array(e1[:-1]) + np.array(e1[1:])), 0.5 * (np.array(e2[:-1]) + np.array(e2[1:]))))).tolist(),
            "note": "weights = the 60 x 60 posterior exp(-chi2/2) m^-2 summed into 12 x 12 blocks of 5 x 5 grid points, normalised to 1; "
                    "x axes are the visible-energy fractions on [m_vis,i^2 / m_tau^2, 1]; m_block_median = m_vis / sqrt(x1 x2) at the block centres"}


def load_data():
    import uproot
    parts = []
    for key in ("data_2016G", "data_2016H"):
        with uproot.open(TAUTAU_NTUPLES / f"{key}.root") as f:
            parts.append(f["ntuple"].arrays(BRANCHES, library="np"))
    return {k: np.concatenate([p[k] for p in parts]) for k in BRANCHES}


def build(seed):
    a = load_data()
    rng = np.random.default_rng(seed)
    t1T, t2T = (a["t1_vsjet"] & 32) > 0, (a["t2_vsjet"] & 32) > 0
    t1L = ((a["t1_vsjet"] & 1) > 0) & ~t1T
    pt = (a["t1_pt"] > 40.0) & (a["t2_pt"] > 40.0)
    os_ = a["os"].astype(bool)
    sr = pt & os_ & t1T & t2T
    n_sr = int(sr.sum())
    n_sr_era = {"G": int((sr & (a["era"] == 0)).sum()), "H": int((sr & (a["era"] == 1)).sum())}
    # ---- the display event -------------------------------------------------------------------------------------------
    near = np.abs(a["m_tt"] - M_Z) < 5.0
    dm1, dm2 = a["t1_dm"], a["t2_dm"]
    mixed = (np.isin(dm1, [0, 1]) & np.isin(dm2, [10, 11])) | (np.isin(dm1, [10, 11]) & np.isin(dm2, [0, 1]))
    metw = (a["met"] > 10.0) & (a["met"] < 60.0)
    cand = np.where(sr & near & mixed & metw)[0]
    relaxed = False
    if len(cand) == 0:
        cand = np.where(sr & near & metw)[0]
        relaxed = True
    k = int(rng.choice(cand))
    ev = {**head(a, k), "t1": tau(a, 1, k), "t2": tau(a, 2, k),
          "met": {"pt": r4(a["met"][k]), "phi": r4(np.arctan2(a["met_y"][k], a["met_x"][k])), "x": r4(a["met_x"][k]), "y": r4(a["met_y"][k]),
                  "covxx": r4(a["met_covxx"][k]), "covxy": r4(a["met_covxy"][k]), "covyy": r4(a["met_covyy"][k])},
          "m_vis": r4(a["m_vis"][k]), "m_col": (r4(a["m_col"][k]) if a["m_col"][k] > 0 else None), "m_tt": r4(a["m_tt"][k]), "pt_tt": r4(a["pt_tt"][k])}
    ev_full = {"t1": {"pt": float(a["t1_pt"][k]), "phi": float(a["t1_phi"][k]), "mass": float(a["t1_mass"][k])},
               "t2": {"pt": float(a["t2_pt"][k]), "phi": float(a["t2_phi"][k]), "mass": float(a["t2_mass"][k])},
               "met": {"x": float(a["met_x"][k]), "y": float(a["met_y"][k]), "covxx": float(a["met_covxx"][k]), "covxy": float(a["met_covxy"][k]), "covyy": float(a["met_covyy"][k])},
               "m_vis": float(a["m_vis"][k])}
    ev["posterior"] = posterior(ev_full)
    ev["m_tt_ntuple"] = float(a["m_tt"][k])
    # cross-check with the channel's own function (same inputs)
    add_ztautau_path()
    from ztautau import mass
    one = lambda v: np.array([v], dtype=np.float64)
    m_chan = float(mass.likelihood_mass((one(ev_full["t1"]["pt"]), one(a["t1_eta"][k]), one(ev_full["t1"]["phi"]), one(ev_full["t1"]["mass"])),
                                        (one(ev_full["t2"]["pt"]), one(a["t2_eta"][k]), one(ev_full["t2"]["phi"]), one(ev_full["t2"]["mass"])),
                                        one(ev_full["met"]["x"]), one(ev_full["met"]["y"]), one(ev_full["met"]["covxx"]), one(ev_full["met"]["covxy"]), one(ev_full["met"]["covyy"]))[0])
    ev["m_tt_recomputed_channel_code"] = m_chan
    # ---- same-sign pairs for the tally ------------------------------------------------------------------------------
    ss_base = pt & ~os_ & t2T & (a["njets"] == 0) & (a["era"] == 0)
    pool_T, pool_L = np.where(ss_base & t1T)[0], np.where(ss_base & t1L)[0]
    pick_T = rng.choice(pool_T, 2, replace=False)
    pick_L = rng.choice(pool_L, 8, replace=False)
    order = [("fail", pick_L[0])] + [t for t in sorted([("pass", i) for i in pick_T] + [("fail", i) for i in pick_L[1:]], key=lambda _: rng.random())]
    pairs = []
    for tag, i in order:
        i = int(i)
        pairs.append({**head(a, i), "t1": tau(a, 1, i), "t2": tau(a, 2, i), "passes": tag == "pass", "os": bool(os_[i]),
                      "met": r4(a["met"][i]), "m_vis": r4(a["m_vis"][i]), "m_tt": r4(a["m_tt"][i])})
    assert not pairs[0]["passes"] and sum(p["passes"] for p in pairs) == 2 and len(pairs) == 10
    # fake factor of the display pair's bin
    ff = load_json(FF_JSON)["mcsub"]["ff"]
    disp = pairs[0]
    ie = 0 if disp["era"] == "G" else 1
    idm = ff["dms"].index(disp["t1"]["dm"])
    inj = int(np.clip(np.searchsorted(np.asarray(ff["njet_bins"]), disp["njets"], side="right") - 1, 0, len(ff["njet_bins"]) - 1))
    ptb = np.asarray(ff["pt_bins"])
    ipt = int(np.clip(np.searchsorted(ptb, disp["t1"]["pt"], side="right") - 1, 0, len(ptb) - 2))
    ceta, cpt2 = ff["closure"]["eta"], ff["closure"]["pt2"]
    ieta = int(np.clip(np.searchsorted(np.asarray(ceta["edges"]), abs(disp["t1"]["eta"]), side="right") - 1, 0, len(ceta["edges"]) - 2))
    ipt2 = int(np.clip(np.searchsorted(np.asarray(cpt2["edges"]), disp["t2"]["pt"], side="right") - 1, 0, len(cpt2["edges"]) - 2))
    f_val, f_err = ff["ff"][ie][idm][inj][ipt], ff["err"][ie][idm][inj][ipt]
    ff_ref = {"era": disp["era"], "dm": disp["t1"]["dm"], "njets": disp["njets"], "njet_bin": ff["njet_bins"][inj], "pt": disp["t1"]["pt"],
              "pt_bin": [float(ptb[ipt]), float(ptb[ipt + 1])], "value": f_val, "err": f_err,
              "num": ff["num"][ie][idm][inj][ipt], "den": ff["den"][ie][idm][inj][ipt], "num_mc": ff["num_mc"][ie][idm][inj][ipt], "den_mc": ff["den_mc"][ie][idm][inj][ipt],
              "closure_eta": ceta["values"][ieta], "closure_pt2": cpt2["values"][ipt2], "value_corrected": f_val * ceta["values"][ieta] * cpt2["values"][ipt2],
              "source": "output/data/fakefactors.json mcsub.ff (era, DM, N_jets, pT bin of the display pair's tau1) and its closure corrections"}
    out = {
        "provenance": provenance("extract_ztautau_events.py", [TAUTAU_NTUPLES / "data_2016G.root", TAUTAU_NTUPLES / "data_2016H.root", FF_JSON],
                                 dataset=DATASET_TAUTAU, version=VERSION_TAUTAU, seed=seed,
                                 rules={"sr": SR_RULE, "ss_ff": SS_RULE}, dm_condition_relaxed=relaxed,
                                 posterior="ztautau/mass.py likelihood_mass re-implemented for one event (60 x 60 grid, MET covariance floor 1 GeV^2, 1/m^2 prior), "
                                           "checked against the ntuple m_tt and against mass.likelihood_mass itself",
                                 anchors=["z-tautau/docs/02-selection.md cutflow (21 160)", "z-tautau/output/RESULTS.md", "z-tautau/docs/04-ditau-mass.md", "z-tautau/docs/05-fake-factors.md"],
                                 display_note="chosen for display with a seeded generator; not a physics result"),
        "m_Z": M_Z,
        "sr": {"selection": SR_RULE, "n_sr": n_sr, "n_sr_per_era": n_sr_era, "n_candidates": int(len(cand)), "dm_condition_relaxed": relaxed, "chosen": [ev]},
        "ss_ff": {"selection": SS_RULE, "pool": {"n_pass": int(len(pool_T)), "n_fail": int(len(pool_L)), "raw_pass_fraction": float(len(pool_T) / (len(pool_T) + len(pool_L))),
                                                "n_T_over_n_L": float(len(pool_T) / len(pool_L))},
                  "pool_all": {"n_pass": int((pt & ~os_ & t2T & t1T).sum()), "n_fail": int((pt & ~os_ & t2T & t1L).sum())},
                  "chosen": pairs, "display": {"n_pass": 2, "n_fail": 8, "ratio": 0.25}},
        "ff_reference": ff_ref,
    }
    return out


def verify(d, ck: Checker):
    ck.check("n_sr 21160", d["sr"]["n_sr"], 21160, 0, "z-tautau/docs/02 cutflow / RESULTS.md")
    ck.check("n_sr per era 10572 / 10588", [d["sr"]["n_sr_per_era"]["G"], d["sr"]["n_sr_per_era"]["H"]], [10572, 10588], 0, "docs/02 cutflow")
    ev = d["sr"]["chosen"][0]
    ck.check_true("exactly one chosen SR event", len(d["sr"]["chosen"]) == 1)
    ck.check("posterior median reproduces the ntuple m_tt within 1.0 GeV", ev["posterior"]["m_median"], ev["m_tt_ntuple"], 1.0, "ztautau/mass.py likelihood_mass")
    ck.check("posterior median == mass.likelihood_mass (channel code) within 1e-3", ev["posterior"]["m_median"], ev["m_tt_recomputed_channel_code"], 1e-3, "ztautau/mass.py")
    ck.check("event |m_tt - m_Z| < 5", abs(ev["m_tt"] - M_Z) < 5.0, True, 0, "display rule")
    ck.check("event 10 < MET < 60", 10 < ev["met"]["pt"] < 60, True, 0, "display rule")
    ck.check_true("event OS, both Tight, pT > 40", ev["t1"]["charge"] * ev["t2"]["charge"] < 0 and ev["t1"]["tight"] and ev["t2"]["tight"] and ev["t1"]["pt"] > 40 and ev["t2"]["pt"] > 40)
    ck.check_true("event one 1-prong and one 3-prong (unless relaxed)", d["sr"]["dm_condition_relaxed"] or {ev["t1"]["prongs"], ev["t2"]["prongs"]} == {1, 3},
                  detail=f"DM {ev['t1']['dm']} / {ev['t2']['dm']}")
    ck.check("posterior weights sum to 1", np.sum(ev["posterior"]["weights"]), 1.0, 1e-9, "consistency")
    ck.check_true("posterior 12 x 12 with 13 + 13 edges", np.asarray(ev["posterior"]["weights"]).shape == (12, 12) and len(ev["posterior"]["x1_edges"]) == 13 and len(ev["posterior"]["x2_edges"]) == 13)
    ck.check("x_lo == first edges", ev["posterior"]["x_lo"], [ev["posterior"]["x1_edges"][0], ev["posterior"]["x2_edges"][0]], 0, "consistency")
    ck.check("m_vis == inv. mass of the two tau_h (from pt, eta, phi, mass)", ev["m_vis"], _mvis(ev), 0.05, "consistency")
    ck.check("MET pt == hypot(x, y)", ev["met"]["pt"], float(np.hypot(ev["met"]["x"], ev["met"]["y"])), 0.01, "consistency")
    ck.check_true("m_vis < m_tt for the event (the neutrinos add mass)", ev["m_vis"] < ev["m_tt"])
    p = d["ss_ff"]["chosen"]
    ck.check_true("10 SS pairs, 2 pass, 8 fail, display pair fails", len(p) == 10 and sum(x["passes"] for x in p) == 2 and not p[0]["passes"])
    ck.check_true("all SS pairs same-sign, tau2 Tight, pT > 40, N_jets 0, era G", all((not x["os"]) and x["t2"]["tight"] and x["t1"]["pt"] > 40 and x["t2"]["pt"] > 40
                                                                                     and x["njets"] == 0 and x["era"] == "G" and x["t1"]["charge"] == x["t2"]["charge"] for x in p))
    ck.check_true("pass <=> tau1 Tight; fail <=> VVVLoose and not Tight", all((x["t1"]["tight"] == x["passes"]) and (x["passes"] or x["t1"]["vvvloose"]) for x in p))
    ck.check_true("10 distinct events", len({(x["run"], x["event"]) for x in p}) == 10)
    ck.check("SS pool pass fraction ~ FF magnitude (0.05 .. 0.3)", 0.05 < d["ss_ff"]["pool"]["n_T_over_n_L"] < 0.3, True, 0, "docs/05")
    ck.check("SS pool sizes: SS_T + SS_L (all eras, all N_jets) == 9554 + 86929", [d["ss_ff"]["pool_all"]["n_pass"], d["ss_ff"]["pool_all"]["n_fail"]], [9554, 86929], 0,
             "docs/05 contamination table / fakefactors.json")
    fr = d["ff_reference"]
    ck.check_true("ff_reference bin matches the display pair", fr["era"] == p[0]["era"] and fr["dm"] == p[0]["t1"]["dm"] and fr["pt_bin"][0] <= p[0]["t1"]["pt"] < fr["pt_bin"][1])
    ck.check("ff_reference value == num-num_mc / den-den_mc", fr["value"], max(fr["num"] - fr["num_mc"], 0) / max(fr["den"] - fr["den_mc"], 1e-9), 1e-9, "fakes.measure", rel=True)
    ck.check_true("ff_reference in (0.02, 0.5)", 0.02 < fr["value"] < 0.5, detail=f"{fr['value']:.4f}")
    return ck


def _mvis(ev):
    p = []
    for t in (ev["t1"], ev["t2"]):
        pt, eta, phi, m = t["pt"], t["eta"], t["phi"], t["mass"]
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        p.append(np.array([np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m ** 2), px, py, pz]))
    s = p[0] + p[1]
    return float(np.sqrt(max(s[0] ** 2 - s[1] ** 2 - s[2] ** 2 - s[3] ** 2, 0.0)))


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_events.json")
    ap.add_argument("--seed", type=int, default=20260916)
    args = ap.parse_args()
    finalize(args, lambda: build(args.seed), verify, "ztautau_events")


if __name__ == "__main__":
    main()
