#!/usr/bin/env python
"""Freeze the Z -> mumu MC corrections: pileup profile/weights, L1 prefiring means, tag-and-probe scale factors,
trigger efficiencies and the Z-peak momentum calibration.

    source fitting/setup.sh && python presentation/data/extract_zmumu_corrections.py [--json PATH] [--check-only]

Reads (read-only): z-mumu/output/v2/tnp/tnp_result.json, tnp/pileup_weights.json, momentum.json, results_v2.json and
presentation/data/zmumu_sr_stack.json (the signal-region event means of every weight, from the new MC pass).
Map means (np.mean over the SF map, as RESULTS_v2.md:128-130 / scripts/v2_6_report.py:141) and signal-region event means
are stored under different names. Writes presentation/data/zmumu_corrections.json.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import DATA_DIR, LUMI_PB, V2, Checker, finalize, load_json, provenance, standard_args  # noqa: E402

TNP = V2 / "tnp" / "tnp_result.json"
PU = V2 / "tnp" / "pileup_weights.json"
MOM = V2 / "momentum.json"
RES = V2 / "results_v2.json"
STACK = DATA_DIR / "zmumu_sr_stack.json"


def _mean(hist):
    h = np.asarray(hist, dtype=float)
    x = np.arange(len(h)) + 0.5
    return float(np.sum(h * x) / h.sum())


def build():
    t, pu, mom, res = load_json(TNP), load_json(PU), load_json(MOM), load_json(RES)
    st = load_json(STACK)
    mf = st["mean_factors"]
    sf_maps = {k: t["sf"][k] for k in ("id", "iso", "antiiso")}
    map_mean = {k: float(np.mean(v)) for k, v in sf_maps.items()}
    out = {
        "provenance": provenance("extract_zmumu_corrections.py", [TNP, PU, MOM, RES, STACK],
                                 anchors=["z-mumu/output/v2/RESULTS_v2.md:123-136", "z-mumu/docs/11-mc-weights.md", "z-mumu/docs/12-tag-and-probe.md",
                                          "z-mumu/docs/14-fit-and-systematics.md"],
                                 naming="map_mean = np.mean over the (pT, |eta|) map (what RESULTS_v2.md quotes); sr_event_mean = weighted mean "
                                        "over the selected signal-region events of the new MC pass (zmumu_sr_stack.json mean_factors)"),
        "lumi_pb": LUMI_PB,
        "pileup": {
            "scale": pu["scale"], "rel_smear": pu["rel_smear"], "sigma_mb": pu["data"]["nominal"]["sigma_mb"], "sigma_mb_unc": 0.046,
            "sigma_mb_bril": 80.0, "bin_edges": list(range(0, 101)), "bin_note": "true pileup (Pileup_nTrueInt) bins [i, i+1), i = 0..99",
            "data_profile": {k: pu["data"][k]["hist"] for k in ("nominal", "up", "down", "csv_raw")},
            "data_profile_unit": "recorded luminosity [pb^-1] per bin (sums to lumi_pb)",
            "data_mean": {k: pu["data"][k]["mean"] for k in ("nominal", "up", "down", "csv_raw")},
            "data_rms": {k: pu["data"][k]["rms"] for k in ("nominal", "up", "down", "csv_raw")},
            "data_lumi_pb": pu["data"]["nominal"]["lumi_pb"], "n_lumisection_classes": pu["data"]["nominal"]["n_ls"],
            "mc_profile": pu["mc_profile"], "mc_profile_note": "genWeight-weighted Pileup_nTrueInt of all generated DY_NLO events (GenSums)",
            "mc_mean": _mean(pu["mc_profile"]), "mc_mean_lower_edge": _mean(pu["mc_profile"]) - 0.5,
            "mean_note": "data_mean and mc_mean use bin centres (i + 0.5, as pileup.py:data_profile); docs/11-mc-weights.md:28 quotes the MC mean "
                         "with the lower-edge convention (21.9 = mc_mean_lower_edge)",
            "data_profile_integral_note": "the smeared nominal/up/down profiles lose ~2e-6 of the luminosity in the Gaussian tails outside the "
                                          "100 bins; csv_raw integrates to lumi_pb exactly",
            "weights": {k: pu["weights"][k] for k in ("nominal", "up", "down")},
            "weights_note": "w(nTrueInt) = data/MC profile ratio, unit mean over the MC profile, clipped at 10 (zmumu/pileup.py:_weights)",
            "npv_match": {k: pu["npv_match"][k] for k in ("scale", "rel_smear", "chi2", "chi2_unmatched", "ndf")},
            "sr_event_mean": {"DYmumu": mf["DYmumu"]["pileup"], "all_mc": mf["all_mc"]["pileup"]},
        },
        "prefiring": {
            "branch": "L1PreFiringWeight_Nom = Muon_Nom x ECAL_Nom (Up/Dn for the systematic)",
            "sr_event_mean": {smp: {"nominal": mf[smp]["prefiring"], "muon": mf[smp]["prefiring_muon"], "ecal": mf[smp]["prefiring_ecal"]}
                              for smp in ("DYmumu", "all_mc")},
            "doc_value": {"nominal": 0.9803, "muon": 0.9818, "ecal": 0.9985, "source": "z-mumu/docs/11-mc-weights.md:31"},
        },
        "scale_factors": {
            "pt_edges": t["pt_edges"], "eta_edges": t["eta_edges"], "shape": "[pT bin][|eta| bin] = (10, 4)",
            "sf": sf_maps, "sf_err": {k: t["sf_err"][k] for k in sf_maps}, "sf_stat": {k: t["sf_stat"][k] for k in sf_maps},
            "sf_syst": {k: t["sf_syst"][k] for k in sf_maps},
            "eff": {f"{k}_{s}": t["eff"][f"{k}_{s}_nominal"] for k in sf_maps for s in ("data", "mc")},
            "eff_err": {f"{k}_{s}": t["eff_err"][f"{k}_{s}_nominal"] for k in sf_maps for s in ("data", "mc")},
            "map_mean": map_mean,
            "map_mean_stat": {k: float(np.mean(t["sf_stat"][k])) for k in sf_maps},
            "map_mean_syst": {k: float(np.mean(t["sf_syst"][k])) for k in sf_maps},
            "sr_event_mean": {smp: {k: mf[smp][k] for k in ("sf_event", "sf_id_pair", "sf_iso_pair", "sf_trigger_event")} for smp in ("DYmumu", "all_mc")},
            "iso_sf_vs_npv": t["meta"]["iso_sf_vs_npv"], "iso_sf_npv_spread": t["meta"]["iso_sf_npv_spread"],
            "definition": "tight ID given loose probe; isolation (pfRelIso04 < 0.15) given tight ID; anti-isolation (0.20 < iso < 1.0) given tight ID; "
                          "T&P fits (Z lineshape + background) in data and DY MC, SF = eff_data / eff_mc",
        },
        "trigger": {
            "pt_edges": t["trig_pt_edges"], "eta_edges": t["eta_edges"], "shape": "[pT bin][|eta| bin] = (14, 4)",
            "eff_data": t["eff"]["trig_data"], "eff_mc": t["eff"]["trig_mc"],
            "eff_err_data": t["eff_err"]["trig_data"], "eff_err_mc": t["eff_err"]["trig_mc"],
            "plateau_data": t["meta"]["trig_plateau_data"], "plateau_mc": t["meta"]["trig_plateau_mc"],
            "n_pairs_data": t["meta"]["n_pairs_data"], "n_pairs_mc": t["meta"]["n_pairs_mc"],
            "event_weight": "[1 - (1 - e1_data)(1 - e2_data)] / [1 - (1 - e1_mc)(1 - e2_mc)] (zmumu/tnp.py:336-341)",
            "path": "HLT_IsoMu24 || HLT_IsoTkMu24, matched to a TrigObj (id 13, filterBits & (2|8), pT > 24) within dR < 0.1",
        },
        "momentum": {
            "eta_edges": mom["eta_edges"], "kappa": mom["kappa"], "kappa_err": mom["kappa_err"], "smear": mom["smear"],
            "fit_window_gev": [80.0, 102.0], "nbins": 88, "m_Z": 91.1876, "Gamma_Z": 2.4952,
            "fits": {j: {s: {k: v for k, v in f.items()} for s, f in mom["fits"][j].items()} for j in mom["fits"]},
            "definition": "kappa_j = 1 + (mu_data - mu_MC)/m_Z and s_j = sqrt(2(sigma_data^2 - sigma_MC^2))/m_Z from Breit-Wigner (x) Gaussian "
                          "fits of the Z peak with both muons in the same |eta| bin (zmumu/momentum.py); MC pT -> pT kappa (1 + s N(0,1))",
        },
    }
    return out


def verify(d, ck: Checker):
    R = "z-mumu/output/v2/RESULTS_v2.md"
    st = load_json(STACK)["mean_factors"]
    m = d["scale_factors"]["map_mean"]
    ck.check("<SF_id> map mean 0.9799", m["id"], 0.9799, 5e-5, R + ":128")
    ck.check("<SF_iso> map mean 1.0059", m["iso"], 1.0059, 5e-5, R + ":129")
    ck.check("<SF_antiiso> map mean 0.8932", m["antiiso"], 0.8932, 5e-5, R + ":130")
    ck.check("map mean stat id/iso/antiiso 0.0015/0.0009/0.0337", [d["scale_factors"]["map_mean_stat"][k] for k in ("id", "iso", "antiiso")], [0.0015, 0.0009, 0.0337], 5e-5, R + ":128-130")
    ck.check("trigger plateau data 0.9071", d["trigger"]["plateau_data"], 0.9071, 5e-5, R + ":123")
    ck.check("trigger plateau MC 0.9227", d["trigger"]["plateau_mc"], 0.9227, 5e-5, R + ":123")
    ck.check("T&P pairs data 20,050,129", d["trigger"]["n_pairs_data"], 20050129, 0, R + ":123")
    ck.check("iso SF vs pileup half-spread 0.0031", d["scale_factors"]["iso_sf_npv_spread"], 0.0031, 5e-5, R + ":124")
    mo = d["momentum"]
    ck.check("kappa 0.99913 / 0.99929 / 0.99876 / 1.00039", mo["kappa"], [0.99913, 0.99929, 0.99876, 1.00039], 5e-6, R + ":136")
    ck.check("kappa_err 0.00050 (floor)", mo["kappa_err"], [0.0005] * 4, 5e-6, R + ":136")
    ck.check("smear 0 / 0.59 / 0.61 / 0.78 %", 100 * np.array(mo["smear"]), [0.0, 0.59, 0.61, 0.78], 5e-3, R + ":136")
    ck.check_true("momentum fits valid", all(mo["fits"][j][s]["valid"] for j in mo["fits"] for s in ("data", "mc")))
    ck.check_true("kappa == 1 + (mu_data - mu_mc)/m_Z", all(abs(mo["kappa"][int(j)] - (1 + (mo["fits"][j]["data"]["mu"] - mo["fits"][j]["mc"]["mu"]) / mo["m_Z"])) < 1e-9 for j in mo["fits"]))
    p = d["pileup"]
    ck.check("pileup profile scale 1.035", p["scale"], 1.035, 1e-9, "output/v2/tnp/pileup_weights.json")
    ck.check("pileup rel_smear 0.09", p["rel_smear"], 0.09, 1e-9, "pileup_weights.json")
    ck.check("pileup sigma_mb 69.2", p["sigma_mb"], 69.2, 1e-9, "zmumu/pileup.py")
    ck.check("lumi_pb 16393.381", d["lumi_pb"], 16393.381, 1e-9, "docs/CONVENTIONS.md")
    ck.check("csv_raw data profile integral == lumi_pb (exact)", sum(p["data_profile"]["csv_raw"]), 16393.381, 1e-7, "pileup_weights.json data.csv_raw", rel=True)
    ck.check("nominal (smeared) data profile integral == lumi_pb within 1e-5 (tails outside 100 bins)", sum(p["data_profile"]["nominal"]), 16393.381, 1e-5,
             "pileup_weights.json data.nominal.hist", rel=True)
    ck.check("data profile lumi_pb", p["data_lumi_pb"], 16393.381, 1e-6, rel=True)
    ck.check("data profile mean 25.52", p["data_mean"]["nominal"], 25.52, 5e-3, "pileup_weights.json")
    ck.check("MC profile mean 21.9 (lower-edge convention of docs/11)", p["mc_mean_lower_edge"], 21.9, 0.05, "z-mumu/docs/11-mc-weights.md:28")
    ck.check("MC profile mean (bin centres) == lower-edge + 0.5", p["mc_mean"], p["mc_mean_lower_edge"] + 0.5, 1e-9)
    ck.check("weights nominal <= 10", max(p["weights"]["nominal"]), 10.0, 0.0) if max(p["weights"]["nominal"]) >= 10 else ck.check_true("weights nominal < 10", True)
    ck.check("npv match chi2 4728.9 / 50", [p["npv_match"]["chi2"], p["npv_match"]["ndf"]], [4728.9, 50], 0.05, "pileup_weights.json")
    ck.check("pileup sr_event_mean == stack mean_factors", [p["sr_event_mean"]["DYmumu"], p["sr_event_mean"]["all_mc"]], [st["DYmumu"]["pileup"], st["all_mc"]["pileup"]], 1e-12)
    pf = d["prefiring"]["sr_event_mean"]["DYmumu"]
    ck.check("<L1 prefiring> DYmumu 0.9803 (muon 0.9818, ECAL 0.9985)", [pf["nominal"], pf["muon"], pf["ecal"]], [0.9803, 0.9818, 0.9985], 5e-4, "z-mumu/docs/11-mc-weights.md:31")
    ck.check("prefiring sr_event_mean == stack", pf["nominal"], st["DYmumu"]["prefiring"], 1e-12)
    se = d["scale_factors"]["sr_event_mean"]["DYmumu"]
    ck.check("sf sr_event_mean == stack", [se["sf_event"], se["sf_id_pair"], se["sf_iso_pair"], se["sf_trigger_event"]],
             [st["DYmumu"][k] for k in ("sf_event", "sf_id_pair", "sf_iso_pair", "sf_trigger_event")], 1e-12)
    ck.check("sf_event == id_pair x iso_pair x trigger (within 0.2%)", se["sf_event"], se["sf_id_pair"] * se["sf_iso_pair"] * se["sf_trigger_event"], 2e-3, rel=True)
    ck.check_true("SF map shapes (10, 4)", all(np.array(d["scale_factors"]["sf"][k]).shape == (10, 4) for k in ("id", "iso", "antiiso")))
    ck.check_true("trigger map shapes (14, 4)", np.array(d["trigger"]["eff_data"]).shape == (14, 4) and np.array(d["trigger"]["eff_mc"]).shape == (14, 4))
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "zmumu_corrections.json")
    args = ap.parse_args()
    finalize(args, build, verify, "zmumu_corrections")


if __name__ == "__main__":
    main()
