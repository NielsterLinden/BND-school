#!/usr/bin/env python
"""Freeze the numbers of the four-channel Z -> tautau measurement (z-tautau v4: tau_h tau_h + mu tau_h + e tau_h + e mu).

    source fitting/setup.sh && python presentation/data/extract_ztautau_v4.py [--json PATH] [--check-only]

Reads (read-only):
  z-tautau/fit/results/ztautau/Plots/<region>_{prefit,postfit}.yaml   TRExFitter per-bin dumps of the nominal four-channel fit
        (git-ignored; reached through a link to the main checkout). Written after cf54006, which put the signal back into the
        plotted total, so data / prediction is the real one.
  z-tautau/fit/results/ztautau_fit_result.json, ztautau_<ch>_fixedid_fit_result.json, ztautau_{emu,taulep}_fit_result.json
  z-tautau/output/results.json (tau_id_sf_pog, fakes, tension), z-tautau/output/data/fakes_{mutau,etau}.json (per-process
        fake factors and application-region fractions)
Writes presentation/data/ztautau_v4.json. Anchors: z-tautau/handoff.md (RESULT block), output/RESULTS.md, docs/10-v4-plan.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import LUMI_PB, ZTAUTAU, Checker, finalize, load_json, provenance, standard_args  # noqa: E402

PLOTS = ZTAUTAU / "fit" / "results" / "ztautau" / "Plots"
FIT_DIR = ZTAUTAU / "fit" / "results"
RES_JSON = ZTAUTAU / "output" / "results.json"
FAKES = {ch: ZTAUTAU / "output" / "data" / f"fakes_{ch}.json" for ch in ("mutau", "etau")}
DMS = ("0", "1", "10", "11")
REGIONS = {"tautau": ["tautau_SR0", "tautau_SR1", "tautau_SR2"],
           "mutau": [f"mutau_SR_dm{d}" for d in DMS], "etau": [f"etau_SR_dm{d}" for d in DMS],
           "emu": ["emu_SR"], "emu_CRtt": ["emu_CRtt"]}
GROUP_OF_TITLE = {"Z/#gamma*#rightarrow#tau#tau (60-120)": "DYtautau", "Z/#gamma*#rightarrow#tau#tau (outside)": "DYtautau_out",
                  "jet#rightarrow#tau_{h} / multijet": "Fakes", "t#bar{t}": "TTbar"}
GROUPS = ("Fakes", "rest", "TTbar", "DYtautau_out", "DYtautau")        # stack order, bottom-up
DATASET = "CMS 2016 Open Data Run2016G+H NanoAODv9: Tau, SingleMuon, SingleElectron, MuonEG"
VERSION = "v4 four channels, DeepTau Tight, tau_h ID SFs free per decay mode"


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def read_yaml(region: str, tag: str) -> dict:
    import yaml
    with open(PLOTS / f"{region}_{tag}.yaml") as fh:
        y = yaml.safe_load(fh)
    edges = [float(v) for v in y["Figure"][0]["BinEdges"]]
    nb = len(edges) - 1
    groups = {g: np.zeros(nb) for g in GROUPS}
    titles = []
    for s in y["Samples"]:
        vals = np.array([_f(v) for v in s["Yield"]])
        assert len(vals) == nb, (region, tag, s["Name"], len(vals), nb)
        groups[GROUP_OF_TITLE.get(s["Name"], "rest")] += vals
        titles.append(s["Name"])
    tot = y["Total"][0]
    return {"edges": edges, "data": [int(round(_f(v))) for v in y["Data"][0]["Yield"]],
            "groups": {g: v.tolist() for g, v in groups.items()},
            "total": [_f(v) for v in tot["Yield"]], "unc_up": [_f(v) for v in tot["UncertaintyUp"]],
            "unc_down": [_f(v) for v in tot["UncertaintyDown"]], "titles": titles}


def summed(parts: list[dict]) -> dict:
    e = parts[0]["edges"]
    assert all(p["edges"] == e for p in parts)
    out = {"edges": e, "data": np.sum([p["data"] for p in parts], axis=0).astype(int).tolist(),
           "groups": {g: np.sum([p["groups"][g] for p in parts], axis=0).tolist() for g in GROUPS},
           "total": np.sum([p["total"] for p in parts], axis=0).tolist(),
           # the uncertainty band of a sum of regions is not the sum of the bands (correlations): not stored
           }
    out["n_data"] = int(sum(out["data"]))
    out["n_pred"] = float(sum(out["total"]))
    out["n"] = {g: float(sum(v)) for g, v in out["groups"].items()}
    return out


def lephad_fakes(ch: str) -> dict:
    fk = load_json(FAKES[ch])
    res = load_json(RES_JSON)["fakes"][ch]
    T = fk["tables"]

    def incl(key, subtract_mc=True):
        t = T[key]
        num, den = np.asarray(t["num"], float), np.asarray(t["den"], float)
        if subtract_mc and "num_mc" in t:
            num = num - np.asarray(t["num_mc"], float)
            den = den - np.asarray(t["den_mc"], float)
        return float(num.sum() / den.sum())

    cnt = fk["fractions"]["counts"]
    n = {p: float(np.asarray(cnt[p], float).sum()) for p in ("qcd", "w", "tt", "other", "data")}
    jet = n["qcd"] + n["w"] + n["tt"]
    cl = res["closure_ss"]
    return {
        "ff_inclusive": {"qcd": incl("qcd"), "w_os": incl("w"), "w_ss": incl("w_ss"), "tt": incl("tt", subtract_mc=False)},
        "ff_note": "inclusive over DM x N_jets x pT: (N_T - MC) / (N_L - MC) summed over the table bins (tt: simulation, no subtraction); "
                   "the fit applies the binned tables",
        "ar_counts": n, "ar_fractions": {p: n[p] / jet for p in ("qcd", "w", "tt")},
        "ar_fraction_note": "application region (tau_h VVVLoose & !Tight, m_T < 40 GeV): multijet = data - simulation, W / tt from simulation; "
                            "fractions of the jet fakes, summed over N_jets and the two m_T bins",
        "osss": {"C": res["osss"]["C"], "stat": res["osss"]["stat"]},
        "closure_ss": {"edges": cl["edges"], "obs": cl["obs_hist"], "pred": cl["pred_hist"], "ratio": cl["ratio"], "stat": cl["stat"]},
        "sr_fakes": res["sr_fakes"], "n_data": dict(res["n_data"]),
        "regions_note": "QCD_T/L: same-sign, isolated lepton (multijet DR); W_T/L: OS, m_T > 70 GeV, no b jet (W DR); WSS_T/L: same, same sign",
    }


def build():
    fit = load_json(FIT_DIR / "ztautau_fit_result.json")
    res = load_json(RES_JSON)
    regions = {tag: {} for tag in ("prefit", "postfit")}
    channels = {tag: {} for tag in ("prefit", "postfit")}
    for tag in ("prefit", "postfit"):
        for ch, rs in REGIONS.items():
            parts = []
            for r in rs:
                regions[tag][r] = read_yaml(r, tag)
                parts.append(regions[tag][r])
            channels[tag][ch] = summed(parts)
    per_channel = {}
    for ch in ("tautau", "mutau", "etau", "emu"):
        f = load_json(FIT_DIR / f"ztautau_{ch}_fixedid_fit_result.json")
        s = f["sigma_60_120_pb"]
        per_channel[ch] = {"mu": f["poi_value"], "mu_up": f["poi_err_up"], "mu_down": f["poi_err_down"],
                           "sigma": s["value"], "sigma_up": s["up"], "sigma_down": s["down"]}
    subs = {}
    for key, job in (("emu", "ztautau_emu"), ("taus", "ztautau_taulep")):
        f = load_json(FIT_DIR / f"{job}_fit_result.json")
        s = f["sigma_60_120_pb"]
        subs[key] = {"mu": f["poi_value"], "mu_up": f["poi_err_up"], "mu_down": f["poi_err_down"],
                     "sigma": s["value"], "sigma_up": s["up"], "sigma_down": s["down"], "tau_id_fixed": f.get("tau_id_fixed")}
    s = fit["sigma_60_120_pb"]
    sources = [PLOTS / f"{r}_{t}.yaml" for t in ("prefit", "postfit") for rs in REGIONS.values() for r in rs]
    sources += [FIT_DIR / "ztautau_fit_result.json", RES_JSON, *FAKES.values()]
    sources += [FIT_DIR / f"ztautau_{ch}_fixedid_fit_result.json" for ch in per_channel]
    return {
        "provenance": provenance("extract_ztautau_v4.py", sources, dataset=DATASET, version=VERSION,
                                 fit_tool="TRExFitter v1.8.0, simultaneous binned profile-likelihood fit of m_tautau in 13 regions "
                                          "(tautau_SR0/1/2, mutau/etau_SR_dm0/1/10/11, emu_SR, emu_CRtt)",
                                 anchors=["z-tautau/handoff.md RESULT block", "z-tautau/output/RESULTS.md", "z-tautau/docs/10-v4-plan.md"],
                                 groups="Fakes = jet->tau_h / multijet (data driven); TTbar; DYtautau = Z/gamma*->tautau 60<m<120 (signal); "
                                        "DYtautau_out = the rest of Z/gamma*->tautau; rest = Z->ee/mumu, low-mass DY, single t, W+jets, WW/WZ/ZZ",
                                 plot_fix="yaml written 17 Sep 13:4x after cf54006 (NOSIG removed: the signal is in Total)"),
        "groups": list(GROUPS),
        "regions": regions,
        "channels": channels,
        "fit": {"mu": fit["poi_value"], "mu_up": fit["poi_err_up"], "mu_down": fit["poi_err_down"],
                "mu_stat": res["fit"]["combined"]["mu_stat"], "gof_p": fit["gof"]["gof_probability"],
                "sigma": s["value"], "sigma_up": s["up"], "sigma_down": s["down"], "sigma_stat": s["stat"],
                "pred": fit["prediction"]["sigma_tautau_60_120_pb"],
                "mu_ttbar": fit["mu_ttbar"],
                "tau_id_sf": {d: {"value": fit["tau_id_sf"][f"DM{d}"]["value"], "err": fit["tau_id_sf"][f"DM{d}"]["err_up"],
                                  "pog": fit["tau_id_sf"][f"DM{d}"]["pog"][0], "pog_err": fit["tau_id_sf"][f"DM{d}"]["pog"][1]} for d in DMS},
                "tau_es": {d: dict(fit["tau_es"][f"DM{d}"]) for d in DMS},
                "lumi_pb": fit["lumi_pb"]},
        "per_channel_fixedid": per_channel,
        "sub_measurements": subs,
        "tension": dict(res["tension_taulep_vs_emu"]),
        "lephad": {ch: lephad_fakes(ch) for ch in ("mutau", "etau")},
        "emu": {"sr_multijet": res["fakes"]["emu"]["sr_multijet"],
                "osss_dr_edges": res["fakes"]["emu"]["osss"]["dr_edges"], "osss_ratio": res["fakes"]["emu"]["osss"]["ratio"],
                "cr_definition": "D_zeta < -40 GeV and p_T^miss > 80 GeV (b inclusive)", "sr_definition": "D_zeta > -20 GeV, no b jet"},
    }


def verify(d, ck: Checker):
    H = "z-tautau/handoff.md RESULT"
    R = "z-tautau/output/RESULTS.md"
    f = d["fit"]
    # since 17 Sep 2026 the uncertainty contains the tau_h ID scale-factor pT dependence (TauIDpT_tautau); without it
    # (z-tautau job ztautau_flatsf) the same fit gives +0.037 -0.036, i.e. 1981 +73 -70 pb
    ck.check("mu_Z 1.019 +0.079 -0.070", [f["mu"], f["mu_up"], f["mu_down"]], [1.019, 0.079, 0.070], 5e-4, H)
    ck.check("sigma(60-120) 1981 +153 -136, stat 9", [f["sigma"], f["sigma_up"], f["sigma_down"], f["sigma_stat"]], [1981, 153, 136, 9], 0.5, H)
    ck.check("prediction 1944.9", f["pred"], 1944.9, 0.05, H)
    ck.check("mu_ttbar 1.11 +- 0.04", [f["mu_ttbar"]["value"], f["mu_ttbar"]["err_up"]], [1.111, 0.038], 5e-4, R)
    ck.check("GoF p 0.15", f["gof_p"], 0.154, 5e-4, R)
    for dm, v, e, pog, pe in (("0", 0.989, 0.040, 0.90, 0.13), ("1", 0.963, 0.034, 0.89, 0.05), ("10", 0.896, 0.035, 0.94, 0.15),
                              ("11", 0.795, 0.050, 0.81, 0.15)):
        t = f["tau_id_sf"][dm]
        ck.check(f"tau ID SF DM{dm} {v} +- {e} (POG {pog} +- {pe})", [t["value"], t["err"], t["pog"], t["pog_err"]], [v, e, pog, pe],
                 5e-3, H)
    ck.check("lumi", f["lumi_pb"], LUMI_PB, 1e-6, "docs/CONVENTIONS.md")
    for ch, mu, s in (("tautau", 1.043, 2029), ("mutau", 1.048, 2038), ("etau", 1.016, 1976), ("emu", 0.960, 1868)):
        p = d["per_channel_fixedid"][ch]
        ck.check(f"{ch} alone (POG SFs) sigma {s}", p["sigma"], s, 0.5, R)
        ck.check(f"{ch} alone mu {mu}", p["mu"], mu, 5e-4, R)
    ck.check("emu alone 0.960, taus alone 1.203", [d["sub_measurements"]["emu"]["mu"], d["sub_measurements"]["taus"]["mu"]], [0.960, 1.203], 5e-4, R)
    pre, post = d["channels"]["prefit"], d["channels"]["postfit"]
    ck.check("data tautau / mutau / etau / emu SR / emu CR", [pre[c]["n_data"] for c in ("tautau", "mutau", "etau", "emu", "emu_CRtt")],
             [21160, 11244 + 23866 + 10795 + 2611, 2639 + 10106 + 4940 + 1211, 79300, 45931], 0, R + " prefit yields table")
    ck.check("tautau fakes 13592", pre["tautau"]["n"]["Fakes"], 13592, 1.0, R)
    # the TRExFitter prefit applies the NormFactors at their nominal (= POG) values; RESULTS.md lists the templates without
    # them: for tau_h tau_h the ratio is ~ SF^2 (0.81), for e mu it is 1
    ck.check("tautau prefit DYtautau / RESULTS.md templates (5174) ~ POG SF^2", pre["tautau"]["n"]["DYtautau"] / (795 + 1630 + 2749), 0.81, 0.05,
             R + " x the nominal TauIDSF NormFactors")
    ck.check("mutau SR fakes 11048", pre["mutau"]["n"]["Fakes"], 11048, 1.5, R + " (mutau fakes: SR fakes 11048)")
    ck.check("etau SR fakes 5012", pre["etau"]["n"]["Fakes"], 5012, 1.5, R)
    ck.check("emu SR multijet 15241", pre["emu"]["n"]["Fakes"], 15241, 1.0, R)
    ck.check("emu SR DYtautau 52985", pre["emu"]["n"]["DYtautau"], 52985, 1.0, R)
    ck.check("emu CR ttbar 37650", pre["emu_CRtt"]["n"]["TTbar"], 37650, 1.0, R)
    ck.check("emu CR ttbar postfit / prefit ~ mu_ttbar", post["emu_CRtt"]["n"]["TTbar"] / pre["emu_CRtt"]["n"]["TTbar"], f["mu_ttbar"]["value"], 0.03,
             "mu_ttbar x the post-fit shifts of the other NPs (EmuTrigger, TopPt, BTag): 1.054, not 1.111", soft=True)
    ck.check_true("postfit total within 2 % of data in every channel",
                  all(abs(post[c]["n_pred"] / post[c]["n_data"] - 1) < 0.02 for c in post),
                  detail=str({c: round(post[c]["n_pred"] / post[c]["n_data"], 4) for c in post}))
    for tag in ("prefit", "postfit"):
        for r, v in d["regions"][tag].items():
            ck.check(f"{tag} {r}: total == sum of groups", sum(v["total"]), sum(sum(g) for g in v["groups"].values()), 1e-3,
                     "signal inside Total (cf54006)", rel=True)
    for ch, C, clo, sr in (("mutau", 1.050, 0.966, 11048), ("etau", 1.144, 1.025, 5012)):
        L = d["lephad"][ch]
        ck.check(f"{ch} C_OS/SS {C}", L["osss"]["C"], C, 5e-4, R)
        ck.check(f"{ch} same-sign closure {clo}", L["closure_ss"]["ratio"], clo, 5e-4, R)
        ck.check(f"{ch} SR fakes {sr}", L["sr_fakes"], sr, 0.5, R)
        ck.check(f"{ch} AR fractions sum to 1", sum(L["ar_fractions"].values()), 1.0, 1e-9, "consistency")
        ck.check(f"{ch} closure hist ratio == obs / pred", sum(L["closure_ss"]["obs"]) / sum(L["closure_ss"]["pred"]), L["closure_ss"]["ratio"], 1e-6,
                 "consistency", rel=True)
    m = d["lephad"]["mutau"]["ff_inclusive"]
    ck.check("mutau FF_W(OS) ~ 0.08", m["w_os"], 0.08, 0.01, "docs/10-v4-plan.md 5.1", soft=True)
    ck.check("mutau FF_W(SS) ~ 0.043", m["w_ss"], 0.043, 0.006, "docs/10-v4-plan.md 5.1", soft=True)
    ck.check_true("mutau FF_W(OS) > FF_W(SS)", m["w_os"] > m["w_ss"], "docs/10 5.1: OS quark-like, SS gluon-like")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_v4.json")
    args = ap.parse_args()
    finalize(args, build, verify, "ztautau_v4")


if __name__ == "__main__":
    main()
