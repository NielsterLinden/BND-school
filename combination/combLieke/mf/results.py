"""Collect every TRExFitter output of work/ into output/result.json (the single source for plots and docs).

Cross sections: sigma = POI x poi.reference_pb. The TRExFitter MINOS interval of the POI is the total
in-fit uncertainty, which includes the acceptance uncertainties of mumu and tautau (config/channels.json).
Data statistics: the stat-only combined fit (no nuisance parameters, no MC-statistics gammas).
Uncertainty groups: TRExFitter's covariance decomposition of the combined fit (<job>_group_errDecomp),
per Category; groups do not add up in quadrature to the total because the NPs are correlated after the fit.
Ranking: the refit-based NPRanking of `mr` when it exists, otherwise the covariance-based one of `mwf`.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import yaml

from . import prediction
from .paths import HERE, OUTPUT, WORK, manifest, references, repo_path

from fitting import run_trex  # noqa: E402  (mf.paths puts the repository root on sys.path)


def _gof(log: Path):
    text = re.sub(r"\x1b\[[0-9;]*m", "", log.read_text(errors="replace")) if log.exists() else ""
    p = re.findall(r"probability\s*=\s*([0-9.eE+-]+)", text)
    status = {k: re.findall(rf"- {k} status\s+(-?\d+)", text) for k in ("minuit", "hess", "minos")}
    edm = re.findall(r"- Edm = ([0-9.eE+-]+)", text)
    return {"gof_probability": float(p[-1]) if p else None, "minuit_status": int(status["minuit"][-1]) if status["minuit"] else None,
            "hesse_status": int(status["hess"][-1]) if status["hess"] else None,
            "minos_status": int(status["minos"][-1]) if status["minos"] else None, "edm": float(edm[-1]) if edm else None,
            "pos_def_forced": bool(re.search(r"made pos-def|forced pos-def", text))}


def _poi(fit: dict, name: str, ref: float) -> dict:
    v, up, down = fit["nps"][name]
    return {"poi": name, "mu": v, "mu_err_up": up, "mu_err_down": down,
            "sigma_pb": v * ref, "err_up_pb": up * ref, "err_down_pb": down * ref}


def _pulls(fit: dict) -> dict:
    return {k: {"value": v[0], "err": 0.5 * (v[1] + v[2])} for k, v in fit["nps"].items() if not k.startswith("mu_Z")}


def _corr(fit: dict, a: str, b: str):
    names = list(fit["nps"])
    if not fit.get("corr") or a not in names or b not in names:
        return None
    # TRExFitter writes the matrix with its rows in reverse parameter order (row k is parameter n-1-k), columns in order
    return fit["corr"]["matrix"][len(names) - 1 - names.index(a)][names.index(b)]


def _grouped(path: Path, ref: float) -> dict | None:
    if not path.exists():
        return None
    return {k: {"impact": v[0], "up": v[1], "down": v[2], "impact_pb": v[0] * ref} for k, v in run_trex.parse_grouped(path).items()}


def _ranking(fits: Path, rankings: Path) -> tuple[list | None, str | None]:
    refit = fits / "NPRanking_mu_Z.txt"
    if refit.exists():
        rows = run_trex.parse_ranking(refit)
        return ([{"name": r["name"], "pull": r["pull"], "err_up": r["err_up"], "err_down": r["err_down"],
                  "post_up": r["dpoi_up_post"], "post_down": r["dpoi_down_post"],
                  "pre_up": r["dpoi_up_pre"], "pre_down": r["dpoi_down_pre"]} for r in rows], "refit (trex-fitter mr)")
    cov = rankings / "Ranking_mu_Z_Breakdown.yaml"
    if cov.exists():
        rows = yaml.safe_load(cov.read_text())
        return ([{"name": r["Name"], "pull": r["NPhat"], "err_up": r["NPerrHi"], "err_down": abs(r["NPerrLo"]),
                  "post_up": r["POIup"], "post_down": r["POIdown"], "pre_up": None, "pre_down": None} for r in rows],
                "post-fit covariance (trex-fitter mwf)")
    return None, None


def _channel_published(key: str, spec: dict, status: dict) -> dict | None:
    """The channel group's own POI, read from the file it published (not from config/channels.json).
    None when the file is not there yet and the manifest marks it "optional" (tautau v4 before its fits are pulled)."""
    src = spec["published"]["source"]
    if key != "ee" and spec["published"].get("optional") and not repo_path(src).exists():
        print(f"NOTE: {src} does not exist yet: no 'published_by_channel' for {key}")
        return None
    if key == "ee":
        v = status["ee_input"]["published_fit"]
        return {"source": src, "mu": v[0], "err_up": v[1], "err_down": v[2]}
    d = json.loads(repo_path(src).read_text())
    if "poi_value" in d:
        return {"source": src, "mu": d["poi_value"], "err_up": d["poi_err_up"], "err_down": d["poi_err_down"]}
    return {"source": src, "mu": d["mu"], "err_up": d["mu_err_up"], "err_down": d["mu_err_down"]}


def _published(entry: dict, scale: float) -> dict:
    k = scale if tuple(entry["window"]) == (66, 116) else 1.0
    err = math.sqrt(sum(e * e for e in entry["errors"].values()))
    return {**entry, "value_60_120_pb": entry["value"] * k, "err_60_120_pb": err * k, "window_scale": k}


def collect(interim: bool = False) -> dict:
    """-> output/result.json; with `interim` -> interim/result.json, leaving output/ (and what reads it) alone."""
    m, refs = manifest(), references()
    ref = m["poi"]["reference_pb"]
    status = json.loads((WORK / "status.json").read_text())
    common = run_trex.parse_fit_txt(WORK / "common/combination/Fits/combination.txt")
    split = run_trex.parse_fit_txt(WORK / "split/combination/Fits/combination.txt")
    out = {"observable": m["observable"], "lumi_pb": m["lumi_pb"], "poi_reference_pb": ref,
           "ee_input": status["ee_input"], "model": {k: v["channels"] for k, v in status["likelihoods"].items() if k == "common"}}

    comb = _poi(common, m["poi"]["name"], ref)
    comb.update(_gof(WORK / "common/logs/multifit_mwf.log"), nll=common["nll"], pulls=_pulls(common))
    fits_dir = WORK / "common/combination/Fits"
    comb["grouped_impacts"] = _grouped(fits_dir / "combination_group_errDecomp_mu_Z.txt", ref)
    full = (comb["grouped_impacts"] or {}).get("FullSyst")
    if full and full["impact"] > 0:
        # the decomposition is HESSE-based and its +/- columns are rescaled to MINOS: their ratio checks the covariance
        comb["hesse_over_minos"] = full["impact"] / (0.5 * (full["up"] + full["down"]))
        if abs(comb["hesse_over_minos"] - 1) > 0.1:
            print(f"WARNING: HESSE/MINOS = {comb['hesse_over_minos']:.2f} for mu_Z in the combined fit: the covariance, and so "
                  f"the uncertainty groups and post-fit NP errors, are unreliable (see CLAUDE.md, FitStrategy)")
    stat = fits_dir / "combination_statOnly.txt"
    if stat.exists():
        s = run_trex.parse_fit_txt(stat)["nps"][m["poi"]["name"]]
        comb["stat_pb"] = 0.5 * (s[1] + s[2]) * ref
        comb["stat_only_fit"] = {"mu": s[0], "err_up": s[1], "err_down": s[2]}
        tot = 0.5 * (comb["err_up_pb"] + comb["err_down_pb"])
        comb["syst_pb"] = math.sqrt(max(tot ** 2 - comb["stat_pb"] ** 2, 0.0))
    comb["ranking"], comb["ranking_method"] = _ranking(fits_dir, WORK / "common/combination/Rankings")
    scan = WORK / "common/combination_scan/LHoodPlots/NLLscan_mu_Z.yaml"
    comb["nll_scan"] = yaml.safe_load(scan.read_text()) if scan.exists() else None
    out["combined"] = comb

    chans = {}
    for key in m["channels"]:
        fit = run_trex.parse_fit_txt(WORK / f"common/{key}/Fits/{key}.txt")
        standalone = _poi(fit, m["poi"]["name"], ref)
        standalone.update(_gof(WORK / f"common/logs/{key}_f.log"), pulls=_pulls(fit))
        joint = _poi(split, f'{m["poi"]["name"]}_{key}', ref)
        pub = _channel_published(key, m["channels"][key], status)
        if pub is not None:
            pub = {**pub, "sigma_pb": pub["mu"] * m["channels"][key]["sigma_reference_pb"]}
        sub = {}
        for label, src in m["channels"][key].get("published_submeasurements", {}).items():
            if not label.startswith("_"):
                d = json.loads(repo_path(src).read_text())
                k = m["channels"][key]["sigma_reference_pb"]
                sub[label] = {"source": src, "mu": d["poi_value"], "sigma_pb": d["poi_value"] * k,
                              "err_up_pb": d["poi_err_up"] * k, "err_down_pb": d["poi_err_down"] * k}
        chans[key] = {"standalone": standalone, "joint": joint,
                      "published_by_channel": pub, "published_submeasurements": sub,
                      "sigma_reference_pb": m["channels"][key]["sigma_reference_pb"]}
    out["channels"] = chans

    split_info = _gof(WORK / "split/logs/multifit_mwf.log")
    q = 2.0 * (common["nll"] - split["nll"])
    ndf = len(m["channels"]) - 1
    out["compatibility"] = {"q": q, "ndf": ndf, "p_value": math.exp(-q / 2.0) if ndf == 2 else None,
                            "definition": "-2 ln(L_common / L_split): one mu_Z against one per channel, all NPs profiled",
                            "split_fit": split_info,
                            "correlations": {f"{a}-{b}": _corr(split, f"mu_Z_{a}", f"mu_Z_{b}")
                                             for a, b in (("ee", "mumu"), ("ee", "tautau"), ("mumu", "tautau"))}}

    variations = {}
    for name, var in m.get("variations", {}).items():
        if name.startswith("_"):
            continue
        path = WORK / f"var_{name}/combination/Fits/combination.txt"
        if not path.exists():
            continue
        v = _poi(run_trex.parse_fit_txt(path), m["poi"]["name"], ref)
        v.update(_gof(WORK / f"var_{name}/logs/multifit_mwf.log"), label=var["label"], note=var.get("note"),
                 shift_pb=v["sigma_pb"] - comb["sigma_pb"])
        variations[name] = v
    out["variations"] = variations

    checks = {}
    for name, chk in m.get("checks", {}).items():
        key = chk.get("channel") if isinstance(chk, dict) else None
        path = WORK / f"check_{name}/{key}/Fits/{key}.txt"
        if name.startswith("_") or not path.exists():
            continue
        v = _poi(run_trex.parse_fit_txt(path), m["poi"]["name"], ref)
        base = chans[key]["standalone"]
        v.update(_gof(WORK / f"check_{name}/logs/{key}_f.log"), label=chk["label"], note=chk.get("note"), channel=key,
                 shift_pb=v["sigma_pb"] - base["sigma_pb"],
                 err_ratio=(v["err_up_pb"] + v["err_down_pb"]) / (base["err_up_pb"] + base["err_down_pb"]))
        checks[name] = v
    out["checks"] = checks

    out["prediction"] = prediction.amcatnlo(m["poi"].get("prediction_flavour", "mumu"))
    scale = prediction.window_ratio()
    if abs(scale - refs["window_scale_66_116_to_60_120"]["value"]) > 1e-4:
        raise ValueError("config/references.json window scale disagrees with the generator sums")
    out["published"] = {"window_scale_66_116_to_60_120": scale,
                        "channel": {k: [_published(e, scale) for e in v] for k, v in refs["channel"].items()},
                        "combined": [_published(e, scale) for e in refs["combined"]]}
    orth = HERE / "checks/orthogonality.json"
    if orth.exists():
        o = json.loads(orth.read_text())
        # measured for ee, mumu and tau_h tau_h; the mu tau_h, e tau_h and e mu channels veto a second lepton by
        # construction and mumu_CRemu, the one overlap with e mu, is dropped (config/channels.json)
        out["orthogonality"] = {"mumu_ee_upper_bound": o["mumu_and_ee"]["overlap_upper_bound"],
                                "tautau_mumu": o["tautau_and_mumu_ee"]["tautau_sr_and_two_mumu_signal_muons"],
                                "tautau_ee_upper_bound": o["tautau_and_mumu_ee"]["tautau_sr_and_ee_selected"],
                                "source": "checks/orthogonality.json"}
    target = HERE / "interim" / "result.json" if interim else OUTPUT / "result.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(out, indent=1) + "\n")
    print(f"-> {target}")
    c = out["combined"]
    print(f"sigma = {c['sigma_pb']:.1f} +{c['err_up_pb']:.1f} -{c['err_down_pb']:.1f} pb  (mu_Z = {c['mu']:.4f}); "
          f"compatibility q = {q:.2f} / {ndf} dof")
    for k, v in chans.items():
        s, j = v["standalone"], v["joint"]
        own = f"{v['published_by_channel']['sigma_pb']:.1f}" if v["published_by_channel"] else "not published yet"
        print(f"  {k:7s} standalone {s['sigma_pb']:7.1f} +{s['err_up_pb']:.1f} -{s['err_down_pb']:.1f}   "
              f"joint {j['sigma_pb']:7.1f} +{j['err_up_pb']:.1f} -{j['err_down_pb']:.1f}   channel's own {own}")
    for k, v in {**variations, **checks}.items():
        print(f"  {k:28s} {v['sigma_pb']:7.1f} +{v['err_up_pb']:.1f} -{v['err_down_pb']:.1f}  ({v['shift_pb']:+.1f})")
    return out
