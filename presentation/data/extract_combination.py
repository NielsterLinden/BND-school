#!/usr/bin/env python
"""Freeze the final combination (four-channel Z -> tautau, 1949 +30 -30 pb) for section 6.

    source setup.sh && python presentation/data/extract_combination.py [--json PATH] [--check-only]

Reads (read-only) combination/combLieke/output/result.json, the file every figure of
combination/combLieke/output/plots/ is drawn from (mf/plots.py), and writes presentation/data/combination_results.json:
  summary     the four rows of summary.png (ee, mumu, tautau standalone fits and the combined fit) + value strings
  prediction  the aMC@NLO line and band (scale + PDF + alpha_s) of summary.png
  published   the CMS and ATLAS rows of combined_vs_published.png (moved to 60-120 GeV as there)
  impacts     the 20 rows of impacts.png, top to bottom: label, post-fit impact up/down in pb, pull and constraint
Row selection, labels and value strings come from mf/plots.py itself (imported), so the clips show what the PNGs show.
Anchors: combination/CLAUDE.md and combination/combLieke/README.md (1949 +30 -30 pb; tautau 1981 +153 -136 pb).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import REPO, Checker, finalize, load_json, provenance, standard_args  # noqa: E402

COMB = REPO / "combination" / "combLieke"
RESULT = COMB / "output" / "result.json"
ORDER = ("ee", "mumu", "tautau")


def _plots():
    sys.path.insert(0, str(COMB))
    from mf import plots   # noqa: E402  (matplotlib, LCG env)
    return plots


def _strip_math(s: str) -> str:
    """matplotlib value string '2094 $^{+122}_{-114}$' / '1952 $\\pm$ 49' -> MathTex '2094\\ {}^{+122}_{-114}'."""
    s = s.replace(" $\\pm$ ", r" \pm ").replace(" $^", r"\ {}^").replace("$", "")
    return s


def build() -> dict:
    plots = _plots()
    res = load_json(RESULT)
    ref = res["poi_reference_pb"]

    rows = []
    for k in ORDER:
        s = res["channels"][k]["standalone"]
        rows.append({"key": k, "sigma_pb": s["sigma_pb"], "up_pb": s["err_up_pb"], "down_pb": s["err_down_pb"]})
    c = res["combined"]
    rows.append({"key": "combined", "sigma_pb": c["sigma_pb"], "up_pb": c["err_up_pb"], "down_pb": c["err_down_pb"],
                 "stat_pb": c.get("stat_pb")})
    for r in rows:
        r["value_tex"] = _strip_math(plots._value_text(r["sigma_pb"], r["up_pb"], r["down_pb"]))

    pred = res["prediction"]
    prediction = {"label": pred["label"], "sigma_pb": pred["sigma_pb"], "up_pb": pred["err_up_pb"],
                  "down_pb": pred["err_down_pb"],
                  "band_lo_pb": pred["sigma_pb"] - pred["err_down_pb"], "band_hi_pb": pred["sigma_pb"] + pred["err_up_pb"]}

    published = []
    for p in res["published"]["combined"]:
        published.append({"label": p["label"], "detail": p["detail"], "window": p["window"],
                          "moved_to_60_120": p["window_scale"] != 1.0, "window_scale": p["window_scale"],
                          "sigma_pb": p["value_60_120_pb"], "err_pb": p["err_60_120_pb"], "citation": p["citation"],
                          "value_tex": _strip_math(plots._value_text(p["value_60_120_pb"], p["err_60_120_pb"], p["err_60_120_pb"]))})

    # impacts.png: the same selection and order as mf/plots.impacts (top 20 by |up| + |down|, no gammas), top to bottom
    ranking = [r for r in c["ranking"] if not r["name"].startswith("gamma")]
    ranking = sorted(ranking, key=lambda r: -(abs(r["post_up"]) + abs(r["post_down"])))[:20]
    impacts = []
    for r in ranking:
        free = r["name"].startswith(plots.FREE_FACTORS)
        row = {"name": r["name"], "label_mpl": plots._np_label(r["name"]),
               "impact_up_pb": r["post_up"] * ref, "impact_down_pb": r["post_down"] * ref,
               "pull": r["pull"], "err_up": r["err_up"], "err_down": r["err_down"], "free": free}
        if free:
            v = c["pulls"][r["name"]]
            row["free_value"], row["free_err"] = v["value"], v["err"]
            row["free_text"] = f"free: {v['value']:.3f} ± {v['err']:.3f}"
        impacts.append(row)
    lim = max(max(abs(r["impact_up_pb"]), abs(r["impact_down_pb"])) for r in impacts)

    return {
        "_doc": "Section 6 (combination): the numbers of combination/combLieke/output/plots/{summary,impacts,"
                "combined_vs_published}.png, frozen from output/result.json with the selection and formatting of mf/plots.py.",
        "provenance": provenance("extract_combination.py", [RESULT, COMB / "mf" / "plots.py", COMB / "config" / "references.json"],
                                 result_inputs=res.get("inputs", {}).get("repository_head"),
                                 observable=res["observable"]),
        "header": {"left": "CMS Open Data", "right": "16.4 fb-1 (13 TeV)", "lumi_pb": res["lumi_pb"]},
        "poi_reference_pb": ref,
        "summary": {"x_range_pb": [1650, 2400], "rows": rows},
        "prediction": prediction,
        "published": published,
        "impacts": {"rows": impacts, "x_lim_pb": 1.35 * lim, "pull_x_range": [-2.6, 2.6],
                    "top3_by_impact": [r["name"] for r in impacts[:3]]},
        "compatibility_p": res["compatibility"]["p_value"],
    }


def verify(d: dict, ck: Checker) -> Checker:
    S = "combination/combLieke/README.md, combination/CLAUDE.md"
    rows = {r["key"]: r for r in d["summary"]["rows"]}
    ck.check("combined 1949.3 +30.2 -29.6 pb", [rows["combined"]["sigma_pb"], rows["combined"]["up_pb"], rows["combined"]["down_pb"]],
             [1949.3, 30.2, 29.6], 0.06, S)
    ck.check("tautau 1981 +153 -136 pb", [rows["tautau"]["sigma_pb"], rows["tautau"]["up_pb"], rows["tautau"]["down_pb"]],
             [1981, 153, 136], 0.5, S)
    ck.check("mumu 1931 +31 -30 pb", [rows["mumu"]["sigma_pb"], rows["mumu"]["up_pb"], rows["mumu"]["down_pb"]], [1931, 31, 30], 0.5, "summary.png")
    ck.check("ee 2094 +122 -114 pb", [rows["ee"]["sigma_pb"], rows["ee"]["up_pb"], rows["ee"]["down_pb"]], [2094, 122, 114], 0.5, "summary.png")
    ck.check_true("value strings as in summary.png",
                  [rows[k]["value_tex"] for k in ("ee", "mumu", "tautau", "combined")]
                  == [r"2094\ {}^{+122}_{-114}", r"1931\ {}^{+31}_{-30}", r"1981\ {}^{+153}_{-136}", r"1949\ {}^{+30}_{-30}"],
                  "summary.png", str([rows[k]["value_tex"] for k in rows]))
    p = d["prediction"]
    ck.check("prediction 1953.9 +56.1 -81.7 pb", [p["sigma_pb"], p["up_pb"], p["down_pb"]], [1953.93, 56.06, 81.66], 0.01,
             "presentation/data/theory_reference.json (1-04)")
    pub = {x["label"]: x for x in d["published"]}
    ck.check("CMS 1952 +- 49 pb", [pub["CMS"]["sigma_pb"], pub["CMS"]["err_pb"]], [1952, 48.63], 0.01, "combined_vs_published.png")
    ck.check("ATLAS 2009 +- 58 pb (66-116 moved to 60-120)", [pub["ATLAS"]["sigma_pb"], pub["ATLAS"]["err_pb"]], [2009.24, 57.88], 0.01,
             "combined_vs_published.png")
    ck.check_true("published strings", [pub["CMS"]["value_tex"], pub["ATLAS"]["value_tex"]] == [r"1952 \pm 49", r"2009 \pm 58"],
                  "combined_vs_published.png")
    im = d["impacts"]["rows"]
    ck.check("20 impact rows", len(im), 20, 0, "impacts.png")
    ck.check_true("impacts.png row order (top 5 / bottom)", [r["name"] for r in im[:5]] + [im[-1]["name"]]
                  == ["Lumi", "ElectronID", "mu_ttbar", "L1Prefiring", "Acc_PDF", "MuonRes"], "impacts.png",
                  str([r["name"] for r in im]))
    ck.check("Lumi impact -22.1 / +22.5 pb", [im[0]["impact_up_pb"], im[0]["impact_down_pb"]], [-22.09, 22.54], 0.01, "impacts.png")
    ck.check("mu_ttbar free 1.113 +- 0.035", [im[2]["free_value"], im[2]["free_err"]], [1.113, 0.035], 0.0006, "impacts.png")
    ck.check("QCDScale (ee shape) pull -1.83", im[11]["pull"], -1.83, 0.005, "impacts.png")
    ck.check("compatibility p = 0.38", d["compatibility_p"], 0.38, 0.005, "CLAUDE.md")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__), "combination_results.json")
    finalize(ap.parse_args(), build, verify, "combination_results")


if __name__ == "__main__":
    main()
