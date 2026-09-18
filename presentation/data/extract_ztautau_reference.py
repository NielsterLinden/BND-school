#!/usr/bin/env python
"""Freeze the reference points of the sigma axis of the last Z -> tautau clip: this measurement, the published CMS and ATLAS
13 TeV measurements of sigma(Z/gamma* -> ll) and the aMC@NLO prediction the signal strength multiplies.

    source fitting/setup.sh && python presentation/data/extract_ztautau_reference.py [--json PATH] [--check-only]

Reads (read-only): presentation/data/ztautau_fit.json (sigma_60_120: value, err_up/down, acc, pred), combination/result.md
("Comparison with published measurements" table: CMS-SMP-20-004 / arXiv:2408.03744 and arXiv:1603.09222, quoted verbatim).
The repository documents no PDF + scale + alpha_s uncertainty on the predicted sigma(60-120) (THEORY_SEARCHED); the
band drawn around the 1944.9 pb line is the relative uncertainty of the published NNLO+NNLL prediction of the same
quantity with the same PDF family (THEORY_PUB: CMS-SMP-20-004, arXiv:2408.03744, Table 5, DYturbo v1.3.2, NNPDF3.1,
13 TeV, 60 < m < 120 GeV: 1940 +15 -21 pb, stat + PDF + alpha_s + muR/muF), frozen here with its reference.
Writes presentation/data/ztautau_reference.json.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _extract_common import DATA_DIR, REPO, Checker, finalize, load_json, provenance, standard_args  # noqa: E402

FIT_JSON = DATA_DIR / "ztautau_fit.json"
RESULT_MD = REPO / "combination" / "result.md"
# published values as they stand in combination/result.md (frozen here; verify() re-reads the file)
CMS = {"value": 1952, "stat": 4, "syst": 18, "lumi": 45, "window": "60-120", "sqrt_s_TeV": 13, "lumi_pb": 206, "label": "CMS",
       "reference": "CMS-SMP-20-004, arXiv:2408.03744"}
ATLAS = {"value": 1981, "stat": 7, "syst": 38, "lumi": 42, "window": "66-116", "sqrt_s_TeV": 13, "lumi_pb": 81, "label": "ATLAS",
         "reference": "arXiv:1603.09222", "note": "narrower mass window, not directly comparable"}
# documents searched for an uncertainty on the predicted sigma(Z/gamma* -> ll, 60-120) (PDF + scale + alpha_s)
THEORY_SEARCHED = ["combination/result.md", "combination/docs/00-overview.md", "combination/docs/01-inputs.md", "combination/docs/02-correlation-model.md",
                   "combination/docs/03-method.md", "combination/docs/04-results.md", "docs/CONVENTIONS.md", "z-tautau/docs/06-cross-section.md",
                   "z-tautau/handoff.md", "z-mumu/docs/06-cross-section.md", "z-mumu/docs/15-combination-inputs.md", "z-mumu/handoff.md"]
# the published prediction whose relative uncertainty is drawn as the band (read off arXiv:2408.03744v2, Table 5,
# "Z -> l+ l-" row, NNPDF3.1 column, on 2026-09-16; the other columns: NNPDF4.0 1970 +11 -14, CT18 1921 +30 -33,
# MSHT20 1935 +23 -27 pb)
THEORY_PUB = {"value": 1940.0, "up": 15.0, "down": 21.0, "pdf": "NNPDF3.1", "order": "NNLO+NNLL QCD (DYturbo v1.3.2)",
              "window": "60-120", "sqrt_s_TeV": 13, "reference": "CMS-SMP-20-004, arXiv:2408.03744, Table 5",
              "components": "statistical, PDF, alpha_s, renormalisation and factorisation scale"}
_PUB = re.compile(r"\|\s*(CMS|ATLAS), 13 TeV, (\d+) pb\^-1, (\d+) < m < (\d+) GeV\s*\|\s*(\d+) \+- (\d+) \(stat\) \+- (\d+) \(syst\) \+- (\d+) \(lumi\) pb\s*\|")


def published_lines(text):
    """{label: (line, dict)} of the 'Comparison with published measurements' rows of combination/result.md."""
    out = {}
    for line in text.splitlines():
        m = _PUB.search(line)
        if m:
            lab, lumi, lo, hi, val, stat, syst, lum = m.groups()
            out[lab] = (line.strip(), {"value": int(val), "stat": int(stat), "syst": int(syst), "lumi": int(lum), "window": f"{lo}-{hi}", "lumi_pb": int(lumi)})
    return out


def total(d):
    return float(np.sqrt(d["stat"] ** 2 + d["syst"] ** 2 + d["lumi"] ** 2))


def build():
    fit = load_json(FIT_JSON)
    s = fit["sigma_60_120"]
    pub = published_lines(RESULT_MD.read_text())
    assert set(pub) == {"CMS", "ATLAS"}, f"published-measurement rows not found in {RESULT_MD}: {list(pub)}"
    cms = {**CMS, "total": total(CMS)}
    atlas = {**ATLAS, "total": total(ATLAS)}
    pred = float(s["pred"])
    theory = {"value": pred,
              "unc_up": pred * THEORY_PUB["up"] / THEORY_PUB["value"], "unc_down": pred * THEORY_PUB["down"] / THEORY_PUB["value"],
              "unc_rel_up": THEORY_PUB["up"] / THEORY_PUB["value"], "unc_rel_down": THEORY_PUB["down"] / THEORY_PUB["value"],
              "unc_source": (f"relative uncertainty of the published {THEORY_PUB['order']} {THEORY_PUB['pdf']} prediction "
                             f"{THEORY_PUB['value']:.0f} +{THEORY_PUB['up']:.0f} -{THEORY_PUB['down']:.0f} pb ({THEORY_PUB['components']}), "
                             f"{THEORY_PUB['reference']}"),
              "published_prediction": THEORY_PUB,
              "label": "aMC@NLO / NNLO",
              "note": ("sigma^pred(Z/gamma* -> tautau, 60 < m_LHE < 120) = 6077.22 pb (FEWZ NNLO, all flavours) x sumw(LHE tautau, 60-120) / sumw of the "
                       "aMC@NLO DYJetsToLL_M-50 sample (z-tautau/docs/06, docs/CONVENTIONS.md section 6). The repository documents no PDF + scale + "
                       "alpha_s uncertainty on this number (only on the acceptance A: 3.7 %, and on the mumu reference 1953.9 pb the 0.47 % definitional spread; "
                       "the documents searched are listed in provenance.theory_searched), so unc_up / unc_down carry the relative uncertainty of the "
                       "published NNLO+NNLL prediction of the same quantity with the same PDF family (unc_source); the two central values agree to 0.25 %"),
              "sigma_m50_pb": float(fit["prediction"]["sigma_tautau_pb"]),
              "reference_mumu_pb": 1953.9}
    out = {
        "provenance": provenance("extract_ztautau_reference.py", [FIT_JSON, RESULT_MD],
                                 dataset="CMS 2016 Open Data, Tau Run2016G+H, NanoAODv9 (records 30532, 30565)",
                                 published_lines={k: v[0] for k, v in pub.items()},
                                 published_note="cms / atlas copied verbatim from combination/result.md, 'Comparison with published measurements'; "
                                                "total = sqrt(stat^2 + syst^2 + lumi^2) computed here",
                                 this_work_note="this_work = ztautau_fit.json sigma_60_120 (value = mu_Z x pred; err_up/down = MINOS syst (+) stat; acc = acceptance "
                                                "uncertainty 3.7 %, quoted separately)",
                                 theory_searched=THEORY_SEARCHED,
                                 anchors=["presentation/data/ztautau_fit.json", "combination/result.md", "z-tautau/handoff.md (2082 +222/-194 +- 76 (acc) pb, NNLO 1945 pb)"]),
        "unit": "pb",
        "this_work": {"value": float(s["value"]), "err_up": float(s["err_up"]), "err_down": float(s["err_down"]), "acc": float(s["acc"]),
                      "stat": float(s["stat"]), "syst": float(s["syst"]), "window": "60-120", "sqrt_s_TeV": 13, "lumi_pb": float(fit["lumi_pb"]),
                      "label": "this work", "source": "ztautau_fit.json sigma_60_120"},
        "cms": cms,
        "atlas": atlas,
        "theory": theory,
    }
    return out


def verify(d, ck: Checker):
    R = "combination/result.md"
    pub = published_lines(RESULT_MD.read_text())
    for lab, key in (("CMS", "cms"), ("ATLAS", "atlas")):
        want = pub[lab][1]
        ck.check(f"{lab} value / stat / syst / lumi == result.md line", [d[key][k] for k in ("value", "stat", "syst", "lumi")], [want[k] for k in ("value", "stat", "syst", "lumi")], 0, R)
        ck.check_true(f"{lab} window and luminosity == result.md line", d[key]["window"] == want["window"] and d[key]["lumi_pb"] == want["lumi_pb"], R,
                      f"{d[key]['window']} GeV, {d[key]['lumi_pb']} pb^-1")
        ck.check(f"{lab} total == sqrt(stat^2 + syst^2 + lumi^2)", d[key]["total"], total(d[key]), 1e-9, "consistency")
    ck.check("CMS 1952 +- 4 +- 18 +- 45 (arXiv:2408.03744)", [d["cms"]["value"], d["cms"]["stat"], d["cms"]["syst"], d["cms"]["lumi"]], [1952, 4, 18, 45], 0, R)
    ck.check("ATLAS 1981 +- 7 +- 38 +- 42 (arXiv:1603.09222)", [d["atlas"]["value"], d["atlas"]["stat"], d["atlas"]["syst"], d["atlas"]["lumi"]], [1981, 7, 38, 42], 0, R)
    ck.check_true("ATLAS window 66-116 flagged as not comparable", d["atlas"]["window"] == "66-116" and "not directly comparable" in d["atlas"].get("note", ""), R)
    fit = load_json(FIT_JSON)["sigma_60_120"]
    ck.check("this_work == ztautau_fit.json sigma_60_120 (value, err_up, err_down, acc)", [d["this_work"][k] for k in ("value", "err_up", "err_down", "acc")],
             [fit[k] for k in ("value", "err_up", "err_down", "acc")], 1e-9, "ztautau_fit.json")
    ck.check("this_work 2082 +222 -194 +- 76 (acc) pb", [d["this_work"]["value"], d["this_work"]["err_up"], d["this_work"]["err_down"], d["this_work"]["acc"]], [2082, 222, 194, 76], 0.5, "z-tautau/handoff.md")
    ck.check("theory value == ztautau_fit.json sigma_60_120.pred (1944.9 pb)", d["theory"]["value"], fit["pred"], 1e-9, "ztautau_fit.json")
    ck.check("theory value 1944.9 pb (docs/06, CONVENTIONS section 6)", d["theory"]["value"], 1944.9, 0.05, "z-tautau/docs/06-cross-section.md")
    ck.check_true("theory uncertainty either documented with a source or null", (d["theory"]["unc_up"] is None and d["theory"]["unc_down"] is None and
                  d["theory"]["unc_source"] == "not documented in the repository") or (d["theory"]["unc_up"] is not None and d["theory"]["unc_source"] not in ("", None)),
                  detail=d["theory"]["unc_source"])
    ck.check("theory band +15 -21 pb (1944.9 x 15/1940, x 21/1940; arXiv:2408.03744 Table 5)", [d["theory"]["unc_up"], d["theory"]["unc_down"]], [15, 21], 0.5,
             "arXiv:2408.03744")
    ck.check_true("published prediction within 0.5 % of the aMC@NLO/FEWZ reference", abs(THEORY_PUB["value"] / d["theory"]["value"] - 1) < 0.005,
                  "arXiv:2408.03744", f"{THEORY_PUB['value']:.0f} vs {d['theory']['value']:.1f} pb")
    ck.check_true("this work within 1 sigma of CMS", abs(d["this_work"]["value"] - d["cms"]["value"]) < d["this_work"]["err_down"] + d["cms"]["total"], "sanity")
    return ck


def main():
    ap = standard_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter), "ztautau_reference.json")
    args = ap.parse_args()
    finalize(args, build, verify, "ztautau_reference")


if __name__ == "__main__":
    main()
