#!/usr/bin/env python
"""Turn the z-ee TRExFitter output into the JSON the combination reads.

    source ../../setup.sh
    python tools/extract_zee.py            # writes inputs/zee_fit_result.json
    python tools/extract_zee.py --print    # and dumps it

**Why this script exists.** z-mumu and z-tautau publish a machine-readable result file
(`fitting/CONVENTIONS.md` section 4); z-ee publishes the raw TRExFitter job instead --
`z-ee/Zee_fit.tar.gz`, 2.7 MB, committed -- plus its histogram inputs in `datasets/z-ee/`. This
script is the adapter: it reads those two and writes `inputs/zee_fit_result.json` in the shape
`comb/inputs.load_ee` expects, so that `run_combination.py` itself stays pure python (no uproot,
no ROOT, no TRExFitter binary -- the rule in `combination/CLAUDE.md`).

Run it again whenever the z-ee group re-publishes `Zee_fit.tar.gz`; `run_combination.py --check`
asserts the numbers it produced and will fail first if the two drift apart.

Three things are computed here rather than read:

1. **The reference cross section.** z-ee's `mu_signal` multiplies a DY_ee template normalised with
   sigma(DY, m > 50) = 6077.22 pb over the *whole* DYJetsToLL_M-50 sample, so on the combination's
   60 < m_LHE < 120 GeV footing its reference is
   6077.22 pb x sumw(LHE ee, 60 < m_LHE < 120) / sumw = 1954.1 pb -- the same construction, on the
   same sample, that z-mumu uses for itself (1953.9 pb) and z-tautau for itself (1944.9 pb).
   The generator weights come from `z-mumu/output/v2/gensums.json` (committed), so the number is
   traceable rather than transcribed from `fitting/CONVENTIONS.md` section 6.

2. **The normalisation part of each systematic template**, integral(Up)/integral(nominal) from the
   TRExFitter histograms inside the tarball. Needed for the `normfix` extraction below.

3. **The `normfix` signal strength.** z-ee's `PDF` and `QCDScale` templates are LHE-weight
   envelopes that were *not* renormalised to a constant yield, so they carry the +-5.9 % scale
   uncertainty of the aMC@NLO cross section itself as a normalisation on the signal -- which
   `fitting/CONVENTIONS.md` section 3 forbids precisely because it is degenerate with the POI.
   The fit pulls `QCDScale` by -1.85 sigma, i.e. it rescales the signal prediction by 0.894, and
   `mu_signal` = 1.049 is measured against *that* rescaled prediction. The cross section that
   corresponds to the measured yield is therefore

       sigma = mu_hat x kappa_theory x sigma^pred,      kappa_theory = prod over PDF, QCDScale

   with kappa the HistFactory interpolation (code 4: exponential outside +-1 sigma). Every other
   nuisance parameter (Lumi, Pileup, L1Prefiring, Electron*) scales the *denominator* L x A x C,
   where a pull belongs inside mu_hat and must not be undone. See docs/01-inputs.md.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMB = HERE.parent
REPO = COMB.parent

TARBALL = REPO / "z-ee/Zee_fit.tar.gz"
ZEE_CONFIG = REPO / "z-ee/fit.config"
GENSUMS = REPO / "z-mumu/output/v2/gensums.json"
HISTO_DIR = REPO / "datasets/z-ee"
OUT = COMB / "inputs/zee_fit_result.json"

JOB = "Zee_fit"
POI = "mu_signal"
LUMI_PB = 16393.381
LUMI_REL_UNC = 0.012
DY_XSEC_PB = 6077.22

#: z-ee `Category`/`SubCategory` string -> the `fitting/CONVENTIONS.md` section 3 Category the
#: combination correlates on. z-ee writes the electron scale factors as Category "Electron" with
#: SubCategory "Electron ID", and TRExFitter reports the SubCategory; "MC stat." is what it calls
#: the per-bin gammas, which the conventions call "Gammas".
CATEGORY_MAP = {
    "Luminosity": "Luminosity",
    "Pileup": "Pileup",
    "L1 prefiring": "L1 prefiring",
    "Background normalisation": "Background normalisation",
    "Signal modelling": "Signal modelling",
    "Electron ID": "Electron efficiency",
    "Electron": "Electron efficiency",
    "MC stat.": "Gammas",
}
#: rows of the grouped-impact file that are totals, not categories
SUMMARY_ROWS = {"FullSyst", "Total syst.", "Stat unc."}

#: the systematics whose templates normalise the *prediction* instead of the response, i.e. the
#: ones z-ee did not renormalise to a constant yield (see the module docstring, point 3)
THEORY_NORM_NPS = ("PDF", "QCDScale")

#: samples of the z-ee fit that are background (everything but the signal and the data)
SIGNAL_SAMPLE = "DY_ee"
DATA_SAMPLE = "Data"
BACKGROUND_SAMPLES = ("DY_tautau", "Diboson", "ttbar", "Wjets")


# ------------------------------------------------------------------------------ tarball readers
def _member(tar: tarfile.TarFile, name: str) -> str:
    f = tar.extractfile(f"{JOB}/{name}")
    if f is None:
        raise RuntimeError(f"{TARBALL.name} has no {JOB}/{name}")
    return f.read().decode()


def parse_fit_txt(text: str) -> tuple[float, float, float, dict[str, list[float]]]:
    """`Fits/<job>.txt`: the POI and every nuisance parameter with its post-fit uncertainty."""
    nps: dict[str, list[float]] = {}
    poi = None
    for line in text.splitlines():
        if line.startswith("CORRELATION_MATRIX"):
            break
        m = re.match(r"^(\S+)\s+([-\d.eE+]+)\s+\+([-\d.eE+]+)\s+-([-\d.eE+]+)\s*$", line)
        if not m:
            continue
        name, val, up, down = m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))
        if name == POI:
            poi = (val, up, down)
        else:
            nps[name] = [val, up, down]
    if poi is None:
        raise RuntimeError(f"no {POI} in Fits/{JOB}.txt")
    return (*poi, nps)


def parse_err_decomp(text: str) -> dict[str, float]:
    """`Fits/<job>_errDecomp_<poi>.txt`: TOT/STAT/SYST/MCSTAT and the per-NP contributions."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].endswith("_ERROR"):
            out[parts[0]] = float(parts[1])
    missing = {"TOT_ERROR", "STAT_ERROR", "SYST_ERROR", "MCSTAT_ERROR"} - set(out)
    if missing:
        raise RuntimeError(f"errDecomp is missing {sorted(missing)}")
    return out


def parse_grouped_impacts(text: str) -> tuple[dict[str, float], dict[str, float]]:
    """`Fits/<job>_group_errDecomp_<poi>.yaml` -> (categories on mu, summary rows)."""
    groups, summary = {}, {}
    name = None
    for line in text.splitlines():
        m = re.match(r"^-?\s*Category:\s*(.+?)\s*$", line)
        if m:
            name = m.group(1).strip().strip('"')
            continue
        m = re.match(r"^\s*Impact:\s*([-\d.eE+]+)\s*$", line)
        if m and name is not None:
            value = abs(float(m.group(1)))
            if name in SUMMARY_ROWS:
                summary[name] = value
            elif name in CATEGORY_MAP:
                groups[CATEGORY_MAP[name]] = groups.get(CATEGORY_MAP[name], 0.0) + value
            else:
                raise RuntimeError(
                    f"z-ee reports a grouped impact for {name!r} which tools/extract_zee.py does "
                    "not map onto a fitting/CONVENTIONS.md Category -- add it to CATEGORY_MAP and "
                    "give it a rho in comb/model.CORRELATION")
            name = None
    if not groups:
        raise RuntimeError("no grouped impacts parsed")
    return groups, summary


def parse_ranking(text: str) -> dict[str, float]:
    """`Rankings/Ranking_<poi>_Breakdown.yaml` -> NP -> symmetrised post-fit impact on mu."""
    out, name = {}, None
    up = None
    for line in text.splitlines():
        m = re.match(r"^-?\s*Name:\s*(.+?)\s*$", line)
        if m:
            name, up = m.group(1).strip(), None
            continue
        m = re.match(r"^\s*POIup:\s*([-\d.eE+]+)\s*$", line)
        if m:
            up = float(m.group(1))
            continue
        m = re.match(r"^\s*POIdown:\s*([-\d.eE+]+)\s*$", line)
        if m and name is not None and up is not None:
            out[name] = 0.5 * (abs(up) + abs(float(m.group(1))))
            name, up = None, None
    return out


def parse_prefit_yields(text: str) -> tuple[dict[str, float], float]:
    """`Plots/SR_prefit.yaml` -> (sample -> prefit yield, data yield). Sample names are TRExFitter
    `Title`s, so they are mapped back onto the config's sample names by the order of the config."""
    samples, order = {}, []
    name = None
    for line in text.splitlines():
        m = re.match(r'^\s*-?\s*Name:\s*"?(.+?)"?\s*$', line)
        if m:
            name = m.group(1)
            continue
        m = re.match(r"^\s*-?\s*Yield:\s*\[(.*)\]\s*$", line)
        if m:
            total = sum(float(x) for x in m.group(1).split(","))
            if name is not None:
                samples[name] = total
                order.append(name)
                name = None
            else:                       # the Total: / Data: blocks carry no Name
                order.append(total)
    data = order[-1] if isinstance(order[-1], float) else None
    if data is None:
        raise RuntimeError("no Data yield in Plots/SR_prefit.yaml")
    return samples, data


# ------------------------------------------------------------------------------ histogram reader
def template_norms(tar: tarfile.TarFile, tmpdir: Path) -> dict[str, dict[str, float]]:
    """integral(Up)/integral(nominal) - 1 for every systematic on the signal sample.

    Read from the TRExFitter histogram file inside the tarball, so it is the template *after*
    TRExFitter's symmetrisation -- exactly what the fit used.
    """
    import uproot                                      # only this script needs it

    path = tmpdir / "Zee_fit_histos.root"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(tar.extractfile(f"{JOB}/Histograms/{JOB}_histos.root").read())
    out = {}
    with uproot.open(path) as f:
        keys = f.keys(cycle=False)
        nominal = f[f"SR/{SIGNAL_SAMPLE}/nominal/SR_{SIGNAL_SAMPLE}"].values().sum()
        prefix = f"SR/{SIGNAL_SAMPLE}/"
        systs = sorted({k[len(prefix):].split("/")[0] for k in keys
                        if k.startswith(prefix) and "/" in k[len(prefix):]} - {"nominal"})
        for syst in systs:
            try:
                up = f[f"{prefix}{syst}/SR_{SIGNAL_SAMPLE}_{syst}_Up"].values().sum()
                dn = f[f"{prefix}{syst}/SR_{SIGNAL_SAMPLE}_{syst}_Down"].values().sum()
            except Exception:
                continue
            out[syst] = {"norm_up": float(up / nominal - 1.0), "norm_down": float(dn / nominal - 1.0)}
    out["_nominal_yield"] = {"value": float(nominal)}
    return out


def histogram_yields() -> dict[str, float]:
    """Sample yields straight from the z-ee fit inputs, as an independent check of the tarball."""
    import uproot

    out = {}
    for name in (SIGNAL_SAMPLE, DATA_SAMPLE, *BACKGROUND_SAMPLES):
        path = HISTO_DIR / f"{name}.root"
        if not path.exists():
            continue
        with uproot.open(path) as f:
            out[name] = float(f["h_mass"].values().sum())
    return out


# ------------------------------------------------------------------------------ the physics bits
def kappa(alpha: float, up: float, down: float) -> float:
    """HistFactory `code 4` normalisation response to a nuisance parameter at `alpha`.

    `up`/`down` are the relative changes of the template integral at +-1 sigma, so the response is
    exponential outside |alpha| > 1 -- what TRExFitter v1.8.0 uses by default -- and the 6th-order
    polynomial of `RooStats::HistFactory::FlexibleInterpVar` (matching value, first and second
    derivative at +-1) inside. Reproduced here rather than imported so that this script needs no
    ROOT.
    """
    hi, lo = 1.0 + up, 1.0 + down
    if alpha >= 1.0:
        return hi ** alpha
    if alpha <= -1.0:
        return lo ** (-alpha)
    log_hi, log_lo = math.log(hi), math.log(lo)
    s0 = 0.5 * (hi + lo)
    a0 = 0.5 * (hi - lo)
    s1 = 0.5 * (hi * log_hi - lo * log_lo)
    a1 = 0.5 * (hi * log_hi + lo * log_lo)
    s2 = 0.5 * (hi * log_hi ** 2 + lo * log_lo ** 2)
    a2 = 0.5 * (hi * log_hi ** 2 - lo * log_lo ** 2)
    a = (15 * a0 - 7 * s1 + a2) / 8.0
    b = (-24 + 24 * s0 - 9 * a1 + s2) / 8.0
    c = (-5 * a0 + 5 * s1 - a2) / 4.0
    d = (12 - 12 * s0 + 7 * a1 - s2) / 4.0
    e = (3 * a0 - 3 * s1 + a2) / 8.0
    f = (-8 + 8 * s0 - 5 * a1 + s2) / 8.0
    x = alpha
    return 1.0 + x * (a + x * (b + x * (c + x * (d + x * (e + x * f)))))


def theory_norm_factor(nps: dict[str, list[float]], norms: dict[str, dict[str, float]]) -> tuple[float, dict]:
    """kappa_theory: how much the *prediction* was rescaled by the un-renormalised theory templates."""
    factor, detail = 1.0, {}
    for name in THEORY_NORM_NPS:
        if name not in nps or name not in norms:
            continue
        alpha = nps[name][0]
        k = kappa(alpha, norms[name]["norm_up"], norms[name]["norm_down"])
        detail[name] = {"pull": alpha, "norm_up": norms[name]["norm_up"],
                        "norm_down": norms[name]["norm_down"], "kappa": k}
        factor *= k
    return factor, detail


def reference_cross_sections() -> dict:
    """sigma^pred(Z/gamma* -> ll, 60 < m_LHE < 120) per flavour, from the generator weight sums."""
    d = json.loads(GENSUMS.read_text())["DY_NLO"]
    sumw = d["sumw"]
    out = {"sumw": sumw, "dy_xsec_pb": DY_XSEC_PB,
           "construction": "6077.22 pb x sumw(LHE flavour, 60 < m_LHE < 120) / sumw of "
                           "DYJetsToLL_M-50 aMC@NLO (recid 35669)",
           "source": str(GENSUMS.relative_to(REPO))}
    for flavour in ("ee", "mumu", "tautau"):
        w = d[f"sumw_lhe_{flavour}_60_120"]
        out[f"sumw_lhe_{flavour}_60_120"] = w
        out[f"sigma_{flavour}_60_120_pb"] = DY_XSEC_PB * w / sumw
        out[f"sigma_{flavour}_m50_pb"] = DY_XSEC_PB * d[f"sumw_lhe_{flavour}"] / sumw
    return out


# ------------------------------------------------------------------------------ main
def build(tmpdir: Path) -> dict:
    with tarfile.open(TARBALL) as tar:
        mu, mu_up, mu_down, nps = parse_fit_txt(_member(tar, f"Fits/{JOB}.txt"))
        errs = parse_err_decomp(_member(tar, f"Fits/{JOB}_errDecomp_{POI}.txt"))
        groups, summary = parse_grouped_impacts(_member(tar, f"Fits/{JOB}_group_errDecomp_{POI}.yaml"))
        ranking = parse_ranking(_member(tar, f"Rankings/Ranking_{POI}_Breakdown.yaml"))
        titles, data_yield = parse_prefit_yields(_member(tar, "Plots/SR_prefit.yaml"))
        norms = template_norms(tar, tmpdir)

    pred = reference_cross_sections()
    sigma_pred = pred["sigma_ee_60_120_pb"]

    nominal_signal = norms.pop("_nominal_yield")["value"]
    prefit = histogram_yields()
    n_bkg = sum(v for k, v in prefit.items() if k in BACKGROUND_SAMPLES)
    n_obs = prefit.get(DATA_SAMPLE, data_yield)
    n_sig = prefit.get(SIGNAL_SAMPLE, nominal_signal)

    k_theory, k_detail = theory_norm_factor(nps, norms)

    # the data statistical uncertainty on mu, analytically: sqrt(N_obs) / (N_obs - N_bkg).
    # TRExFitter's own `Stat unc.` row (STAT_ERROR) is the quadrature remainder of the total after
    # the grouped impacts, and for this fit it is 25x larger -- see docs/01-inputs.md.
    stat_analytic = mu * math.sqrt(n_obs) / (n_obs - n_bkg)
    stat_residual = math.sqrt(max(errs["STAT_ERROR"] ** 2 - stat_analytic ** 2, 0.0))

    common = {"mu_err_up": mu_up, "mu_err_down": mu_down,
              "mu_stat": stat_analytic, "grouped_impact": groups,
              "fit_residual": stat_residual, "ranking": ranking}
    fits = {
        "published": dict(common, mu=mu,
                          label="TRExFitter fit as published by z-ee (mu_signal against the "
                                "un-renormalised signal template)"),
        "normfix": dict(common, mu=mu * k_theory, mu_err_up=mu_up * k_theory,
                        mu_err_down=mu_down * k_theory,
                        mu_stat=stat_analytic * k_theory,
                        fit_residual=stat_residual * k_theory,
                        grouped_impact={k: v * k_theory for k, v in groups.items()},
                        label="the same fit with the normalisation of the un-renormalised PDF and "
                              "QCDScale templates put back into the prediction "
                              "(fitting/CONVENTIONS.md section 3)"),
    }

    return {
        "channel": "Z -> e+ e-",
        "dataset": "CMS Open Data 2016 Electron Run2016G+H (UL NanoAODv9)",
        "generated_by": "combination/tools/extract_zee.py",
        "source": {"tarball": str(TARBALL.relative_to(REPO)),
                   "config": str(ZEE_CONFIG.relative_to(REPO)),
                   "histograms": str(HISTO_DIR.relative_to(REPO)),
                   "gensums": str(GENSUMS.relative_to(REPO))},
        "lumi_pb": LUMI_PB, "lumi_rel_unc": LUMI_REL_UNC,
        "fit_variable": "m_ee", "regions": ["SR"], "poi": POI,
        "nominal_variant": "published",
        "variants": {k: v["label"] for k, v in fits.items()},
        "prediction": pred,
        "sigma_pred_60_120_pb": sigma_pred,
        "fit": fits,
        "error_decomposition": errs,
        "grouped_impact_summary": summary,
        "nuisance_parameters": {k: {"pull": v[0], "err_up": v[1], "err_down": v[2],
                                    **{kk: vv for kk, vv in norms.get(k, {}).items()}}
                                for k, v in nps.items() if not k.startswith("gamma_")},
        "theory_normalisation": {"kappa": k_theory, "nps": k_detail,
                                 "what": "product of the response of the un-renormalised PDF and "
                                         "QCDScale templates at the best fit; sigma = mu_hat x "
                                         "kappa x sigma^pred is the cross section the measured "
                                         "yield corresponds to"},
        "yields_prefit": {**prefit, "background_total": n_bkg},
        "for_combination": {"n_obs": n_obs, "n_bkg_prefit": n_bkg, "n_sig_prefit": n_sig,
                            "mu_counting": (n_obs - n_bkg) / n_sig,
                            "acc_eff": n_sig / (sigma_pred * LUMI_PB),
                            "acc_eff_definition": "N_SR(DY -> ee, prefit, all corrections) / "
                                                  "(sigma^pred(Z/gamma*->ee, 60<m_LHE<120) x L)"},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--print", action="store_true", help="dump the JSON as well")
    args = ap.parse_args()

    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        payload = build(Path(tmp))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=1))

    p, f = payload, payload["fit"]
    print(f"wrote {OUT.relative_to(REPO)}")
    print(f"  sigma^pred(ee, 60-120)      {p['sigma_pred_60_120_pb']:.2f} pb "
          f"(mumu {p['prediction']['sigma_mumu_60_120_pb']:.2f}, "
          f"tautau {p['prediction']['sigma_tautau_60_120_pb']:.2f})")
    print(f"  mu_signal (published)       {f['published']['mu']:.5f} +- {f['published']['mu_err_up']:.5f}"
          f"   -> sigma = {f['published']['mu'] * p['sigma_pred_60_120_pb']:.1f} pb")
    print(f"  kappa_theory                {p['theory_normalisation']['kappa']:.5f}  "
          f"(QCDScale pull {p['nuisance_parameters']['QCDScale']['pull']:+.3f}, "
          f"norm +-{100 * p['nuisance_parameters']['QCDScale']['norm_up']:.2f} %)")
    print(f"  mu_signal (normfix)         {f['normfix']['mu']:.5f}"
          f"   -> sigma = {f['normfix']['mu'] * p['sigma_pred_60_120_pb']:.1f} pb")
    print(f"  counting (data - bkg)/sig   {p['for_combination']['mu_counting']:.5f}"
          f"   -> sigma = {p['for_combination']['mu_counting'] * p['sigma_pred_60_120_pb']:.1f} pb")
    print(f"  data stat: analytic {100 * f['published']['mu_stat'] / f['published']['mu']:.4f} %, "
          f"TRExFitter `Stat unc.` {100 * p['error_decomposition']['STAT_ERROR']:.4f} %")
    if args.print:
        print(json.dumps(payload, indent=1))


if __name__ == "__main__":
    main()
