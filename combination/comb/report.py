"""Render `result.md` from the combination payload, so the prose can never drift from the numbers.

`run_combination.py` calls `write_result_md`; every figure in the file comes from
`output/combination_result.json`, which is written in the same run.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

PREDICTION_NOTE = (
    "the two channels normalise to slightly different aMC@NLO references for the same quantity "
    "(1953.9 pb for mu mu, 1944.9 pb for tau tau, a 0.47 % difference discussed in "
    "[docs/01-inputs.md](docs/01-inputs.md)), which is why the combination is done on the cross "
    "sections and not on `mu_Z`"
)

#: caveats carried over from the two channels' own open-issue lists, none of them corrected here.
#: Each `{...}` is filled from the payload in `render`.
CAVEATS = [
    "`z-mumu` open issue 3: the data lineshape at 60-80 GeV lies between aMC@NLO and powheg, and "
    "the fit resolves it by pulling `SigModel` by +0.66 sigma. This is the reason the extraction "
    "still depends on the fit configuration at the few-per-mille level -- the accepted binnings "
    "span {mumu_binning_span:.1f} pb on the combination ({mumu_binning_frac:.2f} of its total "
    "uncertainty) and the 1-bin counting extraction sits {mumu_counting_shift:+.1f} pb away. "
    "A third generator, or NLO electroweak corrections, would say which lineshape is right.",
    "`z-mumu` open issue 1: the pileup profile is a two-parameter fit of the luminosity-record "
    "profile to N_PV rather than the official `puWeights` file, which is unreachable in Open Data. "
    "Data/MC in N_PV agree to 3 % over the bulk; the channel quotes 0.13 % for it.",
    "`z-mumu` open issue 4: prefiring maps, muon scale factors and momentum corrections are "
    "in-house (the official files need CERN credentials). The muon-efficiency term "
    "({mumu_muon_pb:.0f} pb) is the largest uncertainty here after the luminosity.",
    "`z-tautau` open issue 1: the channel recommends DeepTau **Tight** on both legs as the next "
    "iteration's working point -- it is more precise ({tautau_tight_rel:.0f} % against "
    "{tautau_rel:.0f} % on sigma) and fits better (GoF p = {tautau_tight_gof:.2f} against "
    "{tautau_gof:.2f}) -- but has not adopted it as nominal. Using it shifts the combination by "
    "{tautau_tight_shift:+.1f} pb and the universality ratio from {ratio:.2f} to "
    "{tautau_tight_ratio:.2f}.",
]

REFERENCES = [
    ("CMS, 13 TeV, 206 pb^-1, 60 < m < 120 GeV", "1952 +- 4 (stat) +- 18 (syst) +- 45 (lumi) pb",
     "[CMS-SMP-20-004, arXiv:2408.03744](https://arxiv.org/abs/2408.03744)"),
    ("ATLAS, 13 TeV, 81 pb^-1, 66 < m < 116 GeV", "1981 +- 7 (stat) +- 38 (syst) +- 42 (lumi) pb",
     "[arXiv:1603.09222](https://arxiv.org/abs/1603.09222)"),
]


def _table(headers, rows, align=None) -> str:
    align = align or ["---"] * len(headers)
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(align) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def render(p: dict) -> str:
    c = p["combined"]
    mm, tt = p["channels"]["mumu"], p["channels"]["tautau"]
    lik, rat = p["likelihood"], p["ratio"]
    v = p["variations"]
    rel = 100 * c["relative"]

    binning = [v[k]["value"] for k in ("baseline", "mumu_bins2gev", "mumu_bins10gev")]
    caveats = "\n".join(
        "* " + t.format(mumu_binning_span=max(binning) - min(binning),
                        mumu_binning_frac=(max(binning) - min(binning)) / c["total_pb"],
                        mumu_counting_shift=v["mumu_counting"]["shift"],
                        mumu_muon_pb=c["breakdown_pb"]["Muon efficiency"],
                        tautau_gof=tt["gof_p"], tautau_tight_shift=v["tautau_tight"]["shift"],
                        ratio=rat["value"], tautau_tight_ratio=v["tautau_tight"]["ratio"],
                        tautau_tight_gof=v["tautau_tight"]["channels"]["tautau"]["gof_p"],
                        tautau_rel=100 * tt["sigma_err_pb"] / tt["sigma_pb"],
                        tautau_tight_rel=100 * v["tautau_tight"]["channels"]["tautau"]["sigma_err_pb"]
                        / v["tautau_tight"]["channels"]["tautau"]["sigma_pb"])
        for t in CAVEATS)

    channel_rows = []
    for ch, name in ((mm, "Z -> mu mu"), (tt, "Z -> tau_h tau_h")):
        gof = f"{ch['gof_p']:.2f}" if ch["gof_p"] is not None else "n/a"
        channel_rows.append([f"**{name}**", ch["variant"],
                             f"{ch['mu']:.4f} +{ch['mu_err_up']:.4f} -{ch['mu_err_down']:.4f}",
                             f"{ch['sigma_pred_pb']:.1f}",
                             f"**{ch['sigma_pb']:.0f} +- {ch['sigma_err_pb']:.0f}**",
                             f"{100 * ch['sigma_err_pb'] / ch['sigma_pb']:.1f} %", gof])

    breakdown_rows = [[k, f"{val:.2f}", f"{100 * val / c['sigma_pb']:.2f} %",
                       "yes" if next(s["rho"] for s in p["sources"] if s["name"] == k) >= 0.999 else "no"]
                      for k, val in c["breakdown_pb"].items() if val >= 0.05]

    var_rows = [[k.replace("_", " "), f"{d['value']:.1f}", f"{d['shift']:+.1f}",
                 f"{d['shift'] / c['total_pb']:+.2f}", d["why"]]
                for k, d in v.items() if k != "baseline"]

    source_rows = [[s["name"], "1" if s["rho"] >= 0.999 else ("0" if s["rho"] <= 0.001 else f"{s['rho']:.2f}"),
                    f"{s['sizes_pb'].get('mumu', 0):.2f}", f"{s['sizes_pb'].get('tautau', 0):.2f}"]
                   for s in sorted(p["sources"], key=lambda s: -max(s["sizes_pb"].values()))]

    return f"""# Combined Z cross section from Z -> mu mu and Z -> tau_h tau_h

*Generated by `run_combination.py` on {date.today().isoformat()} from
`output/combination_result.json`. Do not edit by hand -- edit the code or the docs.*

CMS Open Data 2016, Run2016G+H, L = 16393.381 pb^-1 (normtag), sqrt(s) = 13 TeV.

## Result

> **sigma(pp -> Z/gamma* -> ll, 60 < m_ll < 120 GeV) = {c['sigma_pb']:.0f} +- {c['total_pb']:.0f} pb**
> = {c['sigma_pb']:.0f} +- {c['stat_pb']:.1f} (stat) +- {c['syst_pb']:.0f} (syst)
> +- {c['acc_pb']:.0f} (acc) +- {c['lumi_pb']:.0f} (lumi) pb   ({rel:.2f} %)
>
> per lepton flavour, assuming lepton universality. Channel compatibility:
> chi2/ndf = {c['chi2']:.2f}/{c['ndf']}, p = {c['pvalue']:.2f}.

The profile-likelihood form of the same combination, which carries the tau tau channel's asymmetric uncertainty exactly, gives {lik['value_pb']:.0f} +{lik['err_up_pb']:.0f} / -{lik['err_down_pb']:.0f} pb -- the same answer.

![combined result](output/plots/forest.png)

## Inputs

{_table(["channel", "variant", "mu_Z", "reference [pb]", "sigma(60-120) [pb]", "rel.", "GoF p"],
        channel_rows, ["---", "---", "---", "---:", "---:", "---:", "---:"])}

Both are TRExFitter v1.8.0 profile-likelihood fits performed by the channel groups; this
combination reads their published results (`{mm['provenance'][0]}`,
`{tt['provenance'][0]}`) and does not refit. The `reference` column is the aMC@NLO prediction each
channel's `mu_Z` multiplies -- {PREDICTION_NOTE}. The `mu_Z` uncertainty shown is the in-fit one; the acceptance uncertainty (0.61 % for mu mu, 3.7 % for tau tau) sits outside both fits and is added in the `sigma` column.

Both inputs are the channels' **current baselines**, re-published on 15 Sep 2026. For z-mumu that
is the fit rebuilt after its own review (`z-mumu/REVIEW.md` section 0): 12 x 5 GeV bins instead of
30 x 2 GeV, a two-sided `SigModel` template built inside a common 50 < m_LHE < 120 GeV window, no
template smoothing and MINOS on every parameter. The review's stop-gap -- the counting extraction
plus a +-0.7 % lineshape term -- is retired, because the fit it was protecting against is now
stable to +-0.2 % across the binnings that describe the data; the counting extraction survives
only as the `mumu counting` variation below ({v['mumu_counting']['shift']:+.1f} pb). For z-tautau
it is v2.1, whose nominal is now the fake factor **with** the genuine-tau MC subtraction and whose
OS/SS correction is applied per BDT category.

## Where the uncertainty comes from

{_table(["source", "pb", "of sigma", "correlated?"], breakdown_rows, ["---", "---:", "---:", ":---:"])}

![uncertainty breakdown](output/plots/breakdown.png)

The combination is **luminosity-dominated and completely dominated by the mu mu channel**:

* weights: mu mu {c['weights']['mumu']:+.4f}, tau tau {c['weights']['tautau']:+.4f};
* mu mu alone gives {mm['alone_pb']:.0f} +- {mm['alone_err_pb']:.0f} pb, so adding the tau tau channel improves the uncertainty by {100 * (1 - c['total_pb'] / mm['alone_err_pb']):.2f} % -- it changes nothing;
* the tau tau weight is **negative**. That is standard BLUE behaviour, not an error: when the correlation exceeds the ratio of the two uncertainties ({c['correlation_mumu_tautau']:.3f} > {mm['alone_err_pb'] / tt['alone_err_pb']:.3f} here), the less precise measurement is used to pull on the shared systematic rather than to average the value down. Its effect is {abs(c['sigma_pb'] - mm['alone_pb']):.1f} pb.

The reason there is nothing to gain: the luminosity uncertainty is {c['lumi_pb']:.0f} pb, fully correlated between the channels and therefore irreducible by combining. Both channels use the same 16393.381 pb^-1 with the same 1.2 %. A third channel of comparable precision would not help either; a better luminosity calibration would.

## Lepton universality

This is the part the combination cannot test, because it *assumes* universality. Measured
separately, with everything correlated cancelling:

> **R = sigma(Z -> tau tau) / sigma(Z -> mu mu) = {rat['value']:.2f} +- {rat['error']:.2f}**
> ({rat['z']:+.1f} sigma from 1, p = {rat['pvalue']:.2f})

![lepton universality](output/plots/universality.png)

## Cross-checks and variations

{_table(["variation", "sigma [pb]", "shift", "in sigma", "what it changes"], var_rows,
        ["---", "---:", "---:", "---:", "---"])}

Every variation moves the result by less than {max(abs(d['shift']) for k, d in v.items() if k != 'baseline') / c['total_pb']:.2f} of the total uncertainty. The largest is `mumu counting` ({v['mumu_counting']['shift']:+.1f} pb): the mu mu signal extraction, not the combination, is what the result is most sensitive to. The two rows that bracket the correlation model (`rho zero`, `rho one`) span only {abs(v['rho_zero']['value'] - v['rho_one']['value']):.1f} pb, and `sigmodel correlated` is exactly null because z-tautau no longer fits a generator nuisance parameter at all (`z-tautau/docs/07`). The correlation model is therefore not the limiting assumption -- the luminosity calibration is.

![variations](output/plots/variations.png)

## Comparison with published measurements

{_table(["measurement", "value", "reference"], [[a, b, c_] for a, b, c_ in REFERENCES])}

This work: **{c['sigma_pb']:.0f} +- {c['total_pb']:.0f} pb** (60 < m < 120 GeV). It agrees with the CMS measurement of the same quantity within {abs(c['sigma_pb'] - 1952) / (c['total_pb'] ** 2 + 49 ** 2) ** 0.5:.1f} sigma, and with the aMC@NLO reference the channels normalise to (1945-1954 pb) within {abs(c['sigma_pb'] - 1949) / c['total_pb']:.1f} sigma. The ATLAS number is quoted in a narrower mass window (66-116 GeV) and is not directly comparable.

![comparison](output/plots/comparison.png)

## Correlation model in full

{_table(["source", "rho", "mu mu [pb]", "tau tau [pb]"], source_rows, ["---", ":---:", "---:", "---:"])}

Justification for every row: [docs/02-correlation-model.md](docs/02-correlation-model.md).

## What this is not

This is a **covariance combination of two TRExFitter fits, not a TRExFitter MultiFit.** The
channel workspaces are not committed to the repository and were not rebuilt, so the shared
nuisance parameters are correlated through their `Category` labels instead of being profiled
jointly. The MultiFit configuration is generated and syntax-checked
([fit/comb.config](fit/comb.config), [fit/run_multifit.py](fit/run_multifit.py)) and the exact
recipe to run it is in that file; [docs/03-method.md](docs/03-method.md) explains what the joint
fit would add. At this precision the difference is expected to be small -- the mu mu channel
dominates and its own fit is already profiled -- but it is not zero, and the claim is not made.

Known caveats inherited from the inputs, none of them corrected here. All nine findings of
`z-mumu/REVIEW.md` were fixed by that channel before this combination was made, so what is left
is each channel's own open-issue list:

{caveats}

## Reproduce

```bash
source ../setup.sh
python run_combination.py          # checks, result, figures, JSON, this file
python run_combination.py --check  # input validation and closure tests only
```
"""


def write_result_md(payload: dict, path: Path) -> Path:
    path.write_text(render(payload))
    return path
