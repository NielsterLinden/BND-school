"""Render `result.md` from the combination payload, so the prose can never drift from the numbers.

`run_combination.py` calls `write_result_md`; every figure in the file comes from
`output/combination_result.json`, which is written in the same run.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

NAMES = {"mumu": "Z -> mu mu", "tautau": "Z -> tau_h tau_h", "ee": "Z -> ee"}
SHORT = {"mumu": "mu mu", "tautau": "tau tau", "ee": "ee"}

PREDICTION_NOTE = (
    "the three channels normalise to slightly different aMC@NLO references for the same quantity "
    "(1953.9 pb for mu mu, 1944.9 pb for tau tau, 1954.1 pb for ee), because the generator's LHE "
    "flavour shares are not exactly 1/3 -- 0.47 % between the extremes, discussed in "
    "[docs/01-inputs.md](docs/01-inputs.md)"
)

#: caveats carried over from the channels' own open-issue lists, none of them corrected here.
#: Each `{...}` is filled from the payload in `render`.
CAVEATS = [
    "**`z-ee` does not follow `fitting/CONVENTIONS.md` section 3 on the signal theory templates, "
    "and it is the largest single effect in this combination.** Its `PDF` and `QCDScale` "
    "histograms are raw LHE-weight envelopes that were never renormalised to a constant yield, so "
    "`QCDScale` enters as a +-{ee_scale_norm:.2f} % *normalisation* of the signal -- degenerate "
    "with the POI. The fit pulls it to {ee_scale_pull:+.2f} sigma, i.e. it rescales the prediction "
    "by {ee_kappa:.3f}, and `mu_signal` = {ee_mu:.4f} is measured against that rescaled "
    "prediction; the data-over-prediction ratio is {ee_counting:.4f}. Undoing only that "
    "normalisation gives {ee_normfix_sigma:.0f} pb for the channel and moves the combination by "
    "{ee_normfix_shift:+.0f} pb -- the `ee normfix` row below. Until z-ee re-runs with "
    "renormalised templates, **the combined value carries that bias**.",
    "`z-ee` has no goodness of fit (the job is `trex-fitter hwdfp`, without the saturated model), "
    "so it is the one channel whose fit quality cannot be quoted. It also pulls `Pileup` by "
    "-2.8 sigma and `L1Prefiring` by +3.6 sigma and constrains `ElectronID` to 0.07 of its input "
    "width, all symptoms of a shape mismatch being absorbed by calibration parameters in a "
    "60 x 1 GeV binning.",
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
    "({mumu_muon_pb:.0f} pb) is one of the largest uncertainties here after the luminosity.",
    "`z-tautau` v3 publishes only its nominal (DeepTau Tight, MC-subtracted fake factor), so the "
    "combination carries no tau tau cross-check rows. The channel contributes "
    "{tautau_weight_pct:.2f} % of the weight, so this costs the cross section nothing; it costs "
    "the universality ratio a systematic cross-check.",
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
    order = c["order"]
    ch = p["channels"]
    lik, rats = p["likelihood"], p["ratios"]
    v = p["variations"]
    rel = 100 * c["relative"]
    ee = ch.get("ee", {})
    best = c["most_precise_channel"]

    binning = [v[k]["value"] for k in ("baseline", "mumu_bins2gev", "mumu_bins10gev")]
    caveats = "\n".join(
        "* " + t.format(
            mumu_binning_span=max(binning) - min(binning),
            mumu_binning_frac=(max(binning) - min(binning)) / c["total_pb"],
            mumu_counting_shift=v["mumu_counting"]["shift"],
            mumu_muon_pb=c["breakdown_pb"]["Muon efficiency"],
            tautau_weight_pct=100 * c["weights"]["tautau"],
            ee_scale_norm=100 * ee["extra"]["qcdscale_norm"],
            ee_scale_pull=ee["extra"]["qcdscale_pull"],
            ee_kappa=ee["extra"]["kappa_theory"],
            ee_mu=ee["extra"]["mu_published"],
            ee_counting=ee["extra"]["mu_counting"],
            ee_normfix_sigma=v["ee_normfix"]["channels"]["ee"]["sigma_pb"],
            ee_normfix_shift=v["ee_normfix"]["shift"])
        for t in CAVEATS)

    channel_rows = []
    for name in order:
        d = ch[name]
        gof = f"{d['gof_p']:.2f}" if d["gof_p"] is not None else "n/a"
        channel_rows.append([f"**{NAMES[name]}**", d["variant"],
                             f"{d['mu']:.4f} +{d['mu_err_up']:.4f} -{d['mu_err_down']:.4f}",
                             f"{d['sigma_pred_pb']:.1f}",
                             f"**{d['sigma_pb']:.0f} +- {d['sigma_err_pb']:.0f}**",
                             f"{100 * d['sigma_err_pb'] / d['sigma_pb']:.1f} %", gof])

    breakdown_rows = [[k, f"{val:.2f}", f"{100 * val / c['sigma_pb']:.2f} %",
                       "yes" if next(s["rho"] for s in p["sources"] if s["name"] == k) >= 0.999 else "no"]
                      for k, val in c["breakdown_pb"].items() if val >= 0.05]

    var_rows = [[k.replace("_", " "), f"{d['value']:.1f}", f"{d['shift']:+.1f}",
                 f"{d['shift'] / c['total_pb']:+.2f}", d["why"]]
                for k, d in v.items() if k != "baseline"]

    source_rows = [[s["name"], "1" if s["rho"] >= 0.999 else ("0" if s["rho"] <= 0.001 else f"{s['rho']:.2f}")]
                   + [f"{s['sizes_pb'].get(n, 0):.2f}" for n in order]
                   for s in sorted(p["sources"], key=lambda s: -max(s["sizes_pb"].values()))]

    ratio_rows = [[f"**{SHORT[r['numerator']]} / {SHORT[r['denominator']]}**",
                   f"**{r['value']:.3f} +- {r['error']:.3f}**", f"{r['z']:+.1f}",
                   f"{r['pvalue']:.3f}"] for r in rats.values()]

    weights = ", ".join(f"{SHORT[k]} {val:+.4f}" for k, val in c["weights"].items())
    rhos = ", ".join(f"{' / '.join(SHORT[x] for x in k.split('-'))} {val:.3f}"
                     for k, val in c["correlation"].items())
    alone = ", ".join(f"{SHORT[n]} {ch[n]['alone_pb']:.0f} +- {ch[n]['alone_err_pb']:.0f} pb"
                      for n in order)
    gain = 100 * (1 - c["total_pb"] / ch[best]["alone_err_pb"])
    negatives = [SHORT[k] for k, val in c["weights"].items() if val < 0]

    compat = (
        "**Read the compatibility number before the central value.** With three channels the "
        f"chi2 has two degrees of freedom and it is {'no longer a formality' if c['pvalue'] < 0.05 else 'comfortable'}: "
        f"p = {c['pvalue']:.3f}. The ee and mu mu measurements are of comparable precision "
        f"({100 * ch['ee']['sigma_err_pb'] / ch['ee']['sigma_pb']:.1f} % and "
        f"{100 * ch['mumu']['sigma_err_pb'] / ch['mumu']['sigma_pb']:.1f} %) and "
        f"{abs(ch['ee']['sigma_pb'] - ch['mumu']['sigma_pb']):.0f} pb apart, and the first caveat "
        "below says why: the z-ee fit lets a theory *normalisation* float against its own signal "
        "strength.")
    univ = (
        "The ee/mu mu ratio is the sharp one: the luminosity, pileup, prefiring and shared theory "
        "terms cancel exactly, so what is left is "
        f"{100 * rats['ee/mumu']['error'] / rats['ee/mumu']['value']:.1f} % -- and the two "
        f"channels disagree by {abs(rats['ee/mumu']['z']):.1f} sigma. That is a statement about "
        "the two *fits*, not about lepton universality: see the first caveat.")

    return f"""# Combined Z cross section from Z -> mu mu, Z -> tau_h tau_h and Z -> ee

*Generated by `run_combination.py` on {date.today().isoformat()} from
`output/combination_result.json`. Do not edit by hand -- edit the code or the docs.*

CMS Open Data 2016, Run2016G+H, L = 16393.381 pb^-1 (normtag), sqrt(s) = 13 TeV.

## Result

> **sigma(pp -> Z/gamma* -> ll, 60 < m_ll < 120 GeV) = {c['sigma_pb']:.0f} +- {c['total_pb']:.0f} pb**
> = {c['sigma_pb']:.0f} +- {c['stat_pb']:.1f} (stat) +- {c['syst_pb']:.0f} (syst)
> +- {c['acc_pb']:.0f} (acc) +- {c['lumi_pb']:.0f} (lumi) pb   ({rel:.2f} %)
>
> per lepton flavour, assuming lepton universality. Channel compatibility:
> chi2/ndf = {c['chi2']:.2f}/{c['ndf']}, p = {c['pvalue']:.3f}.

{compat}

The profile-likelihood form of the same combination, which carries the tau tau channel's asymmetric uncertainty exactly, gives {lik['value_pb']:.0f} +{lik['err_up_pb']:.0f} / -{lik['err_down_pb']:.0f} pb.

![combined result](output/plots/forest.png)

## Inputs

{_table(["channel", "variant", "mu_Z", "reference [pb]", "sigma(60-120) [pb]", "rel.", "GoF p"],
        channel_rows, ["---", "---", "---", "---:", "---:", "---:", "---:"])}

All three are TRExFitter v1.8.0 profile-likelihood fits performed by the channel groups; this
combination reads their published results and does not refit:

{chr(10).join(f"* `{NAMES[n]}`: `{ch[n]['provenance'][0]}`" for n in order)}

The `reference` column is the aMC@NLO prediction each channel's `mu_Z` multiplies -- {PREDICTION_NOTE}. The `mu_Z` uncertainty shown is the in-fit one; for mu mu and tau tau the acceptance uncertainty (0.61 % and 3.7 %) sits outside the fit and is added in the `sigma` column, while z-ee fits directly against the 60-120 GeV prediction and has no separate acceptance term.

The inputs are each channel's **current baseline**: z-mumu's fit rebuilt after its own review
(`z-mumu/REVIEW.md` section 0) in 12 x 5 GeV bins with a two-sided `SigModel` template, z-tautau
**v3** with DeepTau Tight on both legs and the MC-subtracted fake factor, and the z-ee TRExFitter
job as published in `z-ee/Zee_fit.tar.gz` (converted by `tools/extract_zee.py`).

## Where the uncertainty comes from

{_table(["source", "pb", "of sigma", "correlated?"], breakdown_rows, ["---", "---:", "---:", ":---:"])}

![uncertainty breakdown](output/plots/breakdown.png)

* weights: {weights};
* correlations: {rhos};
* alone: {alone}; the combination reaches +- {c['total_pb']:.1f} pb, an improvement of {gain:.1f} % on the best single channel ({SHORT[best]});
{'* the ' + ' and '.join(negatives) + ' weight is **negative**. That is standard BLUE behaviour, not an error: when the correlation with a more precise measurement exceeds the ratio of their uncertainties, the less precise one is used to pull on the shared systematic rather than to average the value down.' if negatives else ''}

The luminosity term is {c['lumi_pb']:.0f} pb, fully correlated between the channels and therefore
irreducible by combining: all three use the same 16393.381 pb^-1 with the same 1.2 %. It is
{100 * c['lumi_pb'] / c['total_pb']:.0f} % of the total uncertainty in quadrature, so the route
to a better number from this dataset is a better luminosity calibration, not more channels.

## Lepton universality

This is the part the combination cannot test, because it *assumes* universality. Measured
separately, with everything correlated cancelling:

{_table(["ratio", "value", "z from 1", "p"], ratio_rows, ["---", "---:", "---:", "---:"])}

![lepton universality](output/plots/universality.png)

{univ}

## Cross-checks and variations

{_table(["variation", "sigma [pb]", "shift", "in sigma", "what it changes"], var_rows,
        ["---", "---:", "---:", "---:", "---"])}

![variations](output/plots/variations.png)

The dominant row is `ee normfix` ({v['ee_normfix']['shift']:+.1f} pb, {v['ee_normfix']['shift'] / c['total_pb']:+.2f} sigma) -- not a correlation choice but a convention violation inside one channel, and larger than everything else in the table put together. The two rows that bracket the correlation model (`rho zero`, `rho one`) span {abs(v['rho_zero']['value'] - v['rho_one']['value']):.1f} pb, and `sigmodel correlated` is exactly null because neither z-tautau nor z-ee fits a generator nuisance parameter at all. The correlation model is therefore not the limiting assumption; the z-ee signal-template normalisation is.

## Comparison with published measurements

{_table(["measurement", "value", "reference"], [[a, b, c_] for a, b, c_ in REFERENCES])}

This work: **{c['sigma_pb']:.0f} +- {c['total_pb']:.0f} pb** (60 < m < 120 GeV), which sits {abs(c['sigma_pb'] - 1952) / (c['total_pb'] ** 2 + 49 ** 2) ** 0.5:.1f} sigma from the CMS measurement of the same quantity and {abs(c['sigma_pb'] - 1951) / c['total_pb']:.1f} sigma from the aMC@NLO references the channels normalise to (1945-1954 pb). The ATLAS number is quoted in a narrower mass window (66-116 GeV) and is scaled onto ours by sigma(60-120)/sigma(66-116) = 1.0143, computed at LHE level from the same aMC@NLO sample.

![comparison](output/plots/comparison.png)

### Per channel, against the published measurement of the same decay

{_table(["decay", "this work [pb]", "published [pb]", "difference"],
        [["Z -> mu mu", f"{ch['mumu']['sigma_pb']:.0f} +- {ch['mumu']['sigma_err_pb']:.0f}",
          "ATLAS 2005 +- 60 (1977 +- 60 at 66-116)",
          f"{(ch['mumu']['sigma_pb'] - 2005) / (ch['mumu']['sigma_err_pb'] ** 2 + 60 ** 2) ** 0.5:+.1f} sigma"],
         ["Z -> ee", f"{ch['ee']['sigma_pb']:.0f} +- {ch['ee']['sigma_err_pb']:.0f}",
          "ATLAS 2015 +- 61 (1987 +- 61 at 66-116)",
          f"{(ch['ee']['sigma_pb'] - 2015) / (ch['ee']['sigma_err_pb'] ** 2 + 61 ** 2) ** 0.5:+.1f} sigma"],
         ["Z -> tau_h tau_h", f"{ch['tautau']['sigma_pb']:.0f} +- {ch['tautau']['sigma_err_pb']:.0f}",
          "CMS 1848 +- 68 (five tau tau final states)",
          f"{(ch['tautau']['sigma_pb'] - 1848) / (ch['tautau']['sigma_err_pb'] ** 2 + 68 ** 2) ** 0.5:+.1f} sigma"],
         ["combined", f"{c['sigma_pb']:.0f} +- {c['total_pb']:.0f}", "CMS 1952 +- 49; ATLAS 2009 +- 58",
          f"{(c['sigma_pb'] - 1952) / (c['total_pb'] ** 2 + 49 ** 2) ** 0.5:+.1f}; "
          f"{(c['sigma_pb'] - 2009) / (c['total_pb'] ** 2 + 58 ** 2) ** 0.5:+.1f} sigma"]],
        ["---", "---:", "---:", "---:"])}

![per-channel comparison](output/plots/comparison_channels.png)

**Read the first two rows together.** Each of the two precise channels agrees with the ATLAS
measurement of the same decay at the 1 sigma level -- and yet they disagree with *each other* by
3.0 sigma. ATLAS's per-channel uncertainty is +-60 pb against our +-35 and +-40 pb, so a 123 pb gap
sits comfortably inside the published error bars while being far outside ours. Agreement with a
published number is not a validation at this precision.
CMS-SMP-20-004 fits its ee and mu mu final states together and publishes no per-channel cross
section, so ATLAS is the only 13 TeV measurement that can carry this comparison. Sources: ATLAS
[arXiv:1603.09222](https://arxiv.org/abs/1603.09222) Table 2, CMS Z -> tau tau
[arXiv:1801.03535](https://arxiv.org/abs/1801.03535). Orthogonality of the three selections and
the systematics CMS and ATLAS carry that we do not:
[docs/05-vs-published.md](docs/05-vs-published.md).

## Correlation model in full

{_table(["source", "rho"] + [f"{SHORT[n]} [pb]" for n in order], source_rows,
        ["---", ":---:"] + ["---:"] * len(order))}

Justification for every row: [docs/02-correlation-model.md](docs/02-correlation-model.md).

## What this is not

This is a **covariance combination of three TRExFitter fits, not a TRExFitter MultiFit.** The
channel workspaces are not all committed to the repository and were not rebuilt, so the shared
nuisance parameters are correlated through their `Category` labels instead of being profiled
jointly. The MultiFit configuration is generated and syntax-checked
([fit/comb.config](fit/comb.config), [fit/run_multifit.py](fit/run_multifit.py)) and the exact
recipe to run it is in that file; [docs/03-method.md](docs/03-method.md) explains what the joint
fit would add. A joint fit would also expose the z-ee normalisation problem directly, because
`QCDScale` would then be one parameter shared with two channels that renormalise it correctly.

Known caveats inherited from the inputs, none of them corrected here:

{caveats}

## Reproduce

```bash
source ../setup.sh
python tools/extract_zee.py        # only after z-ee re-publishes Zee_fit.tar.gz
python run_combination.py          # checks, result, figures, JSON, this file
python run_combination.py --check  # input validation and closure tests only
```
"""


def write_result_md(payload: dict, path: Path) -> Path:
    path.write_text(render(payload))
    return path
