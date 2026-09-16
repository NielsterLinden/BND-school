"""Load the published per-channel results and turn them into a common structure.

Nothing here recomputes physics: every number is read from a file committed to the repository
and is traceable to the channel that produced it (`provenance` on each `ChannelResult`).

    from comb import inputs
    chans = inputs.load_channels()                     # {"mumu": ..., "tautau": ..., "ee": ...}
    chans = inputs.load_channels(mumu="counting", ee="normfix")
    chans = inputs.load_channels(only=["mumu", "tautau"])

All three channels have now published a result and all three are in the combination:

* **z-mumu** fixed the findings of its own review (`z-mumu/REVIEW.md` section 0): the fit runs in
  12 x 5 GeV bins with a two-sided `SigModel` template built inside a common 50 < m_LHE < 120 GeV
  window, no template smoothing, MINOS on every parameter. That fit is stable to +-0.2 % across
  the binnings that describe the data, so the counting extraction and the +-0.7 % lineshape term
  the review had recommended as a stop-gap are **no longer used**: the shape fit is the baseline
  (`mumu="shapefit"`) and the channel's stability table supplies the alternative configurations.
* **z-tautau v3** makes DeepTau **Tight** on both legs its nominal working point, with the
  MC-subtracted fake factor (`nominal_variant`).
* **z-ee** published its TRExFitter job rather than a result file; `tools/extract_zee.py` turns
  it into `inputs/zee_fit_result.json`, which is what `load_ee` reads.

Three things are deliberately *not* symmetric between the channels:

* the three `mu_Z` do not share a reference prediction (1953.9 pb for mumu, 1944.9 pb for tautau,
  1954.1 pb for ee -- see docs/01-inputs.md), so the combination is done on the cross sections,
  never on `mu_Z`;
* z-mumu fits a `SigModel` nuisance parameter (powheg vs aMC@NLO, both NLO); z-tautau and z-ee do
  not. z-tautau's generator comparison (`SigModel_tautau`, madgraph LO vs aMC@NLO) is reported and
  *not* used as an uncertainty, and `z-tautau/docs/07` forbids correlating the two.
  `ChannelResult.sigmodel` is therefore zero on the tautau and ee sides, which is what makes the
  split in `comb/model.py` a no-op there;
* z-mumu and z-tautau measure a *fiducial* cross section and divide by an acceptance whose
  uncertainty sits outside the fit; z-ee fits directly against the full 60-120 GeV prediction, so
  its acceptance uncertainty is inside the likelihood (`PDF`, `QCDScale`) and `ChannelResult.acc`
  is empty. docs/01-inputs.md says what that costs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# --- files read (all committed) -------------------------------------------------------------
ZMUMU_FIT_JSON = REPO / "z-mumu/fit/results/zmumu_fit_result.json"
ZMUMU_STABILITY = REPO / "z-mumu/fit/results/stability.json"
ZTAUTAU_RESULTS = REPO / "z-tautau/output/results.json"
ZEE_RESULTS = REPO / "combination/inputs/zee_fit_result.json"

#: z-mumu fit configurations: name used here -> `tag` of the entry in `stability.json`.
#: `None` is the nominal fit, which is the whole of `zmumu_fit_result.json`.
#: Of the seven rows of that table only the ones the channel accepts are offered: the nominal
#: (GoF p = 0.79), 6 x 10 GeV (0.50), 30 x 2 GeV (0.16) and the 1-bin counting extraction.
#: `stab_nosig` and `stab_smooth` have p <= 0.01 and are rejected by the channel itself
#: (`z-mumu/REVIEW.md` section 0), so they are not available as combination variants.
ZMUMU_VARIANTS = {"shapefit": None, "bins2gev": "stab_2gev", "bins10gev": "stab_10gev",
                  "counting": "stab_1bin"}

#: z-tautau variants: name used here -> key in the `fit` block of its results file.
#: `None` means "whatever that file calls `nominal_variant`", so the channel stays in charge of
#: which of its fake-factor treatments is nominal. Since v3 that is `mcsub` at the DeepTau Tight
#: working point, and it is the only variant the channel publishes in full.
ZTAUTAU_VARIANTS = {"nominal": None, "mcsub": "mcsub"}

#: z-ee extractions: name used here -> key in the `fit` block of `inputs/zee_fit_result.json`.
#: `published` is the channel's own number; `normfix` is the same fit with the normalisation of
#: its un-renormalised `PDF`/`QCDScale` templates put back into the prediction, which is what
#: `fitting/CONVENTIONS.md` section 3 asks every channel to do before the fit. See docs/01.
ZEE_VARIANTS = {"published": None, "normfix": "normfix"}

#: nuisance parameters that make up the "Signal modelling" category in both channels
#: (fitting/CONVENTIONS.md). `SigModel` is split off because the two channels compare different
#: generator pairs -- see comb/model.py.
SIGNAL_MODELLING_NPS = ["PDF", "AlphaS", "QCDScale", "PS_ISR", "PS_FSR", "SigModel"]


@dataclass
class ChannelResult:
    """One channel's measurement of sigma(Z/gamma* -> ll, 60 < m < 120 GeV), in pb."""

    name: str                       # "mumu" / "tautau" / "ee"
    label: str                      # for plots and tables
    variant: str                    # which extraction / fake-factor variant this is
    mu: float                       # signal strength w.r.t. that channel's own prediction
    mu_err_up: float
    mu_err_down: float
    mu_stat: float
    sigma_pred: float               # sigma^pred(60 < m < 120 GeV) this channel normalises to [pb]
    groups: dict[str, float]        # Category -> impact on mu (absolute, symmetric)
    ranking: dict[str, float]       # NP -> symmetrised impact on mu (used to split categories)
    acc: dict[str, float]           # acceptance component -> relative uncertainty on sigma
    sigmodel: float = 0.0           # impact on mu of `SigModel` alone (split out of the category)
    residual: float = 0.0           # impact on mu that the published grouped impacts do not cover
    gof_p: float | None = None
    provenance: list[str] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    # -- derived ----------------------------------------------------------------------------
    @property
    def sigma(self) -> float:
        """Measured sigma(Z/gamma* -> ll, 60 < m < 120 GeV) in pb."""
        return self.mu * self.sigma_pred

    @property
    def sigma_stat(self) -> float:
        return self.mu_stat * self.sigma_pred

    def group_pb(self, category: str) -> float:
        """In-fit systematic `category` as an absolute uncertainty on sigma [pb].

        In-fit impacts are impacts on mu and mu multiplies a fixed prediction, so they scale with
        sigma_pred, not with the measured sigma. This is the same bookkeeping the channels use
        (z-mumu: sigma_fid_lumi_pb = 0.0116267 x 799.566).
        """
        return self.groups.get(category, 0.0) * self.sigma_pred

    def acc_pb(self, component: str) -> float:
        """Acceptance component as an absolute uncertainty on sigma [pb].

        A enters as sigma = sigma_fid / A, so its uncertainty is relative to the *measured* value.
        """
        return self.acc.get(component, 0.0) * self.sigma

    @property
    def residual_pb(self) -> float:
        """The part of the published total uncertainty the grouped impacts do not account for [pb].

        TRExFitter's grouped impacts are quadrature differences and do not in general add up to
        the total MINOS error. z-mumu and z-tautau over-shoot (their categories sum to more than
        their published total, which this folder keeps, conservatively); z-ee *under*-shoots by
        1.03 % of mu_Z, which its results file carries as `fit_residual`. It is uncorrelated
        between channels either way -- it is a property of one fit's own correlation matrix.
        """
        return self.residual * self.sigma_pred

    @property
    def sigma_err_total(self) -> float:
        terms = [self.sigma_stat, self.residual_pb]
        terms += [self.group_pb(c) for c in self.groups]
        terms += [self.acc_pb(k) for k in self.acc]
        return sum(t * t for t in terms) ** 0.5


# ---------------------------------------------------------------------------- parsers
def _drop_totals(groups: dict[str, float]) -> dict[str, float]:
    """Remove the summary rows and empty categories from a grouped-impact table."""
    return {k: v for k, v in groups.items() if k not in ("FullSyst", "Total") and v > 0.0}


def _ranking_mumu(js: dict) -> dict[str, float]:
    return {r["name"]: 0.5 * (abs(r["dpoi_up_post"]) + abs(r["dpoi_down_post"])) for r in js.get("ranking", [])}


def _ranking_tautau(fit: dict) -> dict[str, float]:
    return {r["name"]: 0.5 * (abs(r["impact_up"]) + abs(r["impact_down"])) for r in fit.get("ranking", [])}


def _stability_entry(tag: str) -> dict:
    """One row of z-mumu's fit-configuration stability table."""
    rows = json.loads(ZMUMU_STABILITY.read_text())
    for row in rows:
        if row["tag"] == tag:
            return row
    raise RuntimeError(f"no entry {tag!r} in {ZMUMU_STABILITY}: have {[r['tag'] for r in rows]}")


# ---------------------------------------------------------------------------- channels
def load_mumu(variant: str = "shapefit") -> ChannelResult:
    """Z -> mu mu. `variant` is a key of `ZMUMU_VARIANTS`.

    "shapefit" is the channel's own post-review baseline: the 12 x 5 GeV profile-likelihood fit of
    m(mu mu) with the two-sided `SigModel` template, MINOS everywhere and no smoothing.

    The other keys are rows of the channel's stability table (`z-mumu/fit/results/stability.json`,
    produced by `scripts/v2_5_fit_variants.py`). That table publishes mu_Z, the MINOS errors, the
    total systematic and the goodness of fit per configuration, but **not** the grouped impacts
    per configuration, so a variant keeps the nominal fit's category composition rescaled to its
    own published total systematic. That is an approximation and it is only ever used for the
    variations in `run_combination.VARIATIONS`, never for the baseline.
    """
    if variant not in ZMUMU_VARIANTS:
        raise ValueError(f"unknown z-mumu variant {variant!r}: use {sorted(ZMUMU_VARIANTS)}")
    js = json.loads(ZMUMU_FIT_JSON.read_text())
    meta = js["meta"]
    sigma_pred = meta["sigma_fid_pred_pb"] / meta["A_60_120"]     # 1953.9 pb
    a = meta["acceptance"]
    acc = {"pdf": a["A_pdf_rel"], "alphas": a["A_alphas_rel"], "scale": a["A_scale_rel"],
           "mcstat": a["A_stat"] / a["A"]}
    prov = [str(ZMUMU_FIT_JSON.relative_to(REPO))]

    ranking = _ranking_mumu(js)
    groups = _drop_totals(js["grouped_impacts_mu"])
    syst_nominal = js["grouped_impacts_mu"]["FullSyst"]
    tag = ZMUMU_VARIANTS[variant]

    if tag is None:
        mu, up, down = js["mu"], js["mu_err_up"], js["mu_err_down"]
        stat = js["mu_stat_only_fit"]
        gof = js.get("gof", {}).get("gof_probability")
        # a shape effect: take it from the post-fit ranking of the nuisance parameter itself
        sigmodel = ranking.get("SigModel", 0.0)
    else:
        row = _stability_entry(tag)
        mu, up, down = row["mu"], row["err_up"], row["err_down"]
        # the data statistical uncertainty is a property of the dataset, not of the binning:
        # sqrt(N_obs) / (N_obs - N_bkg) = 0.031 %, which is what the nominal fit's stat-only fit
        # returns as well (asserted in run_combination.check).
        c = meta["counting"]
        stat = mu * c["n_obs"] ** 0.5 / (c["n_obs"] - c["n_bkg"])
        scale = row["syst_total"] / syst_nominal
        groups = {k: v * scale for k, v in groups.items()}
        gof = row["gof_p"] or None          # the 1-bin fit has no degrees of freedom left
        # in the counting extraction the template's shape cannot act: only the powheg/aMC@NLO
        # difference in the fiducial C factor survives
        sigmodel = (abs(meta["sigmodel"]["C_ratio_powheg_over_nlo"] - 1.0) * mu if tag == "stab_1bin"
                    else ranking.get("SigModel", 0.0) * scale)
        prov.append(f"{ZMUMU_STABILITY.relative_to(REPO)} [{tag}: {row['label']}]")

    return ChannelResult(
        name="mumu", label=r"$Z\to\mu\mu$", variant=variant,
        mu=mu, mu_err_up=up, mu_err_down=down, mu_stat=stat,
        sigma_pred=sigma_pred, groups=groups, ranking=ranking, acc=acc, sigmodel=sigmodel, gof_p=gof,
        provenance=prov,
        extra={"sigma_fid_pb": mu * meta["sigma_fid_pred_pb"],
               "sigma_fid_pred_pb": meta["sigma_fid_pred_pb"],
               "A": meta["A_60_120"], "C": meta["C_factor"], "lumi_pb": meta["lumi_pb"],
               "n_obs": meta["counting"]["n_obs"], "n_bkg": meta["counting"]["n_bkg"],
               "sr_bin_width_gev": meta["sr_bin_width_gev"],
               "shapefit_mu": js["mu"], "counting_sigma_fid_pb": meta["counting"]["sigma_fid_pb"]})


def load_tautau(variant: str = "nominal") -> ChannelResult:
    """Z -> tau_h tau_h. `variant` is a key of `ZTAUTAU_VARIANTS`.

    "nominal" resolves to whatever the channel's results file calls `nominal_variant`, so the
    channel stays in charge of its own baseline. Since v3 that is `mcsub` -- the fake factor with
    the genuine-tau MC subtraction -- measured at the DeepTau **Tight** working point on both legs,
    which v2.1 had recommended as its next iteration and v3 adopted.
    """
    if variant not in ZTAUTAU_VARIANTS:
        raise ValueError(f"unknown z-tautau variant {variant!r}: use {sorted(ZTAUTAU_VARIANTS)}")
    js = json.loads(ZTAUTAU_RESULTS.read_text())
    key = ZTAUTAU_VARIANTS[variant] or js["nominal_variant"]
    if js["fit"].get(key) is None:
        raise ValueError(f"{ZTAUTAU_RESULTS.name} publishes no complete fit for {key!r}: "
                         f"have {sorted(k for k, v in js['fit'].items() if v)}")
    fit = js["fit"][key]
    pred = js["prediction"]
    acc = {"pdf": pred["A_unc"]["A_pdf"], "alphas": pred["A_unc"]["A_alphas"],
           "scale": pred["A_unc"]["A_scale"], "isr": pred["A_unc"]["A_isr"],
           "fsr": pred["A_unc"]["A_fsr"], "mcstat": pred["A_mc_stat"]}
    ranking = _ranking_tautau(fit)
    return ChannelResult(
        name="tautau", label=r"$Z\to\tau_h\tau_h$", variant=key,
        mu=fit["mu"], mu_err_up=fit["mu_err_up"], mu_err_down=fit["mu_err_down"],
        mu_stat=fit["mu_stat"],
        sigma_pred=pred["sigma_tautau_60_120_pb"], groups=_drop_totals(fit["grouped_impact"]),
        ranking=ranking, acc=acc,
        # `SigModel_tautau` is reported, not fitted (z-tautau/docs/07); since v2.1 no `SigModel`
        # appears in the ranking, so the generator part of the category is zero by construction.
        sigmodel=ranking.get("SigModel", 0.0), gof_p=fit.get("gof_probability"),
        provenance=[f"{ZTAUTAU_RESULTS.relative_to(REPO)} [fit/{key}]"],
        extra={"sigma_fid_pb": fit["sigma_fid_pb"]["value"],
               "sigma_fid_pred_pb": pred["sigma_fid_pb"], "A": pred["A"],
               "C": js["for_combination"]["C"], "lumi_pb": js["lumi_pb"],
               "n_obs": js["for_combination"]["n_obs"],
               "n_bkg": js["for_combination"]["n_bkg_prefit"],
               "tau_wp": js.get("tau_wp"),
               "sigmodel_C_LO_over_NLO": js["sigmodel"]["C_LO_over_NLO_fiducial"],
               "mu_expected_asimov": fit["mu_expected_asimov"]})


def load_ee(variant: str = "published") -> ChannelResult:
    """Z -> e+ e-. `variant` is a key of `ZEE_VARIANTS`.

    The z-ee group published the TRExFitter job itself (`z-ee/Zee_fit.tar.gz`) rather than a
    results file; `tools/extract_zee.py` turns it into `inputs/zee_fit_result.json` and this
    function reads that. Two extractions are offered and the difference between them is the
    largest single choice in the whole combination:

    * ``published`` -- `mu_signal` exactly as the channel quotes it, 1.0488 +- 0.0190.
    * ``normfix``   -- the same fit, with the normalisation of the `PDF` and `QCDScale` templates
      put back into the prediction. z-ee built those two from LHE weight envelopes *without*
      renormalising them to a constant yield, so they carry the +-5.9 % scale uncertainty of the
      aMC@NLO cross section as a pure normalisation on the signal -- which
      `fitting/CONVENTIONS.md` section 3 forbids because it is degenerate with the POI. The fit
      pulls `QCDScale` to -1.85 sigma, i.e. it rescales the *prediction* by 0.894, and the
      published `mu_signal` is measured against that rescaled prediction. The cross section the
      measured yield corresponds to is mu_hat x kappa_theory x sigma^pred; docs/01-inputs.md
      derives it and `tools/extract_zee.py` computes kappa_theory.

    Two structural differences from the other two channels, both discussed in docs/01:

    * **no acceptance term.** z-ee fits directly against the 60 < m_LHE < 120 GeV prediction
      instead of quoting a fiducial cross section and dividing by A, so the theory uncertainty on
      the extrapolation is inside the likelihood (`Signal modelling`) rather than outside it.
      `acc` is therefore empty -- not missing, *inside the fit*.
    * **the data statistics are the analytic ones.** TRExFitter's `Stat unc.` row for this fit is
      1.03 % of mu_Z, 25x the Poisson value sqrt(N_obs)/(N_obs - N_bkg) = 0.040 % that 6.3 M
      selected events give. It is the quadrature remainder of the total after the grouped impacts,
      not a measurement of the data statistics, so it is carried as `residual` instead of being
      called "Data statistics". The total uncertainty is unchanged; only the label is.
    """
    if variant not in ZEE_VARIANTS:
        raise ValueError(f"unknown z-ee variant {variant!r}: use {sorted(ZEE_VARIANTS)}")
    js = json.loads(ZEE_RESULTS.read_text())
    key = ZEE_VARIANTS[variant] or js["nominal_variant"]
    fit = js["fit"][key]
    ranking = {k: v for k, v in fit["ranking"].items() if not k.startswith("gamma_")}
    return ChannelResult(
        name="ee", label=r"$Z\to ee$", variant=key,
        mu=fit["mu"], mu_err_up=fit["mu_err_up"], mu_err_down=fit["mu_err_down"],
        mu_stat=fit["mu_stat"], residual=fit["fit_residual"],
        sigma_pred=js["sigma_pred_60_120_pb"], groups=_drop_totals(fit["grouped_impact"]),
        ranking=ranking, acc={},
        sigmodel=ranking.get("SigModel", 0.0),
        # no saturated-model goodness of fit was run for this job (`trex-fitter hwdfp`, no `s`)
        gof_p=None,
        provenance=[f"{ZEE_RESULTS.relative_to(REPO)} [fit/{key}]",
                    js["source"]["tarball"], js["source"]["gensums"]],
        extra={"lumi_pb": js["lumi_pb"], "n_obs": js["for_combination"]["n_obs"],
               "n_bkg": js["for_combination"]["n_bkg_prefit"],
               "n_sig": js["for_combination"]["n_sig_prefit"],
               "acc_eff": js["for_combination"]["acc_eff"],
               "mu_counting": js["for_combination"]["mu_counting"],
               "mu_published": js["fit"]["published"]["mu"],
               "kappa_theory": js["theory_normalisation"]["kappa"],
               "qcdscale_pull": js["nuisance_parameters"]["QCDScale"]["pull"],
               "qcdscale_norm": js["nuisance_parameters"]["QCDScale"]["norm_up"],
               "stat_trex": js["error_decomposition"]["STAT_ERROR"],
               "sigma_counting_pb": js["for_combination"]["mu_counting"] * js["sigma_pred_60_120_pb"]})


def load_channels(mumu: str = "shapefit", tautau: str = "nominal", ee: str = "published",
                  only: list[str] | None = None) -> dict[str, ChannelResult]:
    """The three channels, in the order they enter the combination.

    `only` restricts the combination to a subset -- used by the `no_ee` / `mumu_tautau` variations
    of `run_combination.py`, never for the baseline.
    """
    loaders = {"mumu": lambda: load_mumu(mumu), "tautau": lambda: load_tautau(tautau),
               "ee": lambda: load_ee(ee)}
    names = list(loaders) if only is None else list(only)
    unknown = [n for n in names if n not in loaders]
    if unknown:
        raise ValueError(f"unknown channel(s) {unknown}: use {sorted(loaders)}")
    return {n: loaders[n]() for n in names}
