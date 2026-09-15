"""Load the published per-channel results and turn them into a common structure.

Nothing here recomputes physics: every number is read from a file committed to the repository
and is traceable to the channel that produced it (`provenance` on each `ChannelResult`).

    from comb import inputs
    chans = inputs.load_channels()                     # {"mumu": ..., "tautau": ...}
    chans = inputs.load_channels(mumu="shapefit", tautau="mcsub")

Two things are deliberately *not* symmetric between the channels:

* the z-mumu review of 15 Sep 2026 (`z-mumu/REVIEW.md`, finding F3) showed the 30-bin shape fit is
  not robust and recommends the **counting extraction** with an extra +-0.7 % lineshape term until
  the `SigModel` template is fixed. That is the default here (`mumu="counting"`); the shape fit is
  kept as the `"shapefit"` variant.
* the two channels' `mu_Z` do not share a reference prediction (1954.1 pb for mumu against
  1944.9 pb for tautau, a 0.47 % difference -- see docs/01-inputs.md), so the combination is done
  on the cross sections, never on `mu_Z`.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# --- files read (all committed) -------------------------------------------------------------
ZMUMU_FIT_JSON = REPO / "z-mumu/fit/results/zmumu_fit_result.json"
ZMUMU_RESULTS_V2 = REPO / "z-mumu/output/v2/results_v2.json"
ZMUMU_COUNTING_IMPACTS = REPO / "z-mumu/review/fitcheck/GroupedImpact_rebin30.txt"
ZMUMU_COUNTING_FIT = REPO / "z-mumu/review/fitcheck/zmumu_rebin30.txt"
ZTAUTAU_RESULTS = REPO / "z-tautau/output/results.json"

# z-mumu REVIEW.md section 4, recommendation 3: half the spread of the fit-configuration table
# (mu_Z = 0.990 / 0.996 / 0.9935 / 1.005 / 1.006) assigned as a lineshape-model uncertainty on the
# counting extraction. Channel-specific (it is the mumu mass-template modelling), so uncorrelated.
ZMUMU_LINESHAPE_REL = 0.007

#: nuisance parameters that make up the "Signal modelling" category in both channels
#: (fitting/CONVENTIONS.md). `SigModel` is split off because the two channels compare different
#: generator pairs -- see comb/model.py.
SIGNAL_MODELLING_NPS = ["PDF", "AlphaS", "QCDScale", "PS_ISR", "PS_FSR", "SigModel"]


@dataclass
class ChannelResult:
    """One channel's measurement of sigma(Z/gamma* -> ll, 60 < m < 120 GeV), in pb."""

    name: str                       # "mumu" / "tautau"
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
        (z-mumu: sigma_fid_lumi_pb = 0.0116431 x 799.566).
        """
        return self.groups.get(category, 0.0) * self.sigma_pred

    def acc_pb(self, component: str) -> float:
        """Acceptance component as an absolute uncertainty on sigma [pb].

        A enters as sigma = sigma_fid / A, so its uncertainty is relative to the *measured* value.
        """
        return self.acc.get(component, 0.0) * self.sigma

    @property
    def sigma_err_total(self) -> float:
        terms = [self.sigma_stat]
        terms += [self.group_pb(c) for c in self.groups]
        terms += [self.acc_pb(k) for k in self.acc]
        return sum(t * t for t in terms) ** 0.5


# ---------------------------------------------------------------------------- parsers
def _parse_grouped_impact(path: Path) -> dict[str, float]:
    """TRExFitter `GroupedImpact_<poi>.txt`: '<Category>   <unc>  ( +<up>, -<down> )'."""
    out = {}
    rx = re.compile(r"^(\S.*?)\s{2,}(\S+)\s+\(")
    for line in path.read_text().splitlines():
        m = rx.match(line.strip())
        if not m:
            continue
        try:
            value = float(m.group(2))
        except ValueError:          # TRExFitter writes '-nan' for a group it could not evaluate
            continue
        if value == value:          # drop NaN
            out[m.group(1)] = abs(value)
    return out


def _parse_fit_txt_poi(path: Path, poi: str = "mu_Z") -> tuple[float, float, float]:
    for line in path.read_text().splitlines():
        m = re.match(rf"^{poi}\s+(\S+)\s+\+(\S+)\s+-(\S+)", line.strip())
        if m:
            return float(m.group(1)), float(m.group(2)), float(m.group(3))
    raise RuntimeError(f"{poi} not found in {path}")


def _drop_totals(groups: dict[str, float]) -> dict[str, float]:
    """Remove the summary rows and empty categories from a grouped-impact table."""
    return {k: v for k, v in groups.items() if k not in ("FullSyst", "Total") and v > 0.0}


def _ranking_mumu(js: dict) -> dict[str, float]:
    return {r["name"]: 0.5 * (abs(r["dpoi_up_post"]) + abs(r["dpoi_down_post"])) for r in js.get("ranking", [])}


def _ranking_tautau(fit: dict) -> dict[str, float]:
    return {r["name"]: 0.5 * (abs(r["impact_up"]) + abs(r["impact_down"])) for r in fit.get("ranking", [])}


# ---------------------------------------------------------------------------- channels
def load_mumu(variant: str = "counting") -> ChannelResult:
    """Z -> mu mu. `variant` is "counting" (reviewed recommendation) or "shapefit" (30-bin fit)."""
    js = json.loads(ZMUMU_FIT_JSON.read_text())
    meta = js["meta"]
    sigma_pred = meta["sigma_fid_pred_pb"] / meta["A_60_120"]     # 1954.1 pb
    a = meta["acceptance"]
    acc = {"pdf": a["A_pdf_rel"], "alphas": a["A_alphas_rel"], "scale": a["A_scale_rel"],
           "mcstat": a["A_stat"] / a["A"]}
    prov = [str(ZMUMU_FIT_JSON.relative_to(REPO))]

    if variant == "shapefit":
        groups = _drop_totals(js["grouped_impacts_mu"])
        mu, up, down = js["mu"], js["mu_err_up"], js["mu_err_down"]
        stat = js["mu_stat_only_fit"]
        gof = js.get("gof", {}).get("gof_probability")
    elif variant == "counting":
        groups = _drop_totals(_parse_grouped_impact(ZMUMU_COUNTING_IMPACTS))
        mu, _, _ = _parse_fit_txt_poi(ZMUMU_COUNTING_FIT)
        # the 1-bin Hesse error is degenerate (REVIEW.md section 4); the data statistical
        # uncertainty is sqrt(N_obs) / (N_obs - N_bkg), the same 0.031 % as the shape fit.
        c = meta["counting"]
        stat = mu * c["n_obs"] ** 0.5 / (c["n_obs"] - c["n_bkg"])
        groups["Lineshape model"] = ZMUMU_LINESHAPE_REL * mu
        up = down = sum(v * v for v in groups.values()) ** 0.5
        gof = None
        prov += [str(ZMUMU_COUNTING_IMPACTS.relative_to(REPO)),
                 str(ZMUMU_COUNTING_FIT.relative_to(REPO)), "z-mumu/REVIEW.md (F3, recommendation 3)"]
    else:
        raise ValueError(f"unknown z-mumu variant {variant!r}: use 'counting' or 'shapefit'")

    ranking = _ranking_mumu(js)
    # `SigModel` alone, to be decorrelated from the rest of the Signal modelling category. In the
    # shape fit it is a shape effect (take it from the ranking); in the counting extraction only
    # the powheg/aMC@NLO normalisation difference survives.
    sigmodel = ranking.get("SigModel", 0.0) if variant == "shapefit" else \
        abs(meta["sigmodel_powheg_over_nlo"] - 1.0) * mu

    return ChannelResult(
        name="mumu", label=r"$Z\to\mu\mu$", variant=variant,
        mu=mu, mu_err_up=up, mu_err_down=down, mu_stat=stat,
        sigma_pred=sigma_pred, groups=groups, ranking=ranking, acc=acc, sigmodel=sigmodel, gof_p=gof,
        provenance=prov,
        extra={"sigma_fid_pb": mu * meta["sigma_fid_pred_pb"],
               "sigma_fid_pred_pb": meta["sigma_fid_pred_pb"],
               "A": meta["A_60_120"], "C": meta["C_factor"], "lumi_pb": meta["lumi_pb"],
               "n_obs": meta["counting"]["n_obs"], "n_bkg": meta["counting"]["n_bkg"],
               "shapefit_mu": js["mu"], "counting_sigma_fid_pb": meta["counting"]["sigma_fid_pb"]})


def load_tautau(variant: str = "nominal") -> ChannelResult:
    """Z -> tau_h tau_h. `variant` is the fake-factor variant: "nominal" or "mcsub"."""
    js = json.loads(ZTAUTAU_RESULTS.read_text())
    if variant not in js["fit"]:
        raise ValueError(f"unknown z-tautau variant {variant!r}: use {sorted(js['fit'])}")
    fit = js["fit"][variant]
    pred = js["prediction"]
    acc = {"pdf": pred["A_unc"]["A_pdf"], "alphas": pred["A_unc"]["A_alphas"],
           "scale": pred["A_unc"]["A_scale"], "isr": pred["A_unc"]["A_isr"],
           "fsr": pred["A_unc"]["A_fsr"], "mcstat": pred["A_mc_stat"]}
    return ChannelResult(
        name="tautau", label=r"$Z\to\tau_h\tau_h$", variant=variant,
        mu=fit["mu"], mu_err_up=fit["mu_err_up"], mu_err_down=fit["mu_err_down"],
        mu_stat=fit["mu_stat"],
        sigma_pred=pred["sigma_tautau_60_120_pb"], groups=_drop_totals(fit["grouped_impact"]),
        ranking=_ranking_tautau(fit), acc=acc,
        sigmodel=_ranking_tautau(fit).get("SigModel", 0.0), gof_p=fit.get("gof_probability"),
        provenance=[str(ZTAUTAU_RESULTS.relative_to(REPO))],
        extra={"sigma_fid_pb": fit["sigma_fid_pb"]["value"],
               "sigma_fid_pred_pb": pred["sigma_fid_pb"], "A": pred["A"],
               "C": js["for_combination"]["C"], "lumi_pb": js["lumi_pb"],
               "n_obs": js["for_combination"]["n_obs"],
               "n_bkg": js["for_combination"]["n_bkg_prefit"]})


def load_channels(mumu: str = "counting", tautau: str = "nominal") -> dict[str, ChannelResult]:
    return {"mumu": load_mumu(mumu), "tautau": load_tautau(tautau)}
