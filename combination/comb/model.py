"""The correlation model: which uncertainty is shared between the two channels, and how much.

This is the only place in the combination where a *judgement* is made rather than a number read
from a file, so it is deliberately small, explicit and auditable. Every entry of `CORRELATION` is
justified in `docs/02-correlation-model.md`.

The keys are the `Category` strings that `fitting/CONVENTIONS.md` forces every channel to use, so
a category that appears in a channel's grouped-impact table but not in the map is an error, not a
silent zero: a channel that adds a systematic must state how it correlates.

    from comb import inputs, model
    spec = model.build(inputs.load_channels())
    spec.covariance()            # 2x2 numpy array, pb^2
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .inputs import ChannelResult

#: rho between the two channels for each in-fit systematic category.
#: "split" marks a category that is broken up before the correlation is applied.
CORRELATION: dict[str, float | str] = {
    # --- shared inputs, derived once and used by both channels: fully correlated
    "Luminosity": 1.0,                  # same normtag value, same runs, same 1.2 % calibration
    "Pileup": 1.0,                      # same data profile, same 69.2 mb, same reweighting code
    "L1 prefiring": 1.0,                # same maps, same method
    "Background normalisation": 1.0,    # identical XS_* nuisance parameters on shared MC samples
    # --- theory: correlated except the generator comparison, which is a different pair per channel
    "Signal modelling": "split",
    # --- object calibrations: different objects, independently derived -> uncorrelated
    "Muon efficiency": 0.0,
    "Muon momentum": 0.0,
    "Electron efficiency": 0.0,
    "Tau ID": 0.0,
    "Tau trigger": 0.0,
    "Tau energy scale": 0.0,
    # --- method- and channel-specific
    "Gammas": 0.0,                      # per-bin MC statistics of disjoint selections
    "Fakes": 0.0,                       # fake factor vs fake factor, different regions and objects
    "MET": 0.0,                         # tau tau only
    "Lineshape model": 0.0,             # z-mumu REVIEW.md F3/F4: a mumu mass-template issue
}

#: rho for the acceptance components, which live outside both fits.
ACC_CORRELATION: dict[str, float] = {
    "pdf": 1.0,       # same NNPDF3.1 set, same Hessian members
    "alphas": 1.0,    # same alpha_s members
    "scale": 1.0,     # same 7-point envelope of the same sample (decorrelate with correlate_acc_scale=False)
    "isr": 1.0,       # same parton-shower weights
    "fsr": 1.0,
    "mcstat": 0.0,    # disjoint generated events
}

ACC_LABEL = {"pdf": "Acceptance: PDF", "alphas": r"Acceptance: $\alpha_s$", "scale": "Acceptance: QCD scale",
             "isr": "Acceptance: ISR", "fsr": "Acceptance: FSR", "mcstat": "Acceptance: MC stat."}


@dataclass
class Source:
    """One uncertainty source: its size in each channel and the correlation between channels."""

    name: str
    kind: str                      # "in-fit" | "acceptance" | "statistical"
    rho: float
    sizes: dict[str, float]        # channel name -> absolute uncertainty on sigma [pb]

    def cov(self, order: list[str]) -> np.ndarray:
        s = np.array([self.sizes.get(c, 0.0) for c in order])
        m = self.rho * np.outer(s, s)
        np.fill_diagonal(m, s * s)
        return m


@dataclass
class Spec:
    """Everything the combination needs: the measurements and the decomposed covariance."""

    channels: list[ChannelResult]
    sources: list[Source]

    @property
    def order(self) -> list[str]:
        return [c.name for c in self.channels]

    @property
    def values(self) -> np.ndarray:
        return np.array([c.sigma for c in self.channels])

    def covariance(self, scale: dict[str, float] | None = None) -> np.ndarray:
        """Total covariance. `scale` rescales a channel's multiplicative (acceptance) sources,
        which is what the iterative BLUE of `blue.combine` uses to remove the normalisation bias."""
        order = self.order
        total = np.zeros((len(order), len(order)))
        for src in self.sources:
            if scale and src.kind == "acceptance":
                src = Source(src.name, src.kind, src.rho,
                             {c: v * scale.get(c, 1.0) for c, v in src.sizes.items()})
            total += src.cov(order)
        return total

    def by_name(self, name: str) -> Source:
        return next(s for s in self.sources if s.name == name)


def _split_signal_modelling(ch: ChannelResult) -> tuple[float, float]:
    """(correlated theory part, uncorrelated generator part) of the Signal modelling category, in pb.

    The grouped impact is the quantity to preserve, so the category is split by the fraction that
    `SigModel` carries among the nuisance parameters of the category:
    f = impact(SigModel) / impact(category), capped at 1.
    """
    total = ch.group_pb("Signal modelling")
    if total <= 0:
        return 0.0, 0.0
    gen = min(ch.sigmodel * ch.sigma_pred, total)
    theory = (total * total - gen * gen) ** 0.5
    return theory, gen


def build(channels: dict[str, ChannelResult] | list[ChannelResult], *,
          correlate_sigmodel: bool = False, correlate_acc_scale: bool = True,
          rho_override: float | None = None) -> Spec:
    """Turn the channel results into a list of correlated uncertainty sources.

    correlate_sigmodel   treat the generator comparison as correlated too (a variation: the two
                         channels compare different generator pairs, so the default is False)
    correlate_acc_scale  correlate the QCD-scale part of the acceptance (default True, conservative)
    rho_override         force every *systematic* rho to this value -- used for the rho = 0 and
                         rho = 1 variations. The data-statistics source is never overridden: the
                         two channels read disjoint primary datasets, so its rho is 0 by fact.
    """
    chans = list(channels.values()) if isinstance(channels, dict) else list(channels)
    names = [c.name for c in chans]
    sources: list[Source] = []

    def rho_of(default: float) -> float:
        return default if rho_override is None else rho_override

    # --- statistical: disjoint datasets (SingleMuon and Tau primary datasets, orthogonal triggers)
    sources.append(Source("Data statistics", "statistical", 0.0,
                          {c.name: c.sigma_stat for c in chans}))

    # --- in-fit systematics, category by category
    categories = sorted({k for c in chans for k in c.groups})
    unknown = [k for k in categories if k not in CORRELATION]
    if unknown:
        raise KeyError(f"no correlation assigned to {unknown}; add them to comb/model.CORRELATION "
                       "(see docs/02-correlation-model.md) -- a new systematic must not default to zero")
    for cat in categories:
        rule = CORRELATION[cat]
        if rule == "split":
            theory = {c.name: _split_signal_modelling(c)[0] for c in chans}
            gen = {c.name: _split_signal_modelling(c)[1] for c in chans}
            sources.append(Source("Signal modelling (PDF, scales, PS)", "in-fit", rho_of(1.0), theory))
            sources.append(Source("Signal modelling (generator)", "in-fit",
                                  rho_of(1.0 if correlate_sigmodel else 0.0), gen))
        else:
            sizes = {c.name: c.group_pb(cat) for c in chans}
            if any(v > 0 for v in sizes.values()):
                sources.append(Source(cat, "in-fit", rho_of(float(rule)), sizes))

    # --- acceptance, outside both fits
    for comp in sorted({k for c in chans for k in c.acc}):
        rho = ACC_CORRELATION[comp]
        if comp == "scale" and not correlate_acc_scale:
            rho = 0.0
        sizes = {c.name: c.acc_pb(comp) for c in chans}
        if any(v > 0 for v in sizes.values()):
            sources.append(Source(ACC_LABEL.get(comp, f"Acceptance: {comp}"), "acceptance",
                                  rho_of(rho), sizes))

    assert set(names) == {c.name for c in chans}
    return Spec(channels=chans, sources=sources)
