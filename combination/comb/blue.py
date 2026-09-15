"""Best Linear Unbiased Estimate of the combined cross section.

Standard BLUE (Lyons, Gibaut, Clifford, NIM A270 (1988) 110): with measurements x, covariance C
and u = (1, 1, ...)^T,

    w = C^-1 u / (u^T C^-1 u),   sigma_hat = w^T x,   Var(sigma_hat) = 1 / (u^T C^-1 u)

and chi2 = (x - sigma_hat u)^T C^-1 (x - sigma_hat u) with ndf = n - 1 tests whether the channels
are compatible.

Two details that matter here:

* **iteration.** The acceptance uncertainties are multiplicative (delta sigma / sigma = delta A / A),
  so evaluating them at each channel's own measured value biases the weights towards the channel
  that fluctuated low. `combine` re-evaluates them at the combined value and re-solves until stable
  (Lyons et al.; the "iterative BLUE" of the LHC top-mass combinations). Here the shift is tiny --
  the in-fit uncertainties dominate and are additive -- but it is cheap and it removes the question.
* **uncertainty breakdown.** The combined uncertainty is decomposed source by source as
  delta_s = sqrt(w^T C_s w) with C_s the covariance of that source alone. Those add in quadrature to
  the total, so the combined result gets the same kind of table as the channel fits.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from .model import Spec


@dataclass
class BlueResult:
    value: float
    error: float
    weights: dict[str, float]
    breakdown: dict[str, float]          # source -> contribution to the combined uncertainty [pb]
    chi2: float
    ndf: int
    pvalue: float
    covariance: np.ndarray
    correlation: np.ndarray
    inputs: dict[str, tuple[float, float]]
    iterations: int = 1
    meta: dict = field(default_factory=dict)

    @property
    def rel(self) -> float:
        return self.error / self.value

    def group_breakdown(self) -> dict[str, float]:
        """Breakdown collapsed onto the headline groups quoted in result.md."""
        out = {"statistical": 0.0, "luminosity": 0.0, "acceptance": 0.0, "other systematic": 0.0}
        for name, value in self.breakdown.items():
            if name == "Data statistics":
                key = "statistical"
            elif name == "Luminosity":
                key = "luminosity"
            elif name.startswith("Acceptance"):
                key = "acceptance"
            else:
                key = "other systematic"
            out[key] = (out[key] ** 2 + value ** 2) ** 0.5
        return out


def combine(spec: Spec, *, iterate: bool = True, max_iter: int = 20, tol: float = 1e-9) -> BlueResult:
    order = spec.order
    x = spec.values
    u = np.ones(len(x))
    scale = None
    value = float(x.mean())
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        cov = spec.covariance(scale)
        inv = np.linalg.inv(cov)
        denom = float(u @ inv @ u)
        w = (inv @ u) / denom
        new_value = float(w @ x)
        if not iterate:
            value = new_value
            break
        if abs(new_value - value) < tol * max(abs(new_value), 1.0):
            value = new_value
            break
        value = new_value
        # re-evaluate the multiplicative (acceptance) uncertainties at the combined value
        scale = {c.name: value / c.sigma for c in spec.channels}

    cov = spec.covariance(scale)
    inv = np.linalg.inv(cov)
    denom = float(u @ inv @ u)
    w = (inv @ u) / denom
    value = float(w @ x)
    error = float(np.sqrt(1.0 / denom))

    resid = x - value * u
    chi2 = float(resid @ inv @ resid)
    ndf = len(x) - 1
    pvalue = float(stats.chi2.sf(chi2, ndf)) if ndf > 0 else float("nan")

    breakdown = {}
    for src in spec.sources:
        sizes = src.sizes if scale is None or src.kind != "acceptance" else \
            {c: v * scale.get(c, 1.0) for c, v in src.sizes.items()}
        from .model import Source
        contrib = float(np.sqrt(max(w @ Source(src.name, src.kind, src.rho, sizes).cov(order) @ w, 0.0)))
        if contrib > 0:
            breakdown[src.name] = contrib

    d = np.sqrt(np.diag(cov))
    return BlueResult(
        value=value, error=error, weights=dict(zip(order, (float(v) for v in w))),
        breakdown=dict(sorted(breakdown.items(), key=lambda kv: -kv[1])),
        chi2=chi2, ndf=ndf, pvalue=pvalue, covariance=cov, correlation=cov / np.outer(d, d),
        inputs={c.name: (c.sigma, float(np.sqrt(cov[i, i]))) for i, c in enumerate(spec.channels)},
        iterations=n_iter,
        meta={"order": order, "scale": scale or {c: 1.0 for c in order}})


def single(spec: Spec, channel: str) -> BlueResult:
    """The same machinery restricted to one channel -- used as a closure test and for 'mumu alone'."""
    from .model import Spec as _Spec, Source
    keep = [c for c in spec.channels if c.name == channel]
    if not keep:
        raise KeyError(channel)
    sub = _Spec(channels=keep,
                sources=[Source(s.name, s.kind, s.rho, {channel: s.sizes.get(channel, 0.0)})
                         for s in spec.sources if s.sizes.get(channel, 0.0) > 0])
    return combine(sub, iterate=False)
