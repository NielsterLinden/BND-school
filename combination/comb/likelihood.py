"""Profile-likelihood form of the same combination -- a cross-check of `blue.combine`.

The BLUE algebra assumes Gaussian, symmetric uncertainties. Writing the combination as an explicit
likelihood with one nuisance parameter per shared uncertainty source

    -2 ln L(sigma, theta) = sum_c (x_c - sigma - sum_s d_{c,s} theta_s)^2 / u_c^2 + sum_s theta_s^2

    u_c^2 = sum over the channel's uncorrelated sources (data statistics, object calibrations, ...)

makes two things possible that BLUE cannot do: the tau tau channel's asymmetric uncertainty
(+16.6 / -14.1 % on mu_Z) can be carried through with a bifurcated response, and the result comes
with a -2 Delta ln L curve to plot. With symmetric errors and rho in {0, 1} this construction is
*algebraically identical* to BLUE, which is what `run_combination.py --check` verifies.

A source with 0 < rho < 1 is split into a shared part sqrt(rho) * s and an independent part
sqrt(1 - rho) * s, which reproduces its covariance exactly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from iminuit import Minuit

from .model import Spec


@dataclass
class LikelihoodResult:
    value: float
    error_up: float
    error_down: float
    nuisances: dict[str, float]
    scan: tuple[np.ndarray, np.ndarray] = field(default_factory=lambda: (np.array([]), np.array([])))
    valid: bool = True

    @property
    def error(self) -> float:
        return 0.5 * (self.error_up + self.error_down)


def _decompose(spec: Spec):
    """-> (shared source names, shared response matrix d[channel][source], per-channel width u_c)."""
    order = spec.order
    shared_names, shared = [], {c: [] for c in order}
    indep = {c: 0.0 for c in order}
    for src in spec.sources:
        rho = float(np.clip(src.rho, 0.0, 1.0))
        if rho > 0:
            shared_names.append(src.name)
            for c in order:
                shared[c].append(np.sqrt(rho) * src.sizes.get(c, 0.0))
        if rho < 1:
            for c in order:
                indep[c] += (1.0 - rho) * src.sizes.get(c, 0.0) ** 2
    d = {c: np.array(shared[c]) for c in order}
    u = {c: float(np.sqrt(v)) for c, v in indep.items()}
    return shared_names, d, u


def _asymmetry(spec: Spec) -> dict[str, tuple[float, float]]:
    """(k_up, k_down) per channel: how much wider the channel's +1 sigma / -1 sigma response is
    than the symmetrised one, taken from its published asymmetric uncertainty on mu_Z."""
    out = {}
    for c in spec.channels:
        sym = 0.5 * (c.mu_err_up + c.mu_err_down)
        out[c.name] = (c.mu_err_up / sym, c.mu_err_down / sym) if sym > 0 else (1.0, 1.0)
    return out


def combine(spec: Spec, *, asymmetric: bool = True, scan_points: int = 201,
            scan_range: float = 5.0) -> LikelihoodResult:
    order = spec.order
    x = np.array([c.sigma for c in spec.channels])
    names, d, u = _decompose(spec)
    k = _asymmetry(spec) if asymmetric else {c: (1.0, 1.0) for c in order}
    n_np = len(names)

    def nll(pars):
        sigma, theta = pars[0], np.asarray(pars[1:])
        total = 0.0
        for i, c in enumerate(order):
            shift = d[c] * theta
            if asymmetric:
                shift = np.where(theta >= 0, shift * k[c][0], shift * k[c][1])
            resid = x[i] - sigma - shift.sum()
            width = u[c] * (k[c][0] if resid < 0 else k[c][1]) if asymmetric else u[c]
            total += (resid / width) ** 2 if width > 0 else 0.0
        return total + float((theta ** 2).sum())

    start = [float(x.mean())] + [0.0] * n_np
    m = Minuit(nll, start, name=["sigma"] + [f"theta_{i}" for i in range(n_np)])
    m.errordef = 1.0
    m.strategy = 2
    m.migrad()
    m.hesse()
    best = float(m.values[0])
    try:
        m.minos("sigma")
        err_up = abs(float(m.merrors["sigma"].upper))
        err_down = abs(float(m.merrors["sigma"].lower))
    except Exception:
        err_up = err_down = float(m.errors[0])

    # profile scan
    half = scan_range * max(err_up, err_down)
    grid = np.linspace(best - half, best + half, scan_points)
    nll_min = m.fval
    curve = []
    for value in grid:
        mp = Minuit(nll, [value] + [0.0] * n_np, name=["sigma"] + [f"theta_{i}" for i in range(n_np)])
        mp.errordef = 1.0
        mp.fixed["sigma"] = True
        if n_np:
            mp.migrad()
        curve.append(mp.fval - nll_min)
    return LikelihoodResult(
        value=best, error_up=err_up, error_down=err_down,
        nuisances={n: float(m.values[i + 1]) for i, n in enumerate(names)},
        scan=(grid, np.array(curve)), valid=bool(m.valid))
