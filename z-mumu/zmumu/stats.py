"""Efficiency intervals and uncertainty propagation.

Binomial efficiencies are quoted with Clopper-Pearson ("exact") intervals
rather than the naive sqrt(eps(1-eps)/N), which fails at eps -> 0 or 1 and
would give a zero error on a perfectly efficient bin.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import beta


def clopper_pearson(passed, total, level: float = 0.6827):
    """Central Clopper-Pearson interval. Returns (eff, low_err, high_err).

    Works elementwise on arrays. Bins with `total == 0` return efficiency 0
    with zero errors; callers should mask those rather than plot them.
    """
    passed = np.asarray(passed, dtype=float)
    total = np.asarray(total, dtype=float)
    eff = np.divide(passed, total, out=np.zeros_like(passed), where=total > 0)

    alpha = 1.0 - level
    lo = beta.ppf(alpha / 2, passed, total - passed + 1)
    hi = beta.isf(alpha / 2, passed + 1, total - passed)
    lo = np.nan_to_num(lo, nan=0.0)
    hi = np.where(np.isnan(hi), 1.0, hi)

    lo = np.where(total > 0, lo, 0.0)
    hi = np.where(total > 0, hi, 0.0)
    return eff, eff - lo, hi - eff


def combine_in_quadrature(*terms) -> float:
    """Quadrature sum, for independent relative uncertainties."""
    return float(np.sqrt(sum(np.asarray(t, dtype=float) ** 2 for t in terms)))


def ratio_uncertainty(num, num_err, den, den_err) -> float:
    """Absolute uncertainty on num/den for uncorrelated inputs."""
    if den == 0:
        return float("inf")
    rel = np.sqrt((num_err / num) ** 2 + (den_err / den) ** 2) if num else abs(den_err / den)
    return float(abs(num / den) * rel)


def poisson_err(n) -> np.ndarray:
    """sqrt(N), the Gaussian approximation adequate for the yields here."""
    return np.sqrt(np.asarray(n, dtype=float))


def weighted_mean(values, errors):
    """Inverse-variance weighted mean and its uncertainty."""
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    ok = errors > 0
    if not ok.any():
        return float(np.mean(values)), float("inf")
    w = 1.0 / errors[ok] ** 2
    mean = float(np.sum(w * values[ok]) / np.sum(w))
    return mean, float(np.sqrt(1.0 / np.sum(w)))
