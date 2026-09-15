"""Lepton universality: R = sigma(Z -> tau tau) / sigma(Z -> mu mu).

The ratio is the part of the two measurements that the combination itself cannot test, because the
combination *assumes* universality by fitting one cross section to both channels. It is also the
part where the dominant uncertainty of the mu mu channel disappears: luminosity, pileup, L1
prefiring and the correlated theory terms cancel in the ratio, exactly and by construction.

For a ratio of two measurements with relative uncertainties r_mu, r_tau per source and correlation
rho_s,

    (delta R / R)^2 = sum_s [ r_tau,s^2 + r_mu,s^2 - 2 rho_s r_tau,s r_mu,s ]
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from .model import Spec


@dataclass
class RatioResult:
    value: float
    error: float
    breakdown: dict[str, float] = field(default_factory=dict)   # source -> contribution to dR/R
    numerator: str = "tautau"
    denominator: str = "mumu"

    @property
    def rel(self) -> float:
        return self.error / self.value

    @property
    def z_from_unity(self) -> float:
        return (self.value - 1.0) / self.error

    @property
    def pvalue(self) -> float:
        return float(2 * stats.norm.sf(abs(self.z_from_unity)))


def compute(spec: Spec, numerator: str = "tautau", denominator: str = "mumu") -> RatioResult:
    by_name = {c.name: c for c in spec.channels}
    num, den = by_name[numerator], by_name[denominator]
    value = num.sigma / den.sigma
    var, breakdown = 0.0, {}
    for src in spec.sources:
        r_n = src.sizes.get(numerator, 0.0) / num.sigma
        r_d = src.sizes.get(denominator, 0.0) / den.sigma
        contrib = r_n ** 2 + r_d ** 2 - 2.0 * src.rho * r_n * r_d
        contrib = max(contrib, 0.0)
        if contrib > 0:
            breakdown[src.name] = float(np.sqrt(contrib) * value)
        var += contrib
    return RatioResult(value=value, error=float(np.sqrt(var) * value),
                       breakdown=dict(sorted(breakdown.items(), key=lambda kv: -kv[1])),
                       numerator=numerator, denominator=denominator)
