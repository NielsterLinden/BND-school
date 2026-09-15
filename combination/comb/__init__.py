"""Combination of the BND-school Z cross-section channels.

The package turns the published per-channel TRExFitter results into one measurement of
sigma(pp -> Z/gamma* -> ll, 60 < m < 120 GeV) under lepton universality.

    from comb import inputs, model, blue
    chans = inputs.load_channels()                  # z-mumu + z-tautau
    spec  = model.build(chans)                      # sources + correlations -> covariance
    res   = blue.combine(spec)

Everything is driven by the `Category` strings of `fitting/CONVENTIONS.md`; see
`docs/02-correlation-model.md` for why each category gets the correlation it gets.
"""

from . import inputs, model, blue, likelihood, ratio  # noqa: F401

__all__ = ["inputs", "model", "blue", "likelihood", "ratio"]
