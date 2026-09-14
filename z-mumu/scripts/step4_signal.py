#!/usr/bin/env python
"""Step 4 -- signal extraction.

Two independent estimates of the number of Z -> mu mu decays in the mass
window, so that the answer does not rest on one method:

**Counting.** Observed opposite-sign events minus the data-driven backgrounds
from step 3. This is the nominal result. It makes no assumption about the
shape of the peak, which matters because the Z lineshape is distorted by FSR
and by detector resolution in ways that are hard to model without simulation.

**Fitting.** A Voigt profile (relativistic Breit-Wigner convolved with a
Gaussian resolution) plus an exponential background, fitted over the full
window. The Breit-Wigner width is fixed to the PDG Z width; the Gaussian width
floats and comes out as the effective dimuon mass resolution. This is a
*cross-check*: the spread between the two methods is taken as the signal
extraction systematic in step 5.

The fit sits *below* the count by a few percent, and its chi2/ndf is poor.
Both are expected and neither is a bug: a single Voigt profile has no
radiative tail, so the exponential background bends down to absorb the FSR
shoulder on the low-mass side, stealing events from the signal integral. The
fit is therefore used as a *validation* of the peak position, the resolution,
and the absence of a large unaccounted smooth background -- not as an
alternative yield, and not as the signal-extraction systematic. Quoting
|fit - count| as a systematic would compare two different quantities.

The genuine signal-extraction ambiguity is the FSR treatment (which events
land inside the mass window at all), and that is evaluated in step 5 by
recomputing the yield with bare muons. Replacing the Voigt profile with a
Crystal Ball convolved with a Breit-Wigner would make the fit describe the
tail properly; that is listed as a next step in handoff.md.

Outputs
    output/data/step4_signal.pkl   both yields and the fit parameters
    output/plots/step4_fit.png     fitted peak with residuals

Run:  .venv/bin/python scripts/step4_signal.py   (needs steps 1 and 3)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import voigt_profile

from zmumu import config, hists

Z_MASS_PDG = 91.1876    # GeV
Z_WIDTH_PDG = 2.4952    # GeV


def model(m, norm, peak, sigma, bkg_norm, bkg_slope):
    """Voigt signal + exponential background, both in events per bin."""
    gamma = Z_WIDTH_PDG / 2.0
    signal = norm * voigt_profile(m - peak, sigma, gamma)
    background = bkg_norm * np.exp(-bkg_slope * (m - config.MASS_LO))
    return signal + background


def main():
    step1 = hists.load("step1_selection.pkl")
    step3 = hists.load("step3_backgrounds.pkl")

    h_os = step1["h_os_fsr"]
    centres = h_os.axes[0].centers
    counts = h_os.values()
    width = centres[1] - centres[0]

    # ---- counting -------------------------------------------------------
    n_obs = float(counts.sum())
    n_bkg = step3["n_bkg"]
    n_count = n_obs - n_bkg
    n_count_err = float(np.hypot(np.sqrt(n_obs), step3["n_bkg_err"]))

    # ---- fit ------------------------------------------------------------
    errors = np.sqrt(np.maximum(counts, 1.0))
    p0 = [counts.sum() * width * 2.5, Z_MASS_PDG, 1.5, counts[0], 0.05]
    bounds = ([0, 88.0, 0.2, 0.0, 0.0], [np.inf, 94.0, 10.0, np.inf, 1.0])
    popt, pcov = curve_fit(model, centres, counts, p0=p0, sigma=errors,
                           bounds=bounds, maxfev=40000)
    perr = np.sqrt(np.diag(pcov))

    norm, peak, sigma, bkg_norm, bkg_slope = popt
    # Integral of the signal component across the window, in events.
    grid = np.linspace(config.MASS_LO, config.MASS_HI, 6000)
    sig_density = norm * voigt_profile(grid - peak, sigma, Z_WIDTH_PDG / 2.0)
    n_fit = float(np.trapezoid(sig_density, grid) / width)
    n_fit_err = float(n_fit * perr[0] / norm) if norm > 0 else float("inf")

    residual = (counts - model(centres, *popt)) / errors
    chi2 = float(np.sum(residual**2))
    ndf = len(centres) - len(popt)

    # Integral of the *fitted* background across the window. Compared against
    # the data-driven estimate as a cross-check for a missed smooth component.
    bkg_density = bkg_norm * np.exp(-bkg_slope * (grid - config.MASS_LO))
    n_bkg_fit = float(np.trapezoid(bkg_density, grid) / width)

    result = {
        "n_observed": n_obs,
        "n_background": n_bkg,
        "n_signal_counting": n_count,
        "n_signal_counting_err": n_count_err,
        "n_signal_fit": n_fit,
        "n_signal_fit_err": n_fit_err,
        "fit_peak": float(peak), "fit_peak_err": float(perr[1]),
        "fit_sigma": float(sigma), "fit_sigma_err": float(perr[2]),
        "chi2": chi2, "ndf": ndf,
        "n_background_fit": n_bkg_fit,
        # Kept for the record; NOT used as a systematic (see module docstring).
        "fit_minus_count_rel": abs(n_fit - n_count) / n_count if n_count else 0.0,
    }
    _report(result)
    _plot(centres, counts, errors, popt, residual, result)
    path = hists.save(result, "step4_signal.pkl")
    print(f"[step4] wrote {path}")


def _report(r):
    print("\n[step4] Signal extraction in 60 < m(mumu) < 120 GeV\n")
    print(f"  observed                       : {r['n_observed']:>14,.0f}")
    print(f"  background (step 3)            : {r['n_background']:>14,.1f}")
    print(f"  signal, counting  (NOMINAL)    : {r['n_signal_counting']:>14,.0f}"
          f" +/- {r['n_signal_counting_err']:,.0f}")
    print(f"  signal, Voigt+exp fit (check)  : {r['n_signal_fit']:>14,.0f}"
          f" +/- {r['n_signal_fit_err']:,.0f}")
    print(f"  fit - count                    : {100*r['fit_minus_count_rel']:>13.2f}%"
          "   (expected: Voigt has no FSR tail)")
    print()
    print(f"  fitted peak                    : {r['fit_peak']:.3f} +/- {r['fit_peak_err']:.3f} GeV"
          f"   (PDG {Z_MASS_PDG})")
    print(f"  fitted resolution (Gaussian)   : {r['fit_sigma']:.3f} +/- {r['fit_sigma_err']:.3f} GeV")
    print(f"  chi2/ndf                       : {r['chi2']:.0f}/{r['ndf']} = {r['chi2']/r['ndf']:.2f}"
          "   (poor: single Voigt lacks the radiative tail)")
    print()
    print(f"  background, data-driven        : {r['n_background']:>14,.1f}")
    print(f"  background, from fit           : {r['n_background_fit']:>14,.1f}"
          "   (absorbs the FSR shoulder; upper bound only)")


def _plot(centres, counts, errors, popt, residual, r):
    import matplotlib.pyplot as plt
    import mplhep as hep

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(9, 9), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07})
    fine = np.linspace(config.MASS_LO, config.MASS_HI, 1500)
    ax.errorbar(centres, counts, yerr=errors, fmt="o", color="black",
                markersize=3, label="Data (OS, FSR recovered)")
    ax.plot(fine, model(fine, *popt), color="crimson", linewidth=2,
            label="Voigt + exponential fit")
    ax.plot(fine, popt[3] * np.exp(-popt[4] * (fine - config.MASS_LO)),
            color="royalblue", linestyle="--", linewidth=1.6, label="Background component")
    ax.set_ylabel("Events / 0.5 GeV")
    ax.set_yscale("log")
    ax.legend(fontsize=12)
    ax.text(0.04, 0.55,
            f"peak = {r['fit_peak']:.2f} GeV\n"
            f"$\\sigma$ = {r['fit_sigma']:.2f} GeV\n"
            f"$\\chi^2$/ndf = {r['chi2']/r['ndf']:.2f}",
            transform=ax.transAxes, fontsize=13, va="top")
    hep.cms.label("Open Data", data=True, lumi=round(config.LUMI_PB/1000, 1), year=2016, ax=ax)

    rax.errorbar(centres, residual, yerr=1.0, fmt="o", color="black", markersize=3)
    rax.axhline(0.0, linestyle="--", color="grey")
    rax.set_ylim(-6, 6)
    rax.set_ylabel("pull")
    rax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    hists.save_fig(fig, "step4_fit.png")


if __name__ == "__main__":
    main()
