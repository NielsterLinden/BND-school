#!/usr/bin/env python
"""Step 3 -- background estimation, entirely from data.

No simulation is available for this measurement, so both background classes are
estimated from control regions in the same dataset. They are chosen to be
complementary: between them they cover everything that can fake a
Z -> mu mu candidate.

1. Non-prompt ("fake") muons -- same-sign control region.
   Muons from heavy-flavour decay, pion/kaon decay-in-flight and hadronic
   punch-through carry no charge correlation, so they populate same-sign and
   opposite-sign pairs at nearly equal rates. Real Z decays are almost purely
   opposite-sign. The same-sign yield under the identical selection therefore
   measures the fake contamination directly:

       N_fake(OS) = R_OS/SS * N(SS)

   with R_OS/SS = 1.0 +/- 50% (config). See docs/03-fake-leptons.md.

2. Prompt flavour-symmetric backgrounds -- e-mu control region.
   ttbar, tW, WW and Z -> tautau produce e-mu, mu-mu and ee with the
   combinatorial ratio 2 : 1 : 1. These give *opposite-sign prompt* muon pairs
   and are invisible to the same-sign method. Counting e-mu events under an
   otherwise identical selection gives:

       N_fs(mumu) = 0.5 * k * N(emu),   k = eff_mu / eff_e

   See docs/05-backgrounds.md for why k is taken as 1.0 +/- 50% here.

Not covered by either: WZ and ZZ with both muons prompt and same-flavour.
These are ~0.1% of the Z peak in this selection and are carried as a
systematic rather than subtracted.

Outputs
    output/data/step3_backgrounds.pkl   per-bin and integrated background
    output/plots/step3_*.png            control regions and subtracted spectrum

Run:  .venv/bin/python scripts/step3_backgrounds.py   (needs step 1)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from zmumu import config, hists, stats


def main():
    step1 = hists.load("step1_selection.pkl")

    h_os = step1["h_os_fsr"]
    h_ss = step1["h_ss"]
    h_emu = step1["h_emu"]

    # ---- integrated yields in the mass window ---------------------------
    n_os = float(h_os.values().sum())
    n_ss = float(h_ss.values().sum())
    n_emu = float(h_emu.values().sum())

    # ---- fake background from same-sign ---------------------------------
    n_fake = config.R_OS_SS * n_ss
    # Uncertainty: statistical on the SS count, plus the transfer factor.
    fake_stat = config.R_OS_SS * np.sqrt(max(n_ss, 1.0))
    fake_syst = config.R_OS_SS_UNC * n_fake
    fake_err = float(np.hypot(fake_stat, fake_syst))

    # ---- flavour-symmetric background from e-mu -------------------------
    n_fs = config.FS_TRANSFER * config.FS_EFF_RATIO * n_emu
    fs_stat = config.FS_TRANSFER * config.FS_EFF_RATIO * np.sqrt(max(n_emu, 1.0))
    fs_syst = config.FS_REL_UNC * n_fs
    fs_err = float(np.hypot(fs_stat, fs_syst))

    n_bkg = n_fake + n_fs
    bkg_err = float(np.hypot(fake_err, fs_err))

    # ---- per-bin background shape, for the subtracted spectrum ----------
    h_fake = h_ss.copy() * config.R_OS_SS
    h_fs = h_emu.copy() * (config.FS_TRANSFER * config.FS_EFF_RATIO)
    h_bkg = h_fake + h_fs
    h_sub = h_os.copy()
    h_sub[...] = np.maximum(h_os.values() - h_bkg.values(), 0.0)

    result = {
        "n_os": n_os, "n_ss": n_ss, "n_emu": n_emu,
        "n_fake": n_fake, "n_fake_err": fake_err,
        "n_fs": n_fs, "n_fs_err": fs_err,
        "n_bkg": n_bkg, "n_bkg_err": bkg_err,
        "bkg_fraction": n_bkg / n_os if n_os else 0.0,
        "h_fake": h_fake, "h_fs": h_fs, "h_bkg": h_bkg, "h_subtracted": h_sub,
    }

    _report(result)
    _plot(h_os, h_ss, h_emu, h_fake, h_fs, h_sub)
    path = hists.save(result, "step3_backgrounds.pkl")
    print(f"[step3] wrote {path}")


def _report(r):
    print("\n[step3] Background estimate in 60 < m(mumu) < 120 GeV\n")
    print(f"  observed opposite-sign        : {r['n_os']:>14,.0f}")
    print(f"  same-sign control yield       : {r['n_ss']:>14,.0f}")
    print(f"  e-mu control yield            : {r['n_emu']:>14,.0f}")
    print()
    print(f"  non-prompt (fake) background  : {r['n_fake']:>14,.1f} +/- {r['n_fake_err']:,.1f}")
    print(f"  flavour-symmetric background  : {r['n_fs']:>14,.1f} +/- {r['n_fs_err']:,.1f}")
    print(f"  total background              : {r['n_bkg']:>14,.1f} +/- {r['n_bkg_err']:,.1f}")
    print(f"  background fraction           : {100*r['bkg_fraction']:>13.3f}%")


def _plot(h_os, h_ss, h_emu, h_fake, h_fs, h_sub):
    hists.plot_mass({"Same sign (fakes)": h_ss, r"$e\mu$ (flavour symmetric)": h_emu},
                    "step3_control_regions.png",
                    title="Background control regions")
    hists.plot_stack_with_ratio(
        h_os, {"Non-prompt (from SS)": h_fake, r"Flavour symmetric (from $e\mu$)": h_fs},
        "step3_data_vs_background.png",
        title="Opposite-sign data with data-driven background estimates")
    hists.plot_mass({"Observed (OS)": h_os, "After background subtraction": h_sub},
                    "step3_subtracted.png", logy=True,
                    title="Z peak before and after background subtraction")


if __name__ == "__main__":
    main()
