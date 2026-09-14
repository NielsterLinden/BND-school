#!/usr/bin/env python
"""Step 1 -- event selection.

Reads the full skim once and fills, in a single pass, the three regions the
rest of the measurement needs:

    OS  two opposite-sign selected muons ......... signal region
    SS  two same-sign selected muons ............. fake-muon control region
    EMU one selected muon + one selected electron  flavour-symmetric control

Doing all three in one pass matters: the skim is 31.5 GB, and each extra pass
costs several minutes. The regions share every cut except the pairing, so they
are by construction consistent with one another.

Outputs
    output/data/step1_selection.pkl   histograms + cutflow + raw yields
    output/plots/step1_*.png          mass spectra, cutflow, kinematics

Run:  .venv/bin/python scripts/step1_selection.py
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import OrderedDict
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import awkward as ak
import numpy as np

from zmumu import config, hists, io, objects

# Branches needed for this step. Kept explicit so the I/O cost is visible.
BRANCHES = [
    "run", "luminosityBlock",
    "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
    "Muon_mediumId", "Muon_pfRelIso04_all", "Muon_dxy", "Muon_dz",
    "Muon_isTracker", "Muon_isGlobal",
    "nFsrPhoton", "FsrPhoton_pt", "FsrPhoton_eta", "FsrPhoton_phi",
    "FsrPhoton_relIso03", "FsrPhoton_dROverEt2", "FsrPhoton_muonIdx",
    "nElectron", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass",
    "Electron_charge", "Electron_cutBased", "Electron_pfRelIso03_all",
    "PV_npvsGood",
] + config.TRIGGERS + config.MET_FILTERS

CUTS = [
    "skim (>=2 muons in collection)",
    "golden JSON (certified lumi)",
    "HLT_IsoMu24 or HLT_IsoTkMu24",
    "MET filters + good PV",
    "exactly 2 selected muons",
    "leading muon pT > 26 GeV",
    "opposite sign",
    "60 < m(mumu) < 120 GeV",
]


def blank_results():
    """Histograms and counters produced by one worker."""
    return {
        "cutflow": OrderedDict((c, 0) for c in CUTS),
        "h_os": hists.mass_hist(),
        "h_os_fsr": hists.mass_hist(),
        "h_ss": hists.mass_hist(),
        "h_emu": hists.mass_hist(),
        "h_os_wide": hists.regular(190, 10.0, 200.0, r"$m_{\mu\mu}$ [GeV]"),
        "h_pt_lead": hists.regular(60, 20.0, 200.0, r"leading muon $p_{T}$ [GeV]"),
        "h_pt_sub": hists.regular(60, 20.0, 200.0, r"subleading muon $p_{T}$ [GeV]"),
        "h_eta_lead": hists.regular(48, -2.4, 2.4, r"leading muon $\eta$"),
        "h_eta_sub": hists.regular(48, -2.4, 2.4, r"subleading muon $\eta$"),
        "h_zpt": hists.regular(60, 0.0, 300.0, r"$p_{T}^{\mu\mu}$ [GeV]"),
        "h_zy": hists.regular(48, -2.4, 2.4, r"$y^{\mu\mu}$"),
        "h_npv": hists.regular(60, 0.0, 60.0, "good primary vertices"),
        # (pT, eta) maps of the two signal muons, on the efficiency binning.
        # Step 5 folds these with the tag-and-probe maps to get <eff>.
        "h_mu1_kin": hists.hist2d(config.EFF_PT_BINS, config.EFF_ETA_BINS,
                                  r"leading muon $p_T$ [GeV]", r"leading muon $\eta$"),
        "h_mu2_kin": hists.hist2d(config.EFF_PT_BINS, config.EFF_ETA_BINS,
                                  r"subleading muon $p_T$ [GeV]", r"subleading muon $\eta$"),
        "h_nfsr": hists.regular(5, 0.0, 5.0, "recovered FSR photons per event"),
        "h_dm_fsr": hists.regular(80, 0.0, 20.0, r"$m_{\mu\mu\gamma} - m_{\mu\mu}$ [GeV]"),
        "yields": {"os": 0, "ss": 0, "emu": 0, "os_window": 0, "ss_window": 0, "emu_window": 0},
        "missing_branches": set(),
        "n_files": 0,
    }


def merge(into, other):
    for cut, n in other["cutflow"].items():
        into["cutflow"][cut] += n
    for key, val in other.items():
        if key.startswith("h_"):
            into[key] += val
    for key, val in other["yields"].items():
        into["yields"][key] += val
    into["missing_branches"] |= other["missing_branches"]
    into["n_files"] += other["n_files"]
    return into


def process_file(path: Path) -> dict:
    """Select events in one skim file. Returns histograms and counters."""
    grl = io.load_grl()
    out = blank_results()
    out["n_files"] = 1

    for events, missing in io.read_chunks(path, BRANCHES):
        out["missing_branches"] |= set(missing)
        out["cutflow"][CUTS[0]] += len(events)

        # --- event-level cleaning -----------------------------------------
        certified = io.lumi_mask(events.run, events.luminosityBlock, grl)
        events = events[certified]
        out["cutflow"][CUTS[1]] += len(events)
        if len(events) == 0:
            continue

        events = events[io.trigger_or(events, config.TRIGGERS)]
        out["cutflow"][CUTS[2]] += len(events)
        if len(events) == 0:
            continue

        events = events[io.pass_met_filters(events) & (events.PV_npvsGood >= 1)]
        out["cutflow"][CUTS[3]] += len(events)
        if len(events) == 0:
            continue

        # --- muon selection -----------------------------------------------
        good = objects.good_muon_mask(events)
        n_good = ak.sum(good, axis=1)

        # ---- flavour-symmetric control region (1 muon + 1 electron) ------
        # Filled before the dimuon requirement, since it needs exactly one muon.
        _fill_emu(events, good, n_good, out)

        events = events[n_good == config.N_MUONS_REQUIRED]
        good = good[n_good == config.N_MUONS_REQUIRED]
        out["cutflow"][CUTS[4]] += len(events)
        if len(events) == 0:
            continue

        idx1, idx2 = objects.leading_two(events, good)

        on_plateau = objects.take(events.Muon_pt, idx1) > config.MU_PT_LEAD
        events, idx1, idx2 = events[on_plateau], idx1[on_plateau], idx2[on_plateau]
        out["cutflow"][CUTS[5]] += len(events)
        if len(events) == 0:
            continue

        # --- split by charge ----------------------------------------------
        opposite = (objects.take(events.Muon_charge, idx1)
                    * objects.take(events.Muon_charge, idx2)) < 0
        _fill_same_sign(events[~opposite], idx1[~opposite], idx2[~opposite], out)

        events, idx1, idx2 = events[opposite], idx1[opposite], idx2[opposite]
        out["cutflow"][CUTS[6]] += len(events)
        if len(events) == 0:
            continue

        _fill_signal(events, idx1, idx2, out)

    return out


def _fill_signal(events, idx1, idx2, out):
    """Fill the opposite-sign signal region."""
    mass_bare = objects.dimuon_mass(events, idx1, idx2, with_fsr=False)
    mass_fsr = objects.dimuon_mass(events, idx1, idx2, with_fsr=config.FSR_ENABLED)

    out["h_os_wide"].fill(ak.to_numpy(mass_bare))
    out["h_os"].fill(ak.to_numpy(mass_bare))
    out["h_os_fsr"].fill(ak.to_numpy(mass_fsr))
    out["yields"]["os"] += len(events)

    # How often FSR recovery fires, and how much mass it puts back.
    attached = objects.fsr_photon_mask(events) & (
        (events.FsrPhoton_muonIdx == idx1) | (events.FsrPhoton_muonIdx == idx2)
    )
    n_fsr = ak.to_numpy(ak.sum(attached, axis=1))
    out["h_nfsr"].fill(n_fsr)
    shifted = n_fsr > 0
    if shifted.any():
        dm = ak.to_numpy(mass_fsr)[shifted] - ak.to_numpy(mass_bare)[shifted]
        out["h_dm_fsr"].fill(dm)

    # Kinematics, inside the mass window only.
    in_window = (mass_fsr > config.MASS_LO) & (mass_fsr < config.MASS_HI)
    out["yields"]["os_window"] += int(ak.sum(in_window))
    out["cutflow"][CUTS[7]] += int(ak.sum(in_window))

    w = events[in_window]
    i1, i2 = idx1[in_window], idx2[in_window]
    if len(w) == 0:
        return
    out["h_pt_lead"].fill(ak.to_numpy(objects.take(w.Muon_pt, i1)))
    out["h_pt_sub"].fill(ak.to_numpy(objects.take(w.Muon_pt, i2)))
    out["h_eta_lead"].fill(ak.to_numpy(objects.take(w.Muon_eta, i1)))
    out["h_eta_sub"].fill(ak.to_numpy(objects.take(w.Muon_eta, i2)))
    zpt, zy = objects.dimuon_pt_y(w, i1, i2)
    out["h_zpt"].fill(ak.to_numpy(zpt))
    out["h_zy"].fill(ak.to_numpy(zy))
    out["h_npv"].fill(ak.to_numpy(w.PV_npvsGood))

    pt1 = np.clip(ak.to_numpy(objects.take(w.Muon_pt, i1)),
                  config.EFF_PT_BINS[0], config.EFF_PT_BINS[-1] - 1e-3)
    pt2 = np.clip(ak.to_numpy(objects.take(w.Muon_pt, i2)),
                  config.EFF_PT_BINS[0], config.EFF_PT_BINS[-1] - 1e-3)
    out["h_mu1_kin"].fill(pt1, ak.to_numpy(objects.take(w.Muon_eta, i1)))
    out["h_mu2_kin"].fill(pt2, ak.to_numpy(objects.take(w.Muon_eta, i2)))


def _fill_same_sign(events, idx1, idx2, out):
    """Fill the same-sign control region used for the fake-muon estimate."""
    if len(events) == 0:
        return
    mass = objects.dimuon_mass(events, idx1, idx2, with_fsr=config.FSR_ENABLED)
    out["h_ss"].fill(ak.to_numpy(mass))
    out["yields"]["ss"] += len(events)
    out["yields"]["ss_window"] += int(
        ak.sum((mass > config.MASS_LO) & (mass < config.MASS_HI))
    )


def _fill_emu(events, good_mu, n_good_mu, out):
    """Fill the e-mu control region used for flavour-symmetric backgrounds.

    Exactly one selected muon above the trigger plateau and exactly one
    selected electron, opposite sign. ttbar, tW, WW and Z->tautau all populate
    this region with the same rate as they populate mu-mu, up to the known
    factor of two from flavour combinatorics (docs/05-backgrounds.md).
    """
    if "Electron_pt" not in ak.fields(events):
        return
    good_el = objects.good_electron_mask(events)
    n_el = ak.sum(good_el, axis=1)

    keep = (n_good_mu == 1) & (n_el == 1)
    if not ak.any(keep):
        return
    ev = events[keep]
    mu_mask = good_mu[keep]
    el_mask = good_el[keep]

    mu_pt = ev.Muon_pt[mu_mask][:, 0]
    mu_eta = ev.Muon_eta[mu_mask][:, 0]
    mu_phi = ev.Muon_phi[mu_mask][:, 0]
    mu_m = ev.Muon_mass[mu_mask][:, 0]
    mu_q = ev.Muon_charge[mu_mask][:, 0]

    el_pt = ev.Electron_pt[el_mask][:, 0]
    el_eta = ev.Electron_eta[el_mask][:, 0]
    el_phi = ev.Electron_phi[el_mask][:, 0]
    el_m = ev.Electron_mass[el_mask][:, 0]
    el_q = ev.Electron_charge[el_mask][:, 0]

    plateau = mu_pt > config.MU_PT_LEAD
    opposite = (mu_q * el_q) < 0
    sel = plateau & opposite
    if not ak.any(sel):
        return

    px1, py1, pz1, e1 = objects.p4(mu_pt[sel], mu_eta[sel], mu_phi[sel], mu_m[sel])
    px2, py2, pz2, e2 = objects.p4(el_pt[sel], el_eta[sel], el_phi[sel], el_m[sel])
    mass = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)

    out["h_emu"].fill(ak.to_numpy(mass))
    out["yields"]["emu"] += int(ak.sum(sel))
    out["yields"]["emu_window"] += int(
        ak.sum((mass > config.MASS_LO) & (mass < config.MASS_HI))
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-files", type=int, default=None,
                    help="process only the first N files (for a quick test)")
    ap.add_argument("--workers", type=int, default=config.N_WORKERS)
    args = ap.parse_args()

    files = io.find_files()
    if args.max_files:
        files = files[: args.max_files]
    print(f"[step1] selecting from {len(files)} files with {args.workers} workers")

    t0 = time.time()
    total = blank_results()
    total["n_files"] = 0
    with Pool(args.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(process_file, files), 1):
            total = merge(total, res)
            if i % 20 == 0 or i == len(files):
                print(f"        {i}/{len(files)} files  ({time.time()-t0:.0f}s)")

    total["missing_branches"] = sorted(total["missing_branches"])
    _report(total)
    _plot(total)
    path = hists.save(total, "step1_selection.pkl")
    print(f"[step1] wrote {path}  ({time.time()-t0:.0f}s total)")


def _report(res):
    print("\n  Cutflow")
    print(f"  {'cut':<38s} {'events':>14s} {'abs %':>9s} {'rel %':>8s}")
    first = list(res["cutflow"].values())[0]
    prev = first
    for cut, n in res["cutflow"].items():
        print(f"  {cut:<38s} {n:>14,} {100*n/first:>8.3f}% {100*n/prev if prev else 0:>7.2f}%")
        prev = n
    y = res["yields"]
    print(f"\n  Region yields (60-120 GeV):")
    print(f"    opposite-sign mu-mu : {y['os_window']:>12,}")
    print(f"    same-sign mu-mu     : {y['ss_window']:>12,}")
    print(f"    e-mu                : {y['emu_window']:>12,}")


def _plot(res):
    hists.plot_mass({"Opposite sign (FSR recovered)": res["h_os_fsr"],
                     "Same sign (fake control)": res["h_ss"],
                     r"$e\mu$ (flavour-symmetric control)": res["h_emu"]},
                    "step1_mass_regions.png",
                    title="Dimuon mass in the three analysis regions")
    hists.plot_mass({"Bare muons": res["h_os"],
                     "With FSR recovery": res["h_os_fsr"]},
                    "step1_mass_fsr.png", logy=False,
                    title="Effect of FSR photon recovery on the Z peak")
    hists.plot_mass({"Opposite sign": res["h_os_wide"]},
                    "step1_mass_wide.png",
                    title="Dimuon mass, 10-200 GeV",
                    ylabel="Events / GeV")
    hists.plot_cutflow(list(res["cutflow"].keys()), list(res["cutflow"].values()),
                       "step1_cutflow.png", title="Z$\\to\\mu\\mu$ selection cutflow")

    import matplotlib.pyplot as plt
    import mplhep as hep
    panels = [("h_pt_lead", "leading muon $p_T$"), ("h_pt_sub", "subleading muon $p_T$"),
              ("h_eta_lead", "leading muon $\\eta$"), ("h_eta_sub", "subleading muon $\\eta$"),
              ("h_zpt", "$p_T^{\\mu\\mu}$"), ("h_zy", "$y^{\\mu\\mu}$"),
              ("h_npv", "good PVs"), ("h_nfsr", "FSR photons recovered")]
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    for ax, (key, label) in zip(axes.flat, panels):
        hep.histplot(res[key], ax=ax, histtype="step", color="black", linewidth=1.4)
        ax.set_xlabel(res[key].axes[0].label)
        ax.set_ylabel("Events")
        if key in ("h_pt_lead", "h_pt_sub", "h_zpt", "h_nfsr"):
            ax.set_yscale("log")
    fig.suptitle("Kinematics of selected $Z\\to\\mu\\mu$ candidates (60-120 GeV)", fontsize=16)
    fig.tight_layout()
    hists.save_fig(fig, "step1_kinematics.png")

    fig, ax = plt.subplots(figsize=(9, 7))
    hep.histplot(res["h_dm_fsr"], ax=ax, histtype="step", color="black", linewidth=1.4)
    ax.set_xlabel(r"$m_{\mu\mu\gamma} - m_{\mu\mu}$ [GeV]")
    ax.set_ylabel("Events / 0.25 GeV")
    ax.set_yscale("log")
    ax.set_title("Mass recovered by FSR photons", fontsize=15, loc="left", pad=28)
    hists.save_fig(fig, "step1_fsr_mass_shift.png")


if __name__ == "__main__":
    main()
