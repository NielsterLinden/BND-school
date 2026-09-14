#!/usr/bin/env python
"""Step 2 -- muon efficiencies from tag-and-probe.

The cross section needs the probability that a Z -> mu mu decay inside the
fiducial volume actually lands in the selected sample. That probability is
measured *from the data itself* using tag-and-probe on the Z peak, so it does
not depend on simulation (there is none available here).

Three efficiencies are measured, in the order they are applied:

    eff_ID    P(medium ID + impact parameters | reconstructed muon)
    eff_ISO   P(pfRelIso04 < 0.15 | passes ID)
    eff_TRIG  P(event fires IsoMu24 or IsoTkMu24 | two selected muons)

The tag is a fully selected muon above the trigger plateau. The probe is a
*track-like* muon with no ID or isolation requirement, so the denominator is
not biased by the quantity being measured. Tag and probe must be opposite-sign
and form a mass in 70-110 GeV, which makes the pair overwhelmingly a real Z.
Both orientations of every pair are used, so an event where both muons qualify
as tags contributes two measurements.

What is *not* measured here, and why:

    eff_RECO  P(reconstructed muon | muon track). NanoAOD does not store
              generalTracks, so there is no unbiased denominator available.
              Taken from the CMS Muon POG as an external input
              (config.MU_RECO_EFF), and treated as a systematic.

The trigger efficiency uses the *reference trigger* method rather than
per-muon trigger-object matching, because TrigObj branches were not kept in
the skim. See docs/04-efficiency.md for what that costs.

Outputs
    output/data/step2_efficiency.pkl   efficiency maps + scalar averages
    output/plots/step2_*.png           efficiency vs pT, eta and the 2D maps

Run:  .venv/bin/python scripts/step2_efficiency.py
"""

from __future__ import annotations

import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import awkward as ak
import numpy as np

from zmumu import config, hists, io, objects, stats

BRANCHES = [
    "run", "luminosityBlock",
    "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
    "Muon_mediumId", "Muon_pfRelIso04_all", "Muon_dxy", "Muon_dz",
    "Muon_isTracker", "Muon_isGlobal",
    "PV_npvsGood",
] + config.TRIGGERS + config.REFERENCE_TRIGGERS + config.MET_FILTERS


def blank():
    """Pass/total histogram pairs for each efficiency, in pT, eta and 2D."""
    out = {}
    for eff in ("id", "iso"):
        out[f"{eff}_pass_pt"] = hists.variable(config.EFF_PT_BINS, r"probe $p_T$ [GeV]")
        out[f"{eff}_tot_pt"] = hists.variable(config.EFF_PT_BINS, r"probe $p_T$ [GeV]")
        out[f"{eff}_pass_eta"] = hists.variable(config.EFF_ETA_BINS, r"probe $\eta$")
        out[f"{eff}_tot_eta"] = hists.variable(config.EFF_ETA_BINS, r"probe $\eta$")
        out[f"{eff}_pass_2d"] = hists.hist2d(config.EFF_PT_BINS, config.EFF_ETA_BINS,
                                             r"probe $p_T$ [GeV]", r"probe $\eta$")
        out[f"{eff}_tot_2d"] = hists.hist2d(config.EFF_PT_BINS, config.EFF_ETA_BINS,
                                            r"probe $p_T$ [GeV]", r"probe $\eta$")
    # Tag-and-probe pair mass, for a purity check of the method itself.
    out["h_tp_mass_pass"] = hists.mass_hist()
    out["h_tp_mass_fail"] = hists.mass_hist()
    # Trigger efficiency (event level, reference-trigger method).
    out["trig_pass"] = 0
    out["trig_total"] = 0
    out["trig_pass_pt"] = hists.variable(config.EFF_PT_BINS, r"leading muon $p_T$ [GeV]")
    out["trig_tot_pt"] = hists.variable(config.EFF_PT_BINS, r"leading muon $p_T$ [GeV]")
    out["n_pairs"] = 0
    return out


def merge(a, b):
    for k, v in b.items():
        if isinstance(v, (int, float)):
            a[k] = a.get(k, 0) + v
        else:
            a[k] += v
    return a


def _clip_pt(pt):
    return np.clip(pt, config.EFF_PT_BINS[0], config.EFF_PT_BINS[-1] - 1e-3)


def process_file(path: Path) -> dict:
    grl = io.load_grl()
    out = blank()

    for events, _ in io.read_chunks(path, BRANCHES):
        events = events[io.lumi_mask(events.run, events.luminosityBlock, grl)]
        events = events[io.pass_met_filters(events) & (events.PV_npvsGood >= 1)]
        if len(events) == 0:
            continue

        _tag_and_probe(events[io.trigger_or(events, config.TRIGGERS)], out)
        _trigger_efficiency(events, out)

    return out


def _tag_and_probe(events, out):
    """Fill ID and isolation pass/total histograms from tag-probe pairs."""
    if len(events) == 0:
        return

    # A tag is a fully selected muon on the trigger plateau.
    is_tag = objects.good_muon_mask(events, pt_min=config.TAG_PT_MIN)
    # A probe is any track-like muon in acceptance: no ID, no isolation.
    is_probe = objects.loose_muon_mask(events)

    # Need at least one of each; cheap pre-filter before building pairs.
    keep = (ak.sum(is_tag, axis=1) >= 1) & (ak.sum(is_probe, axis=1) >= 2)
    events, is_tag, is_probe = events[keep], is_tag[keep], is_probe[keep]
    if len(events) == 0:
        return

    # All unordered muon pairs, then both orientations of each.
    idx = ak.local_index(events.Muon_pt, axis=1)
    pairs = ak.combinations(idx, 2, fields=["a", "b"])
    tag = ak.concatenate([pairs.a, pairs.b], axis=1)
    probe = ak.concatenate([pairs.b, pairs.a], axis=1)

    pt, eta, phi = events.Muon_pt, events.Muon_eta, events.Muon_phi
    mass_br, charge = events.Muon_mass, events.Muon_charge

    px1, py1, pz1, e1 = objects.p4(pt[tag], eta[tag], phi[tag], mass_br[tag])
    px2, py2, pz2, e2 = objects.p4(pt[probe], eta[probe], phi[probe], mass_br[probe])
    pair_mass = objects.invariant_mass(px1 + px2, py1 + py2, pz1 + pz2, e1 + e2)

    valid = (
        is_tag[tag]
        & is_probe[probe]
        & ((charge[tag] * charge[probe]) < 0)
        & (pair_mass > config.TP_MASS_LO)
        & (pair_mass < config.TP_MASS_HI)
    )

    probe_pt = ak.to_numpy(ak.flatten(pt[probe][valid]))
    probe_eta = ak.to_numpy(ak.flatten(eta[probe][valid]))
    if probe_pt.size == 0:
        return
    out["n_pairs"] += probe_pt.size

    # --- ID efficiency: denominator is every valid probe ------------------
    passes_id = ak.to_numpy(ak.flatten(
        (events[config.MU_ID_BRANCH][probe]
         & (abs(events.Muon_dxy[probe]) < config.MU_DXY_MAX)
         & (abs(events.Muon_dz[probe]) < config.MU_DZ_MAX))[valid]
    ))
    passes_iso = ak.to_numpy(ak.flatten(
        (events.Muon_pfRelIso04_all[probe] < config.MU_ISO_MAX)[valid]
    ))
    pair_mass_flat = ak.to_numpy(ak.flatten(pair_mass[valid]))

    cpt = _clip_pt(probe_pt)
    out["id_tot_pt"].fill(cpt)
    out["id_tot_eta"].fill(probe_eta)
    out["id_tot_2d"].fill(cpt, probe_eta)
    out["id_pass_pt"].fill(cpt[passes_id])
    out["id_pass_eta"].fill(probe_eta[passes_id])
    out["id_pass_2d"].fill(cpt[passes_id], probe_eta[passes_id])

    # --- isolation efficiency: denominator is probes that passed the ID ---
    out["iso_tot_pt"].fill(cpt[passes_id])
    out["iso_tot_eta"].fill(probe_eta[passes_id])
    out["iso_tot_2d"].fill(cpt[passes_id], probe_eta[passes_id])
    both = passes_id & passes_iso
    out["iso_pass_pt"].fill(cpt[both])
    out["iso_pass_eta"].fill(probe_eta[both])
    out["iso_pass_2d"].fill(cpt[both], probe_eta[both])

    out["h_tp_mass_pass"].fill(pair_mass_flat[both])
    out["h_tp_mass_fail"].fill(pair_mass_flat[~passes_id])


def _trigger_efficiency(events, out):
    """Event-level trigger efficiency with the reference-trigger method.

    Denominator: events with a full offline Z -> mu mu selection that fired one
    of the reference single-muon paths. Numerator: those that additionally
    fired IsoMu24 or IsoTkMu24. Because the reference paths are themselves
    single-muon triggers, this is not perfectly unbiased -- see
    docs/04-efficiency.md for the systematic assigned to it.
    """
    ref = io.trigger_or(events, config.REFERENCE_TRIGGERS)
    events = events[ref]
    if len(events) == 0:
        return

    good = objects.good_muon_mask(events)
    keep = ak.sum(good, axis=1) == config.N_MUONS_REQUIRED
    events, good = events[keep], good[keep]
    if len(events) == 0:
        return

    i1, i2 = objects.leading_two(events, good)
    pt1 = objects.take(events.Muon_pt, i1)
    on_plateau = ak.to_numpy(pt1) > config.MU_PT_LEAD
    opposite = ak.to_numpy(
        objects.take(events.Muon_charge, i1) * objects.take(events.Muon_charge, i2)
    ) < 0
    sel = on_plateau & opposite
    if not sel.any():
        return

    events = events[sel]
    lead_pt = _clip_pt(ak.to_numpy(pt1)[sel])
    fired = ak.to_numpy(io.trigger_or(events, config.TRIGGERS))

    out["trig_total"] += int(sel.sum())
    out["trig_pass"] += int(fired.sum())
    out["trig_tot_pt"].fill(lead_pt)
    out["trig_pass_pt"].fill(lead_pt[fired])


def _efficiency_map(pass_h, tot_h):
    """2D efficiency and its (symmetrised) uncertainty from pass/total."""
    p, t = pass_h.values(), tot_h.values()
    eff, lo, hi = stats.clopper_pearson(p, t)
    return eff, 0.5 * (lo + hi), t


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-files", type=int, default=None)
    ap.add_argument("--workers", type=int, default=config.N_WORKERS)
    args = ap.parse_args()

    files = io.find_files()
    if args.max_files:
        files = files[: args.max_files]
    print(f"[step2] tag-and-probe over {len(files)} files with {args.workers} workers")

    t0 = time.time()
    total = blank()
    with Pool(args.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(process_file, files), 1):
            total = merge(total, res)
            if i % 20 == 0 or i == len(files):
                print(f"        {i}/{len(files)} files  ({time.time()-t0:.0f}s)")

    result = _summarise(total)
    _plot(total, result)
    path = hists.save({"hists": total, "result": result}, "step2_efficiency.pkl")
    print(f"[step2] wrote {path}  ({time.time()-t0:.0f}s total)")


def _summarise(h):
    """Collapse the maps into the scalar numbers step 5 consumes."""
    id_eff, id_err, id_tot = _efficiency_map(h["id_pass_2d"], h["id_tot_2d"])
    iso_eff, iso_err, iso_tot = _efficiency_map(h["iso_pass_2d"], h["iso_tot_2d"])

    # Inclusive efficiencies (all probes pooled) -- quoted as the headline
    # numbers; step 5 uses the 2D maps folded with the signal kinematics.
    id_incl, id_lo, id_hi = stats.clopper_pearson(
        h["id_pass_2d"].values().sum(), h["id_tot_2d"].values().sum())
    iso_incl, iso_lo, iso_hi = stats.clopper_pearson(
        h["iso_pass_2d"].values().sum(), h["iso_tot_2d"].values().sum())
    trig_eff, trig_lo, trig_hi = stats.clopper_pearson(h["trig_pass"], h["trig_total"])

    result = {
        "id_map": id_eff, "id_map_err": id_err, "id_map_n": id_tot,
        "iso_map": iso_eff, "iso_map_err": iso_err, "iso_map_n": iso_tot,
        "id_inclusive": float(id_incl), "id_inclusive_err": float(0.5 * (id_lo + id_hi)),
        "iso_inclusive": float(iso_incl), "iso_inclusive_err": float(0.5 * (iso_lo + iso_hi)),
        "trig_event": float(trig_eff), "trig_event_err": float(0.5 * (trig_lo + trig_hi)),
        "trig_pass": h["trig_pass"], "trig_total": h["trig_total"],
        "n_pairs": h["n_pairs"],
        "reco_eff": config.MU_RECO_EFF, "reco_eff_err": config.MU_RECO_EFF_UNC,
    }

    print(f"\n  Tag-and-probe pairs used : {h['n_pairs']:,}")
    print(f"  eff(ID)   inclusive      : {id_incl:.4f} +/- {0.5*(id_lo+id_hi):.4f}")
    print(f"  eff(ISO)  inclusive      : {iso_incl:.4f} +/- {0.5*(iso_lo+iso_hi):.4f}")
    print(f"  eff(TRIG) event level    : {trig_eff:.4f} +/- {0.5*(trig_lo+trig_hi):.4f}"
          f"   ({h['trig_pass']:,}/{h['trig_total']:,})")
    print(f"  eff(RECO) external input : {config.MU_RECO_EFF:.4f} +/- {config.MU_RECO_EFF_UNC:.4f}")
    return result


def _plot(h, result):
    import matplotlib.pyplot as plt
    import mplhep as hep

    for eff, nice in (("id", "Medium ID"), ("iso", "Isolation")):
        for var, edges, xlabel in (("pt", config.EFF_PT_BINS, r"probe $p_T$ [GeV]"),
                                   ("eta", config.EFF_ETA_BINS, r"probe $\eta$")):
            e, lo, hi = stats.clopper_pearson(h[f"{eff}_pass_{var}"].values(),
                                              h[f"{eff}_tot_{var}"].values())
            hists.plot_efficiency(
                edges, e, lo, hi, f"step2_eff_{eff}_vs_{var}.png", xlabel,
                ylabel=f"{nice} efficiency",
                title=f"{nice} efficiency from tag-and-probe",
                ylim=(0.80, 1.02),
                reference=result[f"{eff}_inclusive"])

    # Trigger efficiency turn-on.
    e, lo, hi = stats.clopper_pearson(h["trig_pass_pt"].values(), h["trig_tot_pt"].values())
    hists.plot_efficiency(config.EFF_PT_BINS, e, lo, hi, "step2_eff_trigger_vs_pt.png",
                          r"leading muon $p_T$ [GeV]",
                          ylabel="Event trigger efficiency",
                          title="IsoMu24 OR IsoTkMu24 efficiency (reference-trigger method)",
                          ylim=(0.80, 1.02), reference=result["trig_event"])

    # 2D maps.
    for eff, nice in (("id", "Medium ID"), ("iso", "Isolation")):
        e, _, _ = _efficiency_map(h[f"{eff}_pass_2d"], h[f"{eff}_tot_2d"])
        fig, ax = plt.subplots(figsize=(10, 7))
        mesh = ax.pcolormesh(config.EFF_PT_BINS, config.EFF_ETA_BINS, e.T,
                             cmap="viridis", vmin=0.80, vmax=1.0)
        fig.colorbar(mesh, ax=ax, label=f"{nice} efficiency")
        ax.set_xscale("log")
        ax.set_xlabel(r"probe $p_T$ [GeV]")
        ax.set_ylabel(r"probe $\eta$")
        hep.cms.label("Open Data", data=True, lumi=round(config.LUMI_PB/1000, 1), year=2016, ax=ax)
        hists._title(ax, f"{nice} efficiency map")
        hists.save_fig(fig, f"step2_eff_{eff}_map.png")

    hists.plot_mass({"probe passes ID+iso": h["h_tp_mass_pass"],
                     "probe fails ID": h["h_tp_mass_fail"]},
                    "step2_tagprobe_mass.png",
                    title="Tag-and-probe pair mass (purity check)")


if __name__ == "__main__":
    main()
