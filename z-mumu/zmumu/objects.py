"""Muon selection, FSR recovery, and dimuon kinematics.

Everything here is pure: functions take an awkward record array of one chunk
and return masks or arrays. No I/O, no configuration reading beyond the
`config` defaults, so each piece can be tested on a single chunk in isolation.
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from . import config


# --------------------------------------------------------------------------
# Four-vector helpers
# --------------------------------------------------------------------------
def p4(pt, eta, phi, mass):
    """(px, py, pz, E) from the NanoAOD (pt, eta, phi, mass) representation."""
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
    return px, py, pz, energy


def invariant_mass(px, py, pz, energy):
    """sqrt(E^2 - p^2), clipped at zero to absorb float rounding."""
    m2 = energy**2 - px**2 - py**2 - pz**2
    return np.sqrt(np.maximum(m2, 0.0))


def delta_phi(a, b):
    """Signed phi difference wrapped into (-pi, pi]."""
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    return np.sqrt((eta1 - eta2) ** 2 + delta_phi(phi1, phi2) ** 2)


# --------------------------------------------------------------------------
# Muon selection
# --------------------------------------------------------------------------
def good_muon_mask(events, pt_min: float = None) -> ak.Array:
    """Per-muon mask for analysis-quality ("tight") muons.

    Medium ID + tight relative isolation + prompt impact parameters, inside the
    muon acceptance. This is the object definition of the fiducial volume.
    """
    pt_min = config.MU_PT_SUBLEAD if pt_min is None else pt_min
    return (
        events[config.MU_ID_BRANCH]
        & (events.Muon_pt > pt_min)
        & (abs(events.Muon_eta) < config.MU_ETA_MAX)
        & (events.Muon_pfRelIso04_all < config.MU_ISO_MAX)
        & (abs(events.Muon_dxy) < config.MU_DXY_MAX)
        & (abs(events.Muon_dz) < config.MU_DZ_MAX)
    )


def loose_muon_mask(events, pt_min: float = None) -> ak.Array:
    """Per-muon mask for the *probe* denominator in tag-and-probe.

    Deliberately as loose as NanoAOD allows -- a track-like muon inside
    acceptance, with no ID or isolation requirement -- so that the measured
    efficiency is not biased by the quantity being measured.
    """
    pt_min = config.PROBE_PT_MIN if pt_min is None else pt_min
    return (
        (events.Muon_isTracker | events.Muon_isGlobal)
        & (events.Muon_pt > pt_min)
        & (abs(events.Muon_eta) < config.MU_ETA_MAX)
    )


def good_electron_mask(events) -> ak.Array:
    """Per-electron mask used only for the flavour-symmetric background."""
    return (
        (events.Electron_pt > config.EL_PT_MIN)
        & (abs(events.Electron_eta) < config.EL_ETA_MAX)
        & (events.Electron_cutBased >= config.EL_ID_MIN)
        & (events.Electron_pfRelIso03_all < config.EL_ISO_MAX)
    )


# --------------------------------------------------------------------------
# Dimuon construction
# --------------------------------------------------------------------------
def take(jagged, idx):
    """Element `idx[i]` of event `i` from a jagged array.

    Awkward reads *every* array-valued index on the outer axis as an event
    selection -- flat, 2D and jagged indices all pick events rather than
    elements. The unambiguous way to pick element `idx[i]` from event `i` is to
    build a positional boolean mask and flatten it, which is what this does.
    Kept as a helper so the mistake cannot be made by accident.
    """
    position = ak.local_index(jagged, axis=1)
    return ak.flatten(jagged[position == np.asarray(idx)])


def leading_two(events, mask):
    """Indices of the two highest-pT muons passing `mask`, pT-ordered.

    Returns (idx1, idx2) as flat integer arrays indexing into the full Muon
    collection. Use `take()` to read muon properties with them. Only defined
    for events with >= 2 muons passing the mask.
    """
    order = ak.argsort(events.Muon_pt, axis=1, ascending=False)
    masked_order = order[mask[order]]
    return (ak.to_numpy(masked_order[:, 0]).astype(np.int64),
            ak.to_numpy(masked_order[:, 1]).astype(np.int64))


def fsr_photon_mask(events) -> ak.Array:
    """Quality mask for FSR photon candidates (see docs/02-fsr-photons.md)."""
    return (
        (events.FsrPhoton_pt > config.FSR_PT_MIN)
        & (events.FsrPhoton_relIso03 < config.FSR_REL_ISO_MAX)
        & (events.FsrPhoton_dROverEt2 < config.FSR_DR_OVER_ET2_MAX)
    )


def dimuon_mass(events, idx1, idx2, with_fsr: bool = False):
    """Invariant mass of the two muons, optionally with FSR photons added.

    `with_fsr` adds every quality FSR photon whose `FsrPhoton_muonIdx` points at
    one of the two muons. NanoAOD has already done the matching: `muonIdx` is
    the index of the muon the photon was associated to during reconstruction,
    so we only have to apply quality cuts and sum four-vectors.
    """
    pt = events.Muon_pt
    eta = events.Muon_eta
    phi = events.Muon_phi
    mass = events.Muon_mass

    px1, py1, pz1, e1 = p4(take(pt, idx1), take(eta, idx1),
                           take(phi, idx1), take(mass, idx1))
    px2, py2, pz2, e2 = p4(take(pt, idx2), take(eta, idx2),
                           take(phi, idx2), take(mass, idx2))
    px, py, pz, energy = px1 + px2, py1 + py2, pz1 + pz2, e1 + e2

    if with_fsr and "FsrPhoton_pt" in ak.fields(events):
        attached = fsr_photon_mask(events) & (
            (events.FsrPhoton_muonIdx == idx1) | (events.FsrPhoton_muonIdx == idx2)
        )
        gpt = events.FsrPhoton_pt[attached]
        geta = events.FsrPhoton_eta[attached]
        gphi = events.FsrPhoton_phi[attached]
        gpx, gpy, gpz, ge = p4(gpt, geta, gphi, 0.0)
        px = px + ak.sum(gpx, axis=1)
        py = py + ak.sum(gpy, axis=1)
        pz = pz + ak.sum(gpz, axis=1)
        energy = energy + ak.sum(ge, axis=1)

    return invariant_mass(px, py, pz, energy)


def dimuon_pt_y(events, idx1, idx2):
    """Transverse momentum and rapidity of the dimuon system."""
    px1, py1, pz1, e1 = p4(
        take(events.Muon_pt, idx1), take(events.Muon_eta, idx1),
        take(events.Muon_phi, idx1), take(events.Muon_mass, idx1),
    )
    px2, py2, pz2, e2 = p4(
        take(events.Muon_pt, idx2), take(events.Muon_eta, idx2),
        take(events.Muon_phi, idx2), take(events.Muon_mass, idx2),
    )
    px, py, pz, energy = px1 + px2, py1 + py2, pz1 + pz2, e1 + e2
    pt = np.sqrt(px**2 + py**2)
    # Rapidity; guard the log against numerical equality of E and |pz|.
    denom = energy - pz
    rapidity = 0.5 * np.log(np.where(denom > 0, (energy + pz) / np.where(denom > 0, denom, 1), 1.0))
    return pt, rapidity
