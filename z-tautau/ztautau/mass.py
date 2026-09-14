"""Di-tau mass reconstruction: visible, collinear, and the MET-likelihood mass (docs/04-ditau-mass.md).

In Z -> tau tau at least four neutrinos escape, so the visible mass m_vis of the two tau_h peaks at
~60 GeV with a broad tail. The missing momentum (MET) carries the neutrino information:

collinear approximation
    Each neutrino system is collinear with its visible tau (the taus are boosted), so the visible
    fraction x_i of tau i's momentum obeys  MET = (1/x1 - 1) pT_vis1 + (1/x2 - 1) pT_vis2  (2 equations,
    2 unknowns) and  m_col = m_vis / sqrt(x1 x2). Exact in the collinear limit, but the 2x2 system is
    singular when the taus are back to back (the typical Z topology) and any MET mismeasurement pushes
    x outside (0, 1]: m_col is undefined in a large fraction of events and has long tails.

MET-likelihood mass (m_tt, the fit variable; same idea as CMS FastMTT / the "SVfit-lite" algorithms)
    Instead of solving the equations, scan (x1, x2) on a grid over the physical range
    x_i in [m_vis,i^2 / m_tau^2, 1], and weigh every point with
      * the MET transfer function  exp(-1/2 r^T V^-1 r),  r = MET_measured - MET_predicted(x1, x2),
        with V the per-event MET covariance matrix from NanoAOD (MET_covXX/XY/YY), and
      * the phase space of the hadronic decay tau -> tau_h nu, which for a two-body decay is flat in x,
      * a prior 1/m^2 on the di-tau mass m = m_vis / sqrt(x1 x2) (falling Drell-Yan spectrum; without it
        the weakly constrained low-x corner pulls the estimate up).
    The estimate is the posterior median of m. It is always defined, uses the MET
    resolution event by event, and its peak sits near m_Z.
"""

from __future__ import annotations

import numpy as np

from . import config
from .objects import inv_mass, p4


def visible_mass(t1, t2):
    """t = (pt, eta, phi, mass) numpy arrays."""
    a, b = p4(*t1), p4(*t2)
    return inv_mass(a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def collinear_mass(t1, t2, metx, mety):
    """(m_col, x1, x2, valid). Invalid events (singular system or x outside (0, 1]) get m_col = -1."""
    p1x, p1y = t1[0] * np.cos(t1[2]), t1[0] * np.sin(t1[2])
    p2x, p2y = t2[0] * np.cos(t2[2]), t2[0] * np.sin(t2[2])
    det = p1x * p2y - p1y * p2x
    with np.errstate(divide="ignore", invalid="ignore"):
        a1 = (metx * p2y - mety * p2x) / det          # a_i = 1/x_i - 1
        a2 = (mety * p1x - metx * p1y) / det
        x1 = 1.0 / (1.0 + a1)
        x2 = 1.0 / (1.0 + a2)
        valid = (np.abs(det) > 1e-6) & (x1 > 0) & (x1 <= 1) & (x2 > 0) & (x2 <= 1)
        m = np.where(valid, visible_mass(t1, t2) / np.sqrt(np.abs(x1 * x2)), -1.0)
    return m, x1, x2, valid


def likelihood_mass(t1, t2, metx, mety, covxx, covxy, covyy, n=None, chunk=None):
    """Posterior-median di-tau mass on an n x n grid of visible-energy fractions (see module doc)."""
    n = n or config.MASS_GRID_N
    chunk = chunk or config.MASS_CHUNK
    n_ev = len(metx)
    out = np.zeros(n_ev)
    mvis = visible_mass(t1, t2)
    u = (np.arange(n) + 0.5) / n                                          # grid in [0, 1]
    for start in range(0, n_ev, chunk):
        s = slice(start, start + chunk)
        pt1, phi1, m1 = t1[0][s], t1[2][s], t1[3][s]
        pt2, phi2, m2 = t2[0][s], t2[2][s], t2[3][s]
        xmin1 = np.clip((m1 / config.M_TAU) ** 2, 1e-4, 0.99)
        xmin2 = np.clip((m2 / config.M_TAU) ** 2, 1e-4, 0.99)
        x1 = xmin1[:, None] + (1 - xmin1)[:, None] * u[None, :]           # (ev, n)
        x2 = xmin2[:, None] + (1 - xmin2)[:, None] * u[None, :]
        a1 = 1.0 / x1 - 1.0
        a2 = 1.0 / x2 - 1.0
        nux = a1[:, :, None] * (pt1 * np.cos(phi1))[:, None, None] + a2[:, None, :] * (pt2 * np.cos(phi2))[:, None, None]
        nuy = a1[:, :, None] * (pt1 * np.sin(phi1))[:, None, None] + a2[:, None, :] * (pt2 * np.sin(phi2))[:, None, None]
        rx = metx[s][:, None, None] - nux
        ry = mety[s][:, None, None] - nuy
        cxx = np.maximum(covxx[s], config.COV_MIN)
        cyy = np.maximum(covyy[s], config.COV_MIN)
        cxy = covxy[s]
        det = cxx * cyy - cxy ** 2
        bad = det <= 0
        cxy = np.where(bad, 0.0, cxy)
        det = np.where(bad, cxx * cyy, det)
        ixx, iyy, ixy = (cyy / det)[:, None, None], (cxx / det)[:, None, None], (-cxy / det)[:, None, None]
        chi2 = ixx * rx ** 2 + iyy * ry ** 2 + 2 * ixy * rx * ry
        chi2 -= chi2.min(axis=(1, 2), keepdims=True)                      # numerical safety
        mgrid = mvis[s][:, None, None] / np.sqrt(x1[:, :, None] * x2[:, None, :])
        like = np.exp(-0.5 * chi2) * mgrid ** (-config.MASS_PRIOR_POW)
        # posterior median of m over the grid
        lf = like.reshape(len(like), -1)
        mf = mgrid.reshape(len(mgrid), -1)
        order = np.argsort(mf, axis=1)
        ms = np.take_along_axis(mf, order, axis=1)
        cdf = np.cumsum(np.take_along_axis(lf, order, axis=1), axis=1)
        cdf /= cdf[:, -1:]
        out[s] = ms[np.arange(len(ms)), np.minimum((cdf < 0.5).sum(axis=1), ms.shape[1] - 1)]
    return out


def all_masses(t1, t2, metx, mety, covxx, covxy, covyy):
    """dict(m_vis, m_col, m_tt, pt_tt, mt_tot) for numpy inputs t = (pt, eta, phi, mass)."""
    m_col, _, _, _ = collinear_mass(t1, t2, metx, mety)
    p1x, p1y = t1[0] * np.cos(t1[2]), t1[0] * np.sin(t1[2])
    p2x, p2y = t2[0] * np.cos(t2[2]), t2[0] * np.sin(t2[2])
    met = np.hypot(metx, mety)

    def mt(ptx, pty, qx, qy):
        pt, qt = np.hypot(ptx, pty), np.hypot(qx, qy)
        return np.sqrt(np.maximum(2 * (pt * qt - ptx * qx - pty * qy), 0.0))

    mt_tot = np.sqrt(mt(p1x, p1y, metx, mety) ** 2 + mt(p2x, p2y, metx, mety) ** 2 + mt(p1x, p1y, p2x, p2y) ** 2)
    return {"m_vis": visible_mass(t1, t2), "m_col": m_col,
            "m_tt": likelihood_mass(t1, t2, metx, mety, covxx, covxy, covyy),
            "pt_tt": np.hypot(p1x + p2x + metx, p1y + p2y + mety), "mt_tot": mt_tot, "met": met}
