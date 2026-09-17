#!/usr/bin/env python
"""v2 step 8 -- generator-level studies of the acceptance for the CMS-SMP-20-004 comparison (docs/16).

    python scripts/v2_8_theory_acceptance.py --files 12 --workers 4    # DY aMC@NLO parents, streamed from EOS
    python scripts/v2_8_theory_acceptance.py --summarise-only

CMS carries a "Resum. + FSR" uncertainty on the acceptance (DYTURBO NNLO+NNLL vs aMC@NLO, and
PYTHIA vs PHOTOS). Neither alternative prediction exists for this analysis, so both are estimated
(docs/16, section 4). This script supplies the numbers those estimates rest on, from the nominal
aMC@NLO sample (the full-sample sums in output/v2/gensums.json have no pT(Z) or bare-muon split):

  * A(pT(Z)): the acceptance in bins of the boson pT, so that a reshaped pT(Z) spectrum -- the
    thing resummation changes -- can be folded in;
  * A with bare (post-FSR, undressed) muons vs dressed muons: the size of the QED FSR effect on
    the fiducial volume, which an FSR-model difference can only change by a fraction of;
  * A in a CMS-like volume (dressed, pT > 25/25 GeV) with its 7-point scale and PS envelopes, to
    see whether CMS's larger scale uncertainty (0.66 %) is a property of their symmetric cuts.

Output: output/v2/cms_parity/theory_parts/*.pkl, output/v2/cms_parity/theory_acceptance.json,
plots output/v2/plots/cms_acc_*.png.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import awkward as ak
import numpy as np

from zmumu import batch, config, gen, hists, objects, samples
from zmumu.acceptance import reweight_acceptance

OUT = config.OUTPUT_DIR / "v2" / "cms_parity"
config.PLOT_DIR = config.OUTPUT_DIR / "v2" / "plots"
PTZ_EDGES = np.array([0, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 1000], dtype=float)
VOLUMES = ("dressed", "bare", "dressed2525", "bare2525")
BRANCHES = (["genWeight", "LHEScaleWeight", "PSWeight"] + gen.LHE_BRANCHES + gen.DRESSED_BRANCHES
            + [f"GenPart_{f}" for f in ("pt", "eta", "phi", "mass", "pdgId", "status", "statusFlags")])
SCALE_7 = [0, 1, 3, 4, 5, 7, 8]


def blank():
    nb = len(PTZ_EDGES) - 1
    out = {"n_files": 0, "n_events": 0, "den": np.zeros(nb), "den_w2": np.zeros(nb),
           "den_scale": np.zeros(9), "den_ps": np.zeros(4), "den_ps_ptz": np.zeros((4, nb)),
           "migr_dressed_not_bare": 0.0, "migr_bare_not_dressed": 0.0}
    for v in VOLUMES:
        out[f"fid_{v}"] = np.zeros(nb)
        out[f"fid_{v}_w2"] = np.zeros(nb)
        out[f"fid_{v}_scale"] = np.zeros(9)
        out[f"fid_{v}_ps"] = np.zeros(4)
    return out


def _regular(arr, n):
    counts = ak.to_numpy(ak.num(arr, axis=1))
    if counts.size == 0 or counts.min() != n or counts.max() != n:
        return None
    return ak.to_numpy(ak.to_regular(arr, axis=1)).astype(float)


def _fid(n, pt1, pt2, m, prod, lead, sub):
    return (n == 2) & (pt1 > lead) & (pt2 > sub) & (prod < 0) & (m > config.MASS_LO) & (m < config.MASS_HI)


def process_file(source):
    import uproot
    out = blank()
    out["n_files"] = 1
    with uproot.open(source) as f:
        tree = f["Events"]
        for ev in tree.iterate([b for b in BRANCHES if b in tree.keys()], step_size="300 MB", library="ak"):
            fill(ev, out)
    return out


def fill(ev, out):
    """Accumulate one chunk of DY aMC@NLO generator information."""
    out["n_events"] += len(ev)
    w = ak.to_numpy(ev.genWeight).astype(float)
    flav, mll, _ = gen.lhe_flavour(ev)
    den = (flav == 13) & (mll > config.MASS_LO) & (mll < config.MASS_HI)
    ev, w = ev[den], w[den]
    if len(ev) == 0:
        return
    # boson pT: last-copy Z if present (94 % of events), else the bare dimuon
    last_z = (ev.GenPart_pdgId == 23) & (((ev.GenPart_statusFlags >> 13) & 1) == 1)
    zpt = ak.to_numpy(ak.fill_none(ak.firsts(ev.GenPart_pt[last_z]), -1.0))
    bare = ((abs(ev.GenPart_pdgId) == 13) & (ev.GenPart_status == 1)
            & (((ev.GenPart_statusFlags >> 8) & 1) == 1)
            & (ev.GenPart_pt > 20) & (abs(ev.GenPart_eta) < config.MU_ETA_MAX))
    nb_, b1, b2, mb, pb = gen.leading_pair(ev.GenPart_pt, ev.GenPart_eta, ev.GenPart_phi, ev.GenPart_mass,
                                           ev.GenPart_pdgId, bare)
    allbare = ((abs(ev.GenPart_pdgId) == 13) & (ev.GenPart_status == 1)
               & (((ev.GenPart_statusFlags >> 8) & 1) == 1))
    rec = ak.zip({"pt": ev.GenPart_pt, "eta": ev.GenPart_eta, "phi": ev.GenPart_phi, "mass": ev.GenPart_mass})[allbare]
    rec = ak.pad_none(rec, 2, axis=1)
    px1, py1, _, _ = objects.p4(*(ak.to_numpy(ak.fill_none(rec[:, 0][k], 0.0)) for k in ("pt", "eta", "phi", "mass")))
    px2, py2, _, _ = objects.p4(*(ak.to_numpy(ak.fill_none(rec[:, 1][k], 0.0)) for k in ("pt", "eta", "phi", "mass")))
    zpt = np.where(zpt >= 0, zpt, np.hypot(px1 + px2, py1 + py2))

    dressed = ((abs(ev.GenDressedLepton_pdgId) == 13) & ~ev.GenDressedLepton_hasTauAnc
               & (ev.GenDressedLepton_pt > 20) & (abs(ev.GenDressedLepton_eta) < config.MU_ETA_MAX))
    nd, d1, d2, md, pd = gen.leading_pair(ev.GenDressedLepton_pt, ev.GenDressedLepton_eta,
                                          ev.GenDressedLepton_phi, ev.GenDressedLepton_mass,
                                          ev.GenDressedLepton_pdgId, dressed)
    fid = {"dressed": _fid(nd, d1, d2, md, pd, config.MU_PT_LEAD, config.MU_PT_SUBLEAD),
           "bare": _fid(nb_, b1, b2, mb, pb, config.MU_PT_LEAD, config.MU_PT_SUBLEAD),
           "dressed2525": _fid(nd, d1, d2, md, pd, 25.0, 25.0),
           "bare2525": _fid(nb_, b1, b2, mb, pb, 25.0, 25.0)}
    scale = _regular(ev.LHEScaleWeight, 9)
    ps = _regular(ev.PSWeight, 4)
    ib = np.clip(np.digitize(zpt, PTZ_EDGES) - 1, 0, len(PTZ_EDGES) - 2)
    nbins = len(PTZ_EDGES) - 1

    def binned(mask, weights):
        return np.bincount(ib[mask], weights=weights[mask], minlength=nbins)

    allm = np.ones(len(w), dtype=bool)
    out["den"] += binned(allm, w)
    out["den_w2"] += binned(allm, w ** 2)
    if scale is not None:
        out["den_scale"] += (w[:, None] * scale).sum(axis=0)
    if ps is not None:
        out["den_ps"] += (w[:, None] * ps).sum(axis=0)
        for k in range(4):
            out["den_ps_ptz"][k] += binned(allm, w * ps[:, k])
    for v, sel in fid.items():
        out[f"fid_{v}"] += binned(sel, w)
        out[f"fid_{v}_w2"] += binned(sel, w ** 2)
        if scale is not None:
            out[f"fid_{v}_scale"] += (w[sel, None] * scale[sel]).sum(axis=0)
        if ps is not None:
            out[f"fid_{v}_ps"] += (w[sel, None] * ps[sel]).sum(axis=0)
    out["migr_dressed_not_bare"] += float(w[fid["dressed"] & ~fid["bare"]].sum())
    out["migr_bare_not_dressed"] += float(w[fid["bare"] & ~fid["dressed"]].sum())


# ----------------------------------------------------------------------------- summary
def summarise(t):
    den = t["den"]
    res = {"n_files": t["n_files"], "n_events": t["n_events"], "ptz_edges": PTZ_EDGES.tolist(), "volumes": {},
           "den_ptz": den.tolist()}
    for v in VOLUMES:
        num = t[f"fid_{v}"]
        A = num.sum() / den.sum()
        A_stat = np.sqrt(max(t[f"fid_{v}_w2"].sum() * (1 - A) ** 2 + (t["den_w2"].sum() - t[f"fid_{v}_w2"].sum()) * A ** 2, 0)) / den.sum()
        sc = (t[f"fid_{v}_scale"] / t["den_scale"]) / A
        sc = sc / sc[4]
        ps = (t[f"fid_{v}_ps"] / t["den_ps"]) / A
        res["volumes"][v] = {"A": float(A), "A_stat": float(A_stat),
                             "scale_members_rel": (sc - 1).tolist(),
                             "scale_7pt_rel": float(np.max(np.abs(sc[SCALE_7] - 1))),
                             "ps_rel": (ps - 1).tolist(),
                             "A_ptz": np.where(den > 0, num / np.maximum(den, 1e-9), 0).tolist()}
    d = res["volumes"]
    res["fsr"] = {"A_bare_over_dressed": d["bare"]["A"] / d["dressed"]["A"],
                  "migr_dressed_not_bare_rel": t["migr_dressed_not_bare"] / t["fid_dressed"].sum(),
                  "migr_bare_not_dressed_rel": t["migr_bare_not_dressed"] / t["fid_dressed"].sum()}

    # --- resummation: reshape the boson-pT spectrum and fold it with A(pT(Z)) -----------------
    c = 0.5 * (PTZ_EDGES[1:] + PTZ_EDGES[:-1])
    A_bin = np.array(d["dressed"]["A_ptz"])
    A_bin_2525 = np.array(d["dressed2525"]["A_ptz"])
    A0, A0_2525 = d["dressed"]["A"], d["dressed2525"]["A"]
    shapes = {}
    for amp in (0.05, 0.10):
        # tilt: +amp at pT = 0 falling linearly to -amp at 30 GeV, flat beyond (rate renormalised)
        shapes[f"tilt_{int(amp*100)}pc_0_30"] = 1 + amp * (1 - 2 * np.clip(c, 0, 30) / 30)
        # Sudakov peak: +amp below 10 GeV, -amp between 10 and 30 GeV
        shapes[f"peak_{int(amp*100)}pc_10"] = np.where(c < 10, 1 + amp, np.where(c < 30, 1 - amp, 1.0))
    # the parton-shower ISR variation of the same sample, as a data-free reference shape
    for k, name in ((0, "ps_isr_up_shape"), (2, "ps_isr_down_shape")):
        ratio = (t["den_ps_ptz"][k] / t["den_ps"][k]) / (den / den.sum())
        shapes[name] = np.where(den > 0, ratio, 1.0)
    res["resummation"] = {}
    for name, s in shapes.items():
        res["resummation"][name] = {"shape": np.asarray(s).tolist(),
                                    "dA_rel": reweight_acceptance(A_bin, den, s) / A0 - 1,
                                    "dA_rel_2525": reweight_acceptance(A_bin_2525, den, s) / A0_2525 - 1}
    return res


def plot(res, t):
    import matplotlib.pyplot as plt
    e = PTZ_EDGES.copy()
    e[-1] = 200
    c = 0.5 * (e[1:] + e[:-1])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.5))
    den = t["den"]
    ax1.stairs(den / den.sum() / np.diff(e), e, color="black", label=r"aMC@NLO, LHE $\mu\mu$ 60-120 GeV")
    for k, lab, col in ((0, "PS ISR up", "#d62728"), (2, "PS ISR down", "#1f77b4")):
        s = t["den_ps_ptz"][k]
        ax1.stairs(s / s.sum() / np.diff(e), e, color=col, label=lab)
    s5 = np.array(res["resummation"]["tilt_5pc_0_30"]["shape"])
    ax1.stairs(den * s5 / (den * s5).sum() / np.diff(e), e, color="#2ca02c", linestyle="--", label="±5 % tilt (0-30 GeV), used for 3c")
    ax1.set_xscale("log")
    ax1.set_xlim(1, 200)
    ax1.set_xlabel(r"$p_T(Z)$ [GeV]")
    ax1.set_ylabel("normalised / GeV")
    ax1.legend(fontsize=10)
    for v, col, lab in (("dressed", "black", "our volume (dressed, 26/20)"), ("bare", "grey", "bare muons, 26/20"),
                        ("dressed2525", "#9467bd", "CMS-like (dressed, 25/25)")):
        ax2.plot(c, res["volumes"][v]["A_ptz"], "o-", color=col, label=lab, markersize=4)
    ax2.set_xscale("log")
    ax2.set_xlim(1, 200)
    ax2.set_xlabel(r"$p_T(Z)$ [GeV]")
    ax2.set_ylabel(r"acceptance A($p_T(Z)$)")
    ax2.legend(fontsize=10)
    for ax in (ax1, ax2):
        hists._decorate(ax)
    fig.suptitle("Acceptance vs boson pT: the input of the resummation estimate", fontsize=15)
    fig.tight_layout()
    hists.save_fig(fig, "cms_acc_vs_ptz.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--files", type=int, default=12, help="DY NLO parent files (of 41); the frozen result used 12")
    ap.add_argument("--prefer", choices=["dcache", "eos"], default="eos")
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--one-file", help=argparse.SUPPRESS)
    ap.add_argument("--key", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.one_file:
        batch.write_part(process_file(args.one_file), args.out)
        return
    parts = OUT / "theory_parts"
    if not args.summarise_only:
        srcs = samples.sources("DY_NLO", prefer=args.prefer)[: args.files]
        alt = {Path(s["name"]).stem: s["fallback"] for s in srcs}
        tasks = [(s["primary"], Path(s["name"]).stem) for s in srcs]
        batch.run_files(tasks, Path(__file__), parts, args.workers, stall_s=5400, retries=1,
                        fallback=lambda src, k: alt.get(k), log=lambda s: print(s, flush=True))
    t = batch.load_parts(parts, blank())
    res = summarise(t)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "theory_acceptance.json", "w") as fh:
        json.dump(res, fh, indent=1)
    plot(res, t)
    v = res["volumes"]
    print(f"[theory] {t['n_files']} files, {t['n_events']:,} events")
    for name in VOLUMES:
        print(f"  A({name}) = {v[name]['A']:.5f} +- {v[name]['A_stat']:.5f}; scale 7-pt {100*v[name]['scale_7pt_rel']:.3f} %;"
              f" PS {np.round(100*np.array(v[name]['ps_rel']), 3)} %")
    print(f"  FSR: A_bare/A_dressed = {res['fsr']['A_bare_over_dressed']:.5f}")
    for name, r in res["resummation"].items():
        print(f"  pT(Z) reshaping {name:22s}: dA/A = {100*r['dA_rel']:+.3f} % (ours), {100*r['dA_rel_2525']:+.3f} % (25/25)")


if __name__ == "__main__":
    main()
