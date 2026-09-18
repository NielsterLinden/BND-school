#!/usr/bin/env python
"""NOTE: written for the v1 review at the DeepTau Medium working point; since v3 the default is Tight, so run
with BND_TAUTAU_WP=Medium to reproduce the numbers quoted in REVIEW.md.

Where does the eta(tau1) non-closure of the fake factor come from? (REVIEW.md section 3.5 / 6)

    python review/eta_study.py        # after source ../fitting/setup.sh; writes review/eta_study.txt, review/review_eta.png
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from ztautau import analysis, config, fakes  # noqa

OUT = open(HERE / "eta_study.txt", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n")

data = analysis.load_data()
reg = analysis.regions(data)
T = (data["t1_vsjet"] & 16) > 0
base = (data["t1_pt"] > 40) & (data["t2_pt"] > 40) & ((data["t2_vsjet"] & 16) > 0) & ~data["os"].astype(bool)   # SS, tau2 Medium
aeta = np.abs(data["t1_eta"]); aeta2 = np.abs(data["t2_eta"])
ETA = np.array([0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1])
def ff_vs(mask_num, mask_den, x=aeta, edges=ETA):
    n = np.histogram(x[mask_num], bins=edges)[0].astype(float); d = np.histogram(x[mask_den], bins=edges)[0].astype(float)
    ff = n / np.maximum(d, 1); err = ff * np.sqrt(1 / np.maximum(n, 1) + 1 / np.maximum(d, 1))
    return ff, err, n, d
def row(label, ff, err, ref=None):
    r = ff / ref if ref is not None else ff
    P(f"{label:34s} " + " ".join(f"{v:6.3f}" for v in r) + ("" if ref is not None else "   +-" + " ".join(f"{e:5.3f}" for e in err)))

P("|eta(tau1)| bins: " + " ".join(f"{a:.1f}-{b:.1f}" for a, b in zip(ETA[:-1], ETA[1:])))
# 1. the FF itself vs |eta|, inclusive, with the loose definition used (VVVLoose & !Medium)
L = ((data["t1_vsjet"] & 1) > 0) & ~T
ff_incl, err_incl, n_incl, d_incl = ff_vs(base & T, base & L)
P("\n=== 1. FF(|eta|) = N(SS, tau1 T) / N(SS, tau1 L), tau2 Medium, inclusive in era/DM/Njets/pT ===")
row("FF inclusive", ff_incl, err_incl)
row("  relative to its mean", ff_incl / np.average(ff_incl, weights=d_incl), err_incl)
P("  N(T) per bin: " + " ".join(f"{v:6.0f}" for v in n_incl))

P("\n=== 2. is the eta shape universal? FF(|eta|) / FF(|eta| inclusive), per era, DM, Njets, pT ===")
mean_shape = ff_incl / np.average(ff_incl, weights=d_incl)
for name, sel in (("Run2016G", data["era"] == 0), ("Run2016H", data["era"] == 1),
                  ("DM0", data["t1_dm"] == 0), ("DM1", data["t1_dm"] == 1), ("DM10", data["t1_dm"] == 10), ("DM11", data["t1_dm"] == 11),
                  ("0 jets", data["njets"] == 0), ("1 jet", data["njets"] == 1), (">=2 jets", data["njets"] >= 2),
                  ("pT 40-45", (data["t1_pt"] < 45)), ("pT 45-50", (data["t1_pt"] >= 45) & (data["t1_pt"] < 50)),
                  ("pT 50-60", (data["t1_pt"] >= 50) & (data["t1_pt"] < 60)), ("pT 60-80", (data["t1_pt"] >= 60) & (data["t1_pt"] < 80)), ("pT > 80", data["t1_pt"] >= 80)):
    ff, err, n, d = ff_vs(base & T & sel, base & L & sel)
    shape = ff / np.average(ff, weights=d)
    row(f"{name} (shape / mean)", shape, err / np.average(ff, weights=d))

P("\n=== 3. dependence on the definition of 'loose': FF shape vs |eta| for tighter lower edges ===")
for lname, bit in (("VVVLoose & !Medium (used)", 1), ("VVLoose & !Medium", 2), ("VLoose & !Medium", 4), ("Loose & !Medium", 8)):
    Lx = ((data["t1_vsjet"] & bit) > 0) & ~T
    ff, err, n, d = ff_vs(base & T, base & Lx)
    row(f"{lname}", ff / np.average(ff, weights=d), err / np.average(ff, weights=d))
    P(f"    mean FF {np.average(ff, weights=d):.3f}, N(L) total {d.sum():.0f}")

P("\n=== 4. what the loose population looks like vs |eta|: fraction of the L taus in each raw-score slice ===")
raw = data["t1_rawvsjet"]
Lm = base & L
q = np.quantile(raw[Lm], [0.25, 0.5, 0.75])
P(f"raw DeepTau VSjet quartiles of the L taus: {np.round(q, 3).tolist()}")
for lo, hi, lab in ((0, q[0], "lowest quartile (least tau-like)"), (q[2], 2, "highest quartile (closest to Medium)")):
    sel = Lm & (raw >= lo) & (raw < hi)
    frac = np.histogram(aeta[sel], bins=ETA)[0] / np.maximum(np.histogram(aeta[Lm], bins=ETA)[0], 1)
    row(f"  {lab}", frac, np.zeros_like(frac))
# the same for the T taus: how tau-like are they vs eta
Tm = base & T
P("  median raw VSjet of T taus per |eta| bin: " + " ".join(f"{np.median(raw[Tm & (aeta >= a) & (aeta < b)]):6.3f}" for a, b in zip(ETA[:-1], ETA[1:])))
P("  median raw VSjet of L taus per |eta| bin: " + " ".join(f"{np.median(raw[Lm & (aeta >= a) & (aeta < b)]):6.3f}" for a, b in zip(ETA[:-1], ETA[1:])))
P("  Medium threshold in raw score (min raw of T): " + f"{raw[Tm].min():.4f}; VVVLoose threshold (min raw of L): {raw[Lm].min():.4f}")
# per-eta fraction of DM among T and L taus (DM composition of the fakes vs eta)
for dm in (0, 1, 10, 11):
    fT = np.histogram(aeta[Tm & (data["t1_dm"] == dm)], bins=ETA)[0] / np.maximum(np.histogram(aeta[Tm], bins=ETA)[0], 1)
    fL = np.histogram(aeta[Lm & (data["t1_dm"] == dm)], bins=ETA)[0] / np.maximum(np.histogram(aeta[Lm], bins=ETA)[0], 1)
    row(f"  DM{dm} fraction among T taus", fT, fT * 0); row(f"  DM{dm} fraction among L taus", fL, fL * 0)

P("\n=== 5. eta(tau2) and the tau2 leg: FF shape vs |eta(tau2)| (tau1 T / tau1 L, unchanged definition) ===")
ff, err, n, d = ff_vs(base & T, base & L, x=aeta2)
row("FF vs |eta(tau2)| / mean", ff / np.average(ff, weights=d), err / np.average(ff, weights=d))

P("\n=== 6. effect of adding |eta(tau1)| to the FF binning on the fake template ===")
table = fakes.measure(data, reg["SS_T"], reg["SS_L"])
osss = fakes.osss_correction(data, reg); C = np.asarray(osss["C"])
wf = fakes.fake_weights(data, reg["AR"], table, C)
# eta-binned: multiply the nominal FF by the per-(era, |eta|) closure ratio measured in SS (a closure correction),
# which is what an |eta| bin in the table would do to first order
ETA3 = np.array([0, 0.8, 1.5, 2.1])
corr = np.ones(len(data["t1_pt"]))
for e in (0, 1):
    me = data["era"] == e
    wss = fakes.fake_weights(data, reg["SS_L"] & me, table, 1.0)
    obs = np.histogram(aeta[reg["SS_T"] & me], bins=ETA3)[0]
    pred = np.histogram(aeta[reg["SS_L"] & me], bins=ETA3, weights=wss[reg["SS_L"] & me])[0]
    r = obs / pred
    P(f"  era {'GH'[e]}: SS closure ratio in |eta| bins {ETA3.tolist()}: " + " ".join(f"{v:.3f}" for v in r))
    idx = np.clip(np.searchsorted(ETA3, aeta, side="right") - 1, 0, len(ETA3) - 2)
    corr[me] = r[idx[me]]
wf_eta = wf * corr
e = np.asarray(config.FIT_BINS)
h0 = np.histogram(np.clip(data["m_tt"][reg["AR"]], 0, 349.9), bins=e, weights=wf[reg["AR"]])[0]
h1 = np.histogram(np.clip(data["m_tt"][reg["AR"]], 0, 349.9), bins=e, weights=wf_eta[reg["AR"]])[0]
P(f"  SR fake yield: nominal {h0.sum():.0f}, with |eta| correction {h1.sum():.0f} ({100*(h1.sum()/h0.sum()-1):+.2f}%)")
P("  m_tt bin        nominal   eta-corr   ratio")
for i in range(len(h0)):
    P(f"  {e[i]:4.0f}-{e[i+1]:4.0f}  {h0[i]:9.0f} {h1[i]:9.0f}   {h1[i]/h0[i]:.3f}")
win = (e[:-1] >= 60) & (e[:-1] < 110)
P(f"  fakes in 60-110 GeV: {h0[win].sum():.0f} -> {h1[win].sum():.0f} ({100*(h1[win].sum()/h0[win].sum()-1):+.2f}%); signal there ~ 4000")
# SS closure in eta after the correction (sanity) and in m_tt
wss0 = fakes.fake_weights(data, reg["SS_L"], table, 1.0); wss1 = wss0 * corr
for var, edges in (("t1_eta", np.linspace(-2.1, 2.1, 15)), ("m_tt", e)):
    x = data[var]
    obs = np.histogram(np.clip(x[reg["SS_T"]], edges[0], edges[-1] - 1e-6), bins=edges)[0]
    p0 = np.histogram(np.clip(x[reg["SS_L"]], edges[0], edges[-1] - 1e-6), bins=edges, weights=wss0[reg["SS_L"]])[0]
    p1 = np.histogram(np.clip(x[reg["SS_L"]], edges[0], edges[-1] - 1e-6), bins=edges, weights=wss1[reg["SS_L"]])[0]
    P(f"  SS closure in {var}: obs/pred nominal   " + " ".join(f"{v:.2f}" for v in obs / np.maximum(p0, 1)))
    P(f"  SS closure in {var}: obs/pred eta-corr  " + " ".join(f"{v:.2f}" for v in obs / np.maximum(p1, 1)))

fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
c = 0.5 * (ETA[1:] + ETA[:-1])
ax = axes[0]
for dm, col in ((0, "k"), (1, "tab:red"), (10, "tab:blue"), (11, "tab:green")):
    sel = data["t1_dm"] == dm
    ff, err, n, d = ff_vs(base & T & sel, base & L & sel)
    ax.errorbar(c, ff, yerr=err, fmt="o-", color=col, label=f"DM{dm}")
ax.set_xlabel(r"$|\eta(\tau_1)|$"); ax.set_ylabel("FF = N(Medium)/N(VVVLoose & !Medium), SS"); ax.legend(); ax.set_title("fake factor vs |eta|, per decay mode")
ax = axes[1]
for lname, bit, col in (("VVVLoose", 1, "k"), ("VVLoose", 2, "tab:red"), ("VLoose", 4, "tab:blue"), ("Loose", 8, "tab:green")):
    Lx = ((data["t1_vsjet"] & bit) > 0) & ~T
    ff, err, n, d = ff_vs(base & T, base & Lx)
    ax.errorbar(c, ff / np.average(ff, weights=d), yerr=err / np.average(ff, weights=d), fmt="o-", color=col, label=f"L = {lname} & !Medium")
ax.axhline(1, color="grey"); ax.set_xlabel(r"$|\eta(\tau_1)|$"); ax.set_ylabel("FF(|eta|) / mean FF"); ax.legend(); ax.set_title("eta shape of the FF vs loose definition")
ax = axes[2]
ax.step(e[:-1], h1 / h0, where="post", color="k"); ax.axhline(1, color="grey"); ax.set_ylim(0.9, 1.1)
ax.set_xlabel(r"$m_{\tau\tau}$ [GeV]"); ax.set_ylabel("fake template: |eta|-corrected / nominal"); ax.set_title("effect of an |eta| closure correction on the fit template")
fig.tight_layout(); fig.savefig(HERE / "review_eta.png", dpi=110)
P(f"\nplot -> {HERE}/review_eta.png")
OUT.close()
