#!/usr/bin/env python
"""NOTE: written for the v1 review at the DeepTau Medium working point; since v3 the default is Tight, so run
with BND_TAUTAU_WP=Medium to reproduce the numbers quoted in REVIEW.md.

Numerical cross-checks for REVIEW.md (run from z-tautau/ after `source ../setup.sh`):

    python review/fake_studies.py

A. DeepTau VSjet working-point scan (both legs): signal, fakes (= data - MC with genuine tau1) and the
   stat-only sensitivity from the m_tt shape.
B. tau pT threshold scan at Medium.
C. k-fold BDT (XGBoost) proof of concept: Z->tautau simulation vs the nominal fake estimate (AR x FF),
   mass-agnostic inputs, evaluated on held-out folds (fold = event number mod 5); S/B per score category,
   stat-only sensitivity of an m_tt fit in score categories, and the same-sign closure of the FF method
   in the BDT score (the validation any ML selection on top of a fake factor must pass).
Outputs: review/fake_studies.txt and review/review_*.png. Nothing in output/ or fit/ is touched.
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from ztautau import analysis, config, fakes, samples  # noqa

OUT = open(HERE / "fake_studies.txt", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n")

t0 = time.time()
data = analysis.load_data()
reg = analysis.regions(data)
P(f"data events in ntuple: {len(data['run']):,}; SR {reg['SR'].sum():,}; AR {reg['AR'].sum():,}; SS_T {reg['SS_T'].sum():,}")
mc = []
for key in analysis.available_mc():
    d, _ = analysis.load(key)
    if not len(d):
        continue
    w = analysis.weights(d, key)
    mc.append((key, d, analysis.regions(d, is_mc=True), w))
sig_key = "DY_NLO"
dsig, _ = analysis.load(sig_key)
wsig = [w for k, _, _, w in mc if k == sig_key][0]
rsig = [r for k, _, r, _ in mc if k == sig_key][0]
is_tt = dsig["gen_lhe_flavour"] == 15
P(f"loaded in {time.time()-t0:.0f}s")

edges = np.asarray(config.FIT_BINS)
def h(x, w, e=edges):
    xc = np.clip(x, e[0], e[-1] - 1e-6); return np.histogram(xc, bins=e, weights=w)[0]
def sens(s, b):
    ok = (s + b) > 0; return 1.0 / np.sqrt(np.sum(s[ok] ** 2 / (s[ok] + b[ok])))

def sr_mask(arr, tbit, ptmin, is_mc):
    base = (arr["t1_pt"] > ptmin) & (arr["t2_pt"] > ptmin) & arr["os"].astype(bool)
    base &= ((arr["t1_vsjet"] & tbit) > 0) & ((arr["t2_vsjet"] & tbit) > 0)
    if is_mc:
        base &= arr["t1_genflav"] != 0
    return base

def counts(tbit, ptmin):
    md = sr_mask(data, tbit, ptmin, False)
    n_data = h(data["m_tt"][md], np.ones(md.sum()))
    n_mc = np.zeros_like(n_data)
    for k, d, r, w in mc:
        m = sr_mask(d, tbit, ptmin, True)
        n_mc += h(d["m_tt"][m], w[m])
    ms = sr_mask(dsig, tbit, ptmin, True) & is_tt
    s = h(dsig["m_tt"][ms], wsig[ms])
    fake = np.maximum(n_data - n_mc, 0)
    win = (edges[:-1] >= 70) & (edges[:-1] < 110)
    return dict(data=n_data.sum(), sig=s.sum(), fake=fake.sum(), s_win=s[win].sum(), f_win=fake[win].sum(),
                sens=sens(s, n_data - s), raw_sig=int(ms.sum()))

P("\n=== A. DeepTau VSjet working point (both legs), pT > 40 GeV; ID SF kept at the Medium values ===")
P(f"{'WP':8s} {'data':>8s} {'signal':>8s} {'fakes':>8s} {'fake frac':>9s} {'S/B(70-110)':>11s} {'sig kept':>8s} {'stat-only d(mu)':>15s} {'raw sig MC':>10s}")
ref = None
for name, bit in (("Medium", 16), ("Tight", 32), ("VTight", 64), ("VVTight", 128)):
    c = counts(bit, 40.0); ref = ref or c
    P(f"{name:8s} {c['data']:8.0f} {c['sig']:8.0f} {c['fake']:8.0f} {c['fake']/c['data']:9.2f} {c['s_win']/max(c['f_win'],1):11.2f} "
      f"{c['sig']/ref['sig']:8.2f} {100*c['sens']:14.2f}% {c['raw_sig']:10d}")

P("\n=== B. tau pT threshold (both legs), Medium ===")
P(f"{'pT>':6s} {'data':>8s} {'signal':>8s} {'fakes':>8s} {'fake frac':>9s} {'S/B(70-110)':>11s} {'sig kept':>8s} {'stat-only d(mu)':>15s} {'<trig SF unc>':>13s}")
from ztautau import corrections
pog = corrections.pog()["trigger_ditau_Medium"]; grid = np.asarray(pog["pt"])
def trig_unc(pt, dm):
    idx = np.clip(np.searchsorted(grid, pt, side="right") - 1, 0, len(grid) - 1)
    out = np.zeros(len(pt))
    for d in config.TAU_DMS:
        sel = dm == d
        out[sel] = np.asarray(pog[str(d)]["sf_err"])[idx[sel]] / np.maximum(np.asarray(pog[str(d)]["sf"])[idx[sel]], 1e-3)
    return out
ref = None
for ptmin in (40.0, 45.0, 50.0, 55.0, 60.0):
    c = counts(16, ptmin); ref = ref or c
    ms = sr_mask(dsig, 16, ptmin, True) & is_tt
    u = np.sqrt(trig_unc(dsig["t1_pt"][ms], dsig["t1_dm"][ms]) ** 2 + trig_unc(dsig["t2_pt"][ms], dsig["t2_dm"][ms]) ** 2)
    P(f"{ptmin:6.0f} {c['data']:8.0f} {c['sig']:8.0f} {c['fake']:8.0f} {c['fake']/c['data']:9.2f} {c['s_win']/max(c['f_win'],1):11.2f} "
      f"{c['sig']/ref['sig']:8.2f} {100*c['sens']:14.2f}% {100*np.average(u, weights=np.maximum(wsig[ms],0)):12.1f}%")

# ------------------------------------------------------------------ C. k-fold BDT
P("\n=== C. k-fold BDT, Z->tautau (SR, simulation) vs fakes (data AR x FF x C, nominal) ===")
table = fakes.measure(data, reg["SS_T"], reg["SS_L"])
osss = fakes.osss_correction(data, reg)
C = np.asarray(osss["C"])
wf = fakes.fake_weights(data, reg["AR"], table, C)
P(f"nominal FF recomputed: fakes in SR = {wf[reg['AR']].sum():.0f} (committed prefit yield 38758); C inclusive {osss['inclusive']['C']:.3f}")

def dphi(a, b):
    return np.abs((a - b + np.pi) % (2 * np.pi) - np.pi)
def features(arr):
    ptt = np.hypot(arr["t1_pt"] * np.cos(arr["t1_phi"]) + arr["t2_pt"] * np.cos(arr["t2_phi"]),
                   arr["t1_pt"] * np.sin(arr["t1_phi"]) + arr["t2_pt"] * np.sin(arr["t2_phi"]))
    metphi = np.arctan2(arr["met_y"], arr["met_x"])
    ttphi = np.arctan2(arr["t1_pt"] * np.sin(arr["t1_phi"]) + arr["t2_pt"] * np.sin(arr["t2_phi"]),
                       arr["t1_pt"] * np.cos(arr["t1_phi"]) + arr["t2_pt"] * np.cos(arr["t2_phi"]))
    metsig = arr["met"] / np.sqrt(np.maximum(arr["met_covxx"] + arr["met_covyy"], 1.0))
    X = np.stack([arr["t1_pt"], arr["t2_pt"], arr["t2_pt"] / arr["t1_pt"], np.abs(arr["t1_eta"]), np.abs(arr["t2_eta"]),
                  arr["dr_tt"], dphi(arr["t1_phi"], arr["t2_phi"]), arr["met"], metsig, ptt, dphi(metphi, ttphi),
                  arr["pt_tt"], arr["njets"], arr["jet1_pt"], arr["t1_dm"], arr["t2_dm"]], axis=1).astype(np.float32)
    return X
FEAT = ["pt1", "pt2", "pt2/pt1", "|eta1|", "|eta2|", "dR", "dphi(t1,t2)", "MET", "MET sig.", "pT(vis tt)", "dphi(MET,tt)",
        "pT(tt+MET)", "Njets", "jet1 pT", "DM1", "DM2"]
import xgboost as xgb
K = 5
sig_m = rsig["SR"] & is_tt
Xs, ws_, fs = features({k: v[sig_m] for k, v in dsig.items()}), np.maximum(wsig[sig_m], 0.0), (dsig["event"][sig_m] % K)
ar = reg["AR"]
Xb, wb, fb = features({k: v[ar] for k, v in data.items()}), wf[ar], (data["event"][ar] % K)
Xssl, wssl, fssl = features({k: v[reg["SS_L"]] for k, v in data.items()}), fakes.fake_weights(data, reg["SS_L"], table, 1.0)[reg["SS_L"]], data["event"][reg["SS_L"]] % K
Xsst, fsst = features({k: v[reg["SS_T"]] for k, v in data.items()}), data["event"][reg["SS_T"]] % K
Xsr, fsr = features({k: v[reg["SR"]] for k, v in data.items()}), data["event"][reg["SR"]] % K
P(f"training samples: signal {len(Xs):,} MC events (sum w {ws_.sum():.0f}), fakes {len(Xb):,} AR events (sum w {wb.sum():.0f}); {K} folds by event number")
ws_eval = wsig[sig_m]                      # signed weights for yields; positive weights for training
models = []
imp = np.zeros(len(FEAT))
params = dict(n_estimators=300, max_depth=4, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=5, n_jobs=8,
              tree_method="hist", eval_metric="auc")
for k in range(K):
    tr_s, tr_b = fs != k, fb != k
    X = np.concatenate([Xs[tr_s], Xb[tr_b]]); y = np.concatenate([np.ones(tr_s.sum()), np.zeros(tr_b.sum())])
    # both classes carry the same total weight, with a mean weight of one per event of the smaller class
    n_ref = float(tr_s.sum())
    w = np.concatenate([ws_[tr_s] / ws_[tr_s].sum() * n_ref, wb[tr_b] / wb[tr_b].sum() * n_ref])
    clf = xgb.XGBClassifier(**params).fit(X, y, sample_weight=w)
    models.append(clf)
    imp += clf.feature_importances_ / K
def score(X, folds):
    out = np.zeros(len(X))
    for k in range(K):
        m = folds == k
        if m.any():
            out[m] = models[k].predict_proba(X[m])[:, 1]
    return out
score_s, score_b = score(Xs, fs), score(Xb, fb)
score_ssl, score_sst, score_sr = score(Xssl, fssl), score(Xsst, fsst), score(Xsr, fsr)
from sklearn.metrics import roc_auc_score
auc = roc_auc_score(np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))]), np.concatenate([score_s, score_b]),
                    sample_weight=np.concatenate([ws_ / ws_.sum(), wb / wb.sum()]))
P(f"held-out AUC (signal vs fakes) = {auc:.3f}")
P("feature importance (gain, mean over folds): " + ", ".join(f"{n} {v:.2f}" for n, v in sorted(zip(FEAT, imp), key=lambda t: -t[1])))

# score categories: quantiles of the signal score
qs = np.quantile(score_s, [0.0, 0.25, 0.5, 0.75, 1.0], method="inverted_cdf") if False else np.quantile(score_s, [0.0, 0.25, 0.5, 0.75])
cat_edges = np.concatenate([[0.0], qs[1:], [1.0001]])
P(f"score categories (signal quartiles): edges {np.round(cat_edges, 3).tolist()}")
P(f"{'cat':>4s} {'signal':>8s} {'fakes':>8s} {'S/B':>6s} {'S/B(70-110)':>11s} {'data SR':>8s}")
mt_s = dsig["m_tt"][sig_m]; mt_b = data["m_tt"][ar]; mt_sr = data["m_tt"][reg["SR"]]
tot_incl = sens(h(mt_s, ws_eval), h(mt_b, wb)); tot_cat = 0.0
win = (edges[:-1] >= 70) & (edges[:-1] < 110)
for i in range(len(cat_edges) - 1):
    cs, cb, csr = (score_s >= cat_edges[i]) & (score_s < cat_edges[i + 1]), (score_b >= cat_edges[i]) & (score_b < cat_edges[i + 1]), (score_sr >= cat_edges[i]) & (score_sr < cat_edges[i + 1])
    hs, hb = h(mt_s[cs], ws_eval[cs]), h(mt_b[cb], wb[cb])
    tot_cat += np.sum(hs ** 2 / np.maximum(hs + hb, 1e-9))
    P(f"{i:4d} {hs.sum():8.0f} {hb.sum():8.0f} {hs.sum()/max(hb.sum(),1):6.2f} {hs[win].sum()/max(hb[win].sum(),1):11.2f} {csr.sum():8d}")
P(f"stat-only d(mu)/mu, backgrounds fixed: inclusive m_tt fit {100*tot_incl:.2f}%  ->  m_tt fit in 4 score categories {100/np.sqrt(tot_cat):.2f}%")
top = score_s >= cat_edges[-2]
P(f"top category alone: signal {ws_eval[top].sum():.0f} ({100*ws_eval[top].sum()/ws_eval.sum():.0f}% of signal), fakes {wb[score_b >= cat_edges[-2]].sum():.0f}")

# same-sign closure of the FF method in the score
sc_edges = np.linspace(0, 1, 11)
obs = np.histogram(score_sst, bins=sc_edges)[0]
pred = np.histogram(score_ssl, bins=sc_edges, weights=wssl)[0]
pred2 = np.histogram(score_ssl, bins=sc_edges, weights=wssl ** 2)[0]
P("same-sign closure in the BDT score (SS_T observed / FF x SS_L predicted):")
for i in range(len(obs)):
    r = obs[i] / max(pred[i], 1e-9); e = r * np.sqrt(1 / max(obs[i], 1) + pred2[i] / max(pred[i], 1e-9) ** 2)
    P(f"   score {sc_edges[i]:.1f}-{sc_edges[i+1]:.1f}: obs {obs[i]:7d}  pred {pred[i]:8.0f}  ratio {r:.3f} +- {e:.3f}")
# OS: data SR vs prediction in the score (nominal, MC not subtracted from AR)
sr_obs = np.histogram(score_sr, bins=sc_edges)[0]
sr_fake = np.histogram(score_b, bins=sc_edges, weights=wb)[0]
sr_mc = np.zeros(len(sr_obs))
for k, d, r, w in mc:
    if k == sig_key:
        continue
    m = r["SR"]
    sc = score(features({kk: v[m] for kk, v in d.items()}), d["event"][m] % K)
    sr_mc += np.histogram(sc, bins=sc_edges, weights=w[m])[0]
sr_sig = np.histogram(score_s, bins=sc_edges, weights=wsig[sig_m])[0]
P("signal region in the BDT score (data / (fakes + signal + other MC)):")
for i in range(len(sr_obs)):
    P(f"   score {sc_edges[i]:.1f}-{sc_edges[i+1]:.1f}: data {sr_obs[i]:7d}  fakes {sr_fake[i]:8.0f}  Z->tautau {sr_sig[i]:7.0f}  other {sr_mc[i]:6.0f}  ratio {sr_obs[i]/max(sr_fake[i]+sr_sig[i]+sr_mc[i],1e-9):.3f}")

fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
ax = axes[0]
ax.hist(score_s, bins=sc_edges, weights=ws_eval / ws_eval.sum(), histtype="step", lw=2, label=r"Z$\to\tau\tau$ (MC, SR)")
ax.hist(score_b, bins=sc_edges, weights=wb / wb.sum(), histtype="step", lw=2, label="fakes (AR x FF)")
ax.set_xlabel("BDT score (held-out fold)"); ax.set_ylabel("normalised"); ax.legend(); ax.set_title(f"k-fold BDT, AUC = {auc:.3f}")
ax = axes[1]
c = 0.5 * (sc_edges[1:] + sc_edges[:-1])
ax.errorbar(c, obs / np.maximum(pred, 1e-9), yerr=np.sqrt(obs) / np.maximum(pred, 1e-9), fmt="o", color="k")
ax.axhline(1, color="grey"); ax.set_ylim(0.5, 1.5); ax.set_xlabel("BDT score"); ax.set_ylabel("SS_T / (FF x SS_L)"); ax.set_title("same-sign closure of the FF in the score")
ax = axes[2]
ax.bar(c, sr_fake, width=0.1, color="#f9a8d4", label="fakes"); ax.bar(c, sr_sig, width=0.1, bottom=sr_fake, color="#fbbf24", label=r"Z$\to\tau\tau$")
ax.bar(c, sr_mc, width=0.1, bottom=sr_fake + sr_sig, color="#94a3b8", label="other MC")
ax.errorbar(c, sr_obs, yerr=np.sqrt(sr_obs), fmt="o", color="k", label="data"); ax.set_yscale("log"); ax.set_xlabel("BDT score"); ax.set_title("signal region (prefit, nominal FF)"); ax.legend()
fig.tight_layout(); fig.savefig(HERE / "review_bdt.png", dpi=110); plt.close(fig)
P(f"\ndone in {time.time()-t0:.0f}s; plots -> {HERE}/review_bdt.png")
OUT.close()
