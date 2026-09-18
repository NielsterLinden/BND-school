"""k-fold BDT separating Z -> tau_h tau_h from jet -> tau_h fakes (docs/09-bdt.md).

Why: with DeepTau Medium on both legs the signal region is 80% fakes. The fakes are estimated from the
application region (AR) of the data, so every fake-related uncertainty scales with B/S. A classifier on
mass-agnostic kinematics (dR, the pT balance of the tau tau + MET system, decay modes, ...) separates the two
at AUC ~ 0.89 and, used as a *category* variable with m_tt still fitted in each category, isolates a third
of the signal at S/B ~ 5 while the fake-dominated category fixes the fake normalisation and shape.

Rules that keep the fake-factor method valid (review/REVIEW_v3.md 4.3):
  * inputs never include the tau_1 isolation (raw DeepTau, isolation sums, leadTkPtOverTauPt) or the
    charge: the classifier must not learn the difference between the AR and the SR, nor OS vs SS;
  * k folds by event number: an event is always scored by the model that never saw it; this holds for
    AR data (the fake template), same-sign data (closure), signal-region data and every simulated sample;
  * the same-sign closure of the FF *in the score* is checked (scripts/step3b_bdt.py) before the score is
    used, and the residual non-closure per category enters the fit (FakeClosure_tautau_c<k>_*).

Training: signal = simulated Z -> tautau in the SR (fiducial part, positive weights), background = AR data
weighted by the nominal fake factor (the fakes as the analysis models them). Both classes carry the same
total weight. Models: XGBoost, config.BDT_PARAMS, saved to config.BDT_DIR/fold<k>.ubj (+ bdt.json).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import config

K = config.BDT_K
FEATURES = list(config.BDT_FEATURES)
_MODELS: list | None = None
_SCORES: dict = {}


def _dphi(a, b):
    return np.abs((a - b + np.pi) % (2 * np.pi) - np.pi)


def features(d, kin=None) -> np.ndarray:
    """(n_events, n_features) float32 matrix; `kin` (analysis.kinematics) overrides the pT / MET columns."""
    kin = kin or d
    pt1, pt2 = kin["t1_pt"].astype(np.float64), kin["t2_pt"].astype(np.float64)
    metx, mety = kin["met_x"].astype(np.float64), kin["met_y"].astype(np.float64)
    phi1, phi2 = d["t1_phi"].astype(np.float64), d["t2_phi"].astype(np.float64)
    vx = pt1 * np.cos(phi1) + pt2 * np.cos(phi2)
    vy = pt1 * np.sin(phi1) + pt2 * np.sin(phi2)
    met = np.hypot(metx, mety)
    cols = {
        "t1_pt": pt1, "t2_pt": pt2, "pt_ratio": pt2 / np.maximum(pt1, 1e-6),
        "t1_abseta": np.abs(d["t1_eta"]), "t2_abseta": np.abs(d["t2_eta"]),
        "dr_tt": d["dr_tt"], "dphi_tt": _dphi(phi1, phi2),
        "met": met, "met_sig": met / np.sqrt(np.maximum(d["met_covxx"] + d["met_covyy"], 1.0)),
        "pt_vis": np.hypot(vx, vy), "dphi_met_tt": _dphi(np.arctan2(mety, metx), np.arctan2(vy, vx)),
        "pt_tt": np.hypot(vx + metx, vy + mety),
        "njets": d["njets"], "jet1_pt": d["jet1_pt"], "t1_dm": d["t1_dm"], "t2_dm": d["t2_dm"],
    }
    return np.stack([np.asarray(cols[f], dtype=np.float32) for f in FEATURES], axis=1)


def folds(d) -> np.ndarray:
    return (d["event"].astype(np.uint64) % K).astype(int)


# ------------------------------------------------------------------------------ training
def train(X_sig, w_sig, f_sig, X_bkg, w_bkg, f_bkg, out_dir: Path | None = None) -> dict:
    """Train the K fold models; returns {'auc': [...], 'importance': {...}, ...} and writes the models."""
    import xgboost as xgb
    from sklearn.metrics import roc_auc_score

    out_dir = Path(out_dir or config.BDT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    w_sig, w_bkg = np.maximum(w_sig, 0.0), np.maximum(w_bkg, 0.0)
    info = {"k": K, "features": FEATURES, "params": config.BDT_PARAMS, "auc_test": [], "auc_train": [],
            "n_sig": int(len(X_sig)), "n_bkg": int(len(X_bkg)), "sumw_sig": float(w_sig.sum()), "sumw_bkg": float(w_bkg.sum())}
    imp = np.zeros(len(FEATURES))
    models = []
    for k in range(K):
        tr_s, tr_b = f_sig != k, f_bkg != k
        n_ref = float(tr_s.sum())          # both classes: the same total weight, mean weight ~1 for the signal
        X = np.concatenate([X_sig[tr_s], X_bkg[tr_b]])
        y = np.concatenate([np.ones(tr_s.sum()), np.zeros(tr_b.sum())])
        w = np.concatenate([w_sig[tr_s] / w_sig[tr_s].sum() * n_ref, w_bkg[tr_b] / w_bkg[tr_b].sum() * n_ref])
        clf = xgb.XGBClassifier(**config.BDT_PARAMS).fit(X, y, sample_weight=w)
        te_s, te_b = ~tr_s, ~tr_b
        Xt = np.concatenate([X_sig[te_s], X_bkg[te_b]])
        yt = np.concatenate([np.ones(te_s.sum()), np.zeros(te_b.sum())])
        wt = np.concatenate([w_sig[te_s] / max(w_sig[te_s].sum(), 1e-9), w_bkg[te_b] / max(w_bkg[te_b].sum(), 1e-9)])
        info["auc_test"].append(float(roc_auc_score(yt, clf.predict_proba(Xt)[:, 1], sample_weight=wt)))
        info["auc_train"].append(float(roc_auc_score(y, clf.predict_proba(X)[:, 1], sample_weight=w)))
        imp += clf.feature_importances_ / K
        clf.get_booster().save_model(str(out_dir / f"fold{k}.ubj"))     # Booster API: the sklearn wrapper's
        models.append(clf.get_booster())                                  # save_model breaks with sklearn 1.9
    info["importance"] = {f: float(v) for f, v in sorted(zip(FEATURES, imp), key=lambda t: -t[1])}
    info["category_edges"] = config.BDT_CATEGORY_EDGES
    (out_dir / "bdt.json").write_text(json.dumps(info, indent=1))
    global _MODELS
    _MODELS = models
    _SCORES.clear()
    return info


def load_models(model_dir: Path | None = None) -> list:
    global _MODELS
    if _MODELS is None:
        import xgboost as xgb
        model_dir = Path(model_dir or config.BDT_DIR)
        info = json.loads((model_dir / "bdt.json").read_text())
        if info["features"] != FEATURES or info["k"] != K:
            raise RuntimeError(f"BDT models in {model_dir} were trained with other settings: retrain (step3b)")
        _MODELS = []
        for k in range(K):
            booster = xgb.Booster()
            booster.load_model(str(model_dir / f"fold{k}.ubj"))
            _MODELS.append(booster)
    return _MODELS


def score_arrays(X: np.ndarray, f: np.ndarray) -> np.ndarray:
    import xgboost as xgb
    models = load_models()
    out = np.zeros(len(X), dtype=np.float32)
    for k in range(K):
        m = f == k
        if m.any():
            out[m] = models[k].predict(xgb.DMatrix(X[m]))
    return out


def score(d, key: str, kin=None, variation: str | None = None, direction: str | None = None) -> np.ndarray:
    """Held-out-fold score of every event of a loaded ntuple (cached per sample and kinematic variation)."""
    ck = (key, variation, direction, len(d["run"]))
    if ck not in _SCORES:
        _SCORES[ck] = score_arrays(features(d, kin), folds(d))
    return _SCORES[ck]


def category(scores) -> np.ndarray:
    """Category index 0..n-1 of every score (config.BDT_CATEGORY_EDGES)."""
    edges = np.asarray(config.BDT_CATEGORY_EDGES)
    return np.clip(np.searchsorted(edges, scores, side="right") - 1, 0, len(edges) - 2)
