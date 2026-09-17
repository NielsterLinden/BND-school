"""The acceptance A = fiducial (dressed) / LHE mu mu and its uncertainties (frozen result, docs/16 §4).

A(60-120) = sumw(fid, dressed) / sumw(LHE mu mu, 60 < m_LHE < 120) turns sigma_fid into sigma(60-120);
A(m>50) = 3 sumw(fid, dressed) / sumw turns it into sigma(m > 50). The uncertainties, all outside the fit:

* `A_pdf_rel`: NNPDF3.1 Hessian members on the ratio of full-sample sums (output/v2/gensums.json);
* `A_alphas_rel`: alpha_s +- 0.002 members, scaled to +- 0.0015;
* `A_scale_rel`: 7-point muR/muF envelope;
* `A_stat`: aMC@NLO statistics of the numerator (absolute, like `A`);
* `A_ptz_rel`: the pT(Z) spectrum reweighted to the measured pT(mumu), folded with A(pT(Z)) of step 8
  (the PYTHIA ISR weights, `A_ps_isr_rel`, give less and are not added on top);
* `A_generator_rel`: powheg vs aMC@NLO;
* `A_ps_fsr_rel`: PYTHIA FSR weights;
* `A_qedfsr_rel`: the QED FSR model, an estimate (option c of docs/16 §5): the bare/dressed acceptance
  difference times the data/MC difference of the fraction of SR events with a recovered FSR photon.

`A_rel_unc` is their sum in quadrature. Until 16 Sep 2026 the fit result carried only PDF, alpha_s,
scale and statistics, evaluated for the m > 50 GeV denominator (scripts/mc_acceptance.py) and applied
to both cross sections; docs/16 §4.2.
"""

from __future__ import annotations

import json

import numpy as np

from . import config, weights

MC_SAMPLES = ["DYmumu", "DYtautau", "DYee", "TTbar", "SingleTop", "WW", "WZ", "ZZ"]
THEORY_JSON = config.OUTPUT_DIR / "v2" / "cms_parity" / "theory_acceptance.json"   # scripts/v2_8_theory_acceptance.py
QEDFSR_OPTION = "c"


def theory(gens, den="lhe_mumu_60_120", num="fid_dressed"):
    """PDF, alpha_s, scale, MC statistics, PS and generator uncertainties of sumw(num)/sumw(den), in percent."""
    g, pw = gens["DY_NLO"], gens["DY_powheg"]
    A = g[f"sumw_{num}"] / g[f"sumw_{den}"]
    pdf = np.array(g[f"pdf_{num}"]) / np.array(g[f"pdf_{den}"])
    pdf = pdf / pdf[0]
    scale = np.array(g[f"scale_{num}"]) / np.array(g[f"scale_{den}"])
    scale = scale / scale[weights.SCALE_NOMINAL]
    ps = (np.array(g[f"ps_{num}"]) / np.array(g[f"ps_{den}"])) / A
    A_pw = pw[f"sumw_{num}"] / pw[f"sumw_{den}"]
    stat_rel = np.sqrt(g[f"sumw2_{num}"]) / g[f"sumw_{num}"]        # dominated by the numerator
    return {"A": A, "denominator": den,
            "pdf_pct": 100 * float(np.sqrt(np.sum((pdf[weights.PDF_MEMBERS] - 1) ** 2))),
            "alphas_pct": 100 * float(0.75 * 0.5 * abs(pdf[102] - pdf[101])),
            "scale_pct": 100 * float(np.max(np.abs(scale[weights.SCALE_7POINT] - 1))),
            "mc_stat_pct": 100 * float(stat_rel),
            "ps_isr_pct": 100 * float(0.5 * (abs(ps[weights.PS_ISR_UP] - 1) + abs(ps[weights.PS_ISR_DOWN] - 1))),
            "ps_fsr_pct": 100 * float(0.5 * (abs(ps[weights.PS_FSR_UP] - 1) + abs(ps[weights.PS_FSR_DOWN] - 1))),
            "generator_pct": 100 * float(abs(A_pw / A - 1)), "A_powheg": A_pw,
            "source": "output/v2/gensums.json (all 71.8 M aMC@NLO events, 2.96 M powheg events)"}


def reweight_acceptance(A_bin, den, shape):
    """Acceptance after multiplying the boson-pT spectrum by `shape` (per bin), total rate unchanged."""
    return float((A_bin * den * shape).sum() / (den * shape).sum())


def ptz_reweighting(hists_all, theory_acc):
    """Data/MC ratio of the reconstructed pT(mumu) in the SR (normalised), mapped onto the generator pT(Z) bins."""
    data = np.asarray(hists_all["Data|SR|zpt|nominal"], dtype=float)
    mc = sum(np.asarray(hists_all[f"{s}|SR|zpt|nominal"], dtype=float) for s in MC_SAMPLES
             if f"{s}|SR|zpt|nominal" in hists_all)
    pw = np.asarray(hists_all["DYmumu_powheg|SR|zpt|nominal"], dtype=float)
    sig = np.asarray(hists_all["DYmumu|SR|zpt|nominal"], dtype=float)
    reco_edges = np.linspace(0, 200, len(data) + 1)
    ratio = (data / data.sum()) / (mc / mc.sum())
    ratio_pw = (pw / pw.sum()) / (sig / sig.sum())
    e = np.array(theory_acc["ptz_edges"])
    c = 0.5 * (e[1:] + e[:-1])
    idx = np.clip(np.digitize(c, reco_edges) - 1, 0, len(ratio) - 1)
    w = np.where(c < 200, ratio[idx], 1.0)
    w_pw = np.where(c < 200, ratio_pw[idx], 1.0)
    return {"reco_edges": reco_edges.tolist(), "ratio_data_mc": ratio.tolist(), "ratio_powheg_amcnlo": ratio_pw.tolist(),
            "gen_weights_data": w.tolist(), "gen_weights_powheg": w_pw.tolist()}


def ptz_shifts(theory_acc, ptz):
    """Relative change of A when the generator pT(Z) spectrum is reweighted (to data, to powheg; ours and 25/25 GeV)."""
    den = np.array(theory_acc["den_ptz"])
    out = {}
    for name, key, vol in (("data", "gen_weights_data", "dressed"), ("powheg", "gen_weights_powheg", "dressed"),
                           ("data_2525", "gen_weights_data", "dressed2525")):
        A_bin = np.array(theory_acc["volumes"][vol]["A_ptz"])
        out[name] = reweight_acceptance(A_bin, den, np.array(ptz[key])) / theory_acc["volumes"][vol]["A"] - 1
    return out


def fsr_photons(hists_all):
    """Fraction of SR events with at least one recovered FSR photon, data vs simulation."""
    d = np.asarray(hists_all["Data|SR|nfsr|nominal"], dtype=float)
    m = sum(np.asarray(hists_all[f"{s}|SR|nfsr|nominal"], dtype=float) for s in MC_SAMPLES if f"{s}|SR|nfsr|nominal" in hists_all)
    fd, fm = d[1:].sum() / d.sum(), m[1:].sum() / m.sum()
    return {"frac_with_fsr_photon_data": fd, "frac_with_fsr_photon_mc": fm, "ratio": fd / fm,
            "n_data": d.tolist(), "n_mc": m.tolist()}


def qedfsr_options(theory_acc, fsr, published_pct):
    """The three options of the QED FSR estimate, in percent: a published, b zero, c (size) x (data/MC)."""
    size = abs(theory_acc["fsr"]["A_bare_over_dressed"] - 1) * 100
    r = abs(fsr["ratio"] - 1)
    return {"a": published_pct, "b": 0.0, "c": float(r * size), "size_pct": size, "data_mc_rel": r}


def _block(t, A, ptz_pct, qedfsr_pct):
    rel = {"A_pdf_rel": t["pdf_pct"] / 100, "A_alphas_rel": t["alphas_pct"] / 100, "A_scale_rel": t["scale_pct"] / 100,
           "A_ptz_rel": ptz_pct / 100, "A_generator_rel": t["generator_pct"] / 100, "A_ps_fsr_rel": t["ps_fsr_pct"] / 100,
           "A_qedfsr_rel": qedfsr_pct / 100}
    stat_rel = t["mc_stat_pct"] / 100
    measured = float(np.sqrt(stat_rel ** 2 + sum(v ** 2 for k, v in rel.items() if k != "A_qedfsr_rel")))
    return {"A": float(A), "A_stat": float(A * stat_rel), **{k: float(v) for k, v in rel.items()},
            "A_ps_isr_rel": t["ps_isr_pct"] / 100, "A_rel_unc_measured": measured,
            "A_rel_unc": float(np.hypot(measured, rel["A_qedfsr_rel"]))}


def frozen(gens, hists_all, theory_acc=None, qedfsr_published_pct=0.12):
    """The acceptance block of the fit result and of zmumu.root.meta.json, plus the inputs it was built from.

    The 60-120 GeV block is at the top level (the combination reads `acceptance.A_*`); the m > 50 GeV
    one is under `m50`. PDF, alpha_s, scale, statistics and PS FSR are evaluated for each
    denominator (m > 50: LHE mu mu with m > 50 GeV); the pT(Z), generator and QED FSR rows act on the
    numerator and are taken from the 60-120 GeV evaluation for both."""
    theory_acc = theory_acc or json.loads(THEORY_JSON.read_text())
    g = gens["DY_NLO"]
    t60 = theory(gens, "lhe_mumu_60_120")
    # the powheg sample is generated in 50 < m_LHE < 120 GeV, so its m > 50 denominator is not comparable:
    # powheg vs aMC@NLO is a numerator effect, taken from the 60-120 GeV evaluation for both
    t50 = {**theory(gens, "lhe_mumu"), "generator_pct": t60["generator_pct"], "A_powheg": None}
    ptz = ptz_reweighting(hists_all, theory_acc)
    dA = ptz_shifts(theory_acc, ptz)
    fsr = fsr_photons(hists_all)
    q = qedfsr_options(theory_acc, fsr, qedfsr_published_pct)
    ptz_pct = max(abs(dA["data"]) * 100, t60["ps_isr_pct"])
    block = _block(t60, g["sumw_fid_dressed"] / g["sumw_lhe_mumu_60_120"], ptz_pct, q[QEDFSR_OPTION])
    block["m50"] = _block(t50, 3.0 * g["sumw_fid_dressed"] / g["sumw"], ptz_pct, q[QEDFSR_OPTION])
    block.update(denominator="LHE mu mu, 60 < m_LHE < 120 GeV (m50: m_LHE > 50 GeV, A = 3 sumw_fid / sumw)",
                 qedfsr_option=QEDFSR_OPTION, qedfsr_options_pct={k: q[k] for k in "abc"},
                 ptz_shift_data_rel=dA["data"],
                 source="zmumu/acceptance.py: output/v2/gensums.json, output/v2/histograms.pkl (pT(mumu), FSR photons), "
                        "output/v2/cms_parity/theory_acceptance.json (A(pT(Z)), bare/dressed); docs/16 §4-5")
    details = {"theory_60_120": t60, "theory_m50": t50, "ptz_reweighting": {**ptz, "dA_rel": dA}, "fsr_photons": fsr,
               "qedfsr": q}
    return block, details
