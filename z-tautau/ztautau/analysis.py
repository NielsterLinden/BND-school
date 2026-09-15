"""Reading the ntuples, region definitions, simulation weights, systematic variations and BDT categories.

Everything downstream of step 2 goes through this module, so the regions and weights are defined once.

Regions (docs/05-fake-factors.md); T = DeepTau VSjet Medium, L = VVVLoose and not Medium,
both taus pT > 40 GeV (after any energy-scale variation):

    SR      OS  tau1 T  tau2 T      signal region (the fit, split into BDT categories: `categories`)
    AR      OS  tau1 L  tau2 T      application region: data x FF -> jet -> tau_h fakes in the SR
    SS_T    SS  tau1 T  tau2 T      fake-factor numerator
    SS_L    SS  tau1 L  tau2 T      fake-factor denominator
    OSAI_T  OS  tau1 T  tau2 L      tau2 anti-isolated sideband, OS  } OS/SS extrapolation check of the FF
    OSAI_L  OS  tau1 L  tau2 L                                       }
    SSAI_T  SS  tau1 T  tau2 L      tau2 anti-isolated sideband, SS  } FF measured here, applied to OSAI_L
    SSAI_L  SS  tau1 L  tau2 L                                       }
    SR_relaxed OS tau1 T tau2 T or L   shape template for W+jets (few simulated events, large weights)

Simulation: only events whose leading tau is *not* a jet (genPartFlav != 0) are kept -- jet fakes of
the leading tau are the fake-factor estimate. Events with a genuine leading tau and a jet as the
subleading tau (W+jets, ttbar) stay in the simulation.

Drell-Yan: the inclusive aMC@NLO sample carries the normalisation and the acceptance; the jet-binned
0J/1J/2J samples, when their ntuples exist, add statistics through `dy_norm` (stitching per LHE_NpNLO
bin -- the parton multiplicity of the NLO matrix element, in which the samples are exclusive -- with the
bin fractions of the inclusive sample, so no extra cross sections enter). The LHE tautau
part is split into the fiducial signal (`DYtautau`, scaled by mu_Z) and the non-fiducial rest
(`DYtautau_nonfid`: m_LHE > 120 GeV mostly, a theory-normalised background; docs/06-cross-section.md).
"""

from __future__ import annotations

import json
from functools import lru_cache

import numpy as np
import uproot

from . import bdt, config, corrections, mass, samples

T_BIT, L_BIT = config.TAU_VSJET_TIGHT_BIT, config.TAU_VSJET_LOOSE_BIT
REGIONS = ["SR", "AR", "SS_T", "SS_L", "OSAI_T", "OSAI_L", "SSAI_T", "SSAI_L"]
# samples whose SR shape is taken from SR_relaxed (tau2 VVVLoose) to tame huge per-event weights
SMOOTHED_SAMPLES = ["WJets"]
# ... and which are subtracted from the fake-factor regions with *uniform* weights (every event carries the
# sample's mean weight, so the total normalisation is kept): the raw weights are 20-200 with 16% negative,
# and a single event with weight -63 in the same-sign top BDT bin faked a 75% non-closure. W+jets is a real
# genuine-tau_1 component of the tau2-anti-isolated sideband (2% of OSAI_T, +1.5% on C_OS/SS), so it must
# not simply be dropped. `subtraction_weights` implements this; SUBTRACT_EXCLUDE drops a sample entirely.
SUBTRACT_EXCLUDE: list[str] = []

# systematic variations of the simulation (name -> kind); see weights() / kinematics()
DMS = config.TAU_DMS
WEIGHT_SYSTS = (["Pileup", "L1Prefiring", "TauFakeEle", "TauFakeMu"]
                + [f"TauID_DM{d}" for d in DMS] + [f"TauTrigger_DM{d}" for d in DMS])
KINEMATIC_SYSTS = [f"TauES_DM{d}" for d in DMS] + ["MET_Unclustered"]
THEORY_SYSTS = ["QCDScale", "PDF", "PS_ISR", "PS_FSR"]
KIN_KEYS = ("t1_pt", "t2_pt", "m_vis", "m_tt", "m_col", "met", "pt_tt", "mt_tot", "met_x", "met_y")


# ------------------------------------------------------------------------------ loading
@lru_cache(maxsize=None)
def _load(key: str):
    path = config.NTUPLE_DIR / f"{key}.root"
    meta = json.loads((config.NTUPLE_DIR / f"{key}.meta.json").read_text())
    with uproot.open(path) as f:
        arrays = f["ntuple"].arrays(library="np") if "ntuple" in [k.split(";")[0] for k in f.keys()] else {}
    return arrays, meta


def has_ntuple(key: str) -> bool:
    return (config.NTUPLE_DIR / f"{key}.meta.json").exists()


def available_mc(keys=None) -> list[str]:
    """MC keys whose ntuple exists. Missing non-optional samples are reported (the result needs them);
    jet-binned DY samples are only used when the inclusive sample carries the LHE_NpNLO sums."""
    keys = keys or samples.NOMINAL_MC_KEYS
    ok = [k for k in keys if has_ntuple(k)]
    missing = sorted(k for k in set(keys) - set(ok) if not samples.SAMPLES[k]["optional"])
    if missing:
        print(f"WARNING: no ntuple for {missing} -- the prediction is incomplete")
    if "sumw_npnlo" not in load(samples.DY_INCLUSIVE)[1]["gensums"]:
        dropped = [k for k in ok if k in samples.DY_STITCHED and k != samples.DY_INCLUSIVE]
        if dropped:
            print(f"WARNING: inclusive DY skim has no LHE_NpNLO sums, not using {dropped}")
        ok = [k for k in ok if k not in dropped]
    return ok


def load(key: str):
    """(arrays, meta) of one ntuple; arrays are shared (do not modify in place)."""
    return _load(key)


def load_data():
    """Both eras concatenated."""
    parts = [load(k)[0] for k in samples.DATA_KEYS]
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def normalisation(key: str) -> float:
    """sigma x L / sum of generator weights of the processed files."""
    s = samples.SAMPLES[key]
    meta = load(key)[1]
    return s["xsec_pb"] * config.LUMI_PB / meta["gensums"]["sumw"]


def dy_stitch_keys() -> list[str]:
    g = load(samples.DY_INCLUSIVE)[1]["gensums"]
    if "sumw_npnlo" not in g:
        return [samples.DY_INCLUSIVE]
    return [k for k in samples.DY_STITCHED if has_ntuple(k) and "sumw_npnlo" in load(k)[1]["gensums"]]


def dy_norm(key: str, d) -> np.ndarray:
    """Per-event normalisation of a Drell-Yan signal-group sample, stitched per LHE_NpNLO bin:

        w_j = sigma x L x f_j / sum_samples sumw_sample(j),   f_j = sumw_incl(j) / sumw_incl

    With the inclusive sample alone this is sigma x L / sumw_incl; every jet-binned sample only enlarges the
    denominator of its bin. The inclusive sample fixes the jet-bin fractions (its FxFx merging), the
    jet-binned samples add statistics.
    """
    keys = dy_stitch_keys()
    if key not in keys or len(keys) < 2 or "lhe_npnlo" not in d:
        return np.full(len(d["run"]), normalisation(key))
    incl = load(samples.DY_INCLUSIVE)[1]["gensums"]
    frac = np.asarray(incl["sumw_npnlo"], dtype=float) / incl["sumw"]
    denom = sum(np.asarray(load(k)[1]["gensums"]["sumw_npnlo"], dtype=float) for k in keys)
    nj = np.clip(d["lhe_npnlo"].astype(int), 0, len(frac) - 1)
    return samples.SAMPLES[key]["xsec_pb"] * config.LUMI_PB * (frac / denom)[nj]


def mc_components(key: str):
    """[(fit sample name, event mask)] of an MC ntuple (LHE flavour split for Drell-Yan, and the tautau part
    split into the fiducial signal and the non-fiducial background)."""
    d, _ = load(key)
    s = samples.SAMPLES[key]
    if not len(d):
        return []
    if s["split_lhe"]:
        out = []
        for code, name in samples.LHE_SPLIT.items():
            m = d["gen_lhe_flavour"] == code
            if name == samples.SIGNAL:
                fid = d["gen_fid"].astype(bool)
                out += [(samples.SIGNAL, m & fid), (samples.SIGNAL_NONFID, m & ~fid)]
            else:
                out.append((name, m))
        return out
    return [(s["fit_sample"], np.ones(len(d["run"]), dtype=bool))]


# ------------------------------------------------------------------------------ kinematics
def _masses_for(d, idx, pt1, m1, pt2, m2, metx, mety):
    t1 = (pt1[idx], d["t1_eta"][idx].astype(np.float64), d["t1_phi"][idx].astype(np.float64), m1[idx])
    t2 = (pt2[idx], d["t2_eta"][idx].astype(np.float64), d["t2_phi"][idx].astype(np.float64), m2[idx])
    return mass.all_masses(t1, t2, metx[idx], mety[idx], d["met_covxx"][idx].astype(np.float64),
                           d["met_covxy"][idx].astype(np.float64), d["met_covyy"][idx].astype(np.float64))


_KIN_CACHE: dict = {}


def kinematics(d, key: str, variation: str | None = None, direction: str | None = None):
    """dict(t1_pt, t2_pt, m_vis, m_tt, m_col, met, pt_tt, mt_tot, met_x, met_y) under a kinematic variation.

    TauES_DM<d>: genuine tau_h of that decay mode scaled by (1 +- sigma_TES), MET adjusted, masses
    recomputed for the affected events. MET_Unclustered: MET +- the NanoAOD unclustered-energy shift.
    """
    base = {k: d[k] for k in KIN_KEYS}
    if variation is None:
        return base
    ck = (key, variation, direction, len(d["run"]))
    if ck in _KIN_CACHE:
        return _KIN_CACHE[ck]
    sign = 1.0 if direction == "Up" else -1.0
    pt1, pt2 = d["t1_pt"].astype(np.float64), d["t2_pt"].astype(np.float64)
    m1, m2 = d["t1_mass"].astype(np.float64), d["t2_mass"].astype(np.float64)
    metx, mety = d["met_x"].astype(np.float64), d["met_y"].astype(np.float64)
    if variation.startswith("TauES_DM"):
        dm = int(variation[len("TauES_DM"):])
        rel = corrections.tes_uncertainty()[dm] / corrections.tes_nominal()[dm]
        aff1 = (d["t1_genflav"] == 5) & (d["t1_dm"] == dm)
        aff2 = (d["t2_genflav"] == 5) & (d["t2_dm"] == dm)
        new1 = np.where(aff1, pt1 * (1 + sign * rel), pt1)
        new2 = np.where(aff2, pt2 * (1 + sign * rel), pt2)
        metx = metx - (new1 - pt1) * np.cos(d["t1_phi"]) - (new2 - pt2) * np.cos(d["t2_phi"])
        mety = mety - (new1 - pt1) * np.sin(d["t1_phi"]) - (new2 - pt2) * np.sin(d["t2_phi"])
        m1 = np.where(aff1, m1 * (1 + sign * rel), m1)
        m2 = np.where(aff2, m2 * (1 + sign * rel), m2)
        pt1, pt2 = new1, new2
        idx = np.where(aff1 | aff2)[0]
    elif variation == "MET_Unclustered":
        metx = metx + sign * d["met_uncl_dx"]
        mety = mety + sign * d["met_uncl_dy"]
        idx = np.arange(len(pt1))
    else:
        raise ValueError(variation)
    out = {k: v.copy() for k, v in base.items()}
    out["t1_pt"], out["t2_pt"] = pt1.astype(np.float32), pt2.astype(np.float32)
    out["met_x"], out["met_y"] = metx.astype(np.float32), mety.astype(np.float32)
    out["met"] = np.hypot(metx, mety).astype(np.float32)
    # only events that can enter a region need the (slow) likelihood mass
    idx = idx[(pt1[idx] > config.TAU_PT_MIN) & (pt2[idx] > config.TAU_PT_MIN)]
    if len(idx):
        ms = _masses_for(d, idx, pt1, m1, pt2, m2, metx, mety)
        for k in ("m_vis", "m_tt", "m_col", "met", "pt_tt", "mt_tot"):
            out[k][idx] = ms[k]
    _KIN_CACHE[ck] = out
    return out


# ------------------------------------------------------------------------------ regions and categories
def regions(d, kin=None, is_mc: bool = False) -> dict:
    kin = kin or d
    base = (kin["t1_pt"] > config.TAU_PT_MIN) & (kin["t2_pt"] > config.TAU_PT_MIN)
    if is_mc:
        base &= d["t1_genflav"] != 0
    t1 = (d["t1_vsjet"] & T_BIT) > 0
    t2 = (d["t2_vsjet"] & T_BIT) > 0
    l1 = ((d["t1_vsjet"] & L_BIT) > 0) & ~t1
    l2 = ((d["t2_vsjet"] & L_BIT) > 0) & ~t2
    os_ = d["os"].astype(bool)
    ss = ~os_
    return {"SR": base & os_ & t1 & t2, "AR": base & os_ & l1 & t2,
            "SS_T": base & ss & t1 & t2, "SS_L": base & ss & l1 & t2,
            "OSAI_T": base & os_ & t1 & l2, "OSAI_L": base & os_ & l1 & l2,
            "SSAI_T": base & ss & t1 & l2, "SSAI_L": base & ss & l1 & l2,
            "SR_relaxed": base & os_ & t1 & (t2 | l2)}


def scores(d, key: str, kin=None, variation: str | None = None, direction: str | None = None) -> np.ndarray:
    """Held-out-fold BDT score of every event (ztautau/bdt.py); recomputed under kinematic variations."""
    return bdt.score(d, key, kin, variation, direction)


def categories(d, key: str, kin=None, variation: str | None = None, direction: str | None = None) -> np.ndarray:
    """BDT category index (0 .. n_cat-1) of every event."""
    return bdt.category(scores(d, key, kin, variation, direction))


# ------------------------------------------------------------------------------ weights
@lru_cache(maxsize=1)
def pileup():
    """Pileup weights from the Drell-Yan generator profile (all UL16 samples share the scenario)."""
    prof = np.asarray(load(samples.DY_INCLUSIVE)[1]["gensums"]["pu_true"])
    return corrections.PileupWeights(prof)


def weights(d, key: str, syst: str | None = None, direction: str | None = None, kin=None):
    """Event weights of an MC ntuple: norm x genWeight x pileup x prefiring x trigger SFs x ID SFs."""
    kin = kin or d
    sign = {"Up": 1, "Down": -1, None: 0}[direction]
    pu_var = {"Up": "up", "Down": "down"}.get(direction) if syst == "Pileup" else "nominal"
    norm = dy_norm(key, d) if key in samples.DY_STITCHED else normalisation(key)
    w = norm * d["genWeight"].astype(np.float64) * pileup()(d["Pileup_nTrueInt"], pu_var)
    pref = {"Up": "L1PreFiringWeight_Up", "Down": "L1PreFiringWeight_Dn"}.get(direction) if syst == "L1Prefiring" else None
    w = w * d[pref or "L1PreFiringWeight_Nom"]
    for i in (1, 2):
        dm = d[f"t{i}_dm"].astype(int)
        flav = d[f"t{i}_genflav"].astype(int)
        abseta = np.abs(d[f"t{i}_eta"])
        trig_dm = int(syst[len("TauTrigger_DM"):]) if syst and syst.startswith("TauTrigger_DM") else None
        w = w * corrections.trigger_sf(dm, kin[f"t{i}_pt"], trig_dm, sign if trig_dm is not None else 0, genflav=flav)
        if syst and syst.startswith("TauID_DM"):
            w = w * corrections.id_sf(dm, flav, abseta, int(syst[len("TauID_DM"):]), sign, "id")
        elif syst == "TauFakeEle":
            w = w * corrections.id_sf(dm, flav, abseta, None, sign, "vse")
        elif syst == "TauFakeMu":
            w = w * corrections.id_sf(dm, flav, abseta, None, sign, "vsmu")
        else:
            w = w * corrections.id_sf(dm, flav, abseta)
    return w


def subtraction_weights(key: str, w):
    """Weights to use when subtracting a simulated sample from the fake-factor regions: the event weights,
    except for SMOOTHED_SAMPLES (W+jets), whose pathological per-event weights are replaced by their mean."""
    if samples.SAMPLES[key]["fit_sample"] in SMOOTHED_SAMPLES:
        return np.full(len(w), float(np.mean(w)))
    return np.asarray(w)


def theory_weights(d, key: str, syst: str):
    """Per-event weight multipliers of the members of a theory variation, shape (n_events, n_members),
    each member renormalised so the *fiducial* generator-level yield of the sample is unchanged: on the
    fiducial signal only the C factor varies (fitting/CONVENTIONS.md section 3); on the non-fiducial
    part the same multipliers vary the non-fiducial / fiducial ratio, i.e. the theory uncertainty of that
    background. None if the sample lacks the weights.

    How the members combine is a property of the *histograms*, never of single events (an event-wise
    envelope or quadrature sum is badly biased): see `combine_theory`.
      QCDScale  6 members (muR, muF) in {0.5, 1, 2}^2 without (0.5, 2), (2, 0.5)  -> envelope per bin
      PDF       100 NNPDF3.1 Hessian eigenvectors + 2 alpha_s members           -> quadrature per bin
      PS_ISR    (ISR x2, ISR x0.5);  PS_FSR (FSR x2, FSR x0.5)                     -> up / down
    """
    g = load(key)[1]["gensums"]
    fid0 = g["sumw_fid"]
    branch, sums_key, idx = {"QCDScale": ("LHEScaleWeight", "scale_fid", [0, 1, 3, 5, 7, 8]),
                             "PDF": ("LHEPdfWeight", "pdf_fid", list(range(103))),
                             "PS_ISR": ("PSWeight", "ps_fid", [0, 2]),
                             "PS_FSR": ("PSWeight", "ps_fid", [1, 3])}[syst]
    if branch not in d or sums_key not in g:
        return None
    v = d[branch].astype(np.float64)[:, idx]
    sums = np.asarray(g[sums_key])[idx]
    return v * (fid0 / sums)[None, :]


def combine_theory(syst: str, nominal, members):
    """(up, down) histogram values from the nominal (nbins) and the member histograms (n_members x nbins)."""
    members = np.asarray(members)
    if syst == "QCDScale":
        return np.maximum(members.max(axis=0), nominal), np.minimum(members.min(axis=0), nominal)
    if syst == "PDF":
        # member 0 is the central PDF, 1-100 the eigenvectors, 101/102 alpha_s down/up
        delta = np.sqrt(((members[1:101] - members[0]) ** 2).sum(axis=0) + (0.5 * (members[102] - members[101])) ** 2)
        return nominal + delta, nominal - delta
    return members[0], members[1]


# ------------------------------------------------------------------------------ generator-level numbers
def signal_prediction(key: str = samples.DY_INCLUSIVE) -> dict:
    """Predicted cross sections of the Drell-Yan sample (pb): total tautau, 60-120 GeV, fiducial; A."""
    g = load(key)[1]["gensums"]
    xs = samples.SAMPLES[key]["xsec_pb"]
    out = {"sigma_tautau_pb": xs * g["sumw_lhe_tautau"] / g["sumw"],
           "sigma_tautau_60_120_pb": xs * g["sumw_lhe_tautau_60_120"] / g["sumw"],
           "sigma_fid_pb": xs * g["sumw_fid"] / g["sumw"],
           "A": g["sumw_fid"] / g["sumw_lhe_tautau_60_120"],
           "A_mc_stat": np.sqrt(g["sumw2_fid"]) / g["sumw_fid"]}
    unc = {}
    for name, key_ in (("scale", "scale"), ("pdf", "pdf"), ("ps", "ps")):
        num = np.asarray(g.get(f"{key_}_fid", []))
        den = np.asarray(g.get(f"{key_}_lhe_tautau_60_120", []))
        if num.size and den.size:
            a = num / den
            a0 = out["A"]
            if name == "scale":
                idx = [0, 1, 3, 5, 7, 8]
                unc["A_scale"] = float(np.max(np.abs(a[idx] - a0)) / a0)
            elif name == "pdf":
                unc["A_pdf"] = float(np.sqrt(np.sum((a[1:101] - a[0]) ** 2)) / a0)
                unc["A_alphas"] = float(np.abs(a[102] - a[101]) / 2 / a0)
            else:
                unc["A_isr"] = float(max(abs(a[0] - a0), abs(a[2] - a0)) / a0)
                unc["A_fsr"] = float(max(abs(a[1] - a0), abs(a[3] - a0)) / a0)
    out["A_unc"] = unc
    return out
