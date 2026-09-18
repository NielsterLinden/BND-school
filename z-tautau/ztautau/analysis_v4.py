"""v4: ntuples, regions, weights and systematic variations of the lepton channels, and the v4 view of the
tau_h tau_h ntuples (docs/10-v4-plan.md sections 3, 5, 7, 8).

Lepton-channel regions (T = tau_h VSjet at the nominal working point, L = VVVLoose and not T; every region
requires the trigger, the pair, no additional lepton, the lepton on its plateau and the tau_h above 30 GeV):

  mu tau_h / e tau_h                        lepton iso   charge   tau_h   m_T(l, MET)
    SR        signal region                  < 0.15       OS       T       < 40
    AR        fake application region        < 0.15       OS       L       < 40
    QCD_T/L   multijet FF determination      < 0.15       SS       T/L     < 40
    W_T/L     W+jets FF determination        < 0.15       OS       T/L     > 70, no b jet
    WSS_T/L   the same in same-sign events (W FF of same-sign events, for the same-sign validation only:
              the jet recoiling against a W is charge-correlated with it, so same-sign W fakes are gluon-like
              and have a lower FF than the opposite-sign ones)
    AI_OS_T/L, AI_SS_T/L  anti-isolated lepton sidebands (0.15 < iso < 0.5) for the OS/SS extrapolation
    SS_SR/SS_AR  same-sign validation of the fake estimate
  e mu (both leptons iso < 0.15, no additional lepton)
    SR        OS, D_zeta > -20, no b jet;  CRtt   OS, D_zeta < -40, MET > 80;  SS  like SR with SS
    SB1_OS/SS both iso < 0.5 and at least one > 0.15;  SB2_OS/SS  both < 0.5, at least one > 0.3 (skim limit 0.5)

Simulation in the lepton channels keeps only events whose tau_h is not a jet (genPartFlav != 0): jet fakes
are the fake-factor estimate. Nothing here applies the DeepTau VSjet scale factor (free in the fit).
"""

from __future__ import annotations

import json
from functools import lru_cache

import numpy as np
import uproot

from . import analysis, config, corrections, mass, pog, samples

LTAU = ("mutau", "etau")
T_BIT, L_BIT = config.TAU_VSJET_TIGHT_BIT, config.TAU_VSJET_LOOSE_BIT
DMS = config.TAU_DMS
LEPTON_PREFIX = {"mutau": "l", "etau": "l", "emu": None}
# systematic variations of the simulation, per channel family
WEIGHT_SYSTS_COMMON = ["Pileup", "L1Prefiring", "TauFakeEle", "TauFakeMu"]
WEIGHT_SYSTS = {
    "mutau": WEIGHT_SYSTS_COMMON + ["MuonID", "MuonIso", "MuonTrigger", "BTag", "TopPt"],
    "etau": WEIGHT_SYSTS_COMMON + ["ElectronReco", "ElectronID", "ElectronTrigger", "ElectronTrigger_lowpt", "BTag", "TopPt"],
    "emu": ["Pileup", "L1Prefiring", "MuonID", "MuonIso", "ElectronReco", "ElectronID", "EmuTrigger", "BTag", "TopPt"],
    "tautau": ["Pileup", "L1Prefiring", "TauFakeEle", "TauFakeMu"] + [f"TauTrigger_DM{d}" for d in DMS],
}
KINEMATIC_SYSTS = {
    "mutau": [f"TauES_DM{d}" for d in DMS] + ["MET_Unclustered", "JES", "MuonScale"],
    "etau": [f"TauES_DM{d}" for d in DMS] + ["MET_Unclustered", "JES", "ElectronScale"],
    "emu": ["MET_Unclustered", "JES", "MuonScale", "ElectronScale"],
    "tautau": [f"TauES_DM{d}" for d in DMS] + ["MET_Unclustered"],
}
THEORY_SYSTS = ["QCDScale", "PDF", "PS_ISR", "PS_FSR"]
TRIG_FLAT_UNC = 0.02            # CMS arXiv:1801.03535 Table 2: 2% per channel on the trigger efficiency
MUON_SCALE_REL, ELECTRON_SCALE_REL = 0.002, (0.005, 0.010)      # estimates (docs/10 section 7): mu; e barrel, endcap
TRIG_INSITU = config.EXTERNAL_DIR / "trigger_insitu_v4.json"
BTAG_EFF = config.DATA_DIR_V4 / "btag_eff.json"


# ------------------------------------------------------------------------------ loading
@lru_cache(maxsize=None)
def _load(key: str, channel: str):
    path = config.NTUPLE_DIR_V4 / f"{key}_{channel}.root"
    meta = json.loads((config.NTUPLE_DIR_V4 / f"{key}_{channel}.meta.json").read_text())
    with uproot.open(path) as f:
        arrays = f["ntuple"].arrays(library="np") if "ntuple" in [k.split(";")[0] for k in f.keys()] else {}
    return arrays, meta


def has_ntuple(key: str, channel: str) -> bool:
    return (config.NTUPLE_DIR_V4 / f"{key}_{channel}.meta.json").exists()


def load(key: str, channel: str):
    if channel == "tautau":
        return analysis.load(key)
    return _load(key, channel)


def load_data(channel: str):
    if channel == "tautau":
        return analysis.load_data()
    stream = {"mutau": "SingleMuon", "etau": "SingleElectron", "emu": "MuonEG", "emu_mu": "SingleMuon", "emu_el": "SingleElectron"}[channel]
    parts = [load(k, channel)[0] for k in samples.DATA_KEYS_V4[stream] if has_ntuple(k, channel)]
    parts = [p for p in parts if len(p)]
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def available_mc(channel: str) -> list[str]:
    if channel == "tautau":
        return analysis.available_mc()
    keys = [k for k in samples.MC_KEYS_LEPTON if has_ntuple(k, channel) and len(load(k, channel)[0])]
    missing = [k for k in samples.MC_KEYS_LEPTON if k not in keys and not samples.SAMPLES[k]["optional"]]
    if missing:
        print(f"WARNING [{channel}]: no ntuple for {missing} -- the prediction is incomplete")
    return keys


def normalisation(key: str, channel: str) -> float:
    s = samples.SAMPLES[key]
    return s["xsec_pb"] * config.LUMI_PB / load(key, channel)[1]["gensums"]["sumw"]


def dy_stitch_keys(channel: str) -> list[str]:
    if channel == "tautau":
        return analysis.dy_stitch_keys()
    return [k for k in samples.DY_STITCHED if has_ntuple(k, channel) and "sumw_npnlo" in load(k, channel)[1]["gensums"]]


def dy_norm(key: str, channel: str, d) -> np.ndarray:
    """Jet-binned stitching as analysis.dy_norm, per channel ntuple."""
    if channel == "tautau":
        return analysis.dy_norm(key, d)
    keys = dy_stitch_keys(channel)
    if key not in keys or len(keys) < 2 or "lhe_npnlo" not in d:
        return np.full(len(d["run"]), normalisation(key, channel))
    incl = load(samples.DY_INCLUSIVE, channel)[1]["gensums"]
    frac = np.asarray(incl["sumw_npnlo"], dtype=float) / incl["sumw"]
    denom = sum(np.asarray(load(k, channel)[1]["gensums"]["sumw_npnlo"], dtype=float) for k in keys)
    nj = np.clip(d["lhe_npnlo"].astype(int), 0, len(frac) - 1)
    return samples.SAMPLES[key]["xsec_pb"] * config.LUMI_PB * (frac / denom)[nj]


def mc_components(key: str, channel: str):
    """[(fit sample, mask)]: Drell-Yan split by LHE flavour, the tau tau part into the 60-120 GeV signal and
    the rest (`DYtautau_out`); everything else one component."""
    d, _ = load(key, channel)
    s = samples.SAMPLES[key]
    if not len(d):
        return []
    if s["split_lhe"]:
        out = []
        for code, name in samples.LHE_SPLIT.items():
            m = d["gen_lhe_flavour"] == code
            if name == samples.SIGNAL:
                win = (d["gen_mll_lhe"] > config.MASS_LO) & (d["gen_mll_lhe"] < config.MASS_HI)
                out += [(config.V4_SIGNAL, m & win), (config.V4_SIGNAL_OUT, m & ~win)]
            else:
                out.append((name, m))
        return out
    return [(s["fit_sample"], np.ones(len(d["run"]), dtype=bool))]


# ------------------------------------------------------------------------------ kinematics
def _leptons_of(d, channel):
    """[(prefix, is_electron)] of the selected light leptons."""
    if channel == "mutau":
        return [("l", False)]
    if channel == "etau":
        return [("l", True)]
    return [("mu", False), ("el", True)]


KIN_KEYS_LTAU = ("l_pt", "t_pt", "m_vis", "m_tt", "met_x", "met_y", "met", "mt_1", "mt_2", "dzeta", "njets", "nbjets", "jet1_pt")
KIN_KEYS_EMU = ("mu_pt", "el_pt", "m_vis", "m_tt", "met_x", "met_y", "met", "mt_1", "mt_2", "dzeta", "njets", "nbjets", "jet1_pt")
_KIN_CACHE: dict = {}


def _derived(out, d, channel):
    """m_T of the leptons, D_zeta and MET magnitude from the (possibly varied) momenta and MET."""
    p1, p2 = (("l", "t") if channel in LTAU else ("mu", "el"))
    pt1, pt2 = out[f"{p1}_pt"].astype(np.float64), out[f"{p2}_pt"].astype(np.float64)
    metx, mety = out["met_x"].astype(np.float64), out["met_y"].astype(np.float64)
    out["met"] = np.hypot(metx, mety).astype(np.float32)
    out["mt_1"] = mass.transverse_mass(pt1, d[f"{p1}_phi"], metx, mety).astype(np.float32)
    out["mt_2"] = mass.transverse_mass(pt2, d[f"{p2}_phi"], metx, mety).astype(np.float32)
    out["dzeta"] = mass.dzeta(pt1, d[f"{p1}_phi"], pt2, d[f"{p2}_phi"], metx, mety).astype(np.float32)


def kinematics(d, key: str, channel: str, variation: str | None = None, direction: str | None = None, fit_only=None):
    """Kinematic quantities under a variation; `fit_only`: boolean mask of the events for which the (slow)
    likelihood mass is recomputed (the others keep the nominal m_tt). tautau: analysis.kinematics with the
    3% energy-scale prior."""
    if channel == "tautau":
        return analysis.kinematics(d, key, variation, direction, tes_rel=config.TES_PRIOR_V4 if variation and variation.startswith("TauES") else None)
    keys = KIN_KEYS_LTAU if channel in LTAU else KIN_KEYS_EMU
    base = {k: d[k] for k in keys}
    if variation is None:
        return base
    ck = (key, channel, variation, direction, len(d["run"]))
    if ck in _KIN_CACHE:
        return _KIN_CACHE[ck]
    sign = 1.0 if direction == "Up" else -1.0
    out = {k: v.copy() for k, v in base.items()}
    p1, p2 = (("l", "t") if channel in LTAU else ("mu", "el"))
    pt = {p: d[f"{p}_pt"].astype(np.float64) for p in (p1, p2)}
    m = {p: d[f"{p}_mass"].astype(np.float64) for p in (p1, p2)}
    metx, mety = d["met_x"].astype(np.float64), d["met_y"].astype(np.float64)
    affected = np.zeros(len(metx), dtype=bool)

    def scale(p, factor, mask):
        nonlocal metx, mety
        new = np.where(mask, pt[p] * factor, pt[p])
        metx = metx - (new - pt[p]) * np.cos(d[f"{p}_phi"])
        mety = mety - (new - pt[p]) * np.sin(d[f"{p}_phi"])
        m[p] = np.where(mask, m[p] * factor, m[p])
        pt[p] = new

    if variation.startswith("TauES_DM"):
        dm = int(variation[len("TauES_DM"):])
        aff = (d["t_genflav"] == 5) & (d["t_dm"] == dm)
        scale("t", 1 + sign * config.TES_PRIOR_V4, aff)
        affected |= aff
    elif variation == "MuonScale":
        p = "l" if channel == "mutau" else "mu"
        scale(p, 1 + sign * MUON_SCALE_REL, np.ones(len(metx), bool))
        affected[:] = True
    elif variation == "ElectronScale":
        p = "l" if channel == "etau" else "el"
        endcap = np.abs(d[f"{p}_sceta"]) > 1.479
        rel = np.where(endcap, ELECTRON_SCALE_REL[1], ELECTRON_SCALE_REL[0])
        scale(p, 1 + sign * rel, np.ones(len(metx), bool))
        affected[:] = True
    elif variation == "MET_Unclustered":
        metx = metx + sign * d["met_uncl_dx"]
        mety = mety + sign * d["met_uncl_dy"]
        affected[:] = True
    elif variation == "JES":
        tag = "_jesUp" if direction == "Up" else "_jesDown"
        for k in ("njets", "nbjets", "jet1_pt"):
            out[k] = d[f"{k}{tag}"]
        metx = metx + d[f"met_dx{tag}"]
        mety = mety + d[f"met_dy{tag}"]
        affected[:] = True
    else:
        raise ValueError(variation)
    out[f"{p1}_pt"], out[f"{p2}_pt"] = pt[p1].astype(np.float32), pt[p2].astype(np.float32)
    out["met_x"], out["met_y"] = metx.astype(np.float32), mety.astype(np.float32)
    _derived(out, d, channel)
    idx = np.where(affected & (fit_only if fit_only is not None else True))[0]
    if len(idx):
        t1 = (pt[p1][idx], d[f"{p1}_eta"][idx].astype(np.float64), d[f"{p1}_phi"][idx].astype(np.float64), m[p1][idx])
        t2 = (pt[p2][idx], d[f"{p2}_eta"][idx].astype(np.float64), d[f"{p2}_phi"][idx].astype(np.float64), m[p2][idx])
        mvis, mtt = masses_parallel(t1, t2, metx[idx], mety[idx], d["met_covxx"][idx].astype(np.float64),
                                    d["met_covxy"][idx].astype(np.float64), d["met_covyy"][idx].astype(np.float64),
                                    leptonic=(True, channel == "emu"))
        out["m_vis"][idx] = mvis
        out["m_tt"][idx] = mtt
    _KIN_CACHE[ck] = out
    return out


# ------------------------------------------------------------------------------ parallel likelihood mass
_POOL = None


def _mass_worker(args):
    t1, t2, metx, mety, cxx, cxy, cyy, leptonic = args
    return mass.visible_mass(t1, t2), mass.likelihood_mass(t1, t2, metx, mety, cxx, cxy, cyy, leptonic=leptonic)


def masses_parallel(t1, t2, metx, mety, cxx, cxy, cyy, leptonic, workers: int | None = None):
    """(m_vis, m_tt) of many events, the likelihood mass evaluated in a persistent process pool (spawn context:
    the parent holds uproot locks, CLAUDE.md pitfall 3). ~1 ms per event per core."""
    global _POOL
    import multiprocessing as mp
    n = len(metx)
    workers = workers or min(config.N_WORKERS, 8)
    if n < 4000 or workers < 2:
        return _mass_worker((t1, t2, metx, mety, cxx, cxy, cyy, leptonic))
    if _POOL is None:
        _POOL = mp.get_context("spawn").Pool(workers)
    edges = np.linspace(0, n, workers * 2 + 1).astype(int)
    jobs = [(tuple(a[lo:hi] for a in t1), tuple(a[lo:hi] for a in t2), metx[lo:hi], mety[lo:hi], cxx[lo:hi], cxy[lo:hi], cyy[lo:hi], leptonic)
            for lo, hi in zip(edges[:-1], edges[1:]) if hi > lo]
    res = _POOL.map(_mass_worker, jobs)
    return np.concatenate([r[0] for r in res]), np.concatenate([r[1] for r in res])


# ------------------------------------------------------------------------------ regions
def _lepton_ok(d, kin, channel):
    if channel == "mutau":
        return (kin["l_pt"] > config.MU_PT_MIN_MUTAU) & (np.abs(d["l_eta"]) < config.MU_ETA_MAX)
    if channel == "etau":
        return (kin["l_pt"] > config.EL_PT_MIN_ETAU) & (np.abs(d["l_eta"]) < config.EL_ETA_MAX)
    mu_ok = (kin["mu_pt"] > config.EMU_MU_PT_MIN) & (np.abs(d["mu_eta"]) < config.MU_ETA_MAX)
    el_ok = (kin["el_pt"] > config.EMU_EL_PT_MIN) & (np.abs(d["el_eta"]) < 2.5)
    # trigger-leg logic: Mu23_Ele12 needs mu > 24, Mu8_Ele23 needs e > 24 (the OR of the paths fired at the skim)
    legs = (kin["mu_pt"] > config.EMU_LEAD_PT_MIN) | (kin["el_pt"] > config.EMU_LEAD_PT_MIN)
    return mu_ok & el_ok & legs


def fit_candidates(d, channel: str, kin=None, is_mc: bool = False):
    """Loose superset of every fitted region (used to limit the recomputation of m_tt under variations)."""
    kin = kin or d
    veto = (d["n_veto_mu"] == 0) & (d["n_veto_el"] == 0)
    if channel in LTAU:
        return veto & d["os"] & (d["l_iso"] < config.MU_ISO if channel == "mutau" else d["l_iso"] < config.EL_ISO) & \
            ((d["t_vsjet"] & T_BIT) > 0) & (d["t_pt"] > config.LTAU_TAU_PT_MIN - 3)
    return veto & d["os"] & (d["mu_iso"] < config.EMU_ISO) & (d["el_iso"] < config.EMU_ISO)


def regions(d, channel: str, kin=None, is_mc: bool = False) -> dict:
    kin = kin or d
    veto = (d["n_veto_mu"] == 0) & (d["n_veto_el"] == 0)
    os_ = d["os"].astype(bool)
    ss = ~os_
    if channel in LTAU:
        iso_max = config.MU_ISO if channel == "mutau" else config.EL_ISO
        base = veto & _lepton_ok(d, kin, channel) & (kin["t_pt"] > config.LTAU_TAU_PT_MIN) & (np.abs(d["t_eta"]) < config.LTAU_TAU_ETA_MAX)
        if is_mc:
            base &= d["t_genflav"] != 0
        iso = d["l_iso"] < iso_max
        ai = (d["l_iso"] > config.LEP_ANTIISO[0]) & (d["l_iso"] < config.LEP_ANTIISO[1])
        t = (d["t_vsjet"] & T_BIT) > 0
        lo = ((d["t_vsjet"] & L_BIT) > 0) & ~t
        mt_lo = kin["mt_1"] < config.LTAU_MT_MAX
        mt_hi = kin["mt_1"] > config.LTAU_MT_WDR_MIN
        bveto = kin["nbjets"] == 0
        return {"SR": base & iso & os_ & t & mt_lo, "AR": base & iso & os_ & lo & mt_lo,
                "QCD_T": base & iso & ss & t & mt_lo, "QCD_L": base & iso & ss & lo & mt_lo,
                "W_T": base & iso & os_ & t & mt_hi & bveto, "W_L": base & iso & os_ & lo & mt_hi & bveto,
                "WSS_T": base & iso & ss & t & mt_hi & bveto, "WSS_L": base & iso & ss & lo & mt_hi & bveto,
                "AI_OS_T": base & ai & os_ & t & mt_lo, "AI_OS_L": base & ai & os_ & lo & mt_lo,
                "AI_SS_T": base & ai & ss & t & mt_lo, "AI_SS_L": base & ai & ss & lo & mt_lo,
                "SS_SR": base & iso & ss & t & mt_lo, "SS_AR": base & iso & ss & lo & mt_lo,
                "HIMT_SR": base & iso & os_ & t & mt_hi, "HIMT_AR": base & iso & os_ & lo & mt_hi}
    base = veto & _lepton_ok(d, kin, channel)
    iso = (d["mu_iso"] < config.EMU_ISO) & (d["el_iso"] < config.EMU_ISO)
    relaxed = (d["mu_iso"] < config.EMU_SB1_ISO) & (d["el_iso"] < config.EMU_SB1_ISO)
    sb1 = relaxed & ((d["mu_iso"] > config.EMU_ISO) | (d["el_iso"] > config.EMU_ISO))
    sb2 = relaxed & ((d["mu_iso"] > config.EMU_SB2_ISO) | (d["el_iso"] > config.EMU_SB2_ISO))
    dz = kin["dzeta"] > config.EMU_DZETA_MIN
    bveto = (kin["nbjets"] == 0) if config.EMU_BVETO else np.ones(len(os_), bool)
    cr = (kin["dzeta"] < config.EMU_CR_DZETA_MAX) & (kin["met"] > config.EMU_CR_MET_MIN)
    return {"SR": base & iso & os_ & dz & bveto, "CRtt": base & iso & os_ & cr, "SS": base & iso & ss & dz & bveto,
            "SB1_OS": base & sb1 & os_ & dz & bveto, "SB1_SS": base & sb1 & ss & dz & bveto,
            "SB2_OS": base & sb2 & os_ & dz & bveto, "SB2_SS": base & sb2 & ss & dz & bveto,
            "CRtt_SS": base & iso & ss & cr}


def dm_key(d, channel: str, kin=None) -> np.ndarray:
    """Per event: the key of the tau_h ID scale-factor template it belongs to ('none' or the decay mode(s) of
    the genuine tau_h legs, e.g. '1', '0_10'); array of str."""
    if channel == "emu":
        return np.full(len(d["run"]), "none", dtype=object)
    if channel in LTAU:
        gen = d["t_genflav"] == 5
        return np.where(gen, d["t_dm"].astype(str), "none").astype(object)
    g1, g2 = d["t1_genflav"] == 5, d["t2_genflav"] == 5
    a = np.where(g1, d["t1_dm"].astype(int), -1)
    b = np.where(g2, d["t2_dm"].astype(int), -1)
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    out = np.full(len(a), "none", dtype=object)
    one = (lo < 0) & (hi >= 0)
    out[one] = hi[one].astype(str)
    both = lo >= 0
    out[both] = np.char.add(np.char.add(lo[both].astype(str), "_"), hi[both].astype(str))
    return out


# ------------------------------------------------------------------------------ weights
@lru_cache(maxsize=1)
def trigger_insitu() -> dict:
    return json.loads(TRIG_INSITU.read_text()) if TRIG_INSITU.exists() else {}


def _lookup2(table, x, y, syst=0):
    """Look up a 2D table {'x_edges', 'y_edges', 'sf', 'err'} (x clipped into range); syst = +-1 shifts by err."""
    xe, ye = np.asarray(table["x_edges"]), np.asarray(table["y_edges"])
    ix = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, len(xe) - 2)
    iy = np.clip(np.searchsorted(ye, np.abs(y) if table.get("abs_y", True) else y, side="right") - 1, 0, len(ye) - 2)
    v = np.asarray(table["sf"])[ix, iy]
    if syst:
        v = v + syst * np.asarray(table["err"])[ix, iy]
    return v


def _insitu_sf(name, pt, eta, syst, flat: float = TRIG_FLAT_UNC):
    """In-situ trigger scale factor of one leg. `syst` = +-1 shifts it by the statistical error of the
    table bin and, on top of that, by the flat `flat` (the paper's 2 % per *channel*, CMS
    arXiv:1801.03535 Table 2). A cross trigger has two legs but only one such flat term: the caller
    passes flat=0 for the second leg (review/REVIEW_v4.md finding 3 -- applying it per leg made the e mu
    trigger prior 5.4 % instead of 2.7 % and moved mu_Z by one standard deviation)."""
    t = trigger_insitu().get(name)
    if t is None:
        return np.ones(len(pt)) * (1 + syst * flat)         # no measurement: SF 1 +- flat
    v = _lookup2(t, pt, eta, syst)
    return v * (1 + syst * flat) if syst else v


@lru_cache(maxsize=1)
def ele27_turnon_rel() -> np.ndarray:
    """Relative uncertainty of the Ele27_WPTight scale factor in the turn-on bin (EL_PT_MIN_ETAU to
    ELE27_PLATEAU_PT), per |eta_SC| bin: the step to the next pT bin of the in-situ table. The bin-averaged
    scale factor mis-models an event at the edge of the bin by that much (review/REVIEW_v4.md finding 7:
    4 % / 2 % / 13 % in the three eta bins the e tau_h channel uses, against the flat 2 % of the plateau)."""
    t = trigger_insitu().get("ele27")
    if t is None:
        return np.full(4, TRIG_FLAT_UNC)
    sf = np.asarray(t["sf"], dtype=float)
    i = int(np.clip(np.searchsorted(np.asarray(t["x_edges"]), config.EL_PT_MIN_ETAU, side="right") - 1, 0, len(sf) - 2))
    return np.abs(sf[i + 1] - sf[i]) / np.maximum(sf[i], 1e-9)


def _ele27_sf(pt, sceta, syst_plateau: int = 0, syst_turnon: int = 0):
    """Ele27_WPTight in-situ scale factor with *two* nuisance parameters: `ElectronTrigger` above
    ELE27_PLATEAU_PT (table statistics + the flat 2 %) and `ElectronTrigger_lowpt` in the turn-on bin
    (table statistics + `ele27_turnon_rel`). The e tau_h electrons start at 29 GeV, on the turn-on."""
    t = trigger_insitu().get("ele27")
    if t is None:
        return np.ones(len(pt)) * (1 + (syst_plateau + syst_turnon) * TRIG_FLAT_UNC)
    v = _lookup2(t, pt, sceta, 0)
    if not (syst_plateau or syst_turnon):
        return v
    low = np.asarray(pt) < config.ELE27_PLATEAU_PT
    syst = np.where(low, syst_turnon, syst_plateau)
    stat = (_lookup2(t, pt, sceta, 1) - v) / np.maximum(v, 1e-9)        # +1 sigma of the bin, signed below
    ye = np.asarray(t["y_edges"])
    iy = np.clip(np.searchsorted(ye, np.abs(sceta), side="right") - 1, 0, len(ye) - 2)
    flat = np.where(low, ele27_turnon_rel()[iy], TRIG_FLAT_UNC)
    return v * (1 + syst * (stat + flat))


@lru_cache(maxsize=1)
def btag_efficiencies() -> dict:
    return json.loads(BTAG_EFF.read_text()) if BTAG_EFF.exists() else {}


def btag_weight(d, kin, syst: int = 0) -> np.ndarray:
    """Event weight of the b-jet veto (jets in the b-tag acceptance): prod over untagged jets of
    (1 - SF eps) / (1 - eps) and over tagged jets of SF, with eps(flavour, pT) from the simulation
    (`btag_efficiencies`). Without an efficiency map: 1 (the BTag NP is then the estimate of docs/10)."""
    eff = btag_efficiencies()
    n = len(d["run"])
    w = np.ones(n)
    if not eff:
        return w
    edges = np.asarray(eff["pt_edges"])
    sname = {0: "central", 1: "up", -1: "down"}[syst]
    for i in range(1, 5):
        pt, eta, fl, tag = d[f"bj{i}_pt"], d[f"bj{i}_eta"], d[f"bj{i}_flav"].astype(int), d[f"bj{i}_tag"].astype(bool)
        has = pt > 0
        if not has.any():
            continue
        sf = np.ones(n)
        sf[has] = pog.btag_sf(fl[has], eta[has], pt[has], syst=sname)
        ib = np.clip(np.searchsorted(edges, pt, side="right") - 1, 0, len(edges) - 2)
        e = np.ones(n) * 0.5
        for f, key in ((5, "b"), (4, "c"), (0, "light")):
            sel = has & (fl == f)
            e[sel] = np.asarray(eff[key])[ib[sel]]
        e = np.clip(e, 1e-4, 0.999)
        wj = np.where(tag, sf, (1 - sf * e) / (1 - e))
        w = w * np.where(has, wj, 1.0)
    return w


def weights(d, key: str, channel: str, syst: str | None = None, direction: str | None = None, kin=None):
    """Event weights of a lepton-channel MC ntuple: norm x genWeight x pileup x prefiring x lepton SFs
    (ID, iso, trigger) x tau_h VSe/VSmu SFs (lepton-faked tau_h) x b-tag veto weight x top pT. tautau: the v3
    weights without the VSjet scale factor."""
    if channel == "tautau":
        return analysis.weights(d, key, syst, direction, kin, vsjet_sf=False)
    kin = kin or d
    sign = {"Up": 1, "Down": -1, None: 0}[direction]
    s = samples.SAMPLES[key]
    pu_var = {"Up": "up", "Down": "down"}.get(direction) if syst == "Pileup" else "nominal"
    norm = dy_norm(key, channel, d) if key in samples.DY_STITCHED else normalisation(key, channel)
    w = norm * d["genWeight"].astype(np.float64) * analysis.pileup()(d["Pileup_nTrueInt"], pu_var)
    pref = {"Up": "L1PreFiringWeight_Up", "Down": "L1PreFiringWeight_Dn"}.get(direction) if syst == "L1Prefiring" else None
    w = w * d[pref or "L1PreFiringWeight_Nom"]
    for prefix, is_el in _leptons_of(d, channel):
        pt, eta = kin[f"{prefix}_pt"], d[f"{prefix}_eta"]
        if is_el:
            sceta = d[f"{prefix}_sceta"]
            w = w * pog.electron_sf("reco", pt, sceta, {1: "sfup", -1: "sfdown"}[sign] if syst == "ElectronReco" and sign else "sf")
            w = w * pog.electron_sf("id", pt, sceta, {1: "sfup", -1: "sfdown"}[sign] if syst == "ElectronID" and sign else "sf")
            if channel == "etau":
                w = w * _ele27_sf(pt, sceta, sign if syst == "ElectronTrigger" else 0,
                                  sign if syst == "ElectronTrigger_lowpt" else 0)
        else:
            sy = {1: "systup", -1: "systdown"}
            w = w * pog.muon_sf("id", pt, eta, sy[sign] if syst == "MuonID" and sign else "nominal")
            w = w * pog.muon_sf("iso", pt, eta, sy[sign] if syst == "MuonIso" and sign else "nominal")
            if channel == "mutau":
                w = w * pog.muon_sf("trig", pt, eta, sy[sign] if syst == "MuonTrigger" and sign else "nominal")
    if channel == "emu":
        # one nuisance parameter for the cross trigger: the statistical error of *both* legs plus the flat
        # 2 % of the paper *once* (finding 3 of review/REVIEW_v4.md); total prior 2.3-2.7 % per event.
        w = w * _insitu_sf("emu_e", kin["el_pt"], d["el_sceta"], sign if syst == "EmuTrigger" else 0)
        w = w * _insitu_sf("emu_mu", kin["mu_pt"], d["mu_eta"], sign if syst == "EmuTrigger" else 0, flat=0.0)
    if channel in LTAU:
        vse_wp = "VVLoose" if channel == "mutau" else "Tight"
        vsmu_wp = "Tight" if channel == "mutau" else "VLoose"
        sy = {1: "up", -1: "down"}
        w = w * pog.tau_vse_sf(vse_wp, d["t_eta"], d["t_genflav"], sy[sign] if syst == "TauFakeEle" and sign else "nom")
        w = w * pog.tau_vsmu_sf(vsmu_wp, d["t_eta"], d["t_genflav"], sy[sign] if syst == "TauFakeMu" and sign else "nom")
    w = w * btag_weight(d, kin, sign if syst == "BTag" else 0)
    if s["group"] == "top" and "gen_top_pt" in d:
        power = {1: 2.0, -1: 0.0}[sign] if syst == "TopPt" and sign else 1.0
        w = w * pog.top_pt_weight(d["gen_top_pt"], d["gen_antitop_pt"], power)
    return w


def theory_weights(d, key: str, channel: str, syst: str):
    """As analysis.theory_weights, renormalised to the sigma(60 < m_LHE < 120) generator sum instead of the
    tau_h tau_h fiducial one (the v4 signal definition, docs/10 section 6)."""
    g = load(key, channel)[1]["gensums"]
    ref0 = g["sumw_lhe_tautau_60_120"]
    branch, sums_key, idx = {"QCDScale": ("LHEScaleWeight", "scale_lhe_tautau_60_120", [0, 1, 3, 5, 7, 8]),
                             "PDF": ("LHEPdfWeight", "pdf_lhe_tautau_60_120", list(range(103))),
                             "PS_ISR": ("PSWeight", "ps_lhe_tautau_60_120", [0, 2]),
                             "PS_FSR": ("PSWeight", "ps_lhe_tautau_60_120", [1, 3])}[syst]
    if branch not in d or sums_key not in g:
        return None
    v = d[branch].astype(np.float64)[:, idx]
    sums = np.asarray(g[sums_key])[idx]
    return v * (ref0 / sums)[None, :]


combine_theory = analysis.combine_theory


def signal_prediction(channel: str = "tautau", key: str = samples.DY_INCLUSIVE) -> dict:
    """sigma(Z/gamma* -> tautau, 60-120) of the sample and its theory uncertainties (the v4 reference)."""
    g = load(key, channel)[1]["gensums"]
    xs = samples.SAMPLES[key]["xsec_pb"]
    return {"sigma_tautau_pb": xs * g["sumw_lhe_tautau"] / g["sumw"],
            "sigma_tautau_60_120_pb": xs * g["sumw_lhe_tautau_60_120"] / g["sumw"],
            "sigma_fid_tautau_pb": xs * g["sumw_fid"] / g["sumw"]}
