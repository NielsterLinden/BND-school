"""MC event weights of the v2 analysis (docs/11-mc-weights.md).

    w = sign(genWeight) x (sigma L / sum genWeight) x w_PU x L1PreFiringWeight_Nom
        x SF_ID(mu1) SF_ID(mu2) x SF_iso(mu1) SF_iso(mu2) x SF_trigger(event)

`Weighter.variations(...)` returns the nominal weight and every Up/Down alternative used for
the fit templates. Theory variations (PDF, QCD scales, parton shower) are renormalised so
that the *fiducial* generator-level yield is unchanged: they vary the C factor only, which is
what a fiducial measurement must carry (`theory_renorm`).
"""

from __future__ import annotations

import numpy as np

from . import config, samples, skim

LUMI_PB = config.LUMI_PB_NORMTAG

PDF_MEMBERS = slice(1, 101)      # NNPDF3.1 Hessian eigenvectors (LHA 325300-325400)
ALPHAS_MEMBERS = (101, 102)      # alpha_s 0.116 / 0.120
SCALE_7POINT = [0, 1, 3, 4, 5, 7, 8]   # drop (muR,muF) = (0.5,2) and (2,0.5)
SCALE_NOMINAL = 4
PS_ISR_UP, PS_FSR_UP, PS_ISR_DOWN, PS_FSR_DOWN = 0, 1, 2, 3


class Weighter:
    """Per-sample normalisation, pileup and prefiring weights; SF application via callbacks."""

    def __init__(self, key: str, pileup, sf=None, base=None):
        self.key = key
        self.sample = samples.SAMPLES[key]
        self.gensums = skim.load_gensums(key, base) if self.sample["is_mc"] else {}
        self.pileup = pileup
        self.sf = sf
        if self.sample["is_mc"]:
            sumw = float(self.gensums.get("sumw", 0.0))
            xsec = self.sample["xsec_pb"]
            self.norm = (xsec * LUMI_PB / sumw) if (xsec and sumw) else 1.0
        else:
            self.norm = 1.0

    # ---- pieces ------------------------------------------------------------
    def base(self, ev):
        """genWeight x normalisation x pileup x prefiring (per event, MC) or prescale (data)."""
        if not self.sample["is_mc"]:
            return np.asarray(ev.skim_prescale, dtype=float)
        w = np.asarray(ev.genWeight, dtype=float) * self.norm
        w = w * self.pileup(np.asarray(ev.Pileup_nTrueInt))
        w = w * np.asarray(ev.L1PreFiringWeight_Nom, dtype=float)
        return w

    def pieces(self, ev):
        """Dict of the individual multiplicative factors (for variations and plots)."""
        if not self.sample["is_mc"]:
            return {}
        n_true = np.asarray(ev.Pileup_nTrueInt)
        return {
            "gen": np.asarray(ev.genWeight, dtype=float) * self.norm,
            "pu": self.pileup(n_true), "pu_up": self.pileup(n_true, "up"), "pu_down": self.pileup(n_true, "down"),
            "pref": np.asarray(ev.L1PreFiringWeight_Nom, dtype=float),
            "pref_up": np.asarray(ev.L1PreFiringWeight_Up, dtype=float),
            "pref_down": np.asarray(ev.L1PreFiringWeight_Dn, dtype=float),
        }

    # ---- theory variations -------------------------------------------------
    def theory_renorm(self, kind: str) -> np.ndarray:
        """Factor per variation member that keeps the fiducial (dressed) yield fixed."""
        g = self.gensums
        fid = float(g["sumw_fid_dressed"])
        var = np.asarray(g[f"{kind}_fid_dressed"], dtype=float)
        return np.where(var != 0, fid / np.where(var != 0, var, 1.0), 1.0)

    def theory_weights(self, ev, kind: str, renorm: bool = True):
        """(n_events, n_members) relative weights of LHEPdfWeight / LHEScaleWeight / PSWeight."""
        import awkward as ak
        branch = {"pdf": "LHEPdfWeight", "scale": "LHEScaleWeight", "ps": "PSWeight"}[kind]
        if branch not in ak.fields(ev):
            return None
        arr = ev[branch]
        counts = ak.to_numpy(ak.num(arr, axis=1))
        if counts.size == 0 or counts.min() != counts.max():
            return None
        v = ak.to_numpy(ak.to_regular(arr, axis=1)).astype(float)
        if renorm and f"{kind}_fid_dressed" in self.gensums:
            v = v * self.theory_renorm(kind)[None, :]
        return v


def pdf_envelope(hist_members, nominal):
    """Hessian: symmetric up/down = nominal +/- sqrt(sum (h_i - h_0)^2) over the 100 eigenvectors."""
    delta = np.sqrt(np.sum((hist_members[PDF_MEMBERS] - nominal[None, :]) ** 2, axis=0))
    return nominal + delta, nominal - delta


def alphas_variation(hist_members, nominal):
    """alpha_s +/- 0.0015 from the +/- 0.002 members, symmetrised."""
    d = 0.75 * 0.5 * (hist_members[ALPHAS_MEMBERS[1]] - hist_members[ALPHAS_MEMBERS[0]])
    return nominal + d, nominal - d


def scale_envelope(hist_members, nominal):
    """7-point envelope of the muR/muF variations, symmetrised around the nominal."""
    members = hist_members[SCALE_7POINT]
    up = np.max(members, axis=0)
    down = np.min(members, axis=0)
    return up, down


def ps_variations(hist_members):
    return {"PS_ISR": (hist_members[PS_ISR_UP], hist_members[PS_ISR_DOWN]),
            "PS_FSR": (hist_members[PS_FSR_UP], hist_members[PS_FSR_DOWN])}
