"""Manim-free palette for the BND-school Z cross-section talk.

Import this from matplotlib scripts (static figures) AND from style/bnd_style.py
(Manim scenes) so both use identical hex values. No manim import here, ever.

THE COLOUR SCHEMA (decided 2026-09-15): six chapter colours, one per section
of the talk, run through everything that is *ours* (particles, plots, fits,
highlights). The CMS detector itself is drawn in the official CMS logo
colours (Izaak Neutelings' CMS_logo.tex, cms-docdb 3045), so the detector
is always recognisably "CMS" and never competes with the chapter code.
Every other colour is a tint (toward white) or a shade (toward black) of one
of these, made with ``tint()`` / ``shade()``.
"""

# -- the six chapters ---------------------------------------------------------
CHAPTER = {
    1: "#702677",   # theory        purple
    2: "#00B3C3",   # detector      cyan
    3: "#059F77",   # Z -> ee       green
    4: "#E3D88B",   # Z -> mumu     gold
    5: "#E6223C",   # Z -> tautau   red
    6: "#2B2D42",   # combination   dark slate
}
CHAPTER_NAME = {1: "theory", 2: "detector", 3: "ee", 4: "mumu", 5: "tautau", 6: "combination"}
PURPLE, CYAN, GREEN, GOLD, RED, SLATE = (CHAPTER[i] for i in range(1, 7))
BLACK, WHITE = "#000000", "#FFFFFF"


def _hex2rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def _rgb2hex(rgb):
    return "#%02X%02X%02X" % tuple(max(0, min(255, round(v))) for v in rgb)


def mix(h1: str, h2: str, a: float) -> str:
    """Linear mix in RGB: a=0 -> h1, a=1 -> h2."""
    c1, c2 = _hex2rgb(h1), _hex2rgb(h2)
    return _rgb2hex(tuple(x + a * (y - x) for x, y in zip(c1, c2)))


def tint(h: str, a: float) -> str:
    """Soft version of a colour: move it toward white by ``a`` in [0, 1]."""
    return mix(h, WHITE, a)


def shade(h: str, a: float) -> str:
    """Move a colour toward black by ``a`` in [0, 1] (strokes on light fills)."""
    return mix(h, BLACK, a)


# -- deck ---------------------------------------------------------------------
BG = WHITE            # white background: reads on a projector, matches PowerPoint
INK = SLATE           # lines, axes, quark lines, symbols: the chapter-6 slate, not pure black
GREY = tint(SLATE, 0.55)        # de-emphasis / extrapolation / not-measured
LIGHT_GREY = tint(SLATE, 0.85)

# -- the three channels (THE colour code of the talk) --------------------------
CHANNEL = {
    "ee":     GREEN,
    "mumu":   GOLD,
    "tautau": RED,
}
# strokes / glyphs in the channel colour: gold is too light on white for thin
# lines and text, so lines and labels use these; fills use CHANNEL[...]
CHANNEL_LINE = {
    "ee":     GREEN,
    "mumu":   shade(GOLD, 0.30),
    "tautau": RED,
}
# a thin dark rim drawn under every channel track (the CMS logo's own
# double-line muon style) so the pale gold track still reads on white
TRACK_OUTLINE = {k: shade(v, 0.55) for k, v in CHANNEL.items()}

COMBINED = SLATE           # the combined result: chapter 6
THEORY = PURPLE            # theory prediction band / line: chapter 1
DETECTOR_ACCENT = CYAN     # anything "detector / method" that is not the detector drawing itself
HIGHLIGHT = RED            # a single bright "look here"; in tautau clips Flash in SLATE instead

# -- samples (fit-input names from docs/CONVENTIONS.md) ---------------------
# Signals in their channel colour; backgrounds are tints of the slate and of the
# two non-channel chapter colours, so they never compete with the signal.
SAMPLE = {
    "Data":      SLATE,
    "DYee":      CHANNEL["ee"],
    "DYmumu":    CHANNEL["mumu"],
    "DYtautau":  CHANNEL["tautau"],
    "TTbar":     tint(SLATE, 0.30),
    "SingleTop": tint(SLATE, 0.55),
    "WW":        tint(CYAN, 0.35),
    "WZ":        tint(PURPLE, 0.45),
    "ZZ":        tint(CYAN, 0.65),
    "Fakes":     tint(SLATE, 0.78),   # data-driven, "not simulation"
}

# -- physics objects / particles ----------------------------------------------
PARTICLE = {
    "e":      CHANNEL_LINE["ee"],
    "mu":     CHANNEL_LINE["mumu"],
    "tau":    CHANNEL_LINE["tautau"],
    "photon": shade(GOLD, 0.15),
    "Z":      INK,
    "W":      INK,
    "quark":  INK,
    "pion":   tint(SLATE, 0.30),
    "hadron": tint(SLATE, 0.35),
    "jet":    tint(SLATE, 0.55),
    "nu":     GREY,        # neutrinos: undetected, drawn dashed
    "met":    SLATE,       # missing transverse momentum arrow
}

# -- CMS detector: the official logo colours ----------------------------------
CMS = {
    "tib":        "#ADEB8F",   # inner strip tracker, light green   (173,235,143)
    "tob":        "#92E569",   # outer strip tracker, green         (146,229,105)
    "ecal":       "#92E569",   # ECAL, same green as the TOB        (146,229,105)
    "band":       "#BEAED4",   # thin lavender band outside the ECAL (190,174,212)
    "hcal":       "#E5E9EE",   # HCAL, light grey                   (229,233,238)
    "magnet":     "#FFDF7F",   # solenoid, yellow                   (255,223,127)
    "muon":       "#85D1FB",   # muon system, blue                  (133,209,251)
    "muon_track": "#F0240B",   # the four red muons of the logo     (240,36,11)
    "bg":         "#FDFAF4",   # logo background, off-white         (253,250,244)
    "font":       "#101177",   # the CMS wordmark, dark blue        (16,17,119)
}
# transverse view: fill = logo colour, stroke = a shade of the same colour.
DETECTOR = {
    "beampipe": {"fill": CMS["tib"],     "stroke": tint(SLATE, 0.35)},
    "tracker":  {"fill": CMS["tib"],     "stroke": shade(CMS["tob"], 0.35)},
    "tob":      {"fill": CMS["tob"],     "stroke": shade(CMS["tob"], 0.35)},
    "ecal":     {"fill": CMS["ecal"],    "stroke": shade(CMS["ecal"], 0.35)},
    "band":     {"fill": CMS["band"],    "stroke": shade(CMS["band"], 0.30)},
    "hcal":     {"fill": CMS["hcal"],    "stroke": shade(CMS["hcal"], 0.35)},
    "solenoid": {"fill": CMS["magnet"],  "stroke": shade(CMS["magnet"], 0.35)},
    "muon":     {"fill": CMS["muon"],    "stroke": shade(CMS["muon"], 0.40)},
    "yoke":     {"fill": CMS["hcal"],    "stroke": shade(CMS["hcal"], 0.30)},
}
# where a hit / energy deposit of a *non-channel* particle is drawn (a channel
# lepton lights its deposits in its own chapter colour instead)
DEPOSIT = {
    "ecal": shade(CMS["ecal"], 0.35),
    "hcal": shade(CMS["band"], 0.20),
    "muon": shade(CMS["muon"], 0.35),
}

# outro only (7_outro): the three countries of the BND school, for the closing
# "mic drop" rays. Not part of the chapter schema; never used for physics content.
FLAG = {
    "NL": ["#AE1C28", "#FFFFFF", "#21468B"],
    "BE": ["#000000", "#FDDA24", "#EF3340"],
    "DE": ["#000000", "#DD0000", "#FFCC00"],
}

FONT = "DejaVu Sans"

# -- section 6: figures redrawn in the colours of the analysis figures ----------
# Deck owner's request (17 Sep 2026): the combination clips rebuild summary.png, impacts.png and
# combined_vs_published.png (combination/combLieke/mf/plots.py) and two e mu mass plots in the colours of
# those figures, not in the chapter colours. Copied from the plotting code; never used elsewhere.
COMB_FIG = {
    "ee":          "#8E44AD",   # mf/plots.py C_EE
    "mumu":        "#2B6CB0",   # C_MUMU
    "tautau":      "#EB811B",   # C_TAUTAU
    "combined":    "#23373B",   # C_COMB
    "prediction":  "#14B03D",   # C_PRED (band drawn at opacity 0.18)
    "published":   "#8C8C94",   # C_GREY (published points)
    "published_text": "#444444",
    "separator":   "#999999",
    "impact_up":   "#2B6CB0",   # theta_hat + delta theta_hat (opacity 0.85)
    "impact_down": "#5FA8D3",   # theta_hat - delta theta_hat (opacity 0.85)
    "pull_1sigma": "#7CCB7C",   # opacity 0.45
    "pull_2sigma": "#FFE066",   # opacity 0.45
}
# z-mumu/zmumu/plotting.py COLOURS (datamc_CRemu_mass_fit.png)
ZMUMU_FIG = {
    "DYmumu": "#F2B134", "DYtautau": "#A04CB0", "DYee": "#4C9EE0", "TTbar": "#D62728", "SingleTop": "#E07B7B",
    "WW": "#2CA6A4", "WZ": "#4DC0BE", "ZZ": "#7FD4D2", "WJets": "#8C6D4F", "MCStat": "#9A9A9A",
}
# z-tautau/ztautau/plotting.py COLORS (step4_*_m_vis.png)
ZTAUTAU_FIG = {
    "DYtautau": "#FFCC66", "Fakes": "#FF99CC", "TTbar": "#9999CC", "Diboson": "#66CC99", "WJets": "#D95F5F",
}
