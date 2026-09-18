"""Section 6: the combination in four clips, one per slide (deck owner's spec, 17 Sep 2026).

    CombSummary         comb_summary         summary.png built from scratch on an empty frame: the axis lines, then the
                                             x labels with "CMS Open Data", then Z->ee, Z->mumu, Z->tautau (label, point, value),
                                             the combined row with the grey divider and 16.4 fb-1, last the aMC@NLO
                                             line, band and legend
    CombOrthogonality   comb_orthogonality   opens on comb_summary's last frame: the plot shrinks x5 into the top-right
                                             corner; a schematic Z->tautau m_tautau stack (four final states) and the
                                             schematic Z->mumu e mu control region (as datamc_CRemu_mass_fit.png) side by
                                             side; the tautau e mu part takes the control-region colours while the rest
                                             goes grey; then the control region goes grey and the tautau plot is back
    CombImpacts         comb_impacts         opens on comb_orthogonality's last frame: the e mu plots leave, the two
                                             panels of impacts.png open side by side at the same time (empty axes), then
                                             each row arrives in both at once, from the smallest impact upward; all rows
                                             but the top three by impact dim, those three are outlined in both panels
    CombPublished       comb_published       opens on comb_impacts' last frame: the impacts fade, the corner summary
                                             returns to its slide-1 size; rows move up and tighten, CMS joins under the
                                             combined row, then again for ATLAS (as combined_vs_published.png)

The corner copy of summary.png is the through-line of slides 2-4: it is the same object in the last frame of
comb_orthogonality, the whole of comb_impacts and the first frame of comb_published, so the three clips cut on
identical frames. s7_outro.MicDropResults then opens on comb_published's last frame.

Numbers: data/combination_results.json (extract_combination.py, from combination/combLieke/output/result.json with
the selection and value strings of mf/plots.py). Colours: the colours of those figures (COMB_FIG, ZMUMU_FIG,
ZTAUTAU_FIG), on the deck owner's request, not the chapter colours. The two e mu plots are schematic (shapes only,
no y numbers). Keeps y > 2.7 and the top-left block (x < -5.85, y > 0.22) empty.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))

import numpy as np  # noqa: E402
from manim import (  # noqa: E402
    BOLD, DOWN, ITALIC, LEFT, RIGHT, UP, Create, DashedLine, Dot, FadeIn, FadeOut, GrowFromCenter, Line,
    MathTex, Polygon, Rectangle, RoundedRectangle, Scene, Square, Transform, VGroup, rate_functions,
)
from style.bnd_style import *  # noqa: E402,F401,F403
from style.bnd_style import _cap_scale  # noqa: E402

D = load_data("combination_results")
ROWS = {r["key"]: r for r in D["summary"]["rows"]}
PRED = D["prediction"]
PUB = {p["label"]: p for p in D["published"]}
IMP = D["impacts"]["rows"]

assert [round(ROWS[k]["sigma_pb"]) for k in ("ee", "mumu", "tautau", "combined")] == [2094, 1931, 1981, 1949]
assert ROWS["combined"]["value_tex"] == r"1949\ {}^{+30}_{-30}" and ROWS["tautau"]["value_tex"] == r"1981\ {}^{+153}_{-136}"
assert round(PRED["sigma_pb"], 1) == 1953.9 and round(PRED["band_lo_pb"]) == 1872 and round(PRED["band_hi_pb"]) == 2010
assert (PUB["CMS"]["value_tex"], PUB["ATLAS"]["value_tex"]) == (r"1952 \pm 49", r"2009 \pm 58")
assert PUB["CMS"]["detail"] == r"206 pb$^{-1}$, ee + $\mu\mu$" and PUB["ATLAS"]["detail"] == r"81 pb$^{-1}$, ee + $\mu\mu$"
assert PUB["ATLAS"]["moved_to_60_120"] and not PUB["CMS"]["moved_to_60_120"]
assert len(IMP) == 20 and [r["name"] for r in IMP[:3]] == ["Lumi", "ElectronID", "mu_ttbar"] and IMP[-1]["name"] == "MuonRes"
assert D["summary"]["x_range_pb"] == [1650, 2400]

EASE = rate_functions.ease_in_out_sine
F = COMB_FIG


def tex(expr: str, h: float = 0.24, color=INK) -> MathTex:
    """MathTex whose digits are ``h`` high."""
    return mathtex(expr, color=color).scale(_cap_scale("tex", h))


def txt(s: str, h: float = 0.2, color=INK, **kw):
    return text(s, color=color, **kw).scale(_cap_scale("text", h))


def right_at(m, x: float, y: float):
    m.move_to([x - m.width / 2, y, 0.0])
    return m


def left_at(m, x: float, y: float):
    m.move_to([x + m.width / 2, y, 0.0])
    return m


def cms_open_data(h: float = 0.25) -> VGroup:
    cms = txt("CMS", h, INK, weight=BOLD)
    od = txt("Open Data", h, INK, slant=ITALIC)
    od.next_to(cms, RIGHT, buff=0.13)
    od.shift(UP * (cms.get_bottom()[1] - od[0].get_bottom()[1]))      # baseline on the "O", not the "p"
    return VGroup(cms, od)


def lumi_label(h: float = 0.21) -> VGroup:
    a = txt("16.4 fb", h)
    sup = tex(r"-1", 0.6 * h).next_to(a, RIGHT, buff=0.03).align_to(a, UP).shift(UP * 0.06 * h / 0.21)
    b = txt("(13 TeV)", h).next_to(sup, RIGHT, buff=0.1)
    b.shift(UP * (a.get_bottom()[1] - b[1].get_bottom()[1]))
    return VGroup(a, sup, b)


def grey_of(hexc: str, lift: float = 0.45) -> str:
    """The grey of the same lightness, lifted toward white: a layer 'switched off' but still readable."""
    c = col(hexc).to_hex()
    r, g, b = (int(c[i:i + 2], 16) for i in (1, 3, 5))
    y = (0.299 * r + 0.587 * g + 0.114 * b) / 255.0
    return tint(mix(BLACK, WHITE, y), lift)


def restack(scene, st: dict, keys, order=None) -> None:
    """Replace whatever is on screen by the builder's objects ``keys``, stacked in ``order``
    (default: the order of ``keys``); same frame, no leftovers."""
    scene.clear()
    order = list(keys) if order is None else order
    scene.add(*[st[k] for k in order if k in keys and k in st])


# ---------------------------------------------------------------------------
# slides 1 and 4: summary.png
# ---------------------------------------------------------------------------

X_LO, X_HI = D["summary"]["x_range_pb"]
XL, XR, YB, YT = -3.65, 4.3, -2.85, 2.25
ROW_ORDER = ("ee", "mumu", "tautau", "comb", "cms", "atlas")
SIZE = {4: 1.0, 5: 0.92, 6: 0.85}
ORDER_SUM = ("grid", "band", "pline", "spines", "xaxis", "header_l", "header_r", "legend",
             "row_ee", "row_mumu", "row_tautau", "sep", "row_comb", "row_cms", "row_atlas")
# The corner copy of summary.png sits as high as the frame allows -- it may enter the title band here
# (deck owner, 17 Sep): it stays on screen from comb_orthogonality through comb_impacts into comb_published,
# so the three slides are one continuous shot.
CORNER_SCALE, CORNER_RIGHT, CORNER_TOP = 0.2, 6.95, 3.88


def sx(v: float) -> float:
    return XL + (float(v) - X_LO) / (X_HI - X_LO) * (XR - XL)


def row_y(n: int, i: int) -> float:
    """Row i (0 = top) of n: matplotlib's ylim (-0.7, n - 0.3) on the axes height."""
    return YB + (n - 1 - i + 0.7) * (YT - YB) / n


def result_row(y, value, up, down, colour, label, value_mob, f=1.0, own=True) -> VGroup:
    cc = col(colour)
    sw = 4.6 if own else 3.2
    tick = Line([XL - 0.09, y, 0], [XL, y, 0], stroke_color=col(INK), stroke_width=2.0)
    right_at(label, XL - 0.22, y)
    x0, x1 = sx(value - down), sx(value + up)
    bar = Line([x0, y, 0], [x1, y, 0], stroke_color=cc, stroke_width=sw)
    ch = 0.1 * f
    caps = VGroup(Line([x0, y - ch, 0], [x0, y + ch, 0], stroke_color=cc, stroke_width=sw),
                  Line([x1, y - ch, 0], [x1, y + ch, 0], stroke_color=cc, stroke_width=sw))
    if own:
        marker = Dot([sx(value), y, 0], radius=0.105 * f, color=cc)
    else:
        marker = Square(side_length=0.17 * f, fill_color=cc, fill_opacity=1.0, stroke_width=0).move_to([sx(value), y, 0])
    left_at(value_mob, XR + 0.3, y)
    return VGroup(tick, label, bar, caps, marker, value_mob)


def summary_state(n: int = 4) -> dict:
    ink = col(INK)
    f = SIZE[n]
    st = {}
    st["grid"] = VGroup(*[DashedLine([sx(v), YB, 0], [sx(v), YT, 0], dash_length=0.03, dashed_ratio=0.45,
                                     stroke_color=col(LIGHT_GREY), stroke_width=1.3) for v in range(1700, 2401, 100)])
    lo, hi = sx(PRED["band_lo_pb"]), sx(PRED["band_hi_pb"])
    st["band"] = Rectangle(width=hi - lo, height=YT - YB, stroke_width=0, fill_color=col(F["prediction"]),
                           fill_opacity=0.18).move_to([(lo + hi) / 2, (YT + YB) / 2, 0])
    st["pline"] = Line([sx(PRED["sigma_pb"]), YB, 0], [sx(PRED["sigma_pb"]), YT, 0], stroke_color=col(F["prediction"]),
                       stroke_width=3.2)
    st["spines"] = VGroup(Line([XL, YB, 0], [XL, YT, 0], stroke_color=ink, stroke_width=2.4),
                          Line([XL, YB, 0], [XR, YB, 0], stroke_color=ink, stroke_width=2.4))
    ticks, labels = VGroup(), VGroup()
    for v in range(1700, 2401, 100):
        ticks.add(Line([sx(v), YB, 0], [sx(v), YB - 0.09, 0], stroke_color=ink, stroke_width=2.0))
        labels.add(txt(f"{v}", 0.19).move_to([sx(v), YB - 0.33, 0]))
    title = tex(r"\sigma(\mathrm{pp}\to Z/\gamma^{*}\to\ell\ell),\ 60<m_{\ell\ell}<120\ \mathrm{GeV}\quad[\mathrm{pb}]", 0.2)
    title.move_to([(XL + XR) / 2, YB - 0.82, 0])
    st["xaxis"] = VGroup(ticks, labels, title)
    hl = cms_open_data(0.25)
    st["header_l"] = hl.move_to([XL + hl.width / 2, YT + 0.1 + hl.height / 2, 0])
    hr = lumi_label(0.2)
    st["header_r"] = hr.move_to([XR - hr.width / 2, YT + 0.1 + hr.height / 2, 0])

    patch = Rectangle(width=0.45, height=0.19, stroke_width=0, fill_color=col(F["prediction"]), fill_opacity=0.18)
    leg_lab = tex(r"\mathrm{aMC@NLO}", 0.16)
    y_leg = YB + (0.25 if n == 4 else 0.17)
    if n == 4:                       # summary.png: lower left, clear of the axis corner
        patch.move_to([XL + 1.35 + 0.225, y_leg, 0])
        leg_lab.next_to(patch, RIGHT, buff=0.18)
    else:                            # combined_vs_published.png: lower right
        right_at(leg_lab, XR - 0.2, y_leg)
        patch.next_to(leg_lab, LEFT, buff=0.18)
    st["legend"] = VGroup(patch, leg_lab)

    specs = {
        "ee": (ROWS["ee"], F["ee"], tex(r"Z\to ee", 0.25 * f, F["ee"])),
        "mumu": (ROWS["mumu"], F["mumu"], tex(r"Z\to\mu\mu", 0.25 * f, F["mumu"])),
        "tautau": (ROWS["tautau"], F["tautau"], tex(r"Z\to\tau\tau", 0.25 * f, F["tautau"])),
        "comb": (ROWS["combined"], F["combined"], txt("Combined", 0.22 * f, F["combined"])),
    }
    keys = ROW_ORDER[:n]
    for i, k in enumerate(keys):
        y = row_y(n, i)
        if k in specs:
            r, c, lab = specs[k]
            st[f"row_{k}"] = result_row(y, r["sigma_pb"], r["up_pb"], r["down_pb"], c, lab,
                                        tex(r["value_tex"], 0.2 * f, c), f)
        else:
            p = PUB["CMS" if k == "cms" else "ATLAS"]
            detail = (r"206\ \mathrm{pb}^{-1},\ \mathrm{ee}+\mu\mu" if k == "cms"
                      else r"81\ \mathrm{pb}^{-1},\ \mathrm{ee}+\mu\mu\ (66\text{--}116\ \mathrm{GeV})")
            name = txt(p["label"], 0.2 * f)
            det = tex(detail, (0.17 if k == "cms" else 0.15) * f)
            lab = VGroup(name, det)
            det.next_to(name, DOWN, buff=0.1 * f)
            det.align_to(name, RIGHT)
            st[f"row_{k}"] = result_row(y, p["sigma_pb"], p["err_pb"], p["err_pb"], F["published"], lab,
                                        tex(p["value_tex"], 0.19 * f, F["published_text"]), f, own=False)
    y_sep = (row_y(n, 2) + row_y(n, 3)) / 2
    st["sep"] = Line([XL, y_sep, 0], [XR, y_sep, 0], stroke_color=col(F["separator"]), stroke_width=1.8)
    return st


def summary_group(n: int = 4, small: bool = False) -> VGroup:
    st = summary_state(n)
    g = VGroup(*[st[k] for k in ORDER_SUM if k in st])
    if small:
        for m in g.family_members_with_points():          # every stroke shrinks with the plot, each by itself
            m.set_stroke(width=m.get_stroke_width() * CORNER_SCALE, family=False)
        g.scale(CORNER_SCALE, about_point=g.get_center())
        g.move_to([CORNER_RIGHT - g.width / 2, CORNER_TOP - g.height / 2, 0])
    return g


def row_in(scene, row: VGroup, run_time: float = 0.95, extra=()) -> None:
    """A result row arrives: label, tick and point appear, the error bar opens from the point, the value appears."""
    tick, lab, bar, caps, marker, val = row
    xm = marker.get_center()
    bar0 = Line(xm + LEFT * 1e-3, xm + RIGHT * 1e-3, stroke_color=bar.get_stroke_color(), stroke_width=bar.get_stroke_width())
    caps0 = VGroup(*[c.copy().move_to([xm[0], c.get_center()[1], 0]) for c in caps])
    scene.play(FadeIn(tick), FadeIn(lab), GrowFromCenter(marker), *extra, run_time=0.4 * run_time)
    scene.add(bar0, caps0)
    scene.bring_to_front(marker)
    scene.play(Transform(bar0, bar), Transform(caps0, caps), FadeIn(val), run_time=0.6 * run_time, rate_func=EASE)
    scene.remove(bar0, caps0)


# ---------------------------------------------------------------------------
# slide 2: the two e mu plots (schematic)
# ---------------------------------------------------------------------------

MAIN_Y, RATIO_Y = (-1.05, 2.0), (-2.6, -1.3)


def _peak(c, mu=90.0, s_lo=15.0, s_hi=20.0):
    s = np.where(c < mu, s_lo, s_hi)
    return np.exp(-0.5 * ((c - mu) / s) ** 2)


TT_EDGES = np.arange(40.0, 161.0, 10.0)
_c = 0.5 * (TT_EDGES[:-1] + TT_EDGES[1:])
# bottom-up; (key, counts, original colour, legend tex). Shapes only: no number of this plot is printed.
TT_LAYERS = [
    ("tt_emu", np.array([300, 330, 360, 390, 420, 440, 460, 480, 500, 510, 520, 520.0]), ZTAUTAU_FIG["TTbar"], r"t\bar{t}\to e\mu"),
    ("sig_emu", 2200 * _peak(_c), shade(ZTAUTAU_FIG["DYtautau"], 0.35), r"Z\to\tau\tau\to e\mu"),
    ("fakes", 2600 * np.exp(-(_c - 40.0) / 60.0), ZTAUTAU_FIG["Fakes"], r"\mathrm{jet}\to\tau_h"),
    ("sig_etau", 1800 * _peak(_c), shade(ZTAUTAU_FIG["DYtautau"], 0.18), r"Z\to\tau\tau\to e\tau_h"),
    ("sig_mutau", 3200 * _peak(_c), ZTAUTAU_FIG["DYtautau"], r"Z\to\tau\tau\to\mu\tau_h"),
    ("sig_tautau", 2000 * _peak(_c), tint(ZTAUTAU_FIG["DYtautau"], 0.45), r"Z\to\tau\tau\to\tau_h\tau_h"),
]
TT_RATIO = np.array([1.02, 0.97, 1.01, 0.99, 1.02, 0.985, 1.01, 1.03, 0.97, 1.0, 1.04, 0.96])
TT_EMU_KEYS = {"tt_emu": ZMUMU_FIG["TTbar"], "sig_emu": ZMUMU_FIG["DYtautau"]}     # the e mu part in the CR's colours

CR_EDGES = np.arange(60.0, 121.0, 5.0)
CR_LAYERS = [      # as datamc_CRemu_mass_fit.png, bottom-up
    ("WJets", np.array([480, 440, 400, 420, 330, 320, 340, 250, 250, 300, 330, 220.0]), ZMUMU_FIG["WJets"], r"W+\mathrm{jets}"),
    ("WW", np.array([480, 470, 420, 480, 470, 480, 480, 440, 440, 440, 430, 350.0]), ZMUMU_FIG["WW"], r"WW"),
    ("WZ", np.full(12, 25.0), ZMUMU_FIG["WZ"], r"WZ"),
    ("ZZ", np.full(12, 20.0), ZMUMU_FIG["ZZ"], r"ZZ"),
    ("SingleTop", np.array([300, 300, 330, 350, 350, 350, 350, 330, 340, 350, 330, 300.0]), ZMUMU_FIG["SingleTop"], r"tW"),
    ("TTbar", np.array([3400, 3550, 3700, 3800, 3850, 3850, 3800, 3700, 3600, 3500, 3400, 3200.0]), ZMUMU_FIG["TTbar"], r"t\bar{t}"),
    ("DYee", np.full(12, 15.0), ZMUMU_FIG["DYee"], r"Z/\gamma^{*}\to ee"),
    ("DYtautau", np.array([4800, 3950, 2650, 1400, 620, 260, 130, 60, 40, 30, 20, 15.0]), ZMUMU_FIG["DYtautau"], r"Z/\gamma^{*}\to\tau\tau"),
    ("DYmumu", np.array([200, 280, 430, 400, 400, 330, 160, 50, 30, 25, 15, 10.0]), ZMUMU_FIG["DYmumu"], r"Z/\gamma^{*}\to\mu\mu"),
]
CR_RATIO = np.array([0.965, 0.985, 1.03, 1.015, 1.02, 1.01, 0.965, 1.03, 1.025, 0.945, 0.915, 0.99])

TT_BOX = dict(x0=-5.3, x1=-1.85, edges=TT_EDGES, ymax=13500.0, xticks=(40, 60, 80, 100, 120, 140, 160),
              xlabels=(40, 80, 120, 160), xtitle=r"m_{\tau\tau}\ [\mathrm{GeV}]", title=r"Z\to\tau\tau",
              title_colour=F["tautau"], legend=(-1.62, 1.9, 0.29), region=None)
CR_BOX = dict(x0=1.15, x1=4.6, edges=CR_EDGES, ymax=13000.0, xticks=tuple(range(60, 121, 10)),
              xlabels=(60, 80, 100, 120), xtitle=r"m_{e\mu}\ [\mathrm{GeV}]", title=r"Z\to\mu\mu",
              title_colour=F["mumu"], legend=(4.78, 1.12, 0.25), region=r"e\mu")


def mass_plot(box: dict, layers, ratio, mode: str) -> VGroup:
    """One schematic data/prediction plot with a ratio panel and a legend beside it.
    mode: 'orig' (the figure's colours), 'emu' (tautau: e mu part in the CR colours, the rest grey), 'grey' (all grey).
    Every mode has the same geometry and structure, so a Transform between two modes only recolours."""
    x0, x1, edges = box["x0"], box["x1"], box["edges"]
    lo, hi = float(edges[0]), float(edges[-1])
    grey_all = mode == "grey"
    ink = col(GREY) if grey_all else col(INK)
    data_c = col(GREY) if mode in ("grey", "emu") else col(INK)

    def px(v):
        return x0 + (float(v) - lo) / (hi - lo) * (x1 - x0)

    def py(v):
        return MAIN_Y[0] + float(v) / box["ymax"] * (MAIN_Y[1] - MAIN_Y[0])

    def pr(r):
        return RATIO_Y[0] + (float(r) - 0.8) / 0.4 * (RATIO_Y[1] - RATIO_Y[0])

    def colour_of(key, orig):
        if mode == "grey":
            return grey_of(orig)
        if mode == "emu":
            return TT_EMU_KEYS[key] if key in TT_EMU_KEYS else grey_of(orig)
        return orig

    cum = np.zeros(len(edges) - 1)
    stack = VGroup()
    swatch_colour = {}
    for key, counts, orig, _ in layers:
        lower, upper = cum.copy(), cum + counts
        top = []
        for i in range(len(counts)):
            top += [[px(edges[i]), py(upper[i]), 0], [px(edges[i + 1]), py(upper[i]), 0]]
        bot = []
        for i in reversed(range(len(counts))):
            bot += [[px(edges[i + 1]), py(lower[i]), 0], [px(edges[i]), py(lower[i]), 0]]
        c = colour_of(key, orig)
        swatch_colour[key] = c
        stack.add(Polygon(*top, *bot, fill_color=col(c), fill_opacity=1.0, stroke_color=darken(c, 0.35), stroke_width=0.8))
        cum = upper
    total = cum
    centres = 0.5 * (edges[:-1] + edges[1:])
    mcstat = VGroup()
    for i in range(len(total)):
        mcstat.add(Rectangle(width=px(edges[i + 1]) - px(edges[i]), height=py(total[i] * 1.025) - py(total[i] * 0.975),
                             stroke_width=0, fill_color=col(GREY), fill_opacity=0.35)
                   .move_to([px(centres[i]), py(total[i]), 0]))
    data = VGroup()
    for i, c in enumerate(centres):
        y = py(total[i] * ratio[i])
        data.add(VGroup(Line([px(c), y - 0.05, 0], [px(c), y + 0.05, 0], stroke_color=data_c, stroke_width=2.0),
                        Dot([px(c), y, 0], radius=0.04, color=data_c)))

    frame_main = Rectangle(width=x1 - x0, height=MAIN_Y[1] - MAIN_Y[0], stroke_color=ink, stroke_width=2.2,
                           fill_opacity=0).move_to([(x0 + x1) / 2, sum(MAIN_Y) / 2, 0])
    frame_ratio = Rectangle(width=x1 - x0, height=RATIO_Y[1] - RATIO_Y[0], stroke_color=ink, stroke_width=2.2,
                            fill_opacity=0).move_to([(x0 + x1) / 2, sum(RATIO_Y) / 2, 0])
    rband = Rectangle(width=x1 - x0, height=pr(1.03) - pr(0.97), stroke_width=0, fill_color=col(GREY),
                      fill_opacity=0.3).move_to([(x0 + x1) / 2, pr(1.0), 0])
    rline = DashedLine([x0, pr(1.0), 0], [x1, pr(1.0), 0], dash_length=0.08, stroke_color=col(GREY), stroke_width=1.8)
    rdots = VGroup(*[VGroup(Line([px(c), pr(r) - 0.045, 0], [px(c), pr(r) + 0.045, 0], stroke_color=data_c, stroke_width=2.0),
                            Dot([px(c), pr(r), 0], radius=0.04, color=data_c)) for c, r in zip(centres, ratio)])
    ticks = VGroup()
    for v in box["xticks"]:
        major = v in box["xlabels"]
        for yb in (MAIN_Y[0], RATIO_Y[0]):
            ticks.add(Line([px(v), yb, 0], [px(v), yb + (0.12 if major else 0.07), 0], stroke_color=ink, stroke_width=1.8))
    xlabels = VGroup(*[txt(f"{v}", 0.16, ink).move_to([px(v), RATIO_Y[0] - 0.24, 0]) for v in box["xlabels"]])
    xtitle = tex(box["xtitle"], 0.19, ink)
    xtitle.move_to([x1 - xtitle.width / 2, RATIO_Y[0] - 0.68, 0])
    ytitle = txt("Events", 0.15, ink).rotate(np.pi / 2)
    ytitle.move_to([x0 - 0.3, MAIN_Y[1] - ytitle.height / 2 - 0.05, 0])
    rtitle = tex(r"\mathrm{Data/pred.}", 0.12, ink).rotate(np.pi / 2).move_to([x0 - 0.3, sum(RATIO_Y) / 2, 0])
    title = tex(box["title"], 0.26, GREY if grey_all else box["title_colour"]).move_to([(x0 + x1) / 2, 2.32, 0])
    region = VGroup()
    if box["region"]:
        region.add(tex(box["region"], 0.17, ink))
        region.move_to([x0 + 0.15 + region.width / 2, MAIN_Y[1] - 0.25, 0])

    lx, ly, pitch = box["legend"]
    legend = VGroup()
    entries = [("Data", None)] + [(lab, key) for key, _, _, lab in reversed(layers)]
    if box is CR_BOX:
        entries.append((r"\mathrm{MC\ stat.}", "mcstat"))
    for j, (lab, key) in enumerate(entries):
        y = ly - j * pitch
        if key is None:
            sw = VGroup(Line([lx + 0.13, y - 0.07, 0], [lx + 0.13, y + 0.07, 0], stroke_color=data_c, stroke_width=2.0),
                        Dot([lx + 0.13, y, 0], radius=0.04, color=data_c))
            label = txt("Data", 0.13, ink)
        elif key == "mcstat":
            sw = Rectangle(width=0.26, height=0.16, stroke_width=0, fill_color=col(GREY), fill_opacity=0.35).move_to([lx + 0.13, y, 0])
            label = tex(lab, 0.13, ink)
        else:
            c = swatch_colour[key]
            sw = Rectangle(width=0.26, height=0.16, fill_color=col(c), fill_opacity=1.0, stroke_color=darken(c, 0.35),
                           stroke_width=0.8).move_to([lx + 0.13, y, 0])
            label = tex(lab, 0.13, ink)
        left_at(label, lx + 0.4, y)
        legend.add(VGroup(sw, label))

    g = VGroup(stack, mcstat, frame_main, frame_ratio, rband, rline, rdots, data, ticks, xlabels, xtitle, ytitle, rtitle,
               title, region, legend)
    g.stack, g.mcstat, g.data, g.rdots, g.legend = stack, mcstat, data, rdots, legend
    g.axes = VGroup(frame_main, frame_ratio, rband, rline, ticks, xlabels, xtitle, ytitle, rtitle, title, region)
    return g


def tautau_plot(mode: str = "orig") -> VGroup:
    return mass_plot(TT_BOX, TT_LAYERS, TT_RATIO, mode)


def cr_plot(mode: str = "orig") -> VGroup:
    return mass_plot(CR_BOX, CR_LAYERS, CR_RATIO, mode)


def orth_state() -> dict:
    return {"summary": summary_group(4, small=True), "tt": tautau_plot("orig"), "cr": cr_plot("grey")}


ORDER_ORTH = ("summary", "tt", "cr")


def plot_in(scene, p: VGroup) -> None:
    """A schematic plot appears: frame and labels, the stack rises from the axis, data and legend."""
    base = p.stack.get_bottom()[1]
    stack0 = p.stack.copy().stretch(0.01, 1, about_point=[0, base, 0])
    scene.play(FadeIn(p.axes), run_time=0.5)
    scene.play(Transform(stack0, p.stack), FadeIn(p.legend), run_time=0.9, rate_func=EASE)
    scene.remove(stack0)
    scene.add(p.stack)
    scene.play(FadeIn(p.mcstat), FadeIn(p.data), FadeIn(p.rdots), run_time=0.4)


# ---------------------------------------------------------------------------
# slide 3: impacts.png
# ---------------------------------------------------------------------------

N_IMP = len(IMP)
LIM = D["impacts"]["x_lim_pb"]
P_LO, P_HI = D["impacts"]["pull_x_range"]
IYT, IYB = 2.00, -2.95     # top left clear of the corner summary that stays on screen from slide 2
PXL, PXR = -3.2, 0.4
IXL, IXR = 0.8, 5.3
TOP3 = (0, 1, 2)     # rows with the three largest impacts (Lumi, ElectronID, mu_ttbar in the final combination)
ORDER_IMP = ("hdr_l", "hdr_r", "p_bands", "p_grid", "i_grid", "i_bars", "i_axes", "i_xaxis", "i_legend", "labels",
             "p_axes", "p_xaxis", "p_points", "boxes", "corner")
assert abs(LIM - 1.35 * 22.54) < 0.05 and (P_LO, P_HI) == (-2.6, 2.6)


def iy(i: int) -> float:
    return IYB + (N_IMP - 1 - i + 0.7) * (IYT - IYB) / N_IMP


def px_(v: float) -> float:
    return PXL + (float(v) - P_LO) / (P_HI - P_LO) * (PXR - PXL)


def ix(v: float) -> float:
    return IXL + (float(v) + LIM) / (2 * LIM) * (IXR - IXL)


ROW_P = (IYT - IYB) / N_IMP


def impact_bar(i: int, value: float, colour, grown: bool = True, opacity: float = 0.85) -> Rectangle:
    x0 = ix(0.0)
    x1 = ix(value) if grown else x0 + (1e-3 if value >= 0 else -1e-3)     # not grown: a sliver on the zero line
    return Rectangle(width=abs(x1 - x0), height=0.72 * ROW_P, stroke_width=0, fill_color=col(colour),
                     fill_opacity=opacity).move_to([(x0 + x1) / 2, iy(i), 0])


def np_label(r: dict, opacity: float = 1.0):
    if r["name"] == "mu_ttbar":
        m = tex(r"\mu_{t\bar{t}}\ (\tau\tau,\ \mathrm{free})", 0.14)
    else:
        assert "$" not in r["label_mpl"], r["label_mpl"]
        m = txt(r["label_mpl"], 0.135)
    right_at(m, PXL - 0.17, iy(IMP.index(r)))
    return m.set_opacity(opacity)


def pull_point(i: int, r: dict, opacity: float = 1.0) -> VGroup:
    y = iy(i)
    if r["free"]:
        return VGroup(txt(r["free_text"], 0.115, F["published_text"]).move_to([px_(0.0), y, 0]).set_opacity(opacity))
    bar = Line([px_(r["pull"] - r["err_down"]), y, 0], [px_(r["pull"] + r["err_up"]), y, 0], stroke_color=col(INK),
               stroke_width=2.4, stroke_opacity=opacity)
    dot = Dot([px_(r["pull"]), y, 0], radius=0.048, color=col(INK), fill_opacity=opacity)
    return VGroup(bar, dot)


def impacts_state(dim: bool = False, boxes: bool = False) -> dict:
    ink = col(INK)
    st = {"corner": summary_group(4, small=True)}     # the combination result stays in the corner, from slide 2
    hl = cms_open_data(0.24)
    st["hdr_l"] = hl.move_to([PXL + hl.width / 2, IYT + 0.12 + hl.height / 2, 0])
    hr = lumi_label(0.2)
    st["hdr_r"] = hr.move_to([IXR - hr.width / 2, IYT + 0.12 + hr.height / 2, 0])

    def op(i, full, faded):
        return faded if (dim and i not in TOP3) else full

    # pulls panel
    yellow = Rectangle(width=px_(2) - px_(-2), height=IYT - IYB, stroke_width=0, fill_color=col(F["pull_2sigma"]),
                       fill_opacity=0.45).move_to([px_(0), (IYT + IYB) / 2, 0])
    green = Rectangle(width=px_(1) - px_(-1), height=IYT - IYB, stroke_width=0, fill_color=col(F["pull_1sigma"]),
                      fill_opacity=0.45).move_to([px_(0), (IYT + IYB) / 2, 0])
    st["p_bands"] = VGroup(yellow, green)
    st["p_grid"] = VGroup(*[DashedLine([px_(v), IYB, 0], [px_(v), IYT, 0], dash_length=0.03, dashed_ratio=0.45,
                                       stroke_color=col(LIGHT_GREY), stroke_width=1.2) for v in (-2, -1, 1, 2)])
    yt = VGroup(*[Line([PXL - 0.07, iy(i), 0], [PXL, iy(i), 0], stroke_color=ink, stroke_width=1.6) for i in range(N_IMP)])
    st["p_axes"] = VGroup(Line([PXL, IYB, 0], [PXL, IYT, 0], stroke_color=ink, stroke_width=2.2),
                          Line([PXL, IYB, 0], [PXR, IYB, 0], stroke_color=ink, stroke_width=2.2),
                          Line([px_(0), IYB, 0], [px_(0), IYT, 0], stroke_color=ink, stroke_width=1.8), yt)
    pt, pl = VGroup(), VGroup()
    for v in (-2, -1, 0, 1, 2):
        pt.add(Line([px_(v), IYB, 0], [px_(v), IYB - 0.08, 0], stroke_color=ink, stroke_width=1.8))
        pl.add(txt(f"{v}" if v >= 0 else f"−{-v}", 0.16).move_to([px_(v), IYB - 0.3, 0]))
    ptitle = tex(r"(\hat{\theta}-\theta_0)/\Delta\theta", 0.19).move_to([(PXL + PXR) / 2, IYB - 0.78, 0])
    st["p_xaxis"] = VGroup(pt, pl, ptitle)
    st["p_points"] = VGroup(*[pull_point(i, r, op(i, 1.0, 0.22)) for i, r in enumerate(IMP)])

    # impacts panel
    st["i_grid"] = VGroup(*[DashedLine([ix(v), IYB, 0], [ix(v), IYT, 0], dash_length=0.03, dashed_ratio=0.45,
                                       stroke_color=col(LIGHT_GREY), stroke_width=1.2) for v in (-30, -20, -10, 10, 20, 30)])
    st["i_bars"] = VGroup(*[VGroup(impact_bar(i, r["impact_up_pb"], F["impact_up"], opacity=op(i, 0.85, 0.18)),
                                   impact_bar(i, r["impact_down_pb"], F["impact_down"], opacity=op(i, 0.85, 0.18)))
                            for i, r in enumerate(IMP)])
    st["i_axes"] = VGroup(Line([IXL, IYB, 0], [IXL, IYT, 0], stroke_color=ink, stroke_width=2.2),
                          Line([IXL, IYB, 0], [IXR, IYB, 0], stroke_color=ink, stroke_width=2.2),
                          Line([ix(0), IYB, 0], [ix(0), IYT, 0], stroke_color=ink, stroke_width=1.8))
    it, il = VGroup(), VGroup()
    for v in (-30, -20, -10, 0, 10, 20, 30):
        it.add(Line([ix(v), IYB, 0], [ix(v), IYB - 0.08, 0], stroke_color=ink, stroke_width=1.8))
        il.add(txt(f"{v}" if v >= 0 else f"−{-v}", 0.16).move_to([ix(v), IYB - 0.3, 0]))
    ititle = tex(r"\Delta\sigma\quad[\mathrm{pb}]", 0.19).move_to([(IXL + IXR) / 2, IYB - 0.78, 0])
    st["i_xaxis"] = VGroup(it, il, ititle)
    leg = VGroup()
    for j, (c, lab) in enumerate(((F["impact_up"], r"\hat{\theta}+\Delta\hat{\theta}"),
                                  (F["impact_down"], r"\hat{\theta}-\Delta\hat{\theta}"))):
        y = iy(N_IMP - 2 + j)
        lab_m = tex(lab, 0.14)
        right_at(lab_m, IXR - 0.12, y)
        patch = Rectangle(width=0.36, height=0.14, stroke_width=0, fill_color=col(c), fill_opacity=0.85)
        patch.next_to(lab_m, LEFT, buff=0.14)
        leg.add(VGroup(patch, lab_m))
    st["i_legend"] = leg
    st["labels"] = VGroup(*[np_label(r, op(i, 1.0, 0.3)) for i, r in enumerate(IMP)])

    bx = VGroup()
    if boxes:
        y_top, y_bot = iy(TOP3[0]) + ROW_P / 2 + 0.03, iy(TOP3[-1]) - ROW_P / 2 - 0.03
        for a, b in ((PXL, PXR), (IXL, IXR)):
            bx.add(RoundedRectangle(width=b - a + 0.1, height=y_top - y_bot, corner_radius=0.07, stroke_color=col(HIGHLIGHT),
                                    stroke_width=3.5, fill_opacity=0).move_to([(a + b) / 2, (y_top + y_bot) / 2, 0]))
    st["boxes"] = bx
    return st


# ---------------------------------------------------------------------------
# scenes
# ---------------------------------------------------------------------------

class CombSummary(Scene):
    """comb_summary: summary.png from scratch."""

    def construct(self):
        white_background(self)
        st = summary_state(4)
        shown = []
        restack(self, st, shown, ORDER_SUM)
        self.wait(0.3)

        left, bottom = st["spines"]
        self.play(Create(left), Create(bottom), run_time=0.9, rate_func=EASE)
        # the x labels and the "CMS Open Data" of the figure arrive with them, not before the clip starts
        self.play(FadeIn(st["xaxis"]), FadeIn(st["grid"]), FadeIn(st["header_l"]), run_time=0.45)
        shown += ["spines", "xaxis", "grid", "header_l"]
        restack(self, st, shown, ORDER_SUM)
        self.wait(0.4)

        for k in ("ee", "mumu", "tautau"):
            row_in(self, st[f"row_{k}"])
            shown.append(f"row_{k}")
            restack(self, st, shown, ORDER_SUM)
            self.wait(0.45)

        self.play(Create(st["sep"]), run_time=0.5)
        shown.append("sep")
        row_in(self, st["row_comb"], extra=(FadeIn(st["header_r"]),))
        shown += ["row_comb", "header_r"]
        restack(self, st, shown, ORDER_SUM)
        self.wait(0.6)

        # the prediction: the line rises, the band opens, the legend appears (all behind the points)
        pl = st["pline"]
        pl0 = Line(pl.get_start(), pl.get_start() + UP * 1e-3, stroke_color=pl.get_stroke_color(), stroke_width=pl.get_stroke_width())
        self.bring_to_back(pl0)
        self.bring_to_back(st["grid"])
        self.play(Transform(pl0, pl), run_time=0.7, rate_func=EASE)
        shown.append("pline")
        restack(self, st, shown, ORDER_SUM)
        band = st["band"]
        band0 = band.copy().stretch(0.002, 0, about_point=pl.get_center())
        self.bring_to_back(band0)
        self.bring_to_back(st["grid"])
        self.play(Transform(band0, band), FadeIn(st["legend"]), run_time=0.8, rate_func=EASE)
        shown += ["band", "legend"]
        restack(self, st, shown, ORDER_SUM)
        assert set(shown) == set(st) - {"row_cms", "row_atlas"}
        self.wait(0.1)


class CombOrthogonality(Scene):
    """comb_orthogonality: opens on comb_summary's last frame."""

    def construct(self):
        white_background(self)
        s1 = summary_state(4)
        restack(self, s1, ORDER_SUM)
        self.wait(0.3)
        end = orth_state()

        g = VGroup(*[s1[k] for k in ORDER_SUM if k in s1])
        self.clear()
        self.add(g)
        self.play(Transform(g, summary_group(4, small=True)), run_time=1.3, rate_func=EASE)
        self.wait(0.3)

        tt, cr = tautau_plot("orig"), cr_plot("orig")
        plot_in(self, tt)
        self.clear()
        self.add(g, tt)
        plot_in(self, cr)
        self.clear()
        self.add(g, tt, cr)
        self.wait(1.0)

        # the e mu events of the tautau measurement are the events of the mumu control region
        self.play(Transform(tt, tautau_plot("emu")), run_time=1.6, rate_func=EASE)
        self.wait(1.2)
        # the control region leaves the combination; the tautau measurement keeps them
        self.play(Transform(cr, cr_plot("grey")), Transform(tt, tautau_plot("orig")), run_time=1.6, rate_func=EASE)
        restack(self, end, ORDER_ORTH)
        self.wait(0.1)


class CombImpacts(Scene):
    """comb_impacts: opens on comb_orthogonality's last frame; both panels at once, row by row from the
    smallest impact upward; the top three outlined."""

    def construct(self):
        white_background(self)
        o = orth_state()
        restack(self, o, ORDER_ORTH)
        self.wait(0.3)

        st = impacts_state()
        shown = ["corner"]
        # the two e mu plots give way; the combination result stays in the corner where slide 2 left it
        self.play(FadeOut(o["tt"]), FadeOut(o["cr"]), run_time=0.7)
        restack(self, st, shown, ORDER_IMP)
        self.wait(0.2)

        # both panels open together: the two frames, then their axes, bands, grids and the legend
        ileft, ibottom, izero = st["i_axes"]
        pleft, pbottom, pzero, pticks = st["p_axes"]
        self.play(Create(ileft), Create(ibottom), Create(pleft), Create(pbottom), FadeIn(pticks),
                  FadeIn(st["hdr_l"]), FadeIn(st["hdr_r"]), run_time=0.8, rate_func=EASE)
        yellow, green = st["p_bands"]
        y0 = yellow.copy().stretch(0.002, 0, about_point=[px_(0), 0, 0])
        g0 = green.copy().stretch(0.002, 0, about_point=[px_(0), 0, 0])
        self.bring_to_back(g0)
        self.bring_to_back(y0)
        self.play(Transform(y0, yellow), Transform(g0, green), Create(izero), Create(pzero),
                  FadeIn(st["i_xaxis"]), FadeIn(st["p_xaxis"]), FadeIn(st["i_grid"]), FadeIn(st["p_grid"]),
                  FadeIn(st["i_legend"]), run_time=0.8, rate_func=EASE)
        shown += ["hdr_l", "hdr_r", "i_axes", "i_xaxis", "i_grid", "i_legend", "p_axes", "p_xaxis", "p_bands", "p_grid"]
        restack(self, st, shown, ORDER_IMP)
        self.wait(0.3)

        # one row of the two panels at a time, from the smallest impact (bottom) to the largest (top)
        bars_done, labels_done, pts_done = VGroup(), VGroup(), VGroup()
        st_bars, st_labels, st_pts = st["i_bars"], st["labels"], st["p_points"]
        for i in reversed(range(N_IMP)):
            r = IMP[i]
            up, dn = st_bars[i]
            up0 = impact_bar(i, r["impact_up_pb"], F["impact_up"], grown=False)
            dn0 = impact_bar(i, r["impact_down_pb"], F["impact_down"], grown=False)
            self.add(up0, dn0)
            anims = [Transform(up0, up), Transform(dn0, dn), FadeIn(st_labels[i])]
            p, tmp = st_pts[i], []
            if r["free"]:
                anims.append(FadeIn(p))
            else:
                bar, dot = p
                c = dot.get_center()
                bar0 = Line(c + LEFT * 1e-3, c + RIGHT * 1e-3, stroke_color=bar.get_stroke_color(),
                            stroke_width=bar.get_stroke_width())
                self.add(bar0)
                tmp = [bar0]
                anims += [GrowFromCenter(dot), Transform(bar0, bar)]
            size = max(abs(r["impact_up_pb"]), abs(r["impact_down_pb"]))
            self.play(*anims, run_time=0.3 + 0.015 * size, rate_func=EASE)
            self.remove(up0, dn0, *tmp)
            bars_done.add(st_bars[i])
            labels_done.add(st_labels[i])
            pts_done.add(p)
            self.clear()
            self.add(*[st[k] for k in ORDER_IMP if k in shown and k not in ("i_bars", "labels", "p_points")])
            # z-order: bars under the axes, labels and pull points on top
            self.add(bars_done)
            self.bring_to_front(st["i_axes"], st["i_legend"])
            self.add(labels_done, pts_done, st["corner"])
            self.wait(0.05)
        shown += ["i_bars", "labels", "p_points"]
        restack(self, st, shown, ORDER_IMP)
        self.wait(0.8)

        # the three largest impacts
        dimmed = impacts_state(dim=True, boxes=True)
        self.play(*[Transform(st[k], dimmed[k]) for k in ("i_bars", "labels", "p_points")], run_time=0.8, rate_func=EASE)
        self.play(*[Create(b) for b in dimmed["boxes"]], run_time=0.8, rate_func=EASE)
        restack(self, dimmed, ORDER_IMP)
        self.wait(1.6)


class CombPublished(Scene):
    """comb_published: opens on comb_impacts' last frame; the impacts leave, the corner result grows back."""

    def construct(self):
        white_background(self)
        imp = impacts_state(dim=True, boxes=True)
        restack(self, imp, ORDER_IMP)
        self.wait(0.3)

        self.play(*[FadeOut(imp[k]) for k in ORDER_IMP if k != "corner"], run_time=0.8)
        g = imp["corner"]
        self.play(Transform(g, summary_group(4)), run_time=1.3, rate_func=EASE)
        s4 = summary_state(4)
        restack(self, s4, ORDER_SUM)
        self.wait(0.5)

        moving = ("row_ee", "row_mumu", "row_tautau", "sep", "row_comb", "legend")
        s5 = summary_state(5)
        self.play(*[Transform(s4[k], s5[k]) for k in moving], run_time=1.0, rate_func=EASE)
        restack(self, s5, [k for k in ORDER_SUM if k != "row_cms"])
        row_in(self, s5["row_cms"])
        restack(self, s5, ORDER_SUM)
        self.wait(0.6)

        s6 = summary_state(6)
        self.play(*[Transform(s5[k], s6[k]) for k in moving + ("row_cms",)], run_time=1.0, rate_func=EASE)
        restack(self, s6, [k for k in ORDER_SUM if k != "row_atlas"])
        row_in(self, s6["row_atlas"])
        restack(self, s6, ORDER_SUM)
        self.wait(0.1)
