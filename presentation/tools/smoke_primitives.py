#!/usr/bin/env python
"""Smoke checks for the plotting / event / story primitives in style/bnd_style.py.

Plain python, no render. Run from presentation/ in the Manim env:

    python tools/smoke_primitives.py            # all checks
    python tools/smoke_primitives.py --no-git   # skip the byte-identity check
                                                # against the committed module

Checks (assert, exit 1 on the first failure):
  DataAxes(y_log)   c2p clamps below the floor, decade ticks sit at c2p(x0, 10**k),
                    the frame centre is the axis centre, titles build, no exponent leaks
  DataAxes(linear)  point arrays of axes + data_bar + step_hist + vref_line identical to
                    the committed module (git show HEAD:…)
  stack_hist        8 layers incl. all-zero ones on the log axis, 8 polygons, equal point counts
  ratio_panel       60 dots, x axis aligned under the main plot
  mini_slice        == place(CMSSlice(), (s, c), DET_CENTER) to 1e-9
  kappa_from_pt     kappa_from_pt(45) == 0.28 exactly, clipping
  event_from_json / pair_from_json / muon_in_jet on hand-made records
  colour_key, pull_plot, value_grid(6x4), slider, spine(7), clock, rain
                    build and stay inside x in [-7, 7], y in [-4, 2.7], outside the top-left block
  add_state / check_order on a fake scene (placeholders and trackers skipped)
"""
from __future__ import annotations

import importlib.util
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from manim import DOWN, Group, Mobject, Square, VGroup, ValueTracker  # noqa: E402
from style.bnd_style import *  # noqa: E402,F401,F403
from style.bnd_style import _CAP_H  # noqa: E402

NO_GIT = "--no-git" in sys.argv
X_MIN, X_MAX, Y_MIN, Y_MAX = -7.0, 7.0, -4.0, 2.7
edges60 = np.arange(60.0, 121.0, 1.0)
rng = np.random.default_rng(20260915)


def ok(label, cond, extra=""):
    if not cond:
        print(f"FAIL  {label} {extra}")
        sys.exit(1)
    print(f"ok    {label} {extra}")


def inside(m, label):
    x0, x1 = m.get_left()[0], m.get_right()[0]
    y0, y1 = m.get_bottom()[1], m.get_top()[1]
    ok(f"{label} bbox inside frame", X_MIN <= x0 and x1 <= X_MAX and Y_MIN <= y0 and y1 <= Y_MAX,
       f"x [{x0:.2f}, {x1:.2f}] y [{y0:.2f}, {y1:.2f}]")
    ok(f"{label} bbox outside the top-left block", not (x0 < CORNER_X_MAX and y1 > CORNER_Y_MIN),
       f"x0 {x0:.2f} y1 {y1:.2f}")


# --- 1. log DataAxes ---------------------------------------------------------
dax = DataAxes([60, 120, 10], [0, 7, 1], 6.4, 2.9, y_log=True,
               x_title=r"m_{\mu\mu}\ [\mathrm{GeV}]", y_title=r"\mathrm{events}\,/\,\mathrm{GeV}")
dax.move_frame_to(np.array([1.6, 0.85, 0.0]))
ok("move_frame_to centres the plot area", np.allclose(dax.frame.get_center(), [1.6, 0.85, 0]))
ok("log axis attributes", dax.y_log and dax.y_base == 1.0 and dax.y_top == 1e7 and dax.y_floor == 1.0)
ok("c2p clamps below the floor", dax.c2p(60, 0.5)[1] == dax.c2p(60, 1)[1] == dax.c2p(60, 0)[1])
ok("c2p exact at the top", abs(dax.c2p(60, 1e7)[1] - dax.ax.c2p(60, 1e7)[1]) < 1e-12)
ok("c2p decade spacing uniform",
   np.allclose(np.diff([dax.c2p(60, 10.0 ** k)[1] for k in range(8)]), 2.9 / 7))
ok("8 decade ticks / labels", len(dax.y_ticks_v) == 8 and len(dax.y_labels) == 8)
for k, tk in enumerate(dax.y_ticks_v):
    assert np.allclose(tk.get_start(), dax.c2p(60, 10.0 ** k)), k
ok("decade ticks at c2p(x0, 10**k)", True)
for k, lab in enumerate(dax.y_labels):
    assert abs(lab[0][0].get_center()[1] - dax.c2p(60, 10.0 ** k)[1]) < 1e-9, k
ok("decade labels' digits centred on the ticks", True)
ok("decade labels are 10^{k}", [str(lab.tex_string) for lab in dax.y_labels] == [f"10^{{{k}}}" for k in range(8)])
ok("frame centre = axis centre", np.allclose(dax.frame.get_center(), dax.c2p(90, 10 ** 3.5)))
ok("x ticks on the axis line", all(abs(t.get_start()[1] - dax.c2p(60, 1)[1]) < 1e-9 for t in dax.x_ticks_v))
ok("x title centred", abs(dax.x_title.get_center()[0] - dax.c2p(90, 1)[0]) < 1e-9)
ok("y title centred", abs(dax.y_title.get_center()[1] - dax.c2p(60, 10 ** 3.5)[1]) < 1e-9)
d2 = DataAxes([60, 120, 10], [0, 7, 1], 6.4, 2.9, y_log=True, y_ticks=[0, 2, 4, 6])
ok("y_ticks as exponents", len(d2.y_labels) == 4 and str(d2.y_labels[-1].tex_string) == "10^{6}")
d3 = DataAxes([60, 120, 10], [0, 7, 1], 6.4, 2.9, y_log=True, y_floor=10.0)
ok("y_floor clamps above the base", d3.c2p(60, 3)[1] == d3.c2p(60, 10)[1] > d3.ax.c2p(60, 1)[1])
bar = data_bar(dax, 90.5, 0.0, 0.5, color=SAMPLE["Data"])
ok("data_bar of zero on log: hairline height", bar.height <= 1e-3 + 1e-9)
sh = step_hist(dax, edges60, np.zeros(60), fill_opacity=0.3)
ok("step_hist of zeros on log builds", sh.has_points())
vl = vref_line(dax, 91.19)
ok("vref_line spans base..top", np.allclose(vl.get_top()[1], dax.c2p(91.19, 1e7)[1], atol=0.05))

# --- 2. linear DataAxes byte-identical to the committed module ----------------
if not NO_GIT:
    old_src = subprocess.run(["git", "show", "HEAD:presentation/style/bnd_style.py"],
                             cwd=ROOT.parent, capture_output=True, text=True, check=True).stdout
    tmp = ROOT / "work" / "smoke_old_bnd_style.py"
    tmp.parent.mkdir(exist_ok=True)
    tmp.write_text(old_src)
    spec = importlib.util.spec_from_file_location("old_bnd_style", tmp)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)

    def linear_scene(mod):
        d = mod.DataAxes([60, 120, 10], [0, 1200, 300], 5.0, 2.4, x_title=r"m_{\ell\ell}", y_title=r"N")
        d.move_to(np.array([0.5, -0.4, 0.0]))
        h = mod.schematic_zpeak(edges60, height=1000.0)
        return VGroup(d, mod.data_bar(d, 91.5, 800.0, 0.5), mod.step_hist(d, edges60, h),
                      mod.step_hist(d, edges60, h, fill_opacity=0.3), mod.vref_line(d, 91.19),
                      mod.data_dot(d, 70, 500), mod.data_errorbar(d, 80, 600, 40, cap=0.5),
                      mod.data_trace(d, [60, 120], [100, 100]))

    a, b = linear_scene(old), linear_scene(sys.modules["style.bnd_style"])
    pa, pb = a.get_all_points(), b.get_all_points()
    ok("linear DataAxes + data_* identical to HEAD", pa.shape == pb.shape and np.array_equal(pa, pb),
       f"{pa.shape}")
    tmp.unlink()

# --- 3. stack_hist -----------------------------------------------------------
peak = schematic_zpeak(edges60, height=3e5)
layers = [("Fakes", np.zeros(60), SAMPLE["Fakes"]), ("WW", 60 * np.ones(60), SAMPLE["WW"]),
          ("WZ", np.zeros(60), SAMPLE["WZ"]), ("ZZ", 100 * np.ones(60), SAMPLE["ZZ"]),
          ("SingleTop", 50 * np.ones(60), SAMPLE["SingleTop"]), ("TTbar", 500 * np.ones(60), SAMPLE["TTbar"]),
          ("DYtautau", 0.02 * peak, SAMPLE["DYtautau"]), ("DYmumu", peak, SAMPLE["DYmumu"])]
st = stack_hist(dax, edges60, layers)
ok("stack_hist: 8 polygons", len(st) == 8 and st.names == [l[0] for l in layers])
ok("stack_hist: equal point counts", len({len(p.points) for p in st}) == 1)
ok("stack_hist: cum shape", st.cum.shape == (9, 60) and np.allclose(st.cum[-1], st.total))
ok("stack_hist: zero layers degenerate", st.layer["Fakes"].height < 1e-9 and st.layer["WZ"].height < 1e-9)
st0 = stack_hist(dax, edges60, [(n, np.zeros(60), c) for n, _, c in layers])
ok("stack_hist: all-zero stack 8 polygons at the floor", len(st0) == 8 and st0.height < 1e-9
   and abs(st0.get_top()[1] - dax.c2p(60, 1)[1]) < 1e-9)
ok("stack_hist: top of stack = c2p(total)", abs(st.layer["DYmumu"].get_top()[1] - dax.c2p(90, st.total.max())[1]) < 1e-9)

# --- 4. ratio_panel ----------------------------------------------------------
data = st.total * (0.94 + 0.02 * rng.standard_normal(60))
rp = ratio_panel(dax, edges60, data, st.total, center=(1.6, -1.55), num_err=np.sqrt(np.abs(data)),
                 show_x_labels=True)
ok("ratio_panel: 60 dots / 60 errs", len(rp.dots) == 60 and len(rp.errs) == 60)
ok("ratio_panel: x axis under the main plot",
   abs(rp.dax.c2p(60, 0.9)[0] - dax.c2p(60, 1)[0]) < 1e-9 and abs(rp.dax.c2p(120, 0.9)[0] - dax.c2p(120, 1)[0]) < 1e-9)
ok("ratio_panel: ref line at y=1", abs(rp.ref.get_center()[1] - rp.dax.c2p(90, 1.0)[1]) < 1e-9)
rp0 = ratio_panel(dax, edges60, data, np.where(np.arange(60) % 2, st.total, 0.0), center=(1.6, -1.55))
ok("ratio_panel: den==0 bins skipped", len(rp0.dots) == 30 and np.isnan(rp0.ratio[0]))
inside(VGroup(dax, st, rp), "log plot + stack + ratio")

# --- 5. mini_slice == place --------------------------------------------------
ms = mini_slice(0.42, (-4.6, -0.7))
pl = place(CMSSlice(), (0.42, (-4.6, -0.7)), DET_CENTER)
pm, pp = ms.get_all_points(), pl.get_all_points()
ok("mini_slice == place(CMSSlice) to 1e-9", pm.shape == pp.shape and np.max(np.abs(pm - pp)) < 1e-9,
   f"max |d| = {np.max(np.abs(pm - pp)):.2e}")
ok("mini_slice keeps .c / .radii", np.allclose(ms.c, [-4.6, -0.7, 0]) and abs(ms.radii["muon"][1] - 0.42 * R_DET) < 1e-12)

# --- 6. kappa_from_pt --------------------------------------------------------
ok("kappa_from_pt(45) == 0.28", kappa_from_pt(45) == 0.28)
ok("kappa_from_pt clips", kappa_from_pt(5) == 0.55 and kappa_from_pt(500) == 0.10 and 0.10 < kappa_from_pt(90) < 0.28)

# --- 7. events ---------------------------------------------------------------
det = CMSSlice()
ev = {"run": 278820, "lumi": 1, "event": 12345, "mass": 91.2,
      "mu": [{"idx": 0, "pt": 52.3, "eta": 0.4, "phi": 0.61, "charge": -1, "iso": 0.02},
             {"idx": 1, "pt": 38.1, "eta": -0.8, "phi": -2.55, "charge": 1, "iso": 0.05}]}
evg = event_from_json(det, ev)
ok("event_from_json: 2 muons, pieces", len(evg) == 2 and all(len(m) == 3 for m in evg)
   and evg[0].trk is evg[0][0] and evg[0].hits is evg[0][1] and evg[0].deposits is evg[0][2])
ok("event_from_json: .muons carried", evg.muons == ev["mu"] and evg.record is ev)
ev2 = {"muons": ev["mu"]}
ok("event_from_json: 'muons' key tolerated", len(event_from_json(det, ev2)) == 2)
ok("event_from_json: charge sign bends", np.sign(evg[0].trk.pts[-1][1] - det.point_at(det.outer_radius, 0.61)[1]) != 0)
pair = {"tag": {"pt": 41.0, "eta": 0.2, "phi": 1.2, "charge": 1, "iso": 0.01},
        "probe": {"pt": 22.0, "eta": -0.5, "phi": 3.9, "charge": 1, "iso": 0.45, "cls": "anti"},
        "passes": False, "jets": [{"pt": 35.0, "eta": -0.5, "phi": 3.95, "dr_probe": 0.05, "contains_probe": True}]}
pg = pair_from_json(det, pair, seed=1)
ok("pair_from_json: probe in jet", len(pg) == 2 and pg.in_jet and pg.probe is pg[1] and len(pg.probe) == 2
   and pg.probe.jet is pg.probe[0] and pg.probe.mu is pg.probe[1] and len(pg.tag) == 3)
pair2 = dict(pair, passes=True, jets=[])
pg2 = pair_from_json(det, pair2)
ok("pair_from_json: passing probe plain", not pg2.in_jet and len(pg2.probe) == 3)
mj = muon_in_jet(det, 0.3, -1, kappa_from_pt(30), seed=2)
ok("muon_in_jet: jet first, muon on top", mj[0] is mj.jet and mj[1] is mj.mu)
inside(VGroup(det, evg, pg), "detector + event + pair")

# --- 8. keys, pulls, grids, slider, spine, clock, rain ------------------------
key = colour_key([(r"\mathrm{Data}", SAMPLE["Data"], "dot"), (r"\mu^{+}\mu^{-}", SAMPLE["DYmumu"]),
                  (r"\mathrm{other}", SAMPLE["TTbar"]), (r"\mathrm{Fakes}", SAMPLE["Fakes"])])
ok("colour_key: 4 rows", len(key.rows) == 4 and len(key.labels) == 4)
try:
    colour_key([("a", INK)] * 5)
    ok("colour_key rejects 5 rows", False)
except ValueError:
    ok("colour_key rejects 5 rows", True)
inside(key, "colour_key")

pp_ = pull_plot(["SigModel", "MuonScale", "MuonRes", "PDF", "QCDScale", "Lumi"],
                [0.66, -0.48, 0.40, 0.56, 0.76, 0.08], [0.14, 0.14, 0.12, None, None, 0.98])
ok("pull_plot: 6 rows", len(pp_.rows) == 6 and abs(pp_.x_of(0.0) - pp_.zero.get_center()[0]) < 1e-9)
ok("pull_plot: dot at the pull", abs(pp_.rows[0][2].get_center()[0] - pp_.x_of(0.66)) < 1e-9)
pp_.shift(np.array([3.0, 0.5, 0.0]))
ok("pull_plot: x_of follows a shift", abs(pp_.x_of(0.0) - pp_.zero.get_center()[0]) < 1e-9)
inside(pp_, "pull_plot")
pps = pull_plot([r"\theta_1", r"\theta_2"], [0.3, -0.5], [0.6, 0.8], tex=True)
ok("pull_plot: tex names", len(pps.rows) == 2)

ff = np.array([[0.093, 0.103, 0.087, 0.118], [0.10, 0.12, 0.11, 0.15], [0.12, 0.14, 0.15, 0.20],
               [0.15, 0.18, 0.20, 0.26], [0.20, 0.25, 0.28, 0.33], [0.26, 0.30, 0.35, 0.39]])
vg = value_grid(ff, [r"20\text{--}26", "26\\text{--}30", "30\\text{--}35", "35\\text{--}40", "40\\text{--}50", ">50"],
                [r"0\text{--}0.9", r"0.9\text{--}1.2", r"1.2\text{--}2.1", r"2.1\text{--}2.4"])
ok("value_grid: 24 cells, vmax from data", len(vg.cells) == 24 and len(vg.texts) == 24 and vg.vmax == 0.39)
inside(vg, "value_grid 6x4")

sl = slider(0.95, 1.03, 0.988, err=0.016, ref=1.0)
ok("slider: marker at value", abs(sl.marker.dot.get_center()[0] - sl.x_of(0.988)) < 1e-9 and sl.marker is sl[-1])
ok("slider: ref at 1.0", abs(sl.ref.get_center()[0] - sl.x_of(1.0)) < 1e-9)
sl.shift(np.array([2.6, 1.9, 0.0])).scale(0.8)
m2 = sl.marker_at(1.0, 0.0)
ok("slider: x_of / marker_at follow move+scale", abs(m2.dot.get_center()[0] - sl.ref.get_center()[0]) < 1e-9
   and abs(m2.dot.get_center()[1] - sl.axis.get_center()[1]) < 1e-9 and m2.bar.get_stroke_width() == 0)
inside(sl, "slider")

NODE_X = (-6.2, -4.13, -2.07, 0, 2.07, 4.13, 6.2)
icons = [Square(side_length=1.6, stroke_color=col(INK)) for _ in NODE_X]
sp = spine(icons, NODE_X, done=(0, 1))
ok("spine: 7 boxes, 6 arrows, icons shrunk", len(sp.boxes) == 7 and len(sp.arrows) == 6
   and all(abs(i.width - 0.8 * 1.24) < 1e-9 for i in sp.icons) and len(sp.nodes) == 7)
ok("spine: node centres", np.allclose(sp.centers[:, 0], NODE_X) and np.allclose(sp.nodes[3].get_center(), [0, 0, 0]))
inside(sp.copy().shift(DOWN * 0.55), "spine(7) at SPINE_Y = -0.55 (s2_pipeline)")

clk = clock((-1.4, 2.3))
ok("clock: hand follows turns", abs(clk.hand.get_end()[1] - (2.3 + 0.72 * 0.32)) < 1e-9)
clk.turns.set_value(0.25)
clk.hand.update()
ok("clock: quarter turn -> 3 o'clock", abs(clk.hand.get_end()[0] - (-1.4 + 0.72 * 0.32)) < 1e-9)
inside(clk, "clock")

det_p = place(CMSSlice(), (0.42, (-4.6, -0.7)), DET_CENTER)
rn = rain(det_p, dax, edges60, peak, n=60, seed=3)
starts = np.array([a.animations[0].mobject.get_center() for a in rn])
ok("rain: 60 successions from inside the placed tracker", len(rn) == 60
   and np.all(np.linalg.norm(starts - np.array([-4.6, -0.7, 0]), axis=1) <= 0.42 * R_DET * LOGO_R["tob"] / LOGO_R["muon"] + 1e-9))

# --- 9. add_state / check_order ---------------------------------------------
class FakeScene:
    def __init__(self):
        self.mobjects = []

    def add(self, *ms):
        self.mobjects += list(ms)


fs = FakeScene()
state = {"det": det, "dax": dax, "stack": st}
add_state(fs, state, ("det", "dax", "stack"))
fs.add(Mobject(), Group(), ValueTracker(2.0))     # wait placeholder, finished Succession, tracker
check_order(fs, state, ("det", "dax", "stack"))
ok("check_order passes with placeholders / trackers", True)
try:
    check_order(fs, state, ("dax", "det", "stack"))
    ok("check_order catches a swapped order", False)
except AssertionError as e:
    ok("check_order catches a swapped order", "['det', 'dax', 'stack']" in str(e))

print("ALL OK")
