#!/usr/bin/env python3
"""
build_comparison_deck.py — a small, reusable library for assembling ONE
scrollable PDF that gathers the plots of one or several analyses/options/methods
side-by-side, so a reviewer just scrolls to compare.

The design goal: every comparison lives on its own slide as a GRID whose rows
and columns are the things you are comparing (e.g. QCD methods x WP options,
or analysis-A vs analysis-B x variables). Because the same (region, variable)
plot for every variant sits on one page, scrolling through the deck walks you
through the whole comparison with nothing to hunt for.

Style: matches Samuel Jankovych's presentation look — near-black background
(#222222, the SAME as the plots' DARK_BG so embedded plots blend seamlessly),
light-blue bold-uppercase titles with a thin grey underline rule + page number
top-right (no filled title bar), green secondary accent, blue-header tables with
a green best-value highlight. All plots are embedded as VECTOR PDFs (via
PyMuPDF `show_pdf_page`) so you can zoom in without pixelation; a PNG next to
the PDF is used as a fallback.

MATH IN TITLES: titles/dividers/table titles are rendered through matplotlib's
mathtext and embedded as vector, so you can write LaTeX inline —
`r"Expected $\\mu_{HH}$ 95% CL upper limit"`, `r"$m(H_{bb})$"`,
`r"$\\kappa_{2V}$"`. PyMuPDF's built-in fonts have no Greek/math glyphs, hence
the detour; it also makes titles match the plots' typography. Text inside
`$...$` is never upper-cased.

PAGE SIZE: default 67.74 x 38.1 cm (= 1920 x 1080 pt, 16:9). Every font size and
offset scales with the page width, so changing the size keeps the layout.
Use `Style(page_w=cm2pt(w), page_h=cm2pt(h))` for a different one.

Requirements: `pymupdf` (fitz) + `matplotlib`. Run in an environment that has
both, e.g. the betterplottingtool venv:
    env -u PYTHONPATH -u LD_LIBRARY_PATH -u PYTHONHOME <venv>/bin/python ...

Two ways to use it:

1. As a library (recommended — full control, see build_wp_summary_pdf.py):
       from build_comparison_deck import Deck
       d = Deck("out.pdf", author="Samuel Jankovych", date="2026-08-04")
       d.cover("My study", "subtitle", ["line 1", "line 2"])
       d.divider(r"1.  $\\mu_{HH}$ results")
       d.tables("Ranking", [{"title": r"Expected $\\mu_{HH}$", "headers": [...],
                             "rows": [...], "best_cols": {1, 2, 3}}])
       d.image("A diagram", "diagram.pdf")
       d.plot_grid(r"SR ggF — $m(H_{bb})$", rows=[("m1", "m1 B/D")],
                   cols=[("TT", "TT")], cell=lambda rk, ck: f"plots/{ck}/{rk}.pdf")
       d.save()

2. As a CLI over a spec file (a python module exposing `build(deck)`):
       <venv>/bin/python build_comparison_deck.py --spec my_spec.py --output out.pdf

IMPORTANT — the plots must exist as .pdf (vector). Make your plotting save both
formats, e.g. matplotlib:
    fig.savefig(path.with_suffix(".png"), dpi=200, facecolor=fig.get_facecolor())
    fig.savefig(path.with_suffix(".pdf"), facecolor=fig.get_facecolor())
and use the SAME face/background colour as the deck (#222222) so they blend.
"""
from __future__ import annotations

import argparse
import importlib.util
import io
import json
import os
import re
from dataclasses import dataclass

import fitz  # PyMuPDF
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CM2PT = 72.0 / 2.54


def cm2pt(cm: float) -> float:
    return cm * CM2PT


# ---- math-capable vector text (matplotlib mathtext -> embedded PDF) --------
_RICH_CACHE: dict = {}


def _rich(text: str, fs: float, color, bold: bool):
    """Render `text` (may contain $math$) to a tight vector PDF page.

    Returns (doc, width, height, baseline_from_top) — all in points. The page is
    sized to the text's exact ink box plus a small pad and the baseline offset
    is exact, so titles never jitter vertically between slides.
    """
    key = (text, round(float(fs), 2), tuple(round(float(c), 4) for c in color), bold)
    if key in _RICH_CACHE:
        return _RICH_CACHE[key]
    weight = "bold" if bold else "normal"
    probe = plt.figure(figsize=(0.1, 0.1), dpi=72)
    t = probe.text(0, 0, text, fontsize=fs, fontweight=weight, va="baseline", ha="left")
    bb = t.get_window_extent(probe.canvas.get_renderer())
    w, asc, desc = float(bb.width), float(bb.y1), float(-bb.y0)
    plt.close(probe)
    pad = 2.0
    W, H = w + 2 * pad, asc + desc + 2 * pad
    fig = plt.figure(figsize=(W / 72.0, H / 72.0), dpi=72)
    fig.text(pad / W, (desc + pad) / H, text, fontsize=fs, fontweight=weight,
             color=color, va="baseline", ha="left")
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", transparent=True)
    plt.close(fig)
    doc = fitz.open(stream=buf.getvalue(), filetype="pdf")
    res = (doc, W, H, H - (desc + pad))
    _RICH_CACHE[key] = res
    return res


def smart_upper(s: str) -> str:
    """Upper-case everything OUTSIDE $...$ (math must keep its case)."""
    return "".join(p if p.startswith("$") else p.upper()
                   for p in re.split(r"(\$[^$]*\$)", s))


@dataclass
class Style:
    # default page: 67.74 x 38.1 cm = 1920 x 1080 pt (16:9)
    page_w: float = cm2pt(67.74)
    page_h: float = cm2pt(38.1)
    margin: float = 32.0          # in 1280-pt reference units; scaled by `k`
    author: str = "Samuel Jankovych"
    # colours (0-1 rgb). page_bg defaults to the plots' DARK_BG so plots blend.
    page_bg: tuple = (0.1333, 0.1333, 0.1333)   # #222222
    title: tuple = (0.169, 0.643, 0.867)        # light blue #2ba4dd
    accent2: tuple = (0.541, 0.776, 0.247)      # green #8bc63f
    rule: tuple = (0.55, 0.55, 0.58)
    text: tuple = (0.90, 0.90, 0.90)
    muted: tuple = (0.60, 0.60, 0.62)
    header_bg: tuple = None       # table header fill; defaults to `title`
    best_bg: tuple = (0.34, 0.49, 0.16)   # green highlight for best table cell
    row_a: tuple = (0.145, 0.145, 0.155)
    row_b: tuple = (0.185, 0.185, 0.20)
    cell_line: tuple = (0.40, 0.40, 0.45)
    missing: tuple = (0.80, 0.30, 0.30)

    def __post_init__(self):
        if self.header_bg is None:
            self.header_bg = self.title


class Deck:
    """Accumulate slides, then .save() to a single PDF.

    All geometry constants below are expressed in a 1280x720 reference frame and
    multiplied by `self.k` (= page_w / 1280), so any page size keeps the layout.
    """

    REF_W = 1280.0

    def __init__(self, output: str, style: Style | None = None,
                 author: str | None = None, date: str = ""):
        self.output = output
        self.s = style or Style()
        if author:
            self.s.author = author
        self.date = date
        self.k = self.s.page_w / self.REF_W
        self.M = self.s.margin * self.k
        self.doc = fitz.open()

    # ---- low-level text/embed helpers ------------------------------------
    def _bg(self, p):
        p.draw_rect(fitz.Rect(0, 0, self.s.page_w, self.s.page_h), color=None,
                    fill=self.s.page_bg)

    def _text(self, p, x, y, s, size, color, bold=True):
        """Plain (non-math) text via a built-in PDF font."""
        p.insert_text((x, y), s, fontsize=size, color=color,
                      fontname=("hebo" if bold else "helv"))

    def _ctext(self, p, cx, y, s, size, color, bold=True):
        font = "hebo" if bold else "helv"
        w = fitz.get_text_length(s, fontname=font, fontsize=size)
        p.insert_text((cx - w / 2, y), s, fontsize=size, color=color, fontname=font)

    def _rtext(self, p, x, y_base, s, size, color, bold=True, align="left", max_w=None):
        """Math-capable vector text; auto-shrinks to `max_w`. Returns (width, fs)."""
        fs = size
        doc, W, H, b = _rich(s, fs, color, bold)
        while max_w and W > max_w and fs > 8:
            fs -= 1
            doc, W, H, b = _rich(s, fs, color, bold)
        x0 = x if align == "left" else (x - W / 2 if align == "center" else x - W)
        y0 = y_base - b
        p.show_pdf_page(fitz.Rect(x0, y0, x0 + W, y0 + H), doc, 0)
        return W, fs

    def _auto(self, p, x, y, s, size, color, bold=True, align="left"):
        """Dispatch: math strings go through mathtext, plain ones stay native."""
        if "$" in s:
            self._rtext(p, x, y, s, size, color, bold=bold, align=align)
        elif align == "center":
            self._ctext(p, x, y, s, size, color, bold)
        else:
            self._text(p, x, y, s, size, color, bold)

    def _header(self, p, title, accent=None):
        """Blue/green bold-uppercase title + grey rule + page number. Returns
        content-start y. Title auto-shrinks to never hit the page number."""
        accent = accent or self.s.title
        k, M = self.k, self.M
        W = self.s.page_w
        self._rtext(p, M, 46 * k, smart_upper(title), 24 * k, accent,
                    max_w=W - 2 * M - 46 * k)
        self._ctext(p, W - M - 12 * k, 44 * k, str(p.number + 1), 15 * k, self.s.muted)
        p.draw_line(fitz.Point(M, 62 * k), fitz.Point(W - M, 62 * k),
                    color=self.s.rule, width=1.4 * k)
        return 78 * k

    def _place(self, p, box, pdf_path):
        """Embed the vector PDF (preferred) or PNG fallback, aspect-preserved."""
        if pdf_path is None:
            pdf_path = ""
        png = pdf_path[:-4] + ".png" if pdf_path.lower().endswith(".pdf") else pdf_path
        try:
            if pdf_path.lower().endswith(".pdf") and os.path.exists(pdf_path):
                src = fitz.open(pdf_path)
                r0 = src[0].rect
                r = min(box.width / r0.width, box.height / r0.height)
                w, h = r0.width * r, r0.height * r
                fit = fitz.Rect(box.x0 + (box.width - w) / 2, box.y0 + (box.height - h) / 2,
                                box.x0 + (box.width - w) / 2 + w, box.y0 + (box.height - h) / 2 + h)
                p.show_pdf_page(fit, src, 0)
                src.close()
                return
            if png and os.path.exists(png):
                p.insert_image(box, filename=png, keep_proportion=True)
                return
        except Exception as e:
            print(f"  [warn] could not embed {pdf_path}: {e}")
        p.draw_rect(box, color=self.s.missing, width=1)
        self._ctext(p, (box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2, "—",
                    14 * self.k, self.s.missing)

    def _new(self):
        p = self.doc.new_page(width=self.s.page_w, height=self.s.page_h)
        self._bg(p)
        return p

    # ---- public slide builders -------------------------------------------
    def cover(self, title, subtitle="", lines=None):
        p = self._new()
        k, M, W = self.k, self.M, self.s.page_w
        self._text(p, M, 150 * k, self.s.author.upper(), 26 * k, self.s.muted)
        p.draw_line(fitz.Point(M, 178 * k), fitz.Point(W - M, 178 * k),
                    color=self.s.rule, width=1.6 * k)
        self._rtext(p, M, 250 * k, smart_upper(title), 46 * k, self.s.title,
                    max_w=W - 2 * M)
        if subtitle:
            self._auto(p, M, 300 * k, subtitle + (("  " + self.date) if self.date else ""),
                       17 * k, self.s.text)
        y = 380 * k
        for ln in (lines or []):
            self._auto(p, M, y, ln, 14 * k, self.s.text, bold=False)
            y += 32 * k
        return self

    def divider(self, text):
        p = self._new()
        k, W, H = self.k, self.s.page_w, self.s.page_h
        self._rtext(p, W / 2, H / 2 - 6 * k, smart_upper(text), 34 * k, self.s.title,
                    align="center", max_w=W - 2 * self.M)
        p.draw_line(fitz.Point(W / 2 - 320 * k, H / 2 + 28 * k),
                    fitz.Point(W / 2 + 320 * k, H / 2 + 28 * k),
                    color=self.s.rule, width=1.6 * k)
        return self

    def image(self, title, pdf, accent=None):
        """One plot filling the page (e.g. a diagram, a ranking dumbbell)."""
        p = self._new()
        top = self._header(p, title, accent)
        self._place(p, fitz.Rect(self.M, top + 6 * self.k, self.s.page_w - self.M,
                                 self.s.page_h - self.M), pdf)
        return self

    def panels(self, title, pdfs, caps, ncols=2, accent=None):
        """N plots in a grid with per-panel captions (e.g. 2 or 4 methods)."""
        p = self._new()
        k, M = self.k, self.M
        top = self._header(p, title, accent)
        n = len(pdfs)
        nrows = (n + ncols - 1) // ncols
        cw = (self.s.page_w - (ncols + 1) * M) / ncols
        ch = (self.s.page_h - top - (nrows + 1) * M - nrows * 20 * k) / nrows
        for i, (pdf, cap) in enumerate(zip(pdfs, caps)):
            r, c = divmod(i, ncols)
            x0 = M + c * (cw + M)
            y0 = top + M + r * (ch + M + 20 * k)
            self._auto(p, x0 + cw / 2, y0 + 14 * k, cap, 14 * k, self.s.accent2,
                       align="center")
            self._place(p, fitz.Rect(x0, y0 + 20 * k, x0 + cw, y0 + 20 * k + ch), pdf)
        return self

    def grid(self, title, rows, cols, cells, accent=None, green_accent=False):
        """rows x cols grid of plots — the scroll-to-compare core.
        rows/cols: list of (key, label). cells: 2D list[str pdf-path] indexed
        [row][col]. Use `green_accent=True` (or accent=...) for a secondary
        colour (e.g. sideband / LNT pages)."""
        p = self._new()
        k, M = self.k, self.M
        acc = accent or (self.s.accent2 if green_accent else self.s.title)
        top = self._header(p, title, acc) + 22 * k
        nrows, ncols = len(rows), len(cols)
        left = 92 * k
        gw = (self.s.page_w - left - M) / ncols
        gh = (self.s.page_h - top - M) / nrows
        for c, (_, clab) in enumerate(cols):
            self._auto(p, left + c * gw + gw / 2, top - 6 * k, clab, 14 * k,
                       self.s.text, align="center")

        for r, (_, rlab) in enumerate(rows):
            self._auto(p, M - 6 * k, top + r * gh + gh / 2, rlab, 12 * k, self.s.muted)
        for r in range(nrows):
            for c in range(ncols):
                box = fitz.Rect(left + c * gw + 3 * k, top + r * gh + 3 * k,
                                left + (c + 1) * gw - 3 * k, top + (r + 1) * gh - 3 * k)
                self._place(p, box, cells[r][c])
        return self

    def plot_grid(self, title, rows, cols, cell, accent=None, green_accent=False):
        """Convenience: build a grid from a cell(row_key, col_key)->pdf_path
        function. rows/cols are lists of (key, label)."""
        cells = [[cell(rk, ck) for ck, _ in cols] for rk, _ in rows]
        return self.grid(title, rows, cols, cells, accent=accent, green_accent=green_accent)

    def tables(self, title, tables):
        """Several stacked tables on one page. Each table is a dict:
        {title, headers:[...], rows:[[label, v1, ...]], best_cols:set(int)}.
        Cells in best_cols are highlighted green at their per-column minimum."""
        p = self._new()
        k = self.k
        self._header(p, title)
        x = self.M + 40 * k
        y = 108 * k
        for t in tables:
            self._auto(p, x, y - 8 * k, t["title"], 15 * k, self.s.text)
            y = self._draw_table(p, x, y, t["headers"], t["rows"],
                                 t.get("col_w"), t.get("best_cols", set())) + 40 * k
        return self

    def _draw_table(self, p, x, y, headers, rows, col_w=None, best_cols=None,
                    rh=26, fs=12):
        k = self.k
        rh, fs = rh * k, fs * k
        best_cols = best_cols or set()
        ncol = len(headers)
        col_w = [w * k for w in (col_w or ([190] + [150] * (ncol - 1)))]
        best_val = {}
        for c in best_cols:
            vs = []
            for r in rows:
                try:
                    vs.append(float(r[c]))
                except (ValueError, IndexError):
                    vs.append(float("inf"))
            best_val[c] = min(vs) if vs else None
        cx = x
        for c, h in enumerate(headers):
            w = col_w[c]
            p.draw_rect(fitz.Rect(cx, y, cx + w, y + rh), color=self.s.cell_line,
                        fill=self.s.header_bg, width=0.5 * k)
            self._auto(p, cx + w / 2, y + rh - 8 * k, str(h), fs, (0.95, 0.95, 0.95),
                       align="center")
            cx += w
        for ri, row in enumerate(rows):
            ry = y + rh + ri * rh
            cx = x
            for c, val in enumerate(row):
                w = col_w[c]
                fill = self.s.row_a if ri % 2 == 0 else self.s.row_b
                if c in best_cols:
                    try:
                        if abs(float(val) - best_val[c]) < 1e-9:
                            fill = self.s.best_bg
                    except ValueError:
                        pass
                p.draw_rect(fitz.Rect(cx, ry, cx + w, ry + rh), color=self.s.cell_line,
                            fill=fill, width=0.5 * k)
                if c == 0:
                    self._auto(p, cx + 8 * k, ry + rh - 8 * k, str(val), fs, self.s.text)
                else:
                    self._auto(p, cx + w / 2, ry + rh - 8 * k, str(val), fs, self.s.text,
                               bold=(c in best_cols), align="center")
                cx += w
        return y + rh + len(rows) * rh

    def save(self):
        os.makedirs(os.path.dirname(os.path.abspath(self.output)), exist_ok=True)
        self.doc.save(self.output, deflate=True)
        print(f"Wrote {self.output} ({self.doc.page_count} pages, "
              f"{self.s.page_w / CM2PT:.2f} x {self.s.page_h / CM2PT:.2f} cm)")
        self.doc.close()


# --------------------------------------------------------------------------- #
#  Small helpers a study spec commonly needs
# --------------------------------------------------------------------------- #
def load_ranking_tables(ranking_json, options, methods, metrics):
    """Turn a WSMaker-style ranking.json into `Deck.tables` input.
    options: list[(key,label)]; methods: list[(json_key, header)];
    metrics: list[(title, kind)] where kind in {'mu','width:<field>'}.
    Titles may contain $math$."""
    r = json.load(open(ranking_json))

    def width(iv):
        try:
            lo, hi = iv[0]
            return hi - lo
        except Exception:
            return float("nan")

    tables = []
    for title, kind in metrics:
        rows = []
        for ok, olab in options:
            row = [olab]
            for mk, _ in methods:
                d = r[ok][mk]
                if kind == "mu":
                    row.append(f"{d['mu_exp']:.1f}")
                elif kind.startswith("width:"):
                    fld = kind.split(":", 1)[1]
                    fmt = "{:.2f}" if "2V" in fld else "{:.1f}"
                    row.append(fmt.format(width(d[fld])))
            rows.append(row)
        tables.append({"title": title, "headers": ["variant"] + [h for _, h in methods],
                       "rows": rows, "best_cols": set(range(1, len(methods) + 1))})
    return tables


def _run_spec(spec_path, output):
    spec = importlib.util.spec_from_file_location("deck_spec", spec_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    d = Deck(output)
    mod.build(d)          # the spec's build(deck) appends slides
    d.save()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", required=True, help="python file exposing build(deck)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    _run_spec(args.spec, args.output)
