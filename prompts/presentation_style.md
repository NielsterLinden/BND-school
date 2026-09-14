---
name: pdf-presentation-style
description: Use whenever building a PDF presentation / slide deck / review document for the user — for ANY study, in any project. Defines the house style (16:9 67.74x38.1 cm slides, near-black #222222 background matching the plots, light-blue uppercase titles with grey rule + page number, green secondary accent, drawn tables with green best-value highlight), the content conventions (cover with the headline result, numbered section dividers, one idea per slide, results first and derivations last), inline LaTeX in titles via matplotlib mathtext, and vector-PDF plot embedding so everything zooms. Ships a reusable PyMuPDF library. Triggers on "make a presentation", "build a slide deck", "review PDF", "summary slides", "presentation of the results", "put the plots into slides", "make it look like my other talks", "beamer/pptx alternative", "build_comparison_deck".
---

# PDF presentation house style

How to build a slide deck the user will recognise as theirs. This is the **style
and content spec** — it is deliberately study-agnostic; nothing here depends on a
particular analysis. For the specific pattern of *gathering many variants of the
same plot so a reviewer scrolls to compare*, also read the **`comparison-deck`**
skill (grids, axes, campaign wiring); this skill is what that one is built on.

Decks are built with **PyMuPDF (fitz)**, not LaTeX/beamer/pptx: the deliverable
is a single self-contained PDF that embeds every plot as **vector**, so the user
can zoom into any figure without pixelation.

## Reusable implementation

`/project/atlas/users/sjankovy/boostHHbbtautau/HHARD_workfolder/build_comparison_deck.py`

Self-contained — needs only `fitz` (pymupdf) and `matplotlib`. **Copy it next to
your project** if you are working outside that workfolder; there is nothing
HHARD-specific in it. It exposes a `Deck` class implementing everything below, so
in practice you should *use it* rather than re-implement the spec:

```python
from build_comparison_deck import Deck, Style, cm2pt
d = Deck("out.pdf", author="Samuel Jankovych", date="2026-08-04")
d.cover(title, subtitle, lines=[...])   # title slide
d.divider(r"1.  $\mu$ results")         # numbered section separator
d.image(title, "plot.pdf")              # one figure, full page
d.panels(title, pdfs, caps, ncols=3)    # N captioned figures
d.grid(title, rows, cols, cells)        # rows x cols comparison grid
d.tables(title, [{...}])                # drawn tables
d.save()
```

Worked example: `HHARD_workfolder/build_wp_summary_pdf.py` (65-page review deck).

Run it in an environment with both deps. If a CVMFS/LCG environment has leaked
into the shell, strip it or the venv python will import the wrong numpy:

```bash
env -u PYTHONPATH -u LD_LIBRARY_PATH -u PYTHONHOME <venv>/bin/python build_<name>_deck.py ...
```

## Slide geometry

- **Page: 67.74 x 38.1 cm = 1920 x 1080 pt, 16:9 landscape.** `cm2pt(cm) = cm *
  72/2.54`. (38.1 cm is exactly 15 in.)
- **Never hard-code sizes against the page.** Express every font size and offset
  in a **1280-pt reference frame** and multiply by `k = page_w / 1280`. That is
  how the same code renders identically at any page size — changing the page must
  not shrink the text relative to the slide. The library does this via `self.k`.
- Margin 32 (ref) = 48 pt at 1920. Content starts at y = 78 (ref) below the title
  rule.

## Palette (exact)

| role | hex | use |
|---|---|---|
| page background | `#222222` | **must equal the plots' background** |
| primary / titles | `#2BA4DD` | slide titles, cover title, dividers, table header |
| secondary accent | `#8AC63F` | panel captions, sideband/cross-check/LNT slides |
| rule | `#8C8C94` | thin line under the title, cover + divider rules |
| body text | `#E6E6E6` | labels, column headers, table cells |
| muted | `#99999E` | author name, page number, row labels |
| best-value fill | `#577D29` | highlighted (winning) table cell |
| table row A / B | `#252528` / `#2F2F33` | zebra striping |
| cell border | `#666673` | table grid lines |
| missing marker | `#CC4C4C` | placeholder box when a figure is absent |

**The single most important rule: the page background must be identical to the
figures' background.** Otherwise every embedded plot shows as a lighter rectangle
floating on the slide. If you change one, change the other. In this user's
matplotlib code that constant is `DARK_BG = "#222222"`.

## Typography and slide furniture

- **Title**: bold, UPPERCASE, top-left, primary blue (or green on
  cross-check/sideband slides), 24 (ref) pt, **auto-shrinking** so it can never
  collide with the page number.
- **Thin grey rule** under the title, full content width. **No filled title bar.**
- **Page number** top-right, muted, small.
- **Cover**: muted uppercase author name → rule → big blue uppercase title →
  subtitle+date → a short block of plain body lines.
- **Divider**: centred blue uppercase section title with a short rule under it,
  numbered (`1. …`, `2. …`). Nothing else on the slide.
- **Panel captions**: centred, green, above each figure.
- **Grid labels**: column labels across the top (body colour), row labels down
  the left (muted).

## Inline LaTeX in titles

PyMuPDF's built-in fonts have **no Greek or math glyphs**, and Unicode subscripts
don't cover most letters — so math is rendered through **matplotlib mathtext**
into a tiny vector PDF and embedded. Write LaTeX directly in any title:

```python
d.divider(r"1.  Fit results  ($\mu_{HH}$, $\kappa_\lambda$, $\kappa_{2V}$)")
d.grid(r"SR VBF   —   $m(HH,\ E_{T}^{miss})$", ...)
d.tables("Ranking", [{"title": r"$\kappa_{2V}$ 68% CL interval width", ...}])
```

Use it for **every** symbol that has a proper form: POIs (`$\mu_{HH}$`,
`$\kappa_\lambda$`, `$\kappa_{2V}$`), observables (`$m(H_{bb})$`, `$p_{T}$`,
`$E_{T}^{miss}$`, `$s_{HH}$`), and operators (`$\chi^{2}$`, `$\geq$`,
`$\approx$`, `$\pm$`). Never write `MUHH`, `KAPPA2V` or `M(HBB)`.

Two implementation details that matter:

- **Upper-casing must skip math.** `smart_upper()` upper-cases only the text
  outside `$...$`, so `SR VBF — $m(H_{bb})$` doesn't become `M(HBB)`.
- **Baseline alignment must be exact**, or titles jitter vertically from slide to
  slide. Measure the text's ink box with a probe figure
  (`text.get_window_extent(renderer)` at dpi=72 so px == pt), build the final
  figure at exactly that size plus a small pad, draw at a known baseline, and
  embed with the baseline offset applied. Cache by (text, size, colour) — a deck
  re-renders the same strings many times.

A side benefit: titles then use the same font as the figures.

## Figures: vector or it doesn't count

- Embed the **`.pdf`** of each figure with `page.show_pdf_page()`, preserving
  aspect ratio and centring it in its box. Fall back to a `.png` of the same
  stem, and draw a red placeholder box if neither exists (never fail the build on
  a missing plot — you want to see the hole).
- Therefore **every plotting script must save both formats on the deck
  background**:

```python
fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight", facecolor=DARK_BG)
fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor=DARK_BG)
```

  If you inherit a codebase that only writes PNGs, patch its savefig sites first —
  that one change is what makes the whole deck zoomable.

## Tables

**Never render a table as markdown or plain text** — it comes out as unaligned
raw characters and the user will reject it. Draw real cells: blue header row,
zebra body rows, thin borders, and (where there is a "best" per column) fill the
winning cell green and bold it. `Deck.tables()` does this; give it
`{"title", "headers", "rows", "best_cols"}`.

## Content conventions

- **One idea per slide.** A slide is a title plus one figure, one grid, or one
  set of tables.
- **The cover carries the headline result**, not just the topic: what was
  compared, the method definitions, the setup in one line each, and the
  conclusion (e.g. "X wins in all 4 methods"). A reviewer who reads only the
  cover should know the answer.
- **Numbered section dividers** for every section, and list the sections on the
  cover so the deck is navigable by scrolling.
- **Order: results first, inputs and derivations last.** What worked:
  1. results / fit numbers (ranking tables, then per-observable comparisons)
  2. cross-checks and closure tests
  3. selection / cut derivations
  4. yields
  5. input distributions
  6. auxiliary derived quantities (transfer factors, calibration constants, …)
  A reader stops when they've seen enough; don't make them scroll past 40 slides
  of inputs to reach the answer.
- **A method/definition diagram early** (right after the cover) pays for itself.
  Keep diagrams *simple* — a labelled shape (e.g. one square split into four
  quadrants for an ABCD method), not a busy figure.
- **Comparison slides**: rows and columns are the things being compared, one
  slide per (region, variable). See `comparison-deck`.
- **Use a secondary-colour title** for a whole family of slides that is
  qualitatively different (sidebands, loose regions, cross-checks) so the reader
  can tell at a glance which family they're in.
- Where a figure family exists in an "overlaid" and a "separated" form, include
  both — overlaid for the comparison, separated when scales differ by orders of
  magnitude and would hide small entries.

## Verify before delivering

Render a few pages to PNG and actually look at them — far cheaper than opening
the whole PDF, and it catches the things that silently look wrong:

```python
doc[i].get_pixmap(dpi=42).save("/tmp/check.png")     # 1920pt page -> ~1120 px
```

Check: cover, one grid slide, one table slide, and any slide with new math.
Look for overlapping labels, clipped titles, figures that don't fill their box,
and — most importantly — a visible rectangle edge around a plot (background
mismatch). Also print the page box to confirm the size:
`doc[0].rect.width / 28.3465` → cm.

## Gotchas

- **✓ / ✗ and other symbol glyphs render as `¤`** in the built-in PDF fonts. Use
  the words "pass"/"fail", or route the text through mathtext.
- A vector-heavy deck is **tens of MB** (a 65-page one is ~29 MB). That's
  expected; `doc.save(..., deflate=True)` is already on. Don't rasterise to
  shrink it — zoomability is the point.
- A 1xN grid of tall figures renders them tiny; use `panels(ncols=3)` for ~5
  figures instead so each stays legible.
- Keep titles short — auto-shrink protects the layout but a 12-pt title on a
  1080-pt page looks broken.
- Build is a few seconds plus ~0.2 s per unique math string; caching matters.
- Save the deck next to its inputs (e.g. `<campaign>/summary/<Name>_summary.pdf`)
  and rebuild it in place, so the link the user has keeps working.