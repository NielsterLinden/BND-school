# CLAUDE.md — presentation

For agents editing the deck. `README.md` says what it is and how to build it; this is about not
breaking it.

## Build

```bash
bash build.sh            # always build this way -- it sets the TeX PATH and runs the checks
bash build.sh --figures  # also regenerate figures/*.pdf
```

`pdflatex` must come from `/cvmfs/sft.cern.ch/lcg/external/texlive/2025/bin/x86_64-linux`. The
system TeX Live 2020 in `/bin` is missing `metropolis`, `pgfopts` and `siunitx` and will fail with
an unhelpful error. `build.sh` exports it; do not call `pdflatex` directly.

## Rules

1. **Numbers come from the combination, never from memory.** Every figure in the text should be
   traceable to `../combination/output/combination_result.json` or `../combination/result.md`. If
   the combination is rerun and a number moves, the deck must be updated — it is not generated.
   The slides that carry numbers: 1, 2, 3, 4, 13, 14, 18, 21, 22, 23.
2. **Figures are light background.** The channel plots are white CMS-style figures; a dark figure
   would be a black rectangle on the slide. Do **not** pull figures from
   `../z-mumu/review/figures/` — those belong to the dark review deck.
3. **Vector where it exists.** `../combination/output/plots/*.pdf` and `figures/*.pdf` are vector;
   the channel plots are only PNG (their groups do not commit PDFs), so those are used as PNG.
4. **`\graphicspath` resolves the short names** — write `\fig{forest.pdf}`, not a path. The search
   order is `figures/`, the combination plots, z-mumu plots, z-tautau plots, z-tautau slides.
5. **Check the build output.** `build.sh` prints undefined references, missing figures and the
   overfull-box count. A missing figure is a silent hole on a slide.

## Traps already hit

* `\tt` is a LaTeX built-in — the \tauh\tauh shortcut is `\thth`.
* A frame containing `verbatim` needs `[fragile]`.
* `$(a)!0.5!(b)$` in TikZ needs `\usetikzlibrary{calc}`.
* Overfull `\vbox` warnings are content running off the bottom of a slide. Under ~5 pt they are
  invisible; above ~10 pt they collide with the page number. Fix by trimming text or shrinking the
  figure, not by `\vspace{-...}`.

## Verifying a change

Render and *look* — cheaper than opening the PDF and it catches what warnings do not:

```bash
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2025/bin/x86_64-linux:$PATH
pdftoppm -r 85 -png bnd_z_combination.pdf /tmp/deck/p     # then read the PNGs
```

Always check the cover, the TikZ strategy diagram (slide 8) and any slide whose figure or numbers
you touched.
