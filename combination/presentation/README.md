# Presentation — the Z cross section from CMS Open Data

The final deliverable of the BND-school project: **[`bnd_z_combination.pdf`](bnd_z_combination.pdf)**
(24 slides, 16:9). It covers the whole analysis — data, the two channel measurements, the
combination method and the results — with the results first and the derivations last.

```bash
bash build.sh              # figures (if missing) + pdflatex twice + sanity checks
bash build.sh --figures    # regenerate figures/*.pdf first
```

## What it needs

| | |
|---|---|
| LaTeX | the **cvmfs TeX Live 2025** (`/cvmfs/sft.cern.ch/lcg/external/texlive/2025`), which `build.sh` puts on `PATH`. The system TeX Live 2020 in `/bin` has beamer but no `metropolis`, `pgfopts` or `siunitx`. |
| combination outputs | `../combination/output/plots/*.pdf` — run `cd ../combination && python run_combination.py` first; `build.sh` checks and says so |
| deck figures | `figures/*.pdf` from `python figures.py` (needs `source ../setup.sh`) |
| channel plots | `../z-mumu/output/v2/plots/*.png`, `../z-tautau/output/plots/*.png`, committed by those groups |

`build.sh` reports undefined references, missing figures and overfull boxes, and prints the page
count — a missing figure shows up as a hole on a slide, so it is worth reading the output.

## Structure

| slides | |
|---|---|
| 1 | cover, carrying the headline result and the section list |
| 2 | the measurement at a glance |
| 3–7 | **the result**: combined cross section, why combining gains nothing, uncertainty breakdown, lepton universality, comparison with CMS and ATLAS |
| 8–9 | strategy diagram and the datasets |
| 10–13 | `Z→μμ` (selection + corrections; which extraction to take) and `Z→τhτh` (selection + fakes; result) |
| 14–17 | the combination: method, correlation model, cross-checks, what it is not |
| 18 | conclusions |
| 19–24 | backup: the 0.47 % reference mismatch, full breakdown, variations, channel numbers, where everything lives |

## Style

beamer + `metropolis`, light background, 9 pt, 16:9. The channels' plots are white-background
CMS-style figures, so the deck is light too and they drop in without a visible rectangle — the
rule from `prompts/presentation_style.md` that matters most, applied to a light palette. Colours
match the figures: <span>`#2B6CB0`</span> for μμ, `#EB811B` for ττ, `#23373B` for the combination,
`#14B03D` for "correlated / reference".

(The dark house style of `prompts/presentation_style.md` is rendered in beamer by
`../z-mumu/review/deck/zmumu_review.tex` if you ever want that variant; its figures are
dark-background and must not be mixed with these.)

## Files

| file | what |
|---|---|
| `bnd_z_combination.tex` | the deck |
| `build.sh` | build + checks |
| `figures.py` | the two deck-specific figures (`mumu_fit_stability`, `precision_budget`) |
| `figures/` | their PDF (vector, used) and PNG (fallback) output |
| `bnd_z_combination.pdf` | the built deck, committed |
