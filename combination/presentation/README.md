# Presentation — the Z cross section from CMS Open Data

The final deliverable of the BND-school project: **[`bnd_z_combination.pdf`](bnd_z_combination.pdf)**
(29 slides, 16:9). It covers the whole analysis — data, the three channel measurements, the
combination method and the results — with the results first and the derivations last.

```bash
bash build.sh              # figures (if missing) + pdflatex twice + sanity checks
bash build.sh --figures    # regenerate figures/*.pdf first
```

## What it needs

| | |
|---|---|
| LaTeX | the **cvmfs TeX Live 2025** (`/cvmfs/sft.cern.ch/lcg/external/texlive/2025`), which `build.sh` puts on `PATH`. The system TeX Live 2020 in `/bin` has beamer but no `metropolis`, `pgfopts` or `siunitx`. |
| combination outputs | `../output/plots/*.pdf` — run `cd .. && python run_combination.py` first; `build.sh` checks and says so |
| deck figures | `figures/*.pdf` from `python figures.py` (needs `source ../../setup.sh`) |
| channel plots | `../../z-mumu/output/v2/plots/*.png`, `../../z-tautau/output/plots/*.png`, committed by those groups; the z-ee plots come out of `../../z-ee/Zee_fit.tar.gz` via `figures.py` |

`build.sh` reports undefined references, missing figures and overfull boxes, and prints the page
count — a missing figure shows up as a hole on a slide, so it is worth reading the output.

## Structure

| slides | |
|---|---|
| 1 | cover, carrying the headline result, the χ² and the section list |
| 2 | the measurement at a glance |
| 3–8 | **the result**: combined cross section, **why z-ee and z-mumu disagree**, why combining now buys something, uncertainty breakdown, lepton universality, per-channel and combined comparison with CMS and ATLAS |
| 9–10 | strategy diagram and the datasets |
| 11–16 | `Z→μμ`, `Z→τhτh` and `Z→ee`, two slides each except ee |
| 17–19 | the combination: method, correlation model, cross-checks, what it is not |
| 20 | conclusions |
| 21–29 | backup: the reference spread, full breakdown, variations, channel numbers, **orthogonality**, **systematics vs CMS/ATLAS**, where everything lives |

## Style

beamer + `metropolis`, light background, 9 pt, 16:9. The channels' plots are white-background
CMS-style figures, so the deck is light too and they drop in without a visible rectangle — the
rule from `prompts/presentation_style.md` that matters most, applied to a light palette. Colours
match the figures: <span>`#2B6CB0`</span> for μμ, `#EB811B` for ττ, `#23373B` for the combination,
`#8E44AD` for ee, `#14B03D` for "correlated / reference", `#C0392B` for warnings.

(The dark house style of `prompts/presentation_style.md` is rendered in beamer by
`../../z-mumu/review/deck/zmumu_review.tex` if you ever want that variant; its figures are
dark-background and must not be mixed with these.)

## Files

| file | what |
|---|---|
| `bnd_z_combination.tex` | the deck |
| `build.sh` | build + checks |
| `figures.py` | the deck-specific figures (`mumu_fit_stability`, `precision_budget`, `zee_normalisation`) and the z-ee TRExFitter plots pulled out of `../../z-ee/Zee_fit.tar.gz` (`zee_SR_postfit`, `zee_pulls`, …) |
| `figures/` | their PDF (vector, used) and PNG (fallback) output |
| `bnd_z_combination.pdf` | the built deck, committed |
