# 05 — Colour schema: six chapter colours + the official CMS colours

Decided 2026-09-15 (user); replaces the flag palette of 2026-09-14. Every colour
in a clip is one of the six chapter colours, a tint/shade of one
(`palette.tint(hex, a)` / `palette.shade(hex, a)`), or — for the detector
drawing only — an official CMS logo colour.

## The six chapters (`CHAPTER` in `style/palette.py`)

| # | chapter | hex | name in the palette |
|---|---------|-----|---------------------|
| 1 | theory | `#702677` purple | `PURPLE`, `THEORY` |
| 2 | detector / methods | `#00B3C3` cyan | `CYAN`, `DETECTOR_ACCENT` |
| 3 | Z→ee | `#059F77` green | `GREEN`, `CHANNEL["ee"]` |
| 4 | Z→μμ | `#E3D88B` gold | `GOLD`, `CHANNEL["mumu"]` |
| 5 | Z→ττ | `#E6223C` red | `RED`, `CHANNEL["tautau"]` |
| 6 | combination | `#2B2D42` dark slate | `SLATE`, `COMBINED`, `INK` |

These six colours are the PowerPoint chapter colours too (the six-bar anchor).

## Semantic assignment

| role | colour | note |
|------|--------|------|
| background | white | |
| ink (axes, quark lines, symbols) | slate | not pure black, so everything sits inside the scheme |
| **Z→ee** lines, glyphs, fills | green | `CHANNEL_LINE["ee"]` = `CHANNEL["ee"]` |
| **Z→μμ** fills | gold | `CHANNEL["mumu"]` |
| **Z→μμ** lines and glyphs | shade(gold, 0.30) `#9F9761` | gold is too light for thin strokes on white |
| **Z→ττ** lines, glyphs, fills | red | |
| channel tracks in the detector | channel colour core on a `TRACK_OUTLINE[ch]` rim (shade 0.55) | the logo's own double-line muon style: `track(..., outline=...)` |
| combined result | slate | |
| theory band / line | purple | |
| highlight / "look here" | red (`HIGHLIGHT`) | in ττ clips Flash in slate instead, red is the channel there |
| photon | shade(gold, 0.15) | |
| pion / hadron / jet / neutrino | tints of slate (greys) | neutrinos dashed, `p_T^miss` a slate dashed arrow |
| fit backgrounds | TTbar, SingleTop, Fakes: slate tints; WW, ZZ: cyan tints; WZ: purple tint | never compete with the signal |

## The CMS detector: official logo colours (`CMS`, `DETECTOR`)

From Izaak Neutelings' `CMS_logo.tex` (cms-docdb 3045). Radii are the logo's
(`LOGO_R`, fraction of the outer muon radius): TIB 0.059, TOB 0.1145, ECAL
0.220, band 0.231, HCAL 0.385, magnet 0.578, muon 0.996.

| subsystem | fill | stroke |
|-----------|------|--------|
| tracker (TIB disc, pixel + strip rings) | `#ADEB8F` light green | shade(TOB, .35) |
| TOB ring | `#92E569` green | shade(., .35) |
| ECAL (60 cells) | `#92E569` green | shade(., .35) |
| band (thin ring outside the ECAL) | `#BEAED4` lavender | shade(., .30) |
| HCAL (36 towers) | `#E5E9EE` light grey | shade(., .35) |
| solenoid | `#FFDF7F` yellow | shade(., .35) |
| muon system (blue annulus; 4 dodecagonal stations on top) | `#85D1FB` blue (stations shade .10) | shade(., .40) |
| return yoke (3 dodecagonal rings between the stations) | `#E5E9EE` light grey | shade(., .30) |
| logo muons | `#F0240B` red with a white rim | only in the logo clip |
| logo background / wordmark | `#FDFAF4` / `#101177` | only in the logo clip |

Deposits (`DEPOSIT`): a *channel lepton* lights its cells in its own chapter
colour with a `TRACK_OUTLINE` rim (electron ECAL cluster green, muon stubs and
MIP dots gold, τ_h cells red); anything else (photon, jet) in the subsystem
shade (`DEPOSIT["ecal"|"hcal"|"muon"]`).

Geometry: `CMSSlice()` is centred at `DET_CENTER = (0, -0.25)` with outer
radius `R_DET = 2.95` (top edge y = 2.70, under the title band). The logo clip
(`scenes/s2_cms_logo.py`) turns a `CMSLogo` into exactly this slice.

## Where it lives

- `style/palette.py` — manim-free source of truth (also for matplotlib).
- `style/bnd_style.py` — re-exports, `col()`, `lighten()`/`darken()`,
  `LOGO_R`, `logo_rings`, `logo_muons`, `CMSLogo`, `CMSSlice`, `track`
  (`outline=`), `neutral_track`, `calo_hit`, `muon_hits`, `tracker_hits`,
  `signature` (`e`, `mu`, `gamma`, `tau_h` with `prongs`/`pi0`, `jet`, `nu`, `met`),
  `shower_tree`, `hadron_blob`, `pion_lines`.
- `scenes/channel_common.py` — the shared choreography of the section 3–5 clips.
