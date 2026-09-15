# BND School 2026 – CMS Open Data: Z boson decays

Shared workspace for our group project at the BND Graduate School (Heimbach, 2026).
We analyse CMS 2016 open data (NanoAOD, 13 TeV) and split into three subgroups by Z decay channel.

## Layout

| Folder | Channel | Hand-in |
|---|---|---|
| `z-ee/` | Z → e⁺e⁻ | `z-ee/handoff.md` |
| `z-mumu/` | Z → μ⁺μ⁻ | `z-mumu/handoff.md` |
| `z-tautau/` | Z → τ⁺τ⁻ | `z-tautau/handoff.md` |
| `combination/` | Combination of the channels — **[`combination/result.md`](combination/result.md)** | `combination/handoff.md` |
| `presentation/` | Final slide deck covering strategy and results | – |
| `reference/` | Shared material (intro slides, links) | – |
| `TRExFitter-v1.8.0/` | TRExFitter source (profile likelihood fits) | – |

## TRExFitter version

**Every agent and every subgroup must use TRExFitter v1.8.0** (the copy in `TRExFitter-v1.8.0/`).
The combination relies on all three channels producing workspaces with the same fitter version, so do not use another release or `master`.

Each subgroup owns its folder: put notebooks, scripts and plots there.
Keep `handoff.md` up to date so the other groups can pick up where you left off.

## Rules of the road

- Work directly on `main`. Pull before you push (`git pull --rebase`).
- Do **not** commit data files (ROOT files, large CSVs). Copy them locally and add the path/DOI to your `handoff.md` instead.
- Strip notebook outputs before committing, but keep the numbers the combination needs (see below) in your `handoff.md`.

## Result so far

The μμ and ττ channels are finished and combined (z-ee is not yet available):

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1940 ± 35 pb** per lepton flavour, assuming lepton
> universality — [`combination/result.md`](combination/result.md), deck in
> [`presentation/`](presentation/).

## What the combination needs

- `n_obs` — data events in the 60–120 GeV window
- `n_bkg` — background in the same window
- `acc_eff` — A·ε from MC
- systematic uncertainties, same names in all three channels

## Useful links

- CERN Open Data Portal: <https://opendata.cern.ch>
- CMS Open Data guide: <https://cms-opendata-guide.web.cern.ch>
- CMS Open Data forum: <https://opendata-forum.cern.ch>
- 2016 NanoAOD workshop lessons: <https://cms-opendata-workshop.github.io>
