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

All three channels are finished and combined:

> **σ(pp → Z/γ* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1980 ± 32 pb** per lepton flavour, assuming lepton
> universality, with **χ²/ndf = 9.67/2 (p = 0.008)** — [`combination/result.md`](combination/result.md),
> deck in [`combination/presentation/`](combination/presentation/).

The compatibility number is part of the result: σ(ee) = 2054 ± 40 pb sits 3.0σ above
σ(μμ) = 1931 ± 35 pb because the z-ee fit lets a ±5.9 % theory normalisation of its own signal
template float against its signal strength (`fitting/CONVENTIONS.md` §3 forbids exactly this).
Diagnosis: [`combination/docs/01-inputs.md`](combination/docs/01-inputs.md); the list of changes to
ask for: **[`ASK_Z_EE.md`](ASK_Z_EE.md)**. Without z-ee the combination is 1931 ± 35 pb with
p = 0.54. Orthogonality of the three selections, the systematics CMS and ATLAS carry that we do
not, and each channel against the published measurement of the same decay:
[`combination/docs/05-vs-published.md`](combination/docs/05-vs-published.md).

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
