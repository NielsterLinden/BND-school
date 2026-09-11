# BND School 2026 – CMS Open Data: Z boson decays

Shared workspace for our group project at the BND Graduate School (Heimbach, 2026).
We analyse CMS 2016 open data (NanoAOD, 13 TeV) and split into three subgroups by Z decay channel.

## Layout

| Folder | Channel | Hand-in |
|---|---|---|
| `z-ee/` | Z → e⁺e⁻ | `z-ee/handoff.md` |
| `z-mumu/` | Z → μ⁺μ⁻ | `z-mumu/handoff.md` |
| `z-tautau/` | Z → τ⁺τ⁻ | `z-tautau/handoff.md` |
| `reference/` | Shared material (intro slides, links) | – |

Each subgroup owns its folder: put notebooks, scripts and plots there.
Keep `handoff.md` up to date so the other groups can pick up where you left off.

## Rules of the road

- Work directly on `main`. Pull before you push (`git pull --rebase`).
- Do **not** commit data files (ROOT files, large CSVs). Copy them locally and add the path/DOI to your `handoff.md` instead.
- Strip notebook outputs before committing if they are large.

## Useful links

- CERN Open Data Portal: <https://opendata.cern.ch>
- CMS Open Data guide: <https://cms-opendata-guide.web.cern.ch>
- CMS Open Data forum: <https://opendata-forum.cern.ch>
- 2016 NanoAOD workshop lessons: <https://cms-opendata-workshop.github.io>
