# BND School 2026 — the Z cross section on CMS Open Data

Group project of the BND Graduate School (Heimbach, 2026): a measurement of the Drell-Yan Z cross section in three
decay channels and their combination, on CMS 2016 Open Data (NanoAODv9, Run2016G+H, 16.4 fb⁻¹, √s = 13 TeV).

**The talk: [`2026.09.18 -- Presentation BND -- OpenData.pptx`](2026.09.18%20--%20Presentation%20BND%20--%20OpenData.pptx)**
(PowerPoint with the animations embedded; stored with Git LFS — `git lfs pull` after cloning, or use GitHub's
download button).

## Result

> **σ(pp → Z/γ\* → ℓℓ, 60 < m_ℓℓ < 120 GeV) = 1949 ⁺³⁰₋₃₀ pb** per lepton flavour, assuming lepton universality
> = 1949.3 ± 0.5 (stat) ⁺³⁰·²₋₂₉.₆ (syst) pb. aMC@NLO (NNLO-normalised): 1954 ⁺⁵⁶₋₈₂ pb; CMS (2024): 1952 ± 49 pb.

One TRExFitter v1.8.0 MultiFit of the three channels, shared nuisance parameters fitted jointly; the channels are
compatible (−2Δln L = 1.9 for 2 dof, p = 0.38).

| channel | σ(60–120 GeV), same likelihood | where |
|---|---|---|
| Z → ee | 2094 ⁺¹²²₋₁₁₄ pb (the channel's own fit: 1841 ± 30 pb) | [`z-ee/`](z-ee/handoff.md) |
| Z → μμ | 1931 ⁺³¹₋₃₀ pb (σ_fid = 790.2 ± 0.2 ± 6.1 ± 9.6 pb) | [`z-mumu/`](z-mumu/README.md) |
| Z → ττ (τhτh + μτh + eτh + eμ) | 1981 ⁺¹⁵³₋₁₃₆ pb | [`z-tautau/`](z-tautau/README.md) |
| combination | **1949 ⁺³⁰₋₃₀ pb** | [`combination/`](combination/combLieke/README.md) |

The value is carried by μμ. The ee inputs have two known limitations (the ECAL barrel–endcap gap is not vetoed, so the
electron-ID uncertainty is ±5.9 % instead of ~1.2 %; no electron trigger scale factor), which is why the combination
separates the ee template shapes from their normalisations — see "The ee channel" in the combination README.
Every final number with its source file: [`docs/FREEZE.md`](docs/FREEZE.md).

## Layout

| folder | what |
|---|---|
| `z-ee/`, `z-mumu/`, `z-tautau/` | the three channel analyses: `README.md` (start here), `handoff.md` (what the combination needs), `docs/`, code, fit configuration and results |
| `combination/` | the MultiFit of the three channels ([`combLieke/README.md`](combination/combLieke/README.md): method, figures, `output/result.json`) |
| `fitting/` | what the channels share: the environment (`setup.sh`), the TRExFitter helpers (`trexhist.py`, `trexconfig.py`, `run_trex.py`) and `build_trexfitter.sh` |
| `datasets/` | the dataset catalogue, golden JSON, official correction files, the z-ee histograms |
| `docs/` | shared documents: [`CONVENTIONS.md`](docs/CONVENTIONS.md) (the fit-input contract), [`FREEZE.md`](docs/FREEZE.md) (final numbers), [`UNCERTAINTY_PARITY.md`](docs/UNCERTAINTY_PARITY.md) (how a channel's uncertainties are compared with a published measurement), `OpenDataIntro.pdf` (the school's introduction slides) |
| `presentation/` | how the 99 animated clips of the talk were made (Manim): [`presentation/README.md`](presentation/README.md) |
| `.claude/` | the Claude Code setup the project was built with: `CLAUDE.md` (repository guide), `agents/`, `prompts/` (the prompts the analyses started from) |

Documentation follows one pattern everywhere: a folder's `README.md` is the way in, details are numbered pages in its
`docs/`, and a `CLAUDE.md` is the guide for an AI agent working in that folder.

## Running it

```bash
source fitting/setup.sh                # LCG_110: ROOT 6.40, python 3.13, uproot, awkward, hist, mplhep, iminuit
bash fitting/build_trexfitter.sh       # once: clones TRExFitter v1.8.0 next to the repository and builds it
python fitting/selftest.py             # the whole chain on a toy
```

TRExFitter is **not** part of this repository; every fit here was made with **v1.8.0**, and the combination relies on
all channels using that version. Each channel's README gives its own chain (`z-mumu/run_v2.py`, `z-tautau/run_all.py`,
`combination/combLieke/run.py`). The analysis is frozen (17 Sep 2026): outputs in the repository are final, and
rerunning needs the NanoAOD skims, which are not in git (`datasets/`, `z-mumu/docs/10-skims.md`).

History: the state before the final cleanup is the tag `pre-cleanup-2026-09-18`; the frozen three-channel
combination (1945 ⁺³¹₋₃₀ pb, τhτh only) is the tag `zmumu-freeze-2026-09-17`.

## Useful links

- CERN Open Data Portal: <https://opendata.cern.ch>
- CMS Open Data guide: <https://cms-opendata-guide.web.cern.ch>
- CMS Open Data forum: <https://opendata-forum.cern.ch>
- 2016 NanoAOD workshop lessons: <https://cms-opendata-workshop.github.io>
