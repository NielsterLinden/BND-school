# 01 — Data and simulation

## Collision data

| era | record | dataset | files | events | recorded L (normtag) |
|---|---|---|---:|---:|---:|
| Run2016G | 30532 | `/Tau/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD` | 45 | 79,578,661 | 7653.261 pb⁻¹ |
| Run2016H | 30565 | `/Tau/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD` | 55 | 76,758,754 | 8740.119 pb⁻¹ |

Location: `/dcache/atlas/sjankovy/BND/collision_data/Tau/` (identical to CERN EOS; a file that fails on dCache
is retried from EOS automatically). Certified lumisections: `datasets/GRL/GRL.txt`, runs 278820–284044.
**L = 16393.381 pb⁻¹ ± 1.2 %** — the normtag value shared by all three channels (`fitting/CONVENTIONS.md`).

The `Tau` primary dataset is the right one for τhτh: it is where the di-τ triggers are routed. Using it
(rather than SingleMuon/SingleElectron) makes this channel's data statistically independent of z-mumu and
z-ee, which the combination needs.

## Simulation (UL16 postVFP NanoAODv9, matches Run2016G+H)

| key | record | sample | role | σ [pb] | files used |
|---|---|---|---|---:|---:|
| `DY_NLO` | 35669 | DYJetsToLL_M-50 amcatnloFXFX | **signal** (LHE ττ, split into fiducial `DYtautau` and `DYtautau_nonfid`), Z→ee/μμ backgrounds; normalisation and acceptance | 6077.22 (FEWZ NNLO, all flavours) | 41/41 |
| `DY_0J`, `DY_1J`, `DY_2J` | 35577, 35595, 35613 | DYJetsToLL_0J/1J/2J amcatnloFXFX | signal statistics, stitched per LHE_NpNLO bin to the inclusive sample (`analysis.dy_norm`) | as `DY_NLO` | 60/60, 76/76, 75/75 |
| `DY_LO` | 35671 | DYJetsToLL_M-50 madgraphMLM | alternative generator (fiducial C factor cross-check) | 6077.22 | 24/61 |
| `DY_lowmass` | 35631 | DYJetsToLL_M-10to50 amcatnloFXFX | Z/γ*→ℓℓ with m < 50 GeV | 18610 | 25/25 |
| `WJets` | 69745 | WJetsToLNu amcatnloFXFX | W+jets with a genuine leading τh | 61526.7 | 28/28 |
| `TTTo2L2Nu` | 67801 | powheg | tt̄ | 88.29 | 20/49 (EOS) |
| `TTToSemiLeptonic` | 67993 | powheg | tt̄ | 365.34 | 14/138 (EOS) |
| `ST_tW_top/antitop` | 64895/64839 | powheg, NoFullyHadronic | single top | 19.47 each | all (EOS) |
| `WW`, `WZ`, `ZZ` | 72696/72754/75593 | inclusive pythia8 | dibosons (all decays) | 118.7 / 47.13 / 16.523 | all |

* The sample normalisation is σ·L / Σ(generator weights of the *processed* files), so processing only
  some files of a large sample (tt̄) is unbiased, just statistically poorer.
* **Jet-binned Drell-Yan stitching.** The 0J/1J/2J aMC@NLO samples have exactly 0, 1 and 2 partons in
  the NLO matrix element (`LHE_NpNLO`; *not* `LHE_Njets`, which also counts the real-emission parton and
  is therefore not exclusive: stitching in it over-counts by 7 %). Every DY event gets the weight
  σ·L · f_j / Σ_samples Σw_sample(j), with f_j = Σw_incl(j)/Σw_incl the jet-bin fraction of the *inclusive*
  sample (its FxFx merging defines the multiplicity mixture) and Σw_sample(j) the generator-weight sums
  per bin from the `GenSums` trees. The inclusive sample alone reproduces σ·L/Σw_incl; the jet-binned
  samples only enlarge the denominators, i.e. add statistics without any extra cross section. The
  acceptance A and the theory sums come from the inclusive sample; the theory variations are renormalised
  per sample to its own fiducial yield (`analysis.theory_weights`).
* **Fiducial / non-fiducial split.** 38 % of the selected Z/γ*→ττ is outside the fiducial volume
  (30 % has m_LHE > 120 GeV: the two 40 GeV visible-pT cuts enrich the γ* continuum by two orders of
  magnitude relative to the peak, `REVIEW.md` 3.1). It is a separate sample `DYtautau_nonfid` with a 5 %
  normalisation uncertainty plus the theory variations, and is *not* scaled by μ_Z.
* **No QCD multijet simulation** is used anywhere: all jet→τh fakes come from the fake-factor method in data
  (`05-fake-factors.md`).
* The inclusive diboson samples are used (not the exclusive decay samples z-mumu uses) because τh can come
  from any W/Z decay chain; one sample per process avoids double counting.
* Samples whose τ pairs come from **jets** are only partly needed: events whose leading τ is a jet are
  removed from all simulation (they are in the fake estimate). What remains of W+jets/tt̄ is the
  "genuine leading τ + jet as second τ" component.

## Why these generator choices

aMC@NLO FxFx is the CMS reference for Drell–Yan and the sample the other channels use for A. It has negative
weights (~16 %), which reduce its effective statistics to ~45 %: the signal template MC statistics (2 % in the
SR, up to 6 % per bin) is a visible uncertainty. The LO sample has no negative weights but lacks NLO accuracy
in the Z pT spectrum, which matters here (see `06-cross-section.md`).
