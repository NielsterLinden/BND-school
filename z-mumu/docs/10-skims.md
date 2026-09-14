# v2 skims

`zmumu/skim.py`, `scripts/v2_1_skim.py`; bulk copy `$BND_SKIM_DIR`
(`/data/atlas/users/nterlind/BND-school-cache/skims_v2`), one output file per input file,
`manifest.json` with event counts and sums of weights, `<file>.root.json` provenance per file.

## Why a new skim

The >= 2-muon skim of v1 dropped `TrigObj_*` (no trigger matching), `L1PreFiringWeight_*`
(no prefiring correction), 47 of the 56 electron branches and the `Runs` tree, and by
construction has no 1-muon events (e-mu region, single-muon fake region). The v2 skim reads
the unskimmed parents (data from dCache with the EOS copy as fallback; ttbar and single top
only from EOS) and keeps everything the analysis needs, with NanoAOD branch names, so any
NanoAOD-aware code reads it.

## Selection

Data: certified lumisections (`datasets/GRL/GRL.txt`), runs 278820-284044, `PV_npvsGood >= 1`,
`HLT_IsoMu24 || HLT_IsoTkMu24`. MC: `PV_npvsGood >= 1` only (the trigger bits are kept so the
trigger efficiency can be studied on unbiased events). Then at least one category:

| bit | category | definition | used for |
|---|---|---|---|
| 1 | A | >= 2 loose muons (global or tracker, pT > 20, \|eta\| < 2.4), >= 1 with pT > 24 | SR, SS, tag-and-probe, fake-factor measurement and application |
| 2 | B | >= 1 loose muon with pT > 24 and >= 1 electron (pT > 20, \|eta\| < 2.5, veto ID) | e-mu region |
| 4 | C | exactly 1 loose muon (pT > 24), MET < 30 GeV, mT(mu, MET) < 30 GeV, >= 1 jet (pT > 30, \|eta\| < 2.4, not the muon) | single-muon + jet fake region; **prescaled 1:10** (`event % 10 == 0`) when it is the only category, `skim_prescale = 10` |

`skim_cat` (bitmask) and `skim_prescale` are stored per event, `era` (0 = G, 1 = H) for data.

## Content

~100 branches: event (`run`, `luminosityBlock`, `event`, PV, rho, MET, PuppiMET), the 8 UL2016
MET filters, 14 HLT bits (analysis, reference and prescaled non-isolated paths; absent ones
are written as False and listed in the provenance), `Muon_*` (28 fields), `Electron_*` (13),
`Jet_*` (10), `FsrPhoton_*` (6), `TrigObj_*` **filtered to muon objects** (the only filtered
collection; all cross-indices stay valid). MC adds `genWeight`, `Pileup_nTrueInt`, the 11
`L1PreFiringWeight_*`, `LHEScaleWeight[9]`, `PSWeight[4]`, `LHEPdfWeight[103]` (DY only),
`LHEPart_*`, `GenDressedLepton_*`, `genPartFlav`, and the derived `gen_lhe_flavour`,
`gen_mll_lhe`, `gen_fid_dressed`, `gen_fid_born`, `gen_m_dressed`, `gen_pt1/2_dressed`.

Per MC file a **`GenSums`** tree holds the sums over *all* generated events (before any
selection): `sumw`, `sumw2`, `sumw_lhe_{ee,mumu,tautau}` and the 60-120 GeV versions,
`sumw_fid_dressed`, `sumw_fid_born`, the PDF/scale/PS weight sums for `all`, `lhe_mumu`,
`lhe_mumu_60_120`, `fid_dressed`, `fid_born`, and the genWeight-weighted `Pileup_nTrueInt`
histogram. This makes A, C, the theory uncertainties and the pileup profile computable from
the skim alone. The `Runs` sums of the parent are copied as well.

Format: ZSTD(5), written with `uproot.mktree/extend` (NanoAOD naming). Sizes: see
`manifest.json` (data ~50 MB per parent file, ~7 GB in total; DY NLO ~170 MB per file).

## Validation

- `n_in` of every file equals the catalogue event count; `sumw` equals `Runs.genEventSumw`.
- The DY `GenSums` reproduce the review's acceptance: A(dressed, 60 < m_LHE < 120) = 0.4092,
  LHE mu mu fraction 0.3338.
- One data file: 156 214 of 2 939 781 events kept (5.3%): 143 339 in A, 3 686 in B, 9 848 in
  C (93 199 before the prescale); the SR selection of that file gives 96 304 events, the
  same-sign region 22, the e-mu region 670 (x3 the v1 rate, as the reviewer found).

## Laptop bundle

`scripts/v2_1_skim.py` writes per-file skims; `hadd`-ed per-sample files for laptops are
produced by `scripts/v2_bundle.py` (see handoff.md for the path). Read with

```python
import uproot
ev = uproot.open("DY_NLO.root")["Events"].arrays(["nMuon", "Muon_pt", "gen_lhe_flavour"])
```
