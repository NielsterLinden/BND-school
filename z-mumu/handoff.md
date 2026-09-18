# Handoff – Z → μ⁺μ⁻

## Team

- Karel de Vries
- Eugene Shalugin
- Niels Ter Linden

## Frozen result (17 Sep 2026) -- **the channel's final result, used by the combination**

> **σ_fid(pp → Z/γ* → μ⁺μ⁻; dressed, pT > 26/20 GeV, |η| < 2.4, 60 < m < 120 GeV) = 790.2 ± 0.2 (stat) ± 6.1 (syst) ± 9.6 (lumi) pb**
> μ_Z = 0.9883 ⁺⁰·⁰¹⁴¹₋₀.₀₁₃₈ w.r.t. the aMC@NLO prediction (799.6 pb); GoF p = 0.79; counting cross-check 794.5 pb; ±0.2 % between binnings.
> **σ(Z/γ* → μμ, 60 < m_Born < 120 GeV) = 1931 ± 0.6 (stat) ± 15.0 (syst) ± 13.5 (acceptance) ± 23.4 (lumi) pb = 1931 ± 30 pb** (A = 0.4092)
> σ(m > 50 GeV) = 2002 ± 32 pb (A = 0.3947)

The analysis is **frozen**: no further reruns or improvements (17 Sep 2026, git tag `zmumu-freeze-2026-09-17`).
It is the v2 analysis of 15 Sep plus what the comparison with CMS-SMP-20-004 (`docs/16`) measured:

* the muon **reconstruction scale factor** from tag-and-probe on the unskimmed NanoAOD, 1.0001 ± 0.0013 per
  muon, applied to every simulated template; `MuonReco` is its measured ±0.27 % per event (was 1 ± 0.4 %/muon
  assigned, ±0.8 % per event);
* the **acceptance uncertainties for the 60–120 GeV denominator** of the quoted cross section (they had been
  evaluated for m > 50 GeV), with the rows the comparison showed to be missing: the boson-pT spectrum reweighted to
  our measured pT(μμ), powheg vs aMC@NLO, PS FSR, and the QED FSR model as an estimate (option c of docs/16 §5);
  total 0.61 → 0.70 %, in `zmumu/acceptance.py`.

Changes against 15 Sep: μ_Z 0.9881 → 0.9883, σ_fid syst 8.5 → 6.1 pb, σ(60–120) 1931 ± 33 → 1931 ± 30 pb. The
15 Sep result is kept as `fit/results/zmumu_v2_15sep_fit_result.json` (and `python scripts/v2_5_fit.py
--no-reco-sf` rebuilds its fit). Against CMS (2024): 1931 ± 30.7 pb vs 1952 ± 49 pb, −21 pb, −0.36σ (docs/16).

Method: own skims of the unskimmed NanoAOD (data + all MC), pileup and L1-prefiring weights,
tight-ID/iso/trigger and reconstruction scale factors from tag-and-probe fits, fake-factor non-prompt estimate,
Z-peak momentum calibration, **TRExFitter v1.8.0** profile-likelihood fit of m(μμ) (12 × 5 GeV
bins, MINOS on every parameter) with 24 nuisance parameters. Everything is in `docs/09`–`16`, `output/v2/RESULTS_v2.md`,
`CLAUDE.md`; slide decks: v2 before the review <https://claude.ai/code/artifact/8072b961-01d5-47ed-9804-ee21f70ac985>
(pre-review numbers), comparison with CMS <https://claude.ai/artifact/VCnvdCpXACmRdkXe9UrJ6A>.
Reproduce: `source ../fitting/setup.sh && python run_v2.py --from 2` (~1 h; `--from 5` for fit, report and comparison only).

### What the combination gets (see `../docs/CONVENTIONS.md`)

| item | where |
|---|---|
| fit-input histograms (TH1D + Sumw2), names `mumu_SR__<sample>[__<NP>Up/Down]`, `mumu_CRemu__...`, reconstruction SF applied | `fit/fitinputs/zmumu.root` (committed, force-added like the ττ inputs; regenerate with `python scripts/v2_5_fit.py --no-fit`; metadata in `zmumu.root.meta.json`, including the `acceptance` block the MultiFit reads) |
| TRExFitter config (POI `mu_Z`, NP names/categories) | `fit/zmumu.config` (committed) |
| workspace | `fit/results/zmumu/RooStats/zmumu_combined_zmumu_model.root` (committed; `mumu_SR` only -- `mumu_CRemu` is a VALIDATION region) |
| fit result (μ, uncertainties, grouped impacts, σ's, acceptance block, tables) | `fit/results/zmumu_fit_result.json`, `output/v2/results_v2.json` |
| skims for laptops (one ROOT file per sample, NanoAOD names, `manifest.json` with SHA-256) | `/project/atlas/users/nterlind/BND-school-skims-lite/` (bulk: `/data/atlas/users/nterlind/BND-school-cache/skims_v2/`) |

The combination (`../combination/combLieke/`, a TRExFitter MultiFit) reads the config, the fit inputs and the
acceptance block of `zmumu.root.meta.json` (`combLieke/config/channels.json`).

Counting-form numbers (for a cross-check by hand):

| variable | value | note |
|---|---:|---|
| `n_obs` | 10,378,567 | tight ID, ≥ 1 trigger-matched muon, FSR-recovered mass |
| `n_bkg` | 68,799 | prompt MC (Z→ττ 11.0k, tt̄ 31.6k, tW 2.9k, WW 3.8k, WZ 9.2k, ZZ 6.3k) + non-prompt 3.9k |
| `C` | 0.7916 | N_sel(all corrections, incl. reconstruction SF) / N_fid(dressed), aMC@NLO |
| `A_60_120` / `A_m50` | 0.4092 / 0.3947 | dressed fiducial ÷ LHE μμ (60–120) / ÷ (6077.22/3) |
| `acc_eff` | 0.3239 (60–120) / 0.3124 (m>50) | A × C |
| `lumi_pb` | 16393.381 | normtag, record 1059 |

Uncertainties (impact on μ_Z / σ), and correlation with the other channels:

| group | relative | correlated | NP names |
|---|---:|---|---|
| luminosity | 1.20% (external; profiled 1.18%) | yes | `Lumi` |
| L1 prefiring | 0.51% | yes | `L1Prefiring` |
| muon efficiency (ID 0.19, iso 0.30, trigger 0.15, reconstruction 0.27 measured) | 0.48% | with μτ_h only | `MuonID`, `MuonIso`, `MuonTrigger`, `MuonReco` |
| signal modelling (PDF, αs, scales, PS, generator) | 0.46% | yes (`SigModel_mumu` in the combination) | `PDF`, `AlphaS`, `QCDScale`, `PS_ISR`, `PS_FSR`, `SigModel` |
| MC statistics (the effective statistical limit: 10× the data statistics) | 0.36% | no | gammas |
| muon momentum scale/resolution | 0.35% | with μτ_h only | `MuonScale`, `MuonRes` |
| pileup | 0.15% | yes | `Pileup` |
| background cross sections | 0.07% | yes | `XS_TTbar`, `XS_SingleTop`, `XS_WW`, `XS_WZ`, `XS_ZZ`, `XS_DYtautau` |
| non-prompt | 0.04% | no | `FakeStat_mumu`, `FakeMethod_mumu` |
| statistical | 0.03% | no | |
| **fit total without luminosity** | **0.77%** | | |
| acceptance (σ(60–120) only, outside the fit): PDF 0.50%, boson pT 0.38%, scales 0.29%, QED FSR 0.10% (estimate), MC stat 0.05%, generator 0.04%, αs 0.03%, PS FSR 0.03% | 0.70% | PDF, αs, scales, PS FSR yes; the rest no | `Acc_PDF`, `Acc_AlphaS`, `Acc_QCDScale`, `Acc_PS_FSR`, `Acc_PTZ_mumu`, `Acc_QEDFSR_mumu`, `Acc_Generator_mumu`, `AccStat_mumu` |

Settled with the other channels: the acceptance denominator is 60 < m_LHE < 120 GeV, NLO (aMC@NLO).

### Cross-channel checks answered (15 Sep 2026)

Asked by the combination: "resolve the review/README.md shape-fit issue or supply a reviewed counting
configuration", "confirm the exact reference σ_fid^pred / A_60_120", "missing fit-input ROOT file",
and for all channels "common truth mass definition, lepton universality, event orthogonality".

- **Shape fit.** Fixed on 15 Sep (`review/README.md` §0, commit 3c95908): 12 × 5 GeV bins, two-sided
  `SigModel` inside 50 < m_LHE < 120 GeV, no smoothing, MINOS on every parameter; GoF p = 0.79, stable
  to ±0.2 % over the accepted binnings. The combination uses this fit (`../combination/combLieke/`).
  The 1-bin counting configuration exists as `fit/stab_1bin.config` + `fit/fitinputs/stab_1bin.root`
  (`python scripts/v2_5_fit_variants.py`), but it is a variation, not the measurement.
- **Exact reference numbers** (`fit/results/zmumu_fit_result.json`, `meta_json` in the fit-input file,
  all from `output/v2/gensums.json`, aMC@NLO `DYJetsToLL_M-50`, 6077.22 pb for the three flavours):

  | quantity | construction | value |
  |---|---|---:|
  | σ_fid^pred (dressed ΔR < 0.1, pT > 26/20, \|η\| < 2.4, 60 < m_dressed < 120) | 6077.22 × Σw_fid / Σw | **799.566 pb** |
  | A_60_120 | Σw_fid / Σw(LHE μμ, 60 < m_LHE < 120) | **0.409209** |
  | σ^pred(Z/γ*→μμ, 60 < m_LHE < 120) = σ_fid^pred / A_60_120 | 6077.22 × Σw(LHE μμ, 60–120) / Σw | **1953.93 pb** |
  | A_m50 | Σw_fid / Σw × 3 | 0.394703 |
  | σ^pred(Z/γ*→μμ, m > 50) | 6077.22 / 3 | 2025.74 pb |

  μ_Z multiplies σ_fid^pred; σ_fid = μ_Z × 799.566 pb does not depend on that choice (it is
  N_sig / (C · L)), only the value of μ_Z does.
- **Common truth mass definition.** The 60–120 GeV denominator of σ(60–120) is the **LHE mass**
  (Born level, pre-FSR, the two outgoing LHE leptons) in this channel and in z-tautau
  (`z-tautau/docs/06`); only the *fiducial* volumes differ (dressed leptons here, `GenVisTau` there),
  and each channel's A converts only its own σ_fid. z-ee (`z-ee/handoff.md`) has no acceptance and no
  truth definition yet; its quoted 1848.6 pb is not a 60–120 GeV Born-level number.
  The ττ reference 1944.9 pb is the **same construction** as ours: 6077.22 × Σw(LHE ττ, 60–120) / Σw
  from the same sample. The 0.47 % difference is the flavour share of the aMC@NLO sample, not a
  definition: Σw(LHE ee) / Σw = 0.33386, μμ 0.33380, ττ 0.33234 (by the same formula the ee reference
  would be 1954.1 pb). ee and μμ agree to 0.01 %; the ττ share is 0.44 % lower in the generator
  (the τ mass is in the matrix element). It is irrelevant for the combination of cross sections and
  a 0.002 % effect on a MultiFit with one shared `mu_Z` (`docs/CONVENTIONS.md` §6).
- **Lepton universality** is *assumed* in the normalisation (one 3-flavour NNLO cross section
  split by the generator's LHE flavour sums) and *tested* by the combination
  (one σ per channel in the MultiFit, compatibility test in `../combination/combLieke/README.md`). Nothing in this channel's σ_fid or σ(60–120) uses another
  flavour's cross section; Z→ττ in the SR is a background with its own `XS_DYtautau`.
- **Event orthogonality.** μμ: `SingleMuon`, `HLT_IsoMu24||IsoTkMu24`, exactly two tight muons
  (no electron or τ veto). ττ: `Tau` dataset, di-τh trigger, **vetoes** muons and electrons
  → disjoint from μμ by construction. ee: `Electron` dataset, `HLT_Ele27`, exactly two medium electrons,
  no muon veto. The only process that can enter both the μμ and the ee selection is ZZ → 4ℓ:
  in the μμ SR the ZZ → 4ℓ sample contributes 1016 events, of which **61 events** (0.0006 % of the
  10.38 M SR events) also have two OS medium electrons with 60 < m_ee < 120 GeV
  (computed from the `ZZ_4L` skim with the SR selection). No overlap treatment is needed;
  an electron veto in the SR would change nothing at the quoted precision.

### Review of v2 (15 Sep 2026) and what was fixed the same day -- `review/README.md`

Slides: `review/deck/zmumu_review.pdf` (pre-fix numbers). Findings and fixes, all in the result above:
the e-mu "excess" was the non-prompt-electron background removed by the prompt-only MC filter -> the e-mu
regions now keep all MC + W+jets (data/MC 0.992) and a same-sign e-mu region was added; the p_T^miss
mismatch is unsmeared MC jets + pileup/unclustered energy (Puppi p_T^miss agrees; not used anywhere);
the 30 x 2 GeV fit was not robust (over-constrained shape NPs, one-sided `SigModel` with a generator-cut
artefact) -> 12 x 5 GeV bins, two-sided `SigModel` built inside the common 50 < m_LHE < 120 window, no
smoothing, MINOS on all parameters, stability table in `RESULTS_v2.md` (+-0.2% between binnings that
describe the data); pileup profile matched to N_PV (`Pileup` pull -1.5 -> +0.06); `MuonReco` 0.4%/muon ->
0.8%/event; lepton-cleaned jets; luminosity quoted as the external 1.2%; `Fakes` template without Sumw2.

### Comparison with CMS-SMP-20-004 (16 Sep 2026, promoted 17 Sep) -- `docs/16-uncertainties-vs-cms.md`

Every row of the CMS systematic table (arXiv:2408.03744, Table 7; e and μ fitted together) is mapped onto
this analysis. Closed by measurements, first in a tagged fit and **part of the frozen result since 17 Sep**:
the muon reconstruction SF by tag-and-probe on the unskimmed NanoAOD, **1.0001 ± 0.0013 per muon**
(was 1 ± 0.004 assigned; `scripts/v2_7_reco_tnp.py`), and the resummation/boson-pT uncertainty of the
acceptance from our own pT(μμ) spectrum (0.38 %, absent before). The acceptance theory uncertainties are
evaluated for the 60–120 GeV denominator. Estimated with three options: the QED FSR model (a 0.12 / b 0 /
c 0.10 %; the result carries c). Result: **σ(60–120) = 1931 ± 0.6 ± 20.2 ± 23.2 pb against CMS
1952 ± 4 ± 18 ± 45 pb: −21 pb, −0.36σ**, whichever option. Plots `output/v2/plots/cms_parity_*`, `reco_*`;
the same procedure for z-ee and z-tautau: `../docs/UNCERTAINTY_PARITY.md`
(agent `.claude/agents/uncertainty-parity-auditor.md`).

### Open issues after v2

1. The pileup profile is a two-parameter (scale, bunch-to-bunch spread) fit of the luminosity-record
   profile to the N_PV distribution, not the official `puWeights` file (unreachable). Data/MC in N_PV
   agree to 3% over the bulk; the impact on σ is 0.13%.
2. The e-μ regions have no electron scale factors and no data-driven non-prompt-electron estimate; the
   same-sign e-μ region is under-predicted by 24% (W+jets MC statistics -- weights of ~36 per event --
   and no QCD sample). Validation only; `--emu-control` would promote it to a control region.
3. The data lineshape at 60-80 GeV lies between aMC@NLO and powheg (`sigmodel_lineshape.png`); the
   fit resolves this with `SigModel` = +0.66 sigma. A third generator or NLO EW corrections would tell
   which is right; the counting extraction (794.7 pb) is 0.6% above the fit for this reason.
4. Official prefiring maps, Muon-POG scale factors and Rochester corrections would replace the
   in-house versions (all three files need CERN credentials). The reconstruction term, the
   second-largest systematic on 15 Sep (0.4%/muon assigned), is measured since 17 Sep (docs/16).
5. The single-muon + jet fake-factor region is biased by the isolated trigger; the
   prescaled non-isolated paths are in the skim if someone wants to fix it.
6. The LO madgraph acceptance (3.2% lower) is a convention question, not a systematic.

## v1: what we worked on (first iteration, data-only)

- Skimmed the CMS 2016 SingleMuon Open Data (Run2016G+H) down to events with ≥2 muons: **259 GB / 324 M events → 31.5 GB / 80.2 M events**, an 8.2× reduction. Kept ~300 of the original 1363 NanoAOD branches. Output lives on dCache (path below).
- Built a complete **fiducial Z → μμ cross-section measurement** from that skim, with **no simulation anywhere** — none is available locally (`/dcache/atlas/sjankovy/BND/mc/` is empty), so efficiencies come from tag-and-probe and backgrounds from control regions in the data itself.
- Six standalone component scripts plus one orchestrator, each with its own documentation page, and a histogram behind every claim.

**Result**

> **σ_fid(pp → Z → μ⁺μ⁻) = 773.2 ± 0.2 (stat) ± 11.9 (syst) pb**

in the fiducial volume defined below. For scale: the inclusive σ(Z→ℓℓ, 60–120 GeV) at 13 TeV is ≈1975 pb, and a typical acceptance of ~0.4 for this phase space puts the expectation near 790 pb — so the number lands where it should.

### Review and cross-checks (Niels, 14 Sep 2026)

Full write-up: [`docs/08-review-and-crosschecks.md`](docs/08-review-and-crosschecks.md). Steps 1–6 are unchanged and reproduce exactly; the cross-checks are two extra scripts (`scripts/xcheck_efficiency.py`, `scripts/mc_acceptance.py`). What they found:

| Finding | Effect on σ_fid |
|---|---:|
| Reference-trigger efficiency is biased high: 0.9977 → **0.9835** with trigger-object tag-and-probe. In DY simulation the reference method is +1.05% off the truth, the TrigObj method +0.15%. | +1.5% |
| Background in the failing ID/iso probes biases both efficiencies low (medium-ID raw efficiency even *falls* with pT). Corrected by same-sign subtraction and a template fit, midpoint used. | −1.9% |
| Efficiency averaged over observed muons instead of per-event 1/ε weighting | +0.5% |
| **Luminosity lacks the normtag**: official G+H value for this GRL is **16393.381 pb⁻¹** (CMS Open Data recid 1059), 0.63% above the number all channels use | −0.6% |
| "A ≡ 1" is not exact: 1.1% of dressed fiducial events migrate out of the reco kinematic cuts (DY simulation) | +1.1% |
| Acceptance for the combination computed from DY NLO/LO, with PDF, α_s and scale uncertainties | — |

> **Revised: σ_fid = 776.9 ± 0.2 (stat) ± 14.8 (syst) pb** (dressed muons; DY NLO prediction 799.6 pb, ratio 0.972)
> **σ(Z/γ* → μμ, m > 50 GeV) = 1968 pb** with A = 0.3947 (NLO), target 2025.7 pb

The input skim was verified complete (80,191,719 of 323,952,013 parent events have ≥ 2 muons) and all 152 parent files match the CERN catalogue checksums.

## Data and simulation used

| Dataset | Type (data / MC) | Record / DOI | Notes |
|---|---|---|---|
| SingleMuon Run2016G | Data | recid 30530 | UL2016 NanoAODv9, 70 files, 149,916,849 events |
| SingleMuon Run2016H | Data | recid 30563 | UL2016 NanoAODv9, 82 files, 174,035,164 events |
| DY_inclusive_NLO | MC | recid 35669 | `DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8`, 41 files, 71,839,442 events. **Not used in steps 1–6.** Used by the review for the acceptance A and closure tests, streamed from EOS over xrootd. |
| DY_inclusive_LO | MC | recid 35671 | `DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8`, 61 files, 82,448,537 events. Same use; the LO/NLO difference in A is 3.2%. |

- **Parent data:** `/dcache/atlas/sjankovy/BND/collision_data/SingleMuon` (259.7 GB, 152 files, 323,952,013 events)
- **Skim used by the analysis:** `/dcache/atlas/kdevries/BND2026/DoubleMuonSkimmed/` (31.5 GB, 152 files, 80,191,719 events)
- **Certification JSON:** `datasets/GRL/GRL.txt`, restricted to runs 278820–284044
- **Integrated luminosity:** **16290.713420 pb⁻¹** ± 1.2%, from `brilcalc lumi -c web --begin 278820 --end 284044` — the same number the z-ee subgroup computed, reused deliberately so the three channels share a normalisation.
  ⚠️ **Review:** this brilcalc call has no `--normtag`. The normtag-corrected value for the same certified lumisections is **16393.381 pb⁻¹** (Run2016G 7653.261 + Run2016H 8740.119, from [CMS Open Data recid 1059](https://opendata.cern.ch/record/1059), cross-checked by summing `pp_2016lumibyls.csv` over `datasets/GRL/GRL.txt`). This affects **all three channels** equally.

⚠️ This is eras **G and H only**, ~16 of the ~36 fb⁻¹ CMS recorded in 2016. Don't compare raw event counts with a full-2016 analysis.

## Selection and method

- **Trigger:** `HLT_IsoMu24 || HLT_IsoTkMu24`. Both are present in all 152 files — but 39 *other* trigger branches are not (the 2016 HLT menu changed mid-data-taking; 24 paths exist only in era G, 15 only in H, and boundaries move within H). `zmumu/io.py` tolerates this; check universality before adding any trigger.
- **Object selection:** `Muon_mediumId`, `pT > 20 GeV`, `|η| < 2.4`, `pfRelIso04_all < 0.15`, `|dxy| < 0.2 cm`, `|dz| < 0.5 cm`.
- **Event selection:** golden JSON → trigger → MET filters + good PV → **exactly two** selected muons → leading muon `pT > 26 GeV` (IsoMu24 plateau) → opposite sign → `60 < m(μμ) < 120 GeV`.
- **Corrections / weights:** none applied to data. FSR photons are recovered from the NanoAOD `FsrPhoton` collection (standard H→4ℓ cuts: `relIso03 < 1.8`, `dROverEt2 < 0.012`, `pT > 2 GeV`) and added back to the dimuon system. Efficiency corrections are applied at the yield level, not per event.

**Why fiducial and not total.** Converting to an inclusive cross section needs an acceptance factor A, which is a ratio of generator-level quantities and can only come from MC or a theory calculation. With no MC available, quoting a total cross section would mean quoting a number whose dominant input was guessed. Instead the fiducial volume is defined so **A ≡ 1 by construction**:

| Requirement | Value |
|---|---|
| muons | exactly 2, opposite sign |
| leading / subleading pT | > 26 / > 20 GeV |
| \|η\| | < 2.4 |
| ID / isolation | `mediumId` / `pfRelIso04_all < 0.15` |
| m(μμ) | 60–120 GeV |
| lepton level | **dressed** (muon + recovered FSR photons) |

Anyone with a DY sample can compute A for exactly this volume and divide. The "dressed" row matters — compare against a dressed-lepton prediction, not Born level.

**Efficiencies (tag-and-probe on the Z peak, 21.5 M pairs):**

| Term | Value | Source |
|---|---:|---|
| reconstruction (per muon) | 0.9960 ± 0.0040 | **external** (CMS Muon POG) — NanoAOD has no `generalTracks`, so not measurable here |
| medium ID | 0.9775 | tag-and-probe |
| isolation | 0.9495 | tag-and-probe |
| trigger (per event) | 0.9977 | reference-trigger method — **see open issues** |
| **total (per event)** | **0.8546** | product, weighted by observed (pT, η) |

**Backgrounds (both data-driven, 0.131% of the observed yield):**

| Component | Method | Yield |
|---|---|---:|
| non-prompt / fake muons | same-sign control region, R_OS/SS = 1.0 ± 50% | 2,375 |
| flavour-symmetric (tt̄, tW, WW, Z→ττ) | eμ control region, N(μμ) = 0.5·k·N(eμ), k = 1.0 ± 50% | 11,736 |
| **total** | | **14,110** |

You need both: the same-sign method is blind to tt̄, which is the *larger* of the two. WZ/ZZ with two prompt same-flavour muons are caught by neither (~0.1%) and are absorbed into the systematic rather than subtracted.

## How to run

```bash
cd /project/atlas/users/kdevries/BND2026/BND-school/z-mumu

# one-time setup — always use the venv, never pip install into the system env
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# the whole measurement, all 152 files, ~80 s on 12 cores
.venv/bin/python run_all.py

# quick smoke test
.venv/bin/python run_all.py --max-files 4

# after changing a systematic in zmumu/config.py (steps 1–2 need not repeat)
.venv/bin/python run_all.py --from 3
```

Each component also runs standalone, e.g. `.venv/bin/python scripts/step2_efficiency.py`.
Read `docs/00-overview.md` first; every physics choice lives in `zmumu/config.py`.

## Results so far

**Cutflow (full dataset):**

| Cut | Events | Abs. |
|---|---:|---:|
| skim (≥2 muons in collection) | 80,191,719 | 100.00% |
| golden JSON | 78,912,289 | 98.41% |
| `HLT_IsoMu24 \|\| HLT_IsoTkMu24` | 47,314,794 | 59.00% |
| MET filters + good PV | 47,292,438 | 58.97% |
| exactly 2 selected muons | 11,647,158 | 14.52% |
| leading muon pT > 26 GeV | 11,477,398 | 14.31% |
| opposite sign | 11,472,359 | 14.31% |
| 60 < m(μμ) < 120 GeV | **10,777,373** | **13.44%** |

- **N(signal) = 10,763,262** after background subtraction.
- **Fit validation:** peak at **90.789 ± 0.026 GeV** (PDG 91.1876), Gaussian resolution **1.448 GeV**. The χ²/ndf is terrible (908/115) **by construction** — a single Voigt profile has no radiative tail, so the exponential bends up to absorb the FSR shoulder at 70–85 GeV. Clearly visible in `output/plots/step4_fit.png`. This is why the fit is used only to validate the peak and resolution, **not** as an alternative yield, and *not* as the signal-extraction systematic (|fit − count| would compare two different quantities).
- **FSR recovery** fires in 4.6% of selected events. It noticeably sharpens the peak but changes the integrated 60–120 GeV yield by only **0.007%** — in a wide window FSR is a lineshape effect, not a yield effect. It would matter far more in an 80–100 GeV window.

**Uncertainty breakdown — the main message of the measurement:**

| Source | Relative |
|---|---:|
| luminosity | 1.200% |
| eff: reconstruction (external) | 0.803% |
| eff: trigger (method bias, *assigned*) | 0.500% |
| muon momentum scale (*assigned*) | 0.200% |
| background estimate | 0.056% |
| **statistical (data)** | **0.031%** |
| eff: isolation / ID / trigger (T&P stat) | < 0.01% each |
| FSR recovery | 0.002% |
| **total** | **1.542%** |

With 10.8 M signal events this measurement is nowhere near statistics-limited. **More data would change nothing; better inputs would.**

**Plots** (all in `output/plots/`): `step1_cutflow`, `step1_mass_regions`, `step1_mass_fsr`, `step1_mass_wide`, `step1_kinematics`, `step1_fsr_mass_shift`, `step2_eff_{id,iso}_vs_{pt,eta}`, `step2_eff_{id,iso}_map`, `step2_eff_trigger_vs_pt`, `step2_tagprobe_mass`, `step3_control_regions`, `step3_data_vs_background`, `step3_subtracted`, `step4_fit`, `step5_systematics`, `step6_summary`.

## v1 inputs for the combination (superseded by the v2 block above)

> ⚠️ **Read this before plugging the numbers in.** The combination notebook computes
> σ = (n_obs − n_bkg) / (A·ε · L) and expects all three channels to land near
> **2025.74 pb** (= 6077.22/3), i.e. the **total** cross section, with `acc_eff` = A·ε **from MC**.
>
> **This measurement is fiducial, not total.** No MC is available locally, so the acceptance
> A could not be computed and the fiducial volume was defined so that **A ≡ 1**. The
> `acc_eff` below is therefore **ε only**. Using it as-is gives **773.2 pb**, not ~2025 pb —
> a factor ~2.6 low, which in a BLUE combination would look like a huge
> lepton-universality violation. It would be an artefact of mixing two different
> definitions, not physics.

**Revised inputs (review, see [`docs/08`](docs/08-review-and-crosschecks.md) §10).** These make the channel combinable: `acc_eff` now includes the acceptance from DY NLO simulation.

| Variable | Value | Note |
|---|---:|---|
| `n_obs` | 10,777,373 | unchanged |
| `n_bkg` | 14,110 | unchanged |
| `acc_eff` | **0.3336** | A (NLO, 0.3947) × ε (0.8544) × migration (0.9891) |
| `lumi_pb` | **16393.381** | normtag-corrected; should change for every channel |
| → σ(Z/γ* → μμ, m > 50) | 1968 pb | 1981 pb if `lumi_pb` stays 16290.713 |

Acceptance uncertainty 0.61% (PDF 0.52%, scale 0.32%, α_s 0.03%, MC stat 0.05%) — **correlated** across channels. A is defined against σ_DY/3 (all LHE μμ events, m > 50 GeV), which is what the notebook's 2025.74 pb target means; the 60–120 GeV denominator proposed below gives a 3.7% larger A for a different quantity. The three channels must agree on one, and on NLO vs LO (3.2% apart).

**Inputs as committed (steps 1–6), in the notebook's variable names:**

| Variable | Value | Note |
|---|---:|---|
| `n_obs` | 10,777,373 | OS μμ, 60–120 GeV |
| `n_bkg` | 14,110 | data-driven (same-sign + eμ) |
| `acc_eff` | **0.8546** | **ε only — A is NOT included** |
| `lumi_pb` | 16290.713420 | same as z-ee ✓ |

**To make this combinable**, someone needs to compute A for the fiducial volume below from a DY
sample — the z-ee group already uses `DY_inclusive_NLO` (recid 35669) and `DY_inclusive_LO`
(recid 35671):

```
A = N(gen-level dressed muons passing the volume) / N(gen-level, 60 < m < 120 GeV)
```

Fiducial volume: exactly 2 opposite-sign muons, pT > 26 / 20 GeV, |η| < 2.4,
60 < m(μμ) < 120 GeV, **dressed** lepton level (muon + recovered FSR photons — use dressed,
not Born, or the comparison is inconsistent at ~1%). Then `acc_eff = A × 0.8546`.

As a sanity check on the eventual A: matching 2025.74 pb would require A ≈ 0.38, which is a
plausible acceptance for these pT/η cuts. **That is a consistency remark, not a measurement —
do not feed 0.38 back in as an input, it assumes the answer.**

**Systematics, for the covariance matrix.** The notebook takes luminosity/acceptance/pileup as
fully correlated between channels and lepton ID/trigger as uncorrelated:

| Source | Relative | Correlated across channels? |
|---|---:|---|
| `lumi` | 1.200% | **yes** — same luminosity, same calibration |
| `lepton_reco` | 0.803% | no — muon-specific (external Muon POG input) |
| `lepton_trigger` | 0.500% | no — muon-specific, and *assigned* not measured |
| `lepton_scale` | 0.200% | no — muon-specific |
| `background` | 0.056% | no — different methods per channel |
| `lepton_id`, `lepton_iso` | < 0.01% | no — muon-specific |
| `stat` | 0.031% | no |
| `acceptance` | **not evaluated** | would be **yes** (same MC/PDF) — missing because A was not computed |

Note that `acceptance` is absent from my list precisely because A was never computed. Once it is,
its uncertainty must be added and treated as correlated with the other channels.

## Open issues / next steps

Roughly in order of how much they'd improve the result.

1. *(Review: addressed — measured with trigger-object tag-and-probe on the parent NanoAOD, which still has `TrigObj_*`; event efficiency 0.9835, see docs/08 §3. Not yet adopted in steps 1–6.)* **Trigger efficiency is the weakest link.** The skim does not contain the `TrigObj_*` branches, so per-muon trigger-object matching is impossible. The fallback is a reference-trigger method using `HLT_Mu50/IsoMu27/Mu27/Mu45_eta2p1` — but those are *not* independent of IsoMu24 (same muons, partly the same L1 seeds), so the measured 0.9977 is biased **high** and covered only by an assigned 0.5% systematic. **Fix:** re-skim including `TrigObj_*` (cheap — they're a small fraction of the file), or measure the efficiency in the orthogonal `DoubleMuon` dataset (recid available in the Open Data release).
2. **Muon reconstruction efficiency is an external input** (0.996 ± 0.004, Muon POG), contributing 0.8%. NanoAOD has no `generalTracks`, so it cannot be verified here at all — it needs MiniAOD. Currently the second-largest systematic.
3. **Replace the Voigt with a Crystal Ball ⊗ Breit-Wigner** so the fit describes the FSR tail. Then |fit − count| *would* be a meaningful cross-check and the χ² would be interpretable.
4. **Fix k = ε_μ/ε_e in the eμ background** (currently 1.0 ± 50%). The clean way is to take it from the `SingleElectron` dataset (recids 30529/30562) via the ratio of Z yields — the z-ee subgroup already has that sample. Harmless now (0.056%), but it must be fixed if the selection is loosened or moved off the Z peak.
5. **Measure R_OS/SS** in a fake-enriched (inverted-isolation) region rather than assuming 1.0 ± 50%.
6. **Rochester muon momentum corrections.** The fitted peak sits ~0.4 GeV below PDG, partly from the uncorrected momentum scale. Currently a flat 0.2% systematic.
7. **L1 prefiring is not applied** — `L1PreFiringWeight` was not kept in the skim. Small for muons at |η| < 2.4, but a genuine omission.
8. *(Review: done — docs/08 §8.)* **Extend to the total cross section** once a DY sample is available: compute A for the fiducial volume above and divide. Nothing in the measurement needs redoing.
9. **Add the remaining 2016 eras** if they are ever released — this is G+H only, ~16 of 36 fb⁻¹.
10. *(Review: done, A = 0.3947 NLO / 0.4074 LO — docs/08 §8.)* **Compute the acceptance A** so this channel can enter the BLUE combination on the same
    footing as z-ee and z-tautau. This is now the blocker for the combination — see the
    section above. Needs a DY sample (recid 35669 / 35671) and a gen-level selection only;
    the measured part of this analysis does not change.

Added by the review, roughly in order of impact:

11. **Adopt the review corrections in steps 1–6** (or decide not to): trigger from TrigObj tag-and-probe, background-corrected ID/iso, per-event 1/ε folding, migration factor, normtag luminosity. Today they live only in the cross-check scripts.
12. **Tell z-ee and the combination about the luminosity** (16393.381 pb⁻¹, normtag) — it shifts every channel by 0.63%.
13. **Replace the ID/iso background bracket by a proper pass/fail fit** per (pT, η) bin. The same-sign vs template-fit spread is the largest efficiency systematic now (0.88%).
14. **Decide NLO vs LO for the acceptance, jointly with the other channels** (3.2% apart, correlated).

### Notes for whoever picks this up

- **Work in the venv** (`.venv/`), never `pip install` into the system or user environment.
- `output/plots/`, `output/RESULTS.md` and `output/results.json` **are** committed (~2 MB). The binary pickles in `output/data/` are git-ignored — regenerate them with `run_all.py`, it takes 80 seconds.
- One awkward-array trap worth knowing, since it cost real debugging time: indexing a jagged array with a *flat* integer array selects whole **events**, not elements within them — and a 2-D numpy index or a jagged index does the same. Use `zmumu.objects.take()`, which builds a positional mask via `ak.local_index`. It exists specifically so this mistake can't be made by accident.
