## 1. Luminosity is 0.63% low in every channel (verified)

The repo uses 16290.713420 pb⁻¹ from `brilcalc lumi -c web --begin 278820 --end 284044 -i GRL.txt -u /fb`.
That is the **online** luminosity: the command has no `--normtag`. Rerun with brilcalc 3.9.4 (brilconda on cvmfs):

| command | recorded |
|---|---|
| without `--normtag` (repo) | 16.290713420 fb⁻¹ |
| with `--normtag normtag_PHYSICS.json` | **16.393380532 fb⁻¹** |

The normtag value equals the open-data luminosity record (recid 1059: 7.653261 + 8.740119 fb⁻¹).
All cross sections in the repo are therefore 0.63% too high. One-line fix in `z-mumu/zmumu/config.py`,
`z-ee/z-ee.ipynb` and `combination/`.

## 2. z-mumu (data-only fiducial, 773.2 pb): issues with estimated impact

| Issue | Direction on σ | Size (estimate) |
|---|---|---|
| luminosity without normtag | too high | 0.63% (exact) |
| **L1 pre-firing not corrected** (events lost at L1 are invisible to tag-and-probe) | too low | **2.0%** (mean L1PreFiringWeight 0.980 for their exact selection; 1.8% muon, 0.15% ECAL) |
| e-μ control region filled from the ≥2-Muon skim, so most 1-muon e-μ events are absent | background too low → σ too high | e-μ yield 23 471 vs 74 803 in an unskimmed selection (×3.2); ≈ +0.25% |
| WZ/ZZ not subtracted | σ too high | ≈ 0.15% (simulation) |
| tag-and-probe by counting in 70–110 GeV, no background subtraction | efficiency too low → σ too high | ID: counting is 0.9% below fitted efficiency per muon (≈2% per event, tight ID); ISO: 0.1% per muon |
| no trigger matching of the tag (TrigObj not in skim) → probes that fired the trigger bias isolation up | σ too low | a few ‰ |
| reference-trigger method for trigger efficiency (acknowledged) | σ too low | ≈ 0.3% |
| "A ≡ 1 by construction" ignores migrations (resolution, FSR beyond recovery, mass-window edges) and the ID/iso rows in the "fiducial volume" table are reco-level | undefined | needs DY simulation (C factor) |

Same definition, our analysis (tight ID, MC-based C factor, full systematics): **σ_fid(dressed ΔR<0.1, pT>26/20, |η|<2.4,
60–120 GeV) = 797.2 pb** (±1.3% syst ±1.2% lumi), i.e. their 773.2 pb is 3.0% lower.

## 3. Acceptance for the z-mumu fiducial volume (the combination blocker)

From DYJetsToLL_M-50 aMC@NLO (all 71.8 M events), denominator = 60 < m_Born < 120 GeV:

| definition | A |
|---|---|
| dressed muons (GenDressedLepton, ΔR<0.1), pT>26/20, \|η\|<2.4, 60–120 GeV | **0.4092** |
| Born muons, same cuts | 0.4167 |
| Born muons, pT>25/25 (CMS-SMP-20-004 definition) | 0.3821 |

Uncertainties on A (dressed 26/20): scales 0.29%, ISR 0.07%,
FSR 0.03%, PDF ≈0.5%. Dressed vs Born differ by 1.8%: the lepton definition matters.
(Do not use "dressed 60–120" as denominator: GenDressedLepton only stores leptons above a pT threshold.)
