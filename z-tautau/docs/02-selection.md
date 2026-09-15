# 02 — Event selection

Every number below lives in `ztautau/config.py`; the code is `ztautau/objects.py` and
`scripts/step2_ntuples.py`.

## Trigger

| era | path | notes |
|---|---|---|
| Run2016G | `HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg` | absent from all H files |
| Run2016H | `HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg` | absent from all G files |
| simulation | OR of both | both exist in UL16 MC; TauPOG trigger SFs are derived for this OR |

Both offline τh must be matched to an HLT τ object: `TrigObj_id == 15`, `filterBits & 2`
("Medium(Comb)Iso" leg), trigger-object pT > 35 GeV, ΔR < 0.5 (97 % of triggered pairs match).

The HLT τ reconstruction is different in the two eras (plain vs combined isolation). This does not
matter for genuine τh (trigger SFs), but it changes which **jets** fire the trigger, so the fake factors are
measured separately per era (`05-fake-factors.md`).

## τh candidates

| requirement | value | why |
|---|---|---|
| pT | > 40 GeV (candidates kept from 38 GeV so the τ energy scale can move them) | above the 35 GeV HLT threshold where trigger SFs are reliable |
| \|η\| | < 2.1 | HLT `eta2p1` |
| \|dz\| | < 0.2 cm | primary vertex |
| decay mode | 0, 1, 10, 11 | the modes with TauPOG ID/ES SFs; the new DMs 5, 6 are not calibrated |
| DeepTau2017v2p1 VSe | ≥ VVLoose | TauPOG recommendation for τhτh |
| DeepTau2017v2p1 VSmu | ≥ VLoose | TauPOG recommendation for τhτh |
| DeepTau2017v2p1 VSjet | ≥ VVVLoose (candidate) / **Medium** (signal region) | "loose" and "tight" of the fake factor |

## Pair selection

1. all candidate pairs with ΔR > 0.5;
2. choose the pair with the most isolated τ (largest raw DeepTau VSjet score), then the more isolated
   second τ, then the larger scalar pT sum (the CMS H→ττ convention);
3. inside the pair, τ1 is the higher-pT τ (the fake factor is applied to τ1).

The pair is chosen **before** the isolation and charge requirements, so the signal, application and
same-sign regions are built from the same pairs.

## Event requirements

* certified lumisection, ≥ 1 good primary vertex, MET filters (goodVertices, globalSuperTightHalo2016,
  HBHENoise, HBHENoiseIso, EcalDeadCellTP, BadPFMuon, BadPFMuonDz, eeBadSc);
* **extra-lepton veto**: no muon (medium ID, pT > 10, |η| < 2.4, |dxy| < 0.045, |dz| < 0.2, relIso < 0.3) and
  no electron (MVA noIso WP90, pT > 10, |η| < 2.5, same IP, relIso < 0.3, conversion veto, ≤ 1 lost hit).
  This keeps the channel orthogonal to eτh, μτh, ee, μμ — essential for the combination;
* both τ matched to the trigger.

## Regions

T = VSjet Medium, L = VVVLoose and not Medium, both τ pT > 40 GeV:

| region | charge | τ1 | τ2 | use |
|---|---|---|---|---|
| **SR** | OS | T | T | the fit |
| AR | OS | L | T | application region of the fake factor |
| SS_T / SS_L | SS | T / L | T | fake-factor measurement |
| OSAI_T / OSAI_L | OS | T / L | L | OS/SS extrapolation check |
| SSAI_T / SSAI_L | SS | T / L | L | fake factor for the OS/SS check |

In simulation only events whose **leading τ is not a jet** (`Tau_genPartFlav != 0`) are used. Each region is
additionally split by the BDT score into three categories (`config.BDT_CATEGORY_EDGES`); the SR categories are
the fit regions `tautau_SR0/1/2`, the same-sign categories give the per-category closure.

## Cutflow (full data set)

| step | Run2016G | Run2016H |
|---|---:|---:|
| NanoAOD | 79,578,661 | 76,758,754 |
| skim: certified + trigger + ≥ 2 candidates (pT > 35) | 2,807,779 | 2,222,148 |
| MET filters | 2,805,181 | 2,219,869 |
| τ pair (ID, pT > 38, ΔR > 0.5) | 2,033,829 | 1,622,206 |
| both legs trigger-matched | 1,993,984 | 1,591,799 |
| extra-lepton veto (= ntuple) | 1,989,586 | 1,588,079 |
| SR (both pT > 40, OS, both Medium) | 24,106 | 23,480 |

Signal region total: **47,586 events**; the expected prefit composition (v2, MC-subtracted fake factor) is
36931 jet→τh fakes (80 %), 4827 fiducial Z→ττ, 3025 non-fiducial Z/γ*→ττ, 575 tt̄, 540 W+jets,
189 Z→ee and 122 dibosons. The signal region is fitted in three BDT categories (`09-bdt.md`).
