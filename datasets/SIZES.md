# CMS Open Data — Z→ττ / Z→ee / Z→μμ: dataset sizes

13 TeV, NanoAOD. MC campaign `RunIISummer20UL16NanoAODv9` (2016 UL post-VFP — the only
NanoAOD MC campaign released). Data `UL2016_MiniAODv2_NanoAODv9`, **Run2016G + Run2016H only**.

Sizes are the full on-disk dataset size as published. Access without downloading:

```bash
cernopendata-client get-file-locations --recid <recid> --protocol xrootd
```

## Summary by group

| Group | Datasets | Events | Files | Size |
|---|---:|---:|---:|---:|
| collision_data | 12 | 1,084,099,678 | 641 | 0.99 TB |
| ztautau_signal | 31 | 957,650,363 | 1383 | 1.29 TB |
| zee | 19 | 450,519,447 | 469 | 0.55 TB |
| zmumu | 29 | 454,370,447 | 544 | 0.55 TB |
| diboson | 19 | 228,154,100 | 527 | 0.33 TB |
| top | 11 | 495,559,419 | 650 | 1.02 TB |
| wjets | 13 | 585,662,941 | 562 | 0.68 TB |
| qcd | 29 | 642,081,522 | 1015 | 1.09 TB |
| **unique total** | **143** | **4,003,744,023** | **4891** | **5.41 TB** |

> Groups overlap: the `DY_inclusive` block is counted inside ztautau_signal, zee *and*
> zmumu, so the group rows sum to more than the unique total.

## Collision data (primary datasets)

| Primary dataset | Events | Files | Size |
|---|---:|---:|---:|
| SingleMuon | 323,952,013 | 152 | 259.7 GB |
| SingleElectron | 282,385,002 | 151 | 257.1 GB |
| Tau | 156,337,415 | 100 | 159.6 GB |
| MuonEG | 63,091,128 | 48 | 67.7 GB |
| DoubleMuon | 94,148,416 | 57 | 88.3 GB |
| DoubleEG | 164,185,704 | 133 | 159.5 GB |
| **all six** | **1,084,099,678** | **641** | **0.99 TB** |

The three you asked about:

| Primary dataset | Era | recid | Events | Files | Size |
|---|---|---|---:|---:|---:|
| Tau | Run2016G | 30532 | 79,578,661 | 45 | 79.7 GB |
| Tau | Run2016H | 30565 | 76,758,754 | 55 | 79.9 GB |
| **Tau total** | | | **156,337,415** | **100** | **159.6 GB** |
| SingleMuon | Run2016G | 30530 | 149,916,849 | 70 | 119.6 GB |
| SingleMuon | Run2016H | 30563 | 174,035,164 | 82 | 140.1 GB |
| **SingleMuon total** | | | **323,952,013** | **152** | **259.7 GB** |
| SingleElectron | Run2016G | 30529 | 153,363,109 | 71 | 137.7 GB |
| SingleElectron | Run2016H | 30562 | 129,021,893 | 80 | 119.4 GB |
| **SingleElectron total** | | | **282,385,002** | **151** | **257.1 GB** |
| **Tau + SingleMuon + SingleElectron** | | | | | **676.3 GB** |

## Monte Carlo — detail


### ztautau_signal — 31 datasets, 1.29 TB


**DY_inclusive_NLO_and_LO** — 4 datasets, 248.0 GB, 225,943,598 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35631](https://opendata.cern.ch/record/35631) | 49,267,069 | 25 | 41.6 GB | `/DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35633](https://opendata.cern.ch/record/35633) | 22,388,550 | 34 | 18.0 GB | `/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35669](https://opendata.cern.ch/record/35669) | 71,839,442 | 41 | 91.2 GB | `/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35671](https://opendata.cern.ch/record/35671) | 82,448,537 | 61 | 97.3 GB | `/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DY_jet_binned_NLO** — 3 datasets, 267.8 GB, 197,984,243 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35577](https://opendata.cern.ch/record/35577) | 73,908,089 | 60 | 84.0 GB | `/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35595](https://opendata.cern.ch/record/35595) | 82,259,479 | 76 | 114.0 GB | `/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35613](https://opendata.cern.ch/record/35613) | 41,816,675 | 75 | 69.7 GB | `/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DY_HT_binned_LO** — 8 datasets, 61.4 GB, 29,716,956 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35651](https://opendata.cern.ch/record/35651) | 8,316,351 | 43 | 13.8 GB | `/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35653](https://opendata.cern.ch/record/35653) | 1,970,857 | 19 | 5.9 GB | `/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35655](https://opendata.cern.ch/record/35655) | 5,653,782 | 26 | 11.5 GB | `/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35657](https://opendata.cern.ch/record/35657) | 696,811 | 22 | 2.3 GB | `/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35659](https://opendata.cern.ch/record/35659) | 2,491,416 | 32 | 6.2 GB | `/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35661](https://opendata.cern.ch/record/35661) | 2,299,853 | 7 | 6.1 GB | `/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35663](https://opendata.cern.ch/record/35663) | 5,893,910 | 35 | 8.8 GB | `/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35665](https://opendata.cern.ch/record/35665) | 2,393,976 | 36 | 6.8 GB | `/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |

**DY_PtZ_binned_NLO** — 6 datasets, 338.6 GB, 215,168,975 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35615](https://opendata.cern.ch/record/35615) | 100,271,951 | 125 | 137.0 GB | `/DYJetsToLL_LHEFilterPtZ-0To50_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35617](https://opendata.cern.ch/record/35617) | 38,786,699 | 76 | 71.2 GB | `/DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35619](https://opendata.cern.ch/record/35619) | 12,318,492 | 47 | 28.7 GB | `/DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35621](https://opendata.cern.ch/record/35621) | 1,980,401 | 18 | 4.9 GB | `/DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35623](https://opendata.cern.ch/record/35623) | 59,792,819 | 90 | 91.7 GB | `/DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35625](https://opendata.cern.ch/record/35625) | 2,018,613 | 17 | 5.2 GB | `/DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DY_dedicated_tautau** — 4 datasets, 81.6 GB, 67,603,242 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35683](https://opendata.cern.ch/record/35683) | 11,279,282 | 73 | 14.1 GB | `/DYJetsToTauTauToMuTauh_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35681](https://opendata.cern.ch/record/35681) | 11,984,279 | 10 | 14.9 GB | `/DYJetsToTauTauToMuTauh_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17_ext1-v1/NANOAODSIM` |
| [35682](https://opendata.cern.ch/record/35682) | 10,243,463 | 7 | 12.8 GB | `/DYJetsToTauTauToMuTauh_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17_ext2-v1/NANOAODSIM` |
| [35677](https://opendata.cern.ch/record/35677) | 34,096,218 | 39 | 39.9 GB | `/DYJetsToTauTau_TauToMuEle_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DY_high_mass** — 5 datasets, 292.8 GB, 220,780,349 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35627](https://opendata.cern.ch/record/35627) | 201,894 | 9 | 0.4 GB | `/DYJetsToLL_M-1000to1500_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35629](https://opendata.cern.ch/record/35629) | 219,887,619 | 242 | 290.9 GB | `/DYJetsToLL_M-100to200_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35635](https://opendata.cern.ch/record/35635) | 202,387 | 7 | 0.4 GB | `/DYJetsToLL_M-1500to2000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35637](https://opendata.cern.ch/record/35637) | 204,033 | 2 | 0.4 GB | `/DYJetsToLL_M-2000to3000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35649](https://opendata.cern.ch/record/35649) | 284,416 | 19 | 0.6 GB | `/DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |

**EWK_Zjj** — 1 datasets, 0.8 GB, 453,000 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35999](https://opendata.cern.ch/record/35999) | 453,000 | 10 | 0.8 GB | `/EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

### zee — 19 datasets, 0.55 TB


**DY_inclusive_NLO_and_LO** — 4 datasets, 248.0 GB, 225,943,598 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35631](https://opendata.cern.ch/record/35631) | 49,267,069 | 25 | 41.6 GB | `/DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35633](https://opendata.cern.ch/record/35633) | 22,388,550 | 34 | 18.0 GB | `/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35669](https://opendata.cern.ch/record/35669) | 71,839,442 | 41 | 91.2 GB | `/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35671](https://opendata.cern.ch/record/35671) | 82,448,537 | 61 | 97.3 GB | `/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DYToEE_high_mass_powheg** — 9 datasets, 4.9 GB, 3,342,500 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35685](https://opendata.cern.ch/record/35685) | 2,500,000 | 5 | 3.5 GB | `/DYToEE_M-120To200_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35687](https://opendata.cern.ch/record/35687) | 2,500 | 1 | 0.0 GB | `/DYToEE_M-1400To2300_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35689](https://opendata.cern.ch/record/35689) | 750,000 | 4 | 1.2 GB | `/DYToEE_M-200To400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35691](https://opendata.cern.ch/record/35691) | 2,500 | 2 | 0.0 GB | `/DYToEE_M-2300To3500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35693](https://opendata.cern.ch/record/35693) | 2,500 | 2 | 0.0 GB | `/DYToEE_M-3500To4500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35695](https://opendata.cern.ch/record/35695) | 75,000 | 1 | 0.1 GB | `/DYToEE_M-400To800_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35697](https://opendata.cern.ch/record/35697) | 2,500 | 1 | 0.0 GB | `/DYToEE_M-4500To6000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35699](https://opendata.cern.ch/record/35699) | 2,500 | 2 | 0.0 GB | `/DYToEE_M-6000ToInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35701](https://opendata.cern.ch/record/35701) | 5,000 | 1 | 0.0 GB | `/DYToEE_M-800To1400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DY_high_mass** — 5 datasets, 292.8 GB, 220,780,349 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35627](https://opendata.cern.ch/record/35627) | 201,894 | 9 | 0.4 GB | `/DYJetsToLL_M-1000to1500_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35629](https://opendata.cern.ch/record/35629) | 219,887,619 | 242 | 290.9 GB | `/DYJetsToLL_M-100to200_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35635](https://opendata.cern.ch/record/35635) | 202,387 | 7 | 0.4 GB | `/DYJetsToLL_M-1500to2000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35637](https://opendata.cern.ch/record/35637) | 204,033 | 2 | 0.4 GB | `/DYJetsToLL_M-2000to3000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35649](https://opendata.cern.ch/record/35649) | 284,416 | 19 | 0.6 GB | `/DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |

**EWK_Zjj** — 1 datasets, 0.8 GB, 453,000 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35999](https://opendata.cern.ch/record/35999) | 453,000 | 10 | 0.8 GB | `/EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

### zmumu — 29 datasets, 0.55 TB


**DY_inclusive_NLO_and_LO** — 4 datasets, 248.0 GB, 225,943,598 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35631](https://opendata.cern.ch/record/35631) | 49,267,069 | 25 | 41.6 GB | `/DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35633](https://opendata.cern.ch/record/35633) | 22,388,550 | 34 | 18.0 GB | `/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35669](https://opendata.cern.ch/record/35669) | 71,839,442 | 41 | 91.2 GB | `/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35671](https://opendata.cern.ch/record/35671) | 82,448,537 | 61 | 97.3 GB | `/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**ZToMuMu_powheg_mass_binned** — 10 datasets, 5.1 GB, 3,851,000 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [75468](https://opendata.cern.ch/record/75468) | 100,000 | 5 | 0.1 GB | `/ZToMuMu_M-120To200_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75470](https://opendata.cern.ch/record/75470) | 100,000 | 1 | 0.2 GB | `/ZToMuMu_M-1400To2300_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75472](https://opendata.cern.ch/record/75472) | 100,000 | 13 | 0.2 GB | `/ZToMuMu_M-200To400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75474](https://opendata.cern.ch/record/75474) | 100,000 | 4 | 0.2 GB | `/ZToMuMu_M-2300To3500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75476](https://opendata.cern.ch/record/75476) | 99,000 | 1 | 0.2 GB | `/ZToMuMu_M-3500To4500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75478](https://opendata.cern.ch/record/75478) | 98,000 | 9 | 0.2 GB | `/ZToMuMu_M-400To800_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75480](https://opendata.cern.ch/record/75480) | 100,000 | 1 | 0.2 GB | `/ZToMuMu_M-4500To6000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75482](https://opendata.cern.ch/record/75482) | 2,955,000 | 28 | 3.4 GB | `/ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75484](https://opendata.cern.ch/record/75484) | 99,000 | 10 | 0.2 GB | `/ZToMuMu_M-6000ToInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75486](https://opendata.cern.ch/record/75486) | 100,000 | 2 | 0.2 GB | `/ZToMuMu_M-800To1400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DYToMuMu_high_mass_powheg** — 9 datasets, 4.2 GB, 3,342,500 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35751](https://opendata.cern.ch/record/35751) | 2,500,000 | 6 | 3.0 GB | `/DYToMuMu_M-120To200_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35753](https://opendata.cern.ch/record/35753) | 2,500 | 3 | 0.0 GB | `/DYToMuMu_M-1400To2300_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35755](https://opendata.cern.ch/record/35755) | 750,000 | 1 | 1.0 GB | `/DYToMuMu_M-200To400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35757](https://opendata.cern.ch/record/35757) | 2,500 | 2 | 0.0 GB | `/DYToMuMu_M-2300To3500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35759](https://opendata.cern.ch/record/35759) | 2,500 | 2 | 0.0 GB | `/DYToMuMu_M-3500To4500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35761](https://opendata.cern.ch/record/35761) | 75,000 | 2 | 0.1 GB | `/DYToMuMu_M-400To800_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35763](https://opendata.cern.ch/record/35763) | 2,500 | 1 | 0.0 GB | `/DYToMuMu_M-4500To6000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35765](https://opendata.cern.ch/record/35765) | 2,500 | 2 | 0.0 GB | `/DYToMuMu_M-6000ToInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [35767](https://opendata.cern.ch/record/35767) | 5,000 | 1 | 0.0 GB | `/DYToMuMu_M-800To1400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**DY_high_mass** — 5 datasets, 292.8 GB, 220,780,349 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35627](https://opendata.cern.ch/record/35627) | 201,894 | 9 | 0.4 GB | `/DYJetsToLL_M-1000to1500_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35629](https://opendata.cern.ch/record/35629) | 219,887,619 | 242 | 290.9 GB | `/DYJetsToLL_M-100to200_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35635](https://opendata.cern.ch/record/35635) | 202,387 | 7 | 0.4 GB | `/DYJetsToLL_M-1500to2000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35637](https://opendata.cern.ch/record/35637) | 204,033 | 2 | 0.4 GB | `/DYJetsToLL_M-2000to3000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [35649](https://opendata.cern.ch/record/35649) | 284,416 | 19 | 0.6 GB | `/DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |

**EWK_Zjj** — 1 datasets, 0.8 GB, 453,000 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [35999](https://opendata.cern.ch/record/35999) | 453,000 | 10 | 0.8 GB | `/EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

### diboson — 19 datasets, 0.33 TB


**WW** — 4 datasets, 82.1 GB, 58,511,798 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [72658](https://opendata.cern.ch/record/72658) | 19,976,139 | 30 | 30.0 GB | `/WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [72676](https://opendata.cern.ch/record/72676) | 2,900,000 | 7 | 4.3 GB | `/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [72684](https://opendata.cern.ch/record/72684) | 19,814,659 | 40 | 29.1 GB | `/WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [72696](https://opendata.cern.ch/record/72696) | 15,821,000 | 41 | 18.6 GB | `/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**WZ** — 6 datasets, 89.4 GB, 66,407,895 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [72738](https://opendata.cern.ch/record/72738) | 3,690,271 | 23 | 5.8 GB | `/WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [72740](https://opendata.cern.ch/record/72740) | 1,229,946 | 21 | 1.7 GB | `/WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [72742](https://opendata.cern.ch/record/72742) | 13,526,954 | 39 | 22.1 GB | `/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [72752](https://opendata.cern.ch/record/72752) | 10,441,724 | 31 | 15.4 GB | `/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [72748](https://opendata.cern.ch/record/72748) | 29,935,000 | 51 | 35.4 GB | `/WZTo3LNu_mllmin0p1_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [72754](https://opendata.cern.ch/record/72754) | 7,584,000 | 16 | 9.0 GB | `/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**ZZ** — 5 datasets, 120.2 GB, 85,392,407 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [75567](https://opendata.cern.ch/record/75567) | 15,928,000 | 15 | 21.5 GB | `/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75571](https://opendata.cern.ch/record/75571) | 2,468,807 | 23 | 3.2 GB | `/ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75573](https://opendata.cern.ch/record/75573) | 13,740,600 | 14 | 22.6 GB | `/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75589](https://opendata.cern.ch/record/75589) | 52,104,000 | 99 | 71.5 GB | `/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [75593](https://opendata.cern.ch/record/75593) | 1,151,000 | 17 | 1.4 GB | `/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**triboson** — 4 datasets, 34.7 GB, 17,842,000 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [72701](https://opendata.cern.ch/record/72701) | 4,159,000 | 24 | 7.9 GB | `/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17_ext1-v1/NANOAODSIM` |
| [72711](https://opendata.cern.ch/record/72711) | 4,595,000 | 7 | 9.2 GB | `/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17_ext1-v1/NANOAODSIM` |
| [72761](https://opendata.cern.ch/record/72761) | 4,554,000 | 9 | 9.2 GB | `/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17_ext1-v1/NANOAODSIM` |
| [75600](https://opendata.cern.ch/record/75600) | 4,534,000 | 20 | 8.3 GB | `/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17_ext1-v1/NANOAODSIM` |

### top — 11 datasets, 1.02 TB


**ttbar** — 4 datasets, 827.3 GB, 384,338,534 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [67731](https://opendata.cern.ch/record/67731) | 89,003,534 | 109 | 200.6 GB | `/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [67801](https://opendata.cern.ch/record/67801) | 43,546,000 | 49 | 91.2 GB | `/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [67841](https://opendata.cern.ch/record/67841) | 107,067,000 | 146 | 227.2 GB | `/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [67993](https://opendata.cern.ch/record/67993) | 144,722,000 | 138 | 308.3 GB | `/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**single_top** — 7 datasets, 187.9 GB, 111,220,885 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [64635](https://opendata.cern.ch/record/64635) | 5,471,000 | 19 | 8.7 GB | `/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [64659](https://opendata.cern.ch/record/64659) | 30,609,000 | 35 | 50.9 GB | `/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [64759](https://opendata.cern.ch/record/64759) | 63,073,000 | 86 | 104.6 GB | `/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [64839](https://opendata.cern.ch/record/64839) | 3,654,510 | 10 | 7.2 GB | `/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [64825](https://opendata.cern.ch/record/64825) | 2,554,000 | 23 | 5.0 GB | `/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [64895](https://opendata.cern.ch/record/64895) | 3,368,375 | 11 | 6.7 GB | `/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [64881](https://opendata.cern.ch/record/64881) | 2,491,000 | 24 | 4.8 GB | `/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |

### wjets — 13 datasets, 0.68 TB


**inclusive** — 2 datasets, 109.4 GB, 109,226,448 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [69745](https://opendata.cern.ch/record/69745) | 28,268,221 | 28 | 30.0 GB | `/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [69747](https://opendata.cern.ch/record/69747) | 80,958,227 | 68 | 79.4 GB | `/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**HT_binned** — 8 datasets, 102.8 GB, 63,561,129 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [69721](https://opendata.cern.ch/record/69721) | 19,753,958 | 20 | 28.1 GB | `/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69723](https://opendata.cern.ch/record/69723) | 2,090,561 | 3 | 5.7 GB | `/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69725](https://opendata.cern.ch/record/69725) | 15,067,621 | 60 | 26.7 GB | `/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69727](https://opendata.cern.ch/record/69727) | 709,514 | 17 | 2.2 GB | `/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [69729](https://opendata.cern.ch/record/69729) | 2,115,509 | 11 | 4.6 GB | `/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69731](https://opendata.cern.ch/record/69731) | 2,251,807 | 13 | 5.4 GB | `/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69733](https://opendata.cern.ch/record/69733) | 19,439,931 | 19 | 24.6 GB | `/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69735](https://opendata.cern.ch/record/69735) | 2,132,228 | 46 | 5.5 GB | `/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**jet_binned_NLO** — 3 datasets, 470.7 GB, 412,875,364 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [69715](https://opendata.cern.ch/record/69715) | 159,756,701 | 97 | 151.5 GB | `/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69717](https://opendata.cern.ch/record/69717) | 167,292,982 | 117 | 196.1 GB | `/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [69719](https://opendata.cern.ch/record/69719) | 85,825,681 | 63 | 123.2 GB | `/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

### qcd — 29 datasets, 1.09 TB


**HT_binned** — 9 datasets, 506.0 GB, 319,422,127 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [63079](https://opendata.cern.ch/record/63079) | 12,254,238 | 24 | 28.4 GB | `/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63087](https://opendata.cern.ch/record/63087) | 73,506,112 | 67 | 80.6 GB | `/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63095](https://opendata.cern.ch/record/63095) | 9,376,965 | 24 | 23.2 GB | `/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63103](https://opendata.cern.ch/record/63103) | 4,867,995 | 9 | 12.7 GB | `/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63111](https://opendata.cern.ch/record/63111) | 43,280,518 | 35 | 57.6 GB | `/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63119](https://opendata.cern.ch/record/63119) | 46,335,846 | 56 | 75.8 GB | `/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63127](https://opendata.cern.ch/record/63127) | 52,661,606 | 74 | 105.1 GB | `/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63132](https://opendata.cern.ch/record/63132) | 35,474,117 | 35 | 32.3 GB | `/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63140](https://opendata.cern.ch/record/63140) | 41,664,730 | 61 | 90.3 GB | `/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**MuEnriched** — 12 datasets, 548.2 GB, 288,891,339 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [63146](https://opendata.cern.ch/record/63146) | 13,905,446 | 19 | 38.3 GB | `/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63154](https://opendata.cern.ch/record/63154) | 19,772,747 | 31 | 33.6 GB | `/QCD_Pt-120To170_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63164](https://opendata.cern.ch/record/63164) | 4,529,975 | 19 | 4.6 GB | `/QCD_Pt-15To20_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [63174](https://opendata.cern.ch/record/63174) | 34,183,334 | 71 | 69.4 GB | `/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63186](https://opendata.cern.ch/record/63186) | 30,853,571 | 23 | 33.4 GB | `/QCD_Pt-20To30_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63192](https://opendata.cern.ch/record/63192) | 29,824,712 | 31 | 68.6 GB | `/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63204](https://opendata.cern.ch/record/63204) | 35,474,172 | 25 | 42.0 GB | `/QCD_Pt-30To50_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63216](https://opendata.cern.ch/record/63216) | 19,771,458 | 22 | 48.4 GB | `/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63222](https://opendata.cern.ch/record/63222) | 21,491,325 | 33 | 28.2 GB | `/QCD_Pt-50To80_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63228](https://opendata.cern.ch/record/63228) | 18,165,741 | 57 | 46.3 GB | `/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63232](https://opendata.cern.ch/record/63232) | 38,913,226 | 83 | 103.0 GB | `/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63238](https://opendata.cern.ch/record/63238) | 22,005,632 | 43 | 32.4 GB | `/QCD_Pt-80To120_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |

**EMEnriched** — 8 datasets, 35.1 GB, 33,768,056 events

| recid | events | files | size | dataset |
|---|---:|---:|---:|---|
| [63152](https://opendata.cern.ch/record/63152) | 5,007,347 | 7 | 6.6 GB | `/QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [63162](https://opendata.cern.ch/record/63162) | 4,026,314 | 54 | 3.2 GB | `/QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM` |
| [63172](https://opendata.cern.ch/record/63172) | 1,861,129 | 4 | 3.0 GB | `/QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [63184](https://opendata.cern.ch/record/63184) | 7,134,788 | 60 | 5.9 GB | `/QCD_Pt-20to30_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [63196](https://opendata.cern.ch/record/63196) | 1,138,742 | 11 | 2.1 GB | `/QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [63202](https://opendata.cern.ch/record/63202) | 4,351,014 | 6 | 3.8 GB | `/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [63220](https://opendata.cern.ch/record/63220) | 5,443,934 | 22 | 5.2 GB | `/QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |
| [63236](https://opendata.cern.ch/record/63236) | 4,804,788 | 9 | 5.3 GB | `/QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM` |

## Is there a Z→ll sample with multiple lepton flavours?

Yes — **`DYJetsToLL_*` is exactly that**. The `LL` means Z/γ*→ℓℓ inclusive over all three
generations, so one sample gives you Z→ee, Z→μμ and Z→ττ together. Verified directly by
reading `LHEPart_pdgId` from recid 35671 (first 50,000 events):

| Decay | Events | Fraction |
|---|---:|---:|
| Z→ee | 16,485 | 32.97% |
| Z→μμ | 16,798 | 33.60% |
| Z→ττ | 16,717 | 33.43% |

So the standard workflow for all three of your channels is **one** sample,
`/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8` (recid 35671, 82.4M events, 97 GB),
split at analysis time on gen-level flavour into the ZTT / ZL / ZJ components. The
`DYToEE_*`, `DYToMuMu_*` and `ZToMuMu_*` powheg sets are *not* alternatives to it — they
are high-mass-binned samples for the Drell-Yan tail, and the `DYJetsToTauTau*` ones are
final-state-filtered subsets for extra statistics in specific channels.

`EWKZ2Jets_ZToLL_M-50` (recid 35999) is the same idea for the electroweak/VBF Zjj production mode.

## Minimal viable starting set

| Block | Datasets | Size |
|---|---:|---:|
| data (Tau, SingleMuon, SingleElectron) | 6 | 676.3 GB |
| DY inclusive (M-50 MLM + M-10to50 MLM) | 2 | 115.3 GB |
| DY dedicated tautau (mu-tauh + e-mu) | 4 | 81.6 GB |
| diboson WW/WZ/ZZ nominal | 4 | 112.7 GB |
| ttbar (2L2Nu + SemiLeptonic) | 2 | 399.5 GB |
| single top tW (top+antitop) | 2 | 13.9 GB |
| W+jets inclusive MLM | 1 | 79.4 GB |
| **total** | | **1479 GB (1.48 TB)** |

## Reading these files remotely (no download)

`eospublic.cern.ch` serves over HTTPS with a certificate chained to **CERN Root
Certification Authority 2**, which is not in the public `certifi` bundle — so `uproot`
(and `requests`) fail with `CERTIFICATE_VERIFY_FAILED` until you point them at a bundle
that includes it. A combined bundle is checked in at `certs/cern-ca-bundle.pem`:

```bash
export SSL_CERT_FILE=$PWD/certs/cern-ca-bundle.pem
export REQUESTS_CA_BUNDLE=$PWD/certs/cern-ca-bundle.pem
```

Then open files directly, e.g.

```python
import uproot
t = uproot.open("https://eospublic.cern.ch//eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/"
                "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/"
                "106X_mcRun2_asymptotic_v17-v1/40000/14B6A8AE-C9FE-D744-80A4-DDE5D008C1CD.root")["Events"]
```

Do not disable certificate verification — the chain is legitimate, it just needs the CERN root.
