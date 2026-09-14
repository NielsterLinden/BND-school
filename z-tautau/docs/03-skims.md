# 03 — Skims and ntuples (and the laptop bundle)

Two reduction levels, both written with uproot (no ROOT needed to read or write them).

## Level 1: NanoAOD-format skims — `scripts/step1_skim.py`, `ztautau/skim.py`

One output file per parent file, same branch names as NanoAOD (`nTau`, `Tau_pt`, …), so any NanoAOD code
reads them.

* **Preselection**: certified lumisection (data), run 278820–284044, `PV_npvsGood ≥ 1`, di-τ trigger OR
  (also in simulation), ≥ 2 τ candidates with pT > 35 GeV, |η| < 2.3, DM ∈ {0,1,10,11}, DeepTau VSjet ≥
  VVVLoose, VSe ≥ VVVLoose, VSmu ≥ VLoose. Deliberately looser than the analysis.
* **Branches** (~150 of ~1350): event info, MET + covariance + unclustered shift, MET filters, the six
  di-τ triggers, Tau (all ID/iso), Muon/Electron (veto), Jet, TrigObj (τ objects only). Simulation adds
  generator weight, pileup, L1 prefiring weights, LHE particles, GenVisTau, LHE scale/PS weights (PDF weights
  for DY NLO only) and `gen_*` truth flags (LHE flavour, m_LHE, fiducial flag).
* **`GenSums` tree** (simulation): sums over *all generated events* of the file before any cut: Σw, Σw²,
  per-flavour and 60–120 GeV sums, fiducial sum, the pileup profile and per-member sums of the
  scale/PS/PDF weights for {all, LHE ττ, LHE ττ 60–120, fiducial}. Normalisation, A and the
  acceptance-normalised theory variations need nothing else.
* **Robustness**: one subprocess per file (`ztautau/batch.py`) with a wall-clock limit; a failed or stalled
  dCache read is retried from EOS. 14 of 100 data files hit dCache I/O errors on 14 Sep 2026 and were read
  from EOS. Rerunning skips finished files. Every file has a `.root.json` provenance sidecar
  (events in/out, missing branches, timing); `skims_v1/manifest.json` aggregates them.

| sample | files | events in | events out | size |
|---|---:|---:|---:|---:|
| data Run2016G | 45 | 79,578,661 | 2,807,779 | 0.84 GB |
| data Run2016H | 55 | 76,758,754 | 2,222,148 | 0.67 GB |
| DY NLO | 41 | 71.8 M | 52,875 | 30 MB |
| everything else | | | | < 50 MB |

Location: `/data/atlas/users/nterlind/BND-school-cache/ztautau/skims_v1/` (≈ 1.6 GB total; ~20 min on 8
workers for data, ~20 min for the simulation).

## Level 2: flat ntuples — `scripts/step2_ntuples.py` (**the laptop bundle**)

One file per sample, tree `ntuple`, one row per event with a selected τhτh pair (all charges and isolations
above VVVLoose, so every region of the analysis can be rebuilt): the two τ (pT with nominal energy scale,
raw pT, η, φ, mass, DM, charge, ID bitmasks, raw DeepTau score, generator match), MET vector and covariance,
unclustered-energy shift, jets, **m_vis, m_col, m_tt**, pT_ττ, m_T^tot, ΔR, and in simulation the weights and
truth. Everything downstream (steps 3–6) reads only these.

Location: `/data/atlas/users/nterlind/BND-school-cache/ztautau/ntuples_v1/` — **≈ 310 MB for data + all
simulation**. To work on a laptop:

```bash
rsync -av stbc-i2.nikhef.nl:/data/atlas/users/nterlind/BND-school-cache/ztautau/ntuples_v1 ~/ztautau/
export BND_TAUTAU_CACHE=~/ztautau          # ntuples_v1/ must be inside
python run_all.py --from 3                  # fake factors, histograms, fits, report
```

(steps 3–4 need python ≥ 3.10 with numpy, uproot, awkward, hist, matplotlib, mplhep; step 5 needs
TRExFitter.)
