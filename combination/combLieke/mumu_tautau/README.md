# Z → μμ + Z → ττ MultiFit

TRExFitter v1.8.0 run on 16 September 2026, using only the μμ and ττ likelihoods.

**σ(pp → Z/γ* → ℓℓ, 60 < m_truth < 120 GeV) = 1933.7 +30.8 / −30.1 pb**, per flavour assuming lepton universality.
μ_comb = 0.96683300 +0.01542174 / -0.01506460; σ = 2000 × μ_comb pb.
The reported interval includes the full in-fit uncertainty and fixes nominal acceptance.
The TRExFitter goodness-of-fit probability is **0.462546**.

## Inputs and model

Uses the current μμ 5 GeV shape templates with two-sided SigModel and ττ Tight working point,
MC-subtracted fake factors (mcsub). The older handover's ττ numerical result is superseded by
this delivery. Reference cross sections are 1953.9313466473218 pb for μμ and
1944.8826556737133 pb for ττ. Constant channel factors 2000/reference normalize both to the
same cross section. Only fiducial DYtautau carries the ττ POI; nonfiducial DY remains background.
Same-name nuisances are shared, with SigModel_mumu channel-specific, following the existing adapter.
Acceptance is already present in the nominal normalization; do not divide by it again.

Both channel workspaces were rebuilt with `hw`, then combined and fitted with `mwf`.
The saved workspace contains mumu_SR (12 bins), tautau_SR0 (7 bins), tautau_SR1 (14 bins),
and tautau_SR2 (14 bins). The μμ validation region is excluded and ττ DropBins is preserved.
No ee fit channel is included; DYee backgrounds remain wherever the ττ model specifies them.
All 18 μμ nominal / 236 variation and 39 ττ nominal / 1186 variation histograms passed the audit.
Input hashes, git revision, and delivery metadata are in `provenance.json`.

## Numerical checks and limits

- Nominal MINIMIZE and HESSE status 0, covariance quality 3; POI is not at a boundary.
- All-parameter MINOS status 160: gamma_stat_tautau_SR2_bin_0 is at its lower boundary with a very large upper error.
- An identical-model cross-check requesting MINOS only for mu_comb gives status 0 and identical central value and profile interval.
- Automatic statistical/systematic error decomposition contains NaN and is not usable; quote only the full in-fit MINOS interval.
- Acceptance uncertainty is deferred; common truth-definition and selected-event orthogonality reviews remain pending.

The largest absolute fitted Gaussian nuisance pull is approximately 1.45 (FakeClosure_tautau_c2_hi).
The full nominal output is retained. `poi_minos_check/` is a diagnostic repeat with the same model
and workspaces and only the MINOS target list changed; it does not remove or fix any nuisance.
Treat this as a preliminary fit, not an acceptance-inclusive final measurement.

## Files

- `fit_result.json`: combined cross section and full in-fit interval.
- `fit_diagnostics.json`: statuses, covariance quality, pulls, boundaries and workspace regions.
- `results/combLieke/Fits/combLieke.root` and `.txt`: original nominal fit output.
- `results/combLieke/ws_combined.root`: combined workspace.
- `results/combLieke/nuisance_pulls.pdf` and `.png`: Gaussian nuisance diagnostics.
- `results/combLieke/CorrMatrix_comb.png`: correlation matrix.
- `results/combLieke/Fits/*errDecomp*` and `Rankings/*Breakdown*`: automatic decomposition outputs;
  these contain nonfinite entries and must not be used as a validated uncertainty breakdown.
- `generated/`, `logs/`, `run.log`: adapted configs and execution logs.

## Reproduce

From this directory:

```bash
source ../../../setup.sh
python combine.py --run --nominal-only --preliminary > run.log 2>&1
root -l -b -q inspect_fit.C > logs/diagnostics.log 2>&1
trex-fitter mwf poi_minos_check/generated/multifit.config > poi_minos_check/run.log 2>&1
(cd poi_minos_check && root -l -b -q ../inspect_fit.C > diagnostics.log 2>&1)
trex-fitter m generated/multifit_plots.config > logs/multifit_plots.log 2>&1
python summarize.py
```

The local runner is an isolated copy of the handover runner with explicit μμ/ττ channel selection
and the repository path adjusted for this folder. The parent three-channel runner and results
are unchanged. Configs use absolute paths for this workspace. External acceptance refits are
not configured in this run; the intended command is the nominal-only command above.

## Cross-section comparison plots

`output/plots/comparison_channels.pdf` and `.png` match the original publication comparison style, with current standalone μμ/ττ results and the saved MultiFit result. `comparison.pdf` and `.png` show the combined result against the published combined measurements. The asymmetric in-fit intervals are preserved; acceptance uncertainty is excluded for this work and published errors remain total. Reproduce with `source ../../../setup.sh` followed by `python plot_comparison.py`. Plot values are recorded in `output/plots/comparison_inputs.json`.
