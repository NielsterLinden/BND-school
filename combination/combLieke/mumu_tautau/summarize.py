"""Summarize saved ROOT diagnostics; never use nonfinite error-decomposition entries."""
import json
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
read = lambda name: json.loads((HERE / name).read_text())
d = read('fit_diagnostics.json')
c = read('poi_minos_check/fit_diagnostics.json')
r = read('fit_result.json')
p = next(p for p in d['parameters'] if p['name'] == 'mu_comb')
pc = next(p for p in c['parameters'] if p['name'] == 'mu_comb')
assert c['status'] == 0 and c['covariance_quality'] == 3
assert all(abs(p[k] - pc[k]) < 1e-9 for k in ['value', 'minos_up', 'minos_down'])
assert not p['at_boundary']
assert {x['name']: x['entries'] for x in d['workspace_regions']} == {
    'mumu_SR': 12, 'tautau_SR0': 7, 'tautau_SR1': 14, 'tautau_SR2': 14}
log = (HERE / 'logs/multifit.log').read_text()
gof = float(re.findall(r'probability = ([0-9.eE+-]+)', log)[-1])
notes = [
    'Nominal MINIMIZE and HESSE status 0, covariance quality 3; POI is not at a boundary.',
    'All-parameter MINOS status 160: gamma_stat_tautau_SR2_bin_0 is at its lower boundary with a very large upper error.',
    'An identical-model cross-check requesting MINOS only for mu_comb gives status 0 and identical central value and profile interval.',
    'Automatic statistical/systematic error decomposition contains NaN and is not usable; quote only the full in-fit MINOS interval.',
    'Acceptance uncertainty is deferred; common truth-definition and selected-event orthogonality reviews remain pending.'
]
d.update(reported_gof_pvalue=gof, poi=p, poi_minos_crosscheck={k:c[k] for k in ['status','covariance_quality','edm']},
         parameters_at_boundary=[x['name'] for x in d['parameters'] if x['at_boundary']],
         error_decomposition_valid=False, valid_for_final_measurement=False, notes=notes)
(HERE/'fit_diagnostics.json').write_text(json.dumps(d,indent=2)+'\n')
r.update(mu_comb=p['value'],sigma_pb=p['value']*r['reference_pb'],error_up_pb=p['minos_up']*r['reference_pb'],
         error_down_pb=-p['minos_down']*r['reference_pb'],channels=['mumu','tautau'],
         validation='Minimum and covariance checks passed; POI MINOS verified independently. See fit_diagnostics.json for nuisance boundary and invalid decomposition.',
         reported_gof_pvalue=gof,diagnostics_file=str(HERE/'fit_diagnostics.json'),error_decomposition_valid=False)
(HERE/'fit_result.json').write_text(json.dumps(r,indent=2)+'\n')
params=sorted([v for v in d['parameters'] if v['name'].startswith('alpha_')],key=lambda v:abs(v['value']),reverse=True)
fig,ax=plt.subplots(figsize=(9, max(6,len(params)*0.23)))
ax.axvspan(-2,2,color='#fff0b3');ax.axvspan(-1,1,color='#d8ebcb');ax.axvline(0,color='gray',lw=0.8)
ax.errorbar([v['value'] for v in params],range(len(params)),
            xerr=[[-v['minos_down'] for v in params],[v['minos_up'] for v in params]],fmt='o',ms=3,color='#243b53',capsize=2)
ax.set_yticks(range(len(params)),[v['name'].removeprefix('alpha_') for v in params],fontsize=8)
ax.invert_yaxis();ax.set_xlim(-3,3);ax.set_xlabel('Nuisance pull and MINOS interval');ax.set_title('Z → μμ + Z → ττ combined fit')
fig.tight_layout();fig.savefig(HERE/'results/combLieke/nuisance_pulls.png',dpi=160);fig.savefig(HERE/'results/combLieke/nuisance_pulls.pdf');plt.close(fig)
(HERE/'README.md').write_text(f'''# Z → μμ + Z → ττ MultiFit

TRExFitter v1.8.0 run on 16 September 2026, using only the μμ and ττ likelihoods.

**σ(pp → Z/γ* → ℓℓ, 60 < m_truth < 120 GeV) = {r['sigma_pb']:.1f} +{r['error_up_pb']:.1f} / −{r['error_down_pb']:.1f} pb**, per flavour assuming lepton universality.
μ_comb = {r['mu_comb']:.8f} +{p['minos_up']:.8f} / {p['minos_down']:.8f}; σ = 2000 × μ_comb pb.
The reported interval includes the full in-fit uncertainty and fixes nominal acceptance.
The TRExFitter goodness-of-fit probability is **{gof:.6f}**.

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

''' + '\n'.join('- '+n for n in notes) + '''

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
''')
print(f"sigma = {r['sigma_pb']:.6f} +{r['error_up_pb']:.6f} / -{r['error_down_pb']:.6f} pb; p = {gof}")
