## Latest ee intake — reference and repaired templates integrated

The combination manifest now uses sigma_reference_pb = 1954.1032472811223 from
z-ee/z-ee.ipynb at commit 85c3a52: inclusive DY cross section times the generator-weighted
LHE ee fraction in 60 <= m < 120 GeV, divided by the full Runs-tree weight sum.
The generated ee config includes constant reference_ee = 2000 / 1954.1032472811223.
All six current datasets/z-ee ROOT paths are used directly; sizes and SHA256 hashes
are recorded in inputs.json. All 6 nominal and 44 variation histograms pass the audit,
with no negative-bin warnings. The delivery floors nonpositive contents to 1e-6 by
default and preserves variances; its effect remains part of the measurement review.
The updated archive reports mu_signal = 1.05087 +/- 0.0202394.

The missing-reference and negative-bin-presence blockers are resolved. ee measurement
review (trigger treatment, inclusive-DY theory envelopes, flooring effects) and common
observable/event orthogonality remain pending. Acceptance uncertainty stays deferred
until after nominal MultiFit. No fit has been run. Earlier intake statements below
are historical and superseded by this update.

## Acceptance workflow decision — 16 September 2026

User instruction supersedes the earlier requirement for external acceptance refits:
first make the nominal MultiFit work, then evaluate acceptance following the channel analyses.
`inputs.json` now selects `acceptance.method: postfit_deferred`. `combine.py --run`
runs the nominal fit without requiring signed acceptance modes. Its `fit_result.json`
records `acceptance_status: pending_postfit`, with in-fit errors only. It does not
claim a final acceptance-inclusive result or assign a zero acceptance error.

μμ (`scripts/v2_5_fit.py`) extrapolates sigma_fid/A and assigns sigma × A_rel_unc
outside TRExFitter. ττ (`scripts/step6_report.py`) reports sigma × sqrt(sum(A_unc²)
+ A_mc_stat²) separately; this report includes MC stat omitted from its step5 theory-only
total. ee has not supplied a corresponding truth-window acceptance prescription.
The combination's nominal reference normalization already incorporates the nominal
acceptance conversion: do not divide the common cross section by acceptance again.
Propagating channel acceptance errors to the combined result is deferred until the
MultiFit is working; channel sensitivities and correlations must then be accounted for.

The existing `external_refits` implementation remains optional, not the default or a
prerequisite for nominal fitting. The sections below describing mandatory signed
acceptance inputs apply only if that optional mode is explicitly selected.
Nominal normalization and other fit-input requirements still apply.

## Delivery update — 16 September 2026 (supersedes older input-status statements below)

The μμ ROOT file and sidecar are now present. Sidecar sigma_fid_pred/A_60_120
matches the configured 1953.9313466473218 pb reference; 5 GeV bins and SigModel are recorded.
The ee config now contains 11 systematic blocks, and all six ROOT files are present
with ASCII names DY_tautau.root and Wjets.root. The adapter retains these inputs and
maps ee LUMI to NuisanceParameter Lumi, shared with μμ/ττ. Other same-name nuisances
retain the existing correlation policy; the ee inclusive-DY PDF/scale envelope construction
still needs physics review before treating it as the same response as the other channels.

Read-only uproot validation through repository setup.sh checked μμ 18 nominal/236 variation,
ee 6 nominal/44 variation, and ττ 39 nominal/1186 variation histograms: all names, binning,
finite values and explicit MC Sumw2 pass. ee has 18 templates with negative bins;
DY_tautau and Wjets retain negative bins after Rebin=2 (nominal: 3 and 13 respectively).
No input histogram was modified. Review the background model/TRExFitter handling before fitting.

The updated ee archive has mu_signal=1.04882 ± 0.0189931 and both HistoSys/OverallSys;
it is no longer statistical-only. Handoff records corrected luminosity 16393.381 pb^-1
and PU/prefiring/electron reco/ID SFs. Its quoted 2124.6 pb is numerically consistent
with mu_signal × (6077.22/3); this does not establish the agreed truth 60–120 GeV reference.
No supplied weighted truth-window denominator/acceptance was found. Trigger efficiency
treatment also remains to be reviewed. Older statistical-only archive values are historical.

Remaining blockers: ee truth-window reference and extraction review, negative-bin treatment,
common observable/event overlap validation, and signed external acceptance modes/correlations.
The local TRExFitter executable is not built. No real fit was run.
The framework now audits histogram contents on each preparation and records systematic/NP
mappings in status.json. If uproot is unavailable, preparation records a blocker and asks
for setup.sh rather than claiming histogram validation succeeded.

# Combination handover — 15 September 2026

## Goal and scope

Build a simultaneous TRExFitter **MultiFit of Z→μμ, Z→ee and Z→ττ**, then quote a common
per-flavour cross section assuming lepton universality. Work only inside:

```
/project/atlas/users/lgijsen/BND-school/combination/combLieke/
```

The environment sometimes resolves this as `/project/atlas/Users/lgijsen/...`.
All assistant-created/modified files have stayed in `combLieke`. Channel files and the existing
combination outside this directory were read only. The user explicitly checked this restriction.
The user performed several git pulls; reread current channel inputs before resuming.

**No real TRExFitter fit has been run. No combined numerical result has been produced here.**
At the last test run, **11 unit/integration checks passed**. These include simulated acceptance
refits, not a validation against a running TRExFitter binary. Last input audit still found the
missing μμ ROOT file and incomplete ee normalization/systematics and acceptance responses.

## Important decisions

### Common observable and normalization

- Target: σ(pp→Z/γ*→ℓℓ), per lepton flavour, 60 < m_truth(ℓℓ) < 120 GeV, at 13 TeV.
- This is not the sum of three differently defined fiducial cross sections, nor the total Z
  production cross section before branching fractions.
- Agree the exact common Born/dressed truth denominator and verify event orthogonality.
  `common_observable_and_overlap_validated` remains false.
- Shared parameter: `mu_comb`, with σ_comb = mu_comb × 2000 pb. The 2000 pb is a numerical
  scale, not a theory constraint. Range [0,3] is used in all adapted channel configs.
- Channel references differ. A fixed `reference_<channel>` NormFactor on its POI samples is
  `2000 / sigma_reference_pb`, so a common mu_comb means a common cross section.
- μμ reference: **1953.9313466473218 pb**, exactly sigma_fid_pred/A_60_120 from its result JSON.
- ττ reference: **1944.8826556737133 pb**.
- ee reference: **unknown**. Do not use 6077.22/3 blindly: that refers to M>50 GeV, not the
  agreed truth 60–120 window. A reconstructed 60–120 mass cut does not establish this reference.
- Different fiducial cuts are allowed if each channel's extrapolation is defined consistently.
  Numerical rescaling does not fix incompatible physical observable definitions.

### Acceptance uncertainties: user chose external refits (option 2)

We originally blocked fits unless acceptance uncertainties were inside each likelihood.
The user explicitly chose a separate external evaluation instead, and the code was changed:

1. Fit the nominal combination at fixed nominal acceptances.
2. For each independent acceptance source/mode, vary affected channels together with signed
   responses according to the chosen correlation model, rebuild, and refit.
3. Report signed cross-section shifts and a separate acceptance uncertainty.
4. Optionally report an **approximate** quadrature total with the nominal in-fit error.

For channel i the varied signal factor is:

```
(2000 / sigma_reference_i) * (1 + relative_acceptance_shift_i)
```

Existing C-only theory nuisance parameters stay in the likelihood and are reprofiled.
At fixed yield, raising acceptance by 10% lowers σ by a factor 1/1.1.
Per source, take positive/negative envelopes of its two signed cross-section shifts, then add
independent source errors in quadrature separately for the up/down directions.

**Approximation:** A and C can respond to the same theory variation. This external prescription
does not implement their full joint correlation. `approximate_total` assumes independence
between external acceptance and in-fit errors; this is explicitly labelled in the output.
Do not double count acceptance inside the likelihood or add a total on top of its components.

`acceptance.sources` is still empty; no signs/correlations have been invented from magnitudes.
The README documents the full JSON format with a synthetic example. All modes must explicitly
list mumu, ee, tautau (zero if unaffected). One source encodes correlated channel responses;
separate sources are independent. Partial correlations need a documented mode decomposition.

## Implementation and commands

Files:

- `combine.py`: adapts source configs, audits inputs, builds workspaces, runs MultiFit and refits.
- `acceptance.py`: validates external modes and aggregates signed responses/errors.
- `ee_archive.py`: reads archived fit text/XML directly without extracting/executing members.
- `inputs.json`: live manifest, paths, references, readiness and external acceptance model.
- `README.md`: detailed method, channel comparison, required inputs and command guide.
- `mumu_inputs.md`, `tautau_inputs.md`, `ee_inputs.md`: per-channel intake notes.
- `tautau_inputs.json`, `ee_inputs.json`: recorded reference metadata, not executable
  acceptance variation models.
- `test_combine.py`, `test_acceptance.py`: checks.

From the repository root:

```bash
# Prepare configs and list requirements; standard-library Python only.
python3 combination/combLieke/combine.py

# Tests.
python3 -m unittest discover -s combination/combLieke -p 'test_*.py'

# Once nominal inputs and TRExFitter are available: preliminary fixed-acceptance fit.
source setup.sh
python3 combination/combLieke/combine.py --run --nominal-only

# Once external acceptance inputs are also complete: full nominal + acceptance refits.
python3 combination/combLieke/combine.py --run
```

TRExFitter **v1.8.0** is required by the repository. Source `setup.sh` for the environment;
follow repository build instructions if necessary. The local submodule was empty at initial
inspection, and ROOT/uproot were not importable in the default Python. Recheck current state.
Do not install into the system/user Python environment (repository instruction).

Preparation outputs:

- `generated/{mumu,ee,tautau}.config`, `generated/multifit.config`.
- `generated/ee_archive.json` when the manifest names the archive.
- `status.json`: separate nominal readiness and external-acceptance readiness.

Fit outputs:

- Nominal `results/`, `logs/`, `fit_result.json` (in-fit-only uncertainty).
- Each variation in `variations/<source>/<up|down>/`, isolated configs/results/logs.
- `result.json` only after a successful full run: nominal, signed shifts, acceptance errors,
  varied results and optional approximate total.

The runner executes `hw` per channel and `mwf` for MultiFit, rebuilding every time to avoid
old POIs/normalizations. It preserves `UseMinos: all` where present (important for μμ).
For DropBins jobs it copies the generated `*_allBinsFitRegions_combined_*_model.root` to
its conventional workspace name inside that run directory. It never substitutes the original
channel workspace for an adapted one. It does not silently drop incomplete channels.
Review fit logs, convergence, covariance, pulls, boundaries and histogram validity before
quoting results; Python tests and a parsed fit text file are not sufficient scientific validation.

## μμ inputs and status

Sources:

- `z-mumu/summary_deck/zmumu_summary.pdf` — actual filename has no extra hyphen.
  Final numbered slide 19 and review-fixes slide 17.
- `z-mumu/REVIEW.md` **section 0**, which supersedes historical findings below it.
- `z-mumu/fit/results/zmumu_fit_result.json` for exact metadata.

Result: σ_fid = **790.1 ± 0.2 (stat) ± 8.5 (syst) ± 9.6 (lumi) pb**;
mu_Z = 0.988 ± 0.016; recommended A_60_120 = 0.4092; GoF p=0.79.

**Correction to our early conversation:** the old shape-fit problem was fixed. Current result
uses 12×5 GeV bins, two-sided SigModel and MINOS on all parameters; stable to about ±0.2%.
Use the current shape config, not the historical counting fallback. Do not add the historical
extra 0.7% lineshape term: the revised SigModel covers it. `measurement_ready` is true.

**Missing physical file, checked again after the user asked:**

```
z-mumu/fit/fitinputs/zmumu.root
```

The expected config is present, but the last filesystem search found no ROOT files/archives
under z-mumu and no *mumu*.root elsewhere in this repository. Documentation saying what the
channel delivers does not mean the binary exists locally. The handoff calls it uncommitted.
The user pasted `docs/15-combination-inputs.md`; it confirms the required file contract:

- TH1 templates with Sumw2, `mumu_SR__<sample>[__<NP>Up|Down]`, `mumu_CRemu__...`.
- `meta_json` inside and/or `zmumu.root.meta.json` alongside: L, C, sigma_fid_pred,
  A_60_120, A_m50, acceptance uncertainties and counting values.
- `fit/zmumu.config`.
- Optional original workspace `fit/results/zmumu/RooStats/zmumu_combined_zmumu_model.root`.
- `fit/results/zmumu_fit_result.json` is available but cannot replace the histogram likelihood.

Preferred next action: obtain the matching ROOT file and sidecar from the μμ group or point
`histogram_dir` to its actual location. The handoff gives regeneration from existing analysis
outputs: `cd z-mumu; python scripts/v2_5_fit.py --no-fit` (requires environment and upstream
outputs; do not assume they exist). Original channel writes would be outside the user's scope;
ask or arrange an isolated output if regeneration is needed. Signed acceptance responses and
cross-channel correlations are still needed.

## ττ inputs and status

Sources: `z-tautau/output/RESULTS.md`, `output/results.json`, `fit/ztautau.config`,
`fit/fitinputs/ztautau.root.meta.json`, `handoff.md`, `fit/results/ztautau_fit_result.json`.

Current v2.1 nominal: **genuine-tau-subtracted fake factors (`mcsub`)**, three BDT categories
`tautau_SR0`, `tautau_SR1`, `tautau_SR2`. `nosub` is a cross-check.

- mu_Z = **1.205 +0.145 -0.125**.
- σ_fid = **5.421 ± 0.091 (stat) ± 0.601 (syst) pb**.
- σ_60_120 = **2343 ± 39 (stat) ± 260 (syst) ± 86 (acceptance) pb**.
- Exact central σ = 2342.6306075855446 pb; in-fit +282.3463946547757 / −243.13756031639358 pb.
- GoF p = 0.21943; L=16393.381 pb^-1, 1.2% lumi error already in fit.
- A=0.0023140314215664886; C=0.06542356849711156.

Only fiducial `DYtautau` is scaled by the POI. `DYtautau_nonfid` is background, preserved by
our adapter and external variations. This follows the channel's extraction prescription,
not a universal model tying all nonfiducial DY to σ. Preserve category DropBins (SR0 has bins
1–7 dropped); reported prefit sums need not equal sums over actually fitted bins.

**Available now:** nominal and nosub ROOT files, metadata, conventional-name workspace and
fit-result JSON, committed in aa5be6b. Nominal ROOT size was 782677 bytes. Earlier missing-file
statements were corrected; git-aware file listing once omitted a locally present ROOT file.
Use actual `Path.is_file()` checks. Binary histogram contents have not been validated in ROOT.

Latest check at repository 7a4e540: nominal fit, regions, prediction and combination metadata
were unchanged from the recorded intake. The newest top-level commit also updated the older
combination outside combLieke; we did not edit or substitute that implementation.

Correlation clarification: **SigModel_tautau is reported, not fitted** in the current config.
Do not insert it or correlate it with μμ SigModel. Handoff recommends shared Lumi, Pileup,
L1Prefiring, QCDScale, PDF, PS_ISR/FSR and matching background cross-section sources. Tau/fake
parameters are channel-specific. It recommends correlated acceptance but supplies no signed
three-channel response model. Sample/region-specific details still need review.

Unsigned relative acceptance magnitudes (exact in intake JSON):

- scales 3.3771%, PDF 0.4877%, alpha_s 0.5864%, ISR 0.9320%, FSR 0.2467%, MC stat 0.7121%.
- Fit-result JSON `acceptance_rel_unc=0.03593866580978255` includes theory only.
- Including acceptance MC stat gives **0.03663730863528309**, consistent with 3.7% / 86 pb.
- Record this bookkeeping difference; use one consistent component model, not both totals.
- None of these magnitudes have been silently converted to signed acceptance sources.

## ee inputs and status

Available original nominal templates in `datasets/z-ee/`:

```
Data.root  DY_ee.root  DY_ττ.root  Diboson.root  ttbar.root  W+jets.root
```

Current source `z-ee/fit.config`: h_mass in SR, Rebin=2, POI mu_signal on DY_ee, no systematic
blocks. Our adapter fixes the original other-user HistoPath, renames SR→ee_SR and
mu_signal→mu_comb, and changes DY_ee from BACKGROUND to SIGNAL.

The user pulled **`z-ee/Zee_fit.tar.gz`** after it was initially absent. Now integrated via
`result_archive` in inputs.json and the read-only archive reader. It contains processed
histograms, XML, a RooStats workspace, fit results, errors, plots and yield tables.

Verified archived result:

- mu_signal = **0.912567 ± 0.000662336**.
- Data-statistical contribution = **0.000411409**.
- MC-statistical contribution = **0.000519068**.
- Other systematics = **0**; 60 bin-statistical gamma parameters.
- XML enables MC-stat errors on all five MC samples, Lumi=1 and LumiRelErr=0.
  Lumi=1 here is an event-template bookkeeping setting, not physical integrated luminosity.
- `FullSyst` in the grouped output means gammas only, not complete experimental systematics.

Use the original nominal histograms for rebuilding, not the archived postfit templates.
Archived workspace uses mu_signal/SR and is a reference check, not plug-compatible with our
adapted normalization/POI. Archive XML/text have been inspected, but direct binary histogram
contents and Sumw2 have not been independently read (uproot unavailable in default Python).

**Still required for final combination:** per-flavour truth 60–120 reference, acceptance,
nominal electron ID/reco/trigger and pileup corrections, consistent luminosity and their
uncertainties, background/theory systematics, and external acceptance responses/correlations.
The handoff previously used luminosity 16290.713420 pb^-1 vs shared normtag 16393.381 pb^-1;
check whether the next ee delivery corrects this. Changing LumiLabel does not rescale events.
`measurement_ready` remains false; `sigma_reference_pb` remains null.

The user asked whether to fit now or first add ee systematics/SFs. Advice: **finish ee corrections
and systematics for the final result**; its statistical-only model would otherwise receive
excessive weight and could bias the central result. A preliminary machinery test is possible
once nominal references/templates are available, but must be labelled incomplete. The current
runner intentionally still blocks an incomplete nominal model rather than treating it as final.

## Tomorrow's practical checklist

1. Check current git state and reread updated handoffs/configs. Do not overwrite user edits.
2. Locate/obtain the actual missing μμ ROOT templates (documentation alone is insufficient).
3. Obtain ee truth-window normalization and corrected/systematic inputs; update manifest only
   using supplied values with provenance.
4. Agree common observable and orthogonal fitted data selections; different datasets alone
   do not establish absence of overlap. Review any shared MC-stat contributions too.
5. Populate signed external acceptance source responses and correlation descriptions. Resolve
   tau acceptance MC-stat bookkeeping; do not guess missing directions.
6. Check/build TRExFitter environment, prepare/audit, then run nominal and external refits.
7. Validate actual fits and uncertainty accounting before quoting a combined cross section.

## Environment/tool limitation

Shell commands in the sandbox repeatedly failed before starting:
`bwrap: Creating new namespace failed: No space left on device`.
Commands were run with explicit outside-sandbox approvals. This is an execution-environment
issue, not evidence of a bad analysis file. No dependencies were installed. All subsequent
file writes and generated artifacts remained inside combLieke.
