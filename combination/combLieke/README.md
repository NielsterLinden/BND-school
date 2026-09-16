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

# Three-channel TRExFitter combination (Lieke)

This directory prepares a simultaneous **μμ + ee + ττ MultiFit** using the existing channel
likelihoods. It keeps generated configs, workspaces, logs and results here. The older
`../run_combination.py` is a separate two-channel covariance combination.

**Status (16 September 2026):** μμ ROOT templates/metadata and the ee systematic
histograms/config have arrived and are integrated. All referenced nominal/variation
histograms pass presence, binning, finiteness and MC Sumw2 checks in LCG_110.
The ee truth 60–120 GeV reference/acceptance is still missing; review its trigger treatment,
inclusive-DY theory envelopes and negative DY_tautau/Wjets bins surviving Rebin=2.
Common-observable/overlap validation remains pending. Acceptance uncertainties are deferred until after the nominal fit.
No combined fit has been run; the local TRExFitter executable is not built.

## Run

From the repository root:

```bash
source setup.sh  # provides uproot for histogram content validation
python3 combination/combLieke/combine.py
python3 -m unittest discover -s combination/combLieke -p 'test_*.py'
```

Preparation writes the following files. Without uproot it still generates configs, but records
unvalidated histogram contents as a readiness blocker:

- `generated/{mumu,ee,tautau}.config`: adapted copies of the channel configs;
- `generated/multifit.config`: all three Fit blocks, `Combine: TRUE`, `POIName: mu_comb`;
- `status.json`: paths, reference cross sections and outstanding requirements.

Edit `inputs.json` when the groups deliver their inputs. Config and histogram directory paths
are relative to the repository root (absolute paths also work). `measurement_ready` records
that the channel's corrections and extraction have been reviewed;
`acceptance_in_likelihood` must be **false** for this external-acceptance workflow; true
blocks the run to prevent double counting. These are scientific declarations, not automatic validation switches:
changing them does not add corrections or uncertainties. Confirm the common observable and
check event overlap before setting `common_observable_and_overlap_validated`.

Once the input requirements are satisfied:

```bash
source setup.sh
python3 combination/combLieke/combine.py --run
```

This requires the repository's TRExFitter **v1.8.0** build. It runs `hw` separately for each
adapted channel, then `mwf` for the nominal MultiFit. The default stops there and marks acceptance as pending.
Only optional `external_refits` mode repeats these steps for each acceptance source/direction
in isolated `variations/<source>/<up|down>/` directories. Channel fits need not be performed first.
Workspaces are always rebuilt so old normalization/correlation settings cannot survive.
A missing or incomplete channel stops the run; the code never silently substitutes a
2-channel result. Each command's full output goes to `logs/`.

`fit_result.json` reports the nominal cross section and asymmetric **in-fit-only** errors in pb.
After all acceptance refits succeed, `result.json` reports the nominal result, signed source
shifts, separate acceptance errors and (if enabled) an explicitly approximate quadrature total. Review the
fit log, convergence/covariance quality, nuisance pulls and POI boundaries before quoting it;
a successful process exit or parsed text file alone is not proof of a valid fit.
The parser supports the simple HIST configurations currently in this repository; elaborate
include/replacement/multiple-file configurations need an explicit adapter.

## What the analyses currently measure

| Channel | Implementation and extraction | Combination issue |
|---|---|---|
| μμ | `z-mumu/scripts/v2_5_fit.py`, `fit/zmumu.config`: fit the 60–120 GeV mass shape with a `mu_Z` NormFactor on `DYmumu`; efficiencies, backgrounds, detector and theory variations. `sigma_fid = mu_Z * sigma_fid_pred`, then divide by `A_60_120`. The eμ region is VALIDATION and is not a fitted control region. | The final summary slide and `z-mumu/REVIEW.md` §0 confirm that the shape-fit issues were fixed (12 × 5 GeV, two-sided SigModel, MINOS on all parameters). Use the updated shape fit; the older counting recommendation and extra 0.7% lineshape term are superseded. Matching ROOT templates are still required. |
| ee | `fit.config` now includes 11 systematics; corrected templates use luminosity 16393.381 pb⁻¹ and PU/electron SFs. Updated archive: mu_signal=1.04882 ± 0.0189931. | Adapter uses new DY_tautau/Wjets filenames, maps LUMI to shared Lumi, retains all variations and renames SR/POI. Truth-window reference, trigger/theory review and negative-background-bin treatment remain pending. |
| ττ | `z-tautau/scripts/step5_fit.py`, `fit/ztautau.config`: fit reconstructed ditau mass in three BDT categories, using genuine-tau-subtracted fake factors as nominal; only fiducial DYtautau is scaled by the POI and DYtautau_nonfid is background. `mu_Z` scales `DYtautau`; `sigma_60_120 = mu_Z * 1944.8826556737133 pb`. | Acceptance includes hadronic τ branching and visible-τ cuts. The shape is essential for separating fakes. Acceptance errors are currently quoted outside the fit. |

Source normalization metadata: `z-mumu/fit/results/zmumu_fit_result.json`
(`sigma_fid_pred_pb / A_60_120`) and `z-tautau/output/results.json`
(the nominal `sigma_60_120_pb.prediction`). Treat metadata as belonging to a specific set of
templates: recheck it when templates change. Older overview files saying ee has no result
are less current than its handoff; its result remains preliminary.

## What to request from each channel

1. **Nominal binned data and every signal/background template**, including fitted control
   regions, with bin edges and MC/fake sum of squared weights. Supply ROOT files and their
   complete TRExFitter config. Different observables and binning between channels are fine;
   within each region they must match. Do not send only fitted cross sections/error bars:
   MultiFit needs likelihood workspaces, built here from templates.
2. **Systematic model:** up/down histograms or signed relative normalization shifts, sample
   and region scope, constraint type, nuisance name and definition. Explain one-sided
   variations, smoothing, pruning, fake-factor correlations and MC-stat treatment. Preserve
   Sumw2 and inspect TRExFitter's histogram validation output; no missing templates or
   invalid bins should be accepted because a source config uses `HistoChecks: NOCRASH`.
3. **Normalization metadata:** luminosity, MC cross sections, generator sum of weights,
   per-flavour `sigma_reference_pb` for the agreed truth mass window, acceptance A,
   correction/migration factor C, fiducial definition and branching conventions. The nominal
   selected signal must satisfy `N_sig = L * sigma_reference * A * C` with the channel's
   stated definitions (C may include out-of-fiducial migrations).
4. **External acceptance variations:** signed relative `A_up/A_nom - 1` and
   `A_down/A_nom - 1` for each independent source/mode and each channel. Record signs and
   cross-channel correlations. Keep the current C-only theory templates inside the likelihood.
   Do not insert acceptance nuisances into those workspaces as well. The external procedure
   approximates the relationship between acceptance and in-fit theory errors; see below.
5. **Data/selection overlap check:** runs, certification, triggers, vetoes and event identifiers
   `(run, lumi, event)` or a demonstrated orthogonal assignment across all fitted regions.
   Different primary dataset names alone do not prove disjoint events. Check MC statistical
   correlations too if the same generated events contribute to multiple selections.
6. **Channel validation:** nominal yields, independent fit result, uncertainties, goodness of
   fit and known issues. Already-built workspaces can be useful cross-checks, but this runner
   rebuilds them after applying the common POI and nuisance mapping.

For ee specifically, provide corrected templates at **16393.381 pb⁻¹** (if the same certified
G+H data are used), electron ID/reconstruction/trigger corrections and uncertainties, pileup,
L1 prefiring, luminosity, background and signal-modelling systematics, acceptance and truth
mass definition. Updating `LumiLabel` does not renormalize events. Histograms are already in
events: do not add Job `Lumi` and scale them a second time.

## Meaning of the common cross section

Assume lepton universality and define one per-flavour
`σ(pp → Z/γ* → ℓℓ, 60 < m_truth(ℓℓ) < 120 GeV)` at 13 TeV. Agree exactly on Born/dressed
truth definitions first. This is not the sum of three fiducial cross sections (their cuts
differ), nor the total Z production cross section before branching fractions.

For channel i and bin b the model is

```
E[n_ib] = mu_comb * (sigma0 / sigma_ref_i) * s_ib(theta) + backgrounds_ib(theta)
sigma_comb = mu_comb * sigma0,  sigma0 = 2000 pb
L_comb = product of channel Poisson likelihoods, with shared nuisance constraints
```

The fixed NormFactor `reference_<channel>` supplies `sigma0/sigma_ref_i` on the signal;
`mu_comb` has range [0,3] in all channels. This also scales systematic templates through the
sample normalization. 2000 pb is just a parameter scale, not a theory constraint. A different
positive scale leaves the physical result unchanged, provided the POI range still contains it.
At `mu_comb = sigma_ref_i/sigma0`, each channel recovers its original nominal signal.

The existing references differ (about 1953.93 and 1944.88 pb); simply giving both a shared
`mu_Z` would combine different cross sections. The fixed factors solve the numerical
parameterization mismatch; they do **not** fix incompatible physical definitions. For ee,
6077.22/3 refers to M>50 GeV, not automatically 60–120 GeV: obtain its weighted truth-window
fraction/acceptance before filling `sigma_reference_pb`.

DY feed-through backgrounds currently retain the channel treatment: e.g. DYtautau in μμ has
its own background uncertainty rather than the shared POI. Review that approximation under
universality, especially for ee. A fully universal model can tie the matching-mass DY
components to the common cross section after separating other mass ranges and avoiding a
duplicate background cross-section nuisance. It cannot safely be inferred from sample names.

## External acceptance workflow (option 2)

Nominal fitting no longer requires acceptance variations inside the likelihood. To run only
that first stage when nominal inputs are available, even before acceptance responses arrive:

```bash
python3 combination/combLieke/combine.py --run --nominal-only
```

This writes `fit_result.json` labelled **in-fit only; acceptance fixed**. It does not produce
an acceptance-inclusive `result.json`. The full `--run` requires complete external inputs.
Preparation without `--run` lists nominal and external-acceptance requirements separately.

Fill the `acceptance` section in `inputs.json`. This **synthetic example** shows the format;
the 1% values are illustrative, not measured channel inputs:

```json
{
  "method": "external_refits",
  "inputs_complete": true,
  "combine_with_infit_quadrature": true,
  "sources": [{
    "name": "ExampleSharedSource",
    "description": "Illustration only: one mode with fully correlated channel responses",
    "relative_shifts": {
      "mumu": {"up": 0.01, "down": -0.01},
      "ee": {"up": 0.01, "down": -0.01},
      "tautau": {"up": 0.01, "down": -0.01}
    }
  }]
}
```

Every source must explicitly list all three channels: use zeros for unaffected channels.
`up` and `down` label the underlying variation, not the direction of the cross-section shift.
An up variation can decrease acceptance in one channel while increasing it in another.

- **Correlated across channels:** put all signed responses in one source; vary them together.
- **Uncorrelated:** use separate sources, each affecting just one channel.
- **Partially correlated:** supply independent modes from a documented covariance
  decomposition, with signed responses in each channel. The code does not infer correlations.
- Different entries in `sources` are treated as independent. Do not also add a total acceptance
  uncertainty when its components are already supplied. A total such as ττ's approximately
  3.7% is insufficient to define a source-by-source correlated model on its own.

For a variation, the fixed signal factor becomes

```
(sigma0 / sigma_ref_i) * (1 + relative_acceptance_shift_i)
```

The reference cross section and nominal templates stay unchanged. All signal variations get
this same sample factor; original in-fit nuisances remain free and are profiled again. With
one channel and fixed yield, +10% acceptance gives `sigma_var = sigma_nom / 1.1`.
This implements a scalar acceptance-only response, holding C's model fixed; it does not
represent acceptance-induced bin-shape changes.

The two refits produce signed shifts `delta_up` and `delta_down` relative to the nominal
combined cross section. For each source:

```
positive_error = max(0, delta_up, delta_down)
negative_error = max(0, -delta_up, -delta_down)
```

Square and sum these errors separately over independent sources. The output retains both
signed shifts, each varied fit result, and the correlation descriptions. If
`combine_with_infit_quadrature` is true, `approximate_total` combines this external error and
the **nominal** in-fit error in quadrature. The varying fits' error bars are diagnostic and
are not added as extra uncertainties. Set the flag false to report only the separate errors.

**Approximation:** A and C may share PDF/scale/shower origins. Reprofiling C-only nuisances
while changing A does not implement their joint physical correlation. The approximate total
assumes external acceptance and in-fit errors are independent, as recorded in `result.json`.
Assess that approximation before quoting a precision result. Luminosity and detector errors
already inside the likelihood must not be supplied again as acceptance errors.

## Correlations

- Equal nuisance names are shared. Start from `fitting/CONVENTIONS.md` §6 and review physical
  definitions, not only spelling: luminosity, pileup, prefiring and matching background/theory
  sources can be common. `Category` only groups reports; it does not establish correlations.
- For any source named `SigModel`, the adapter sets `NuisanceParameter: SigModel_<channel>`
  while retaining histogram suffixes. The current ττ config does not fit its generator
  comparison (`SigModel_tautau` is reported only); the adapter does not add it.
- Electron, muon and tau efficiencies apply only where the corresponding objects are used.
  Fake and MC statistical parameters remain channel/region-specific.
- Shared envelope names (PDF or QCDScale) are a modelling assumption. Confirm the variation
  construction before using full correlation; partial correlations require an explicit
  nuisance decomposition, not renaming a whole source.

## Validation and limitations

Unit checks exercise quoted labels, normalization closure, all-three-channel generation,
SigModel decorrelation, external response validation, inverse acceptance scaling, correlated
and anticorrelated channel responses, isolated refits, asymmetric error aggregation and
refusal to run incomplete inputs. Controlled refits in tests are simulated, not TRExFitter runs. The local TRExFitter submodule
is empty and the required μμ ROOT input is absent: **TRExFitter runtime/schema validation
and a real fit have not been performed**. Generated files are reviewable preparation artifacts.
The repository's v1.8.0 conventions and helpers are the implementation reference. Framework
background: [TRExFitter tutorial](https://indico.cern.ch/event/1283925/contributions/5394308/attachments/2643277/4574715/TrexFitter.pdf).

## Channel inputs received

- [Muon summary intake](mumu_inputs.md).
- [Updated tau results intake](tautau_inputs.md), with exact metadata in `tautau_inputs.json`.

- [Electron histogram and archive intake](ee_inputs.md), with exact fit metadata in `ee_inputs.json`.
