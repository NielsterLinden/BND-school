# Updated tau input

Sources: `z-tautau/output/RESULTS.md`, exact values in `output/results.json`,
`fit/ztautau.config` and `fit/fitinputs/ztautau.root.meta.json`.
Machine-readable intake: `tautau_inputs.json` (reference information; not an executable variation model).

## Nominal measurement

- Genuine-tau subtraction in fake factors (`mcsub`) is now nominal; `nosub` is a cross-check.
- Three BDT categories: `tautau_SR0`, `tautau_SR1`, `tautau_SR2`.
- Signal strength: 1.205 +0.145 -0.125.
- Fiducial cross section: 5.421 ± 0.091 (stat) ± 0.601 (syst) pb.
- Inclusive 60–120 GeV cross section: 2343 ± 39 (stat) ± 260 (syst) ± 86 (acceptance) pb.
- Goodness of fit: p = 0.21943.
- Reference cross section: 1944.8826556737133 pb (unchanged).
- Acceptance: 0.0023140314215664886; C = 0.06542356849711156.
- Luminosity: 16393.381 pb^-1, with 1.2% uncertainty already in the fit.

The signal POI multiplies only fiducial `DYtautau`. `DYtautau_nonfid` is now a separate
background and remains so in the adapter, including during external acceptance refits.
This follows the channel's fiducial extraction/extrapolation prescription; it is not a fully
universal model tying all non-fiducial DY production to the common cross section.
The original category `DropBins` settings and all nuisance scopes are retained.
Prefit totals in the report are not necessarily totals of fitted bins after DropBins.

## Acceptance magnitudes received

| Source | Relative magnitude (rounded) |
|---|---:|
| QCD scales | 3.3771% |
| PDF | 0.4877% |
| alpha_s | 0.5864% |
| ISR | 0.9320% |
| FSR | 0.2467% |
| MC statistics | 0.7121% |

Exact values are copied in `tautau_inputs.json`. These are magnitudes, not signed up/down
responses. We do not infer signs, symmetrize them, or set correlations with ee/mumu silently.
They therefore have not been inserted in `acceptance.sources`. The approximately 86 pb
acceptance error is a channel cross-check; do not add it directly to the combination error,
or add it again on top of its individual components. Grouped impacts also overlap and must
not be added as independent extra errors.

## Outstanding

Both the metadata sidecar and `fit/fitinputs/ztautau.root` are now present locally.
The runner passes the tau ROOT-file existence check; histogram contents have not yet been
validated through TRExFitter. For external refits,
supply signed acceptance responses and agree correlations (or explicitly adopt a documented
symmetric-magnitude approximation). Common truth-definition and event-overlap checks remain.

## Recheck at repository commit 7a4e540

The recorded nominal fit, regions, prediction and combination metadata are unchanged.
Commit aa5be6b commits the ROOT templates, conventional-name RooStats workspace and
`fit/results/ztautau_fit_result.json`. We still rebuild workspaces for our common POI and
external variations; the channel workspace is a reference, not a directly substituted input.
The runner now supplies the conventional-name workspace alias when TRExFitter produces the
`allBinsFitRegions` filename for DropBins, in each isolated run directory.

The handoff says `SigModel_tautau` is reported, not fitted; the current config has no such
Systematic block. No tau generator nuisance is added by the adapter. The handoff recommends
correlated acceptance, but does not supply signed source responses across the three channels.
That recommendation is recorded; it does not resolve the source-level correlation model.

Acceptance bookkeeping: `fit/results/ztautau_fit_result.json` reports
`acceptance_rel_unc = 0.03593866580978255`, which equals quadrature of the theory components
without acceptance MC statistics. Including `A_mc_stat` gives approximately 0.03663730
(the 3.7% and 86 pb quoted in RESULTS.md). Use one consistent component model for external
refits, with MC acceptance statistics handled explicitly; do not mix these two totals.
