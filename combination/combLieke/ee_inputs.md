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

# Electron input from histograms and archived fit

Sources: `datasets/z-ee/`, `z-ee/fit.config`, `z-ee/Zee_fit.tar.gz`.
`ee_inputs.json` records the archive digest, fit result, uncertainty decomposition and sample model.
The runner reads the archive on preparation and writes `generated/ee_archive.json`, also
referenced in `status.json`. It does not extract archive files or execute anything in them.

## Available inputs

Six nominal ROOT files: Data, DY_ee, DY_ττ, Diboson, ttbar, W+jets. The source config
reads `h_mass` and rebins by 2. The combination uses these original templates, changes SR to
ee_SR and mu_signal to mu_comb, and marks DY_ee as signal. No ROOT contents are rewritten.
Binary histogram bin contents/Sumw2 have not been independently inspected in this environment;
the archived HistFactory XML confirms MC-stat errors were enabled for all five MC samples.

The archive includes processed histograms, pre/postfit plots and yields, the RooStats
workspace, HistFactory XML and fit/error-decomposition text files. It supplies a reference fit:

| Quantity | Value |
|---|---:|
| mu_signal | 0.912567 ± 0.000662336 |
| Data-statistical error on mu_signal | 0.000411409 |
| MC-statistical error on mu_signal | 0.000519068 |
| Other systematic error | 0 |

There are 60 per-bin gamma parameters in the archived text output. MC statistics are already
part of the quoted total; do not add them twice. `FullSyst` in the grouped file contains only
gammas here and does not indicate a complete experimental systematic model.
The XML has Lumi=1 and LumiRelErr=0: templates are already event yields, with no luminosity
uncertainty in this likelihood. These are not the physical integrated luminosity values.

## Use in MultiFit

The original workspace has POI mu_signal, region SR, and normalization limits [-100,100].
It cannot be directly substituted into our common mu_comb workspace without adaptation.
The runner rebuilds ee from the original histograms/config, preserving MC statistics and
using the common parameter range; the archive is an independent reference check.
Archive postfit templates must not become nominal inputs (that would reuse the data).

The archive does not contain an agreed per-flavour truth 60–120 GeV reference cross section,
acceptance, or electron correction/systematic inputs. Therefore `sigma_reference_pb` remains
null and `measurement_ready` remains false. In particular mu_signal times 6077.22/3 relates
to the M>50 reference, not automatically the required truth 60–120 window, even though the
reconstructed fit mass is 60–120 GeV. No normalization is guessed from the fitted yield.

Still needed: that truth reference and acceptance, confirmation/correction of luminosity and
nominal electron/pileup efficiencies, detector/background/theory systematics, signed external
acceptance responses and cross-channel correlation/overlap decisions. The raw statistical
model is available; it is not yet a complete cross-section measurement for the final combination.
