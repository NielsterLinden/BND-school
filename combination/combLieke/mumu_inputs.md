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

# Muon input received from the summary deck

Sources: `z-mumu/summary_deck/zmumu_summary.pdf`, final numbered slide 19,
review-fixes slide 17; `z-mumu/REVIEW.md` §0 confirms the updated status.
The actual PDF filename is `zmumu_summary.pdf` (without the extra hyphen).

## Recorded values and decisions

- Fiducial measurement: **790.1 ± 0.2 (stat) ± 8.5 (syst) ± 9.6 (lumi) pb**.
- Fiducial prediction: **799.6 pb**, rounded on the slide.
- Signal strength: **0.988 ± 0.016**, rounded on the slide.
- Recommended acceptance: **A(60–120) = 0.4092**, NLO.
- Fit goodness of fit: **p = 0.79**; binning stability approximately **±0.2%**.
- Config: `z-mumu/fit/zmumu.config`; templates: `z-mumu/fit/fitinputs/zmumu.root`.
- Keep the exact reference **1953.9313466473218 pb** already obtained from
  `sigma_fid_pred_pb / A_60_120` in `z-mumu/fit/results/zmumu_fit_result.json`.
  Do not replace it with a ratio of rounded slide numbers.

The current shape fit is now marked measurement-ready, based on the group's updated review.
This corrects our earlier reading of the historical review recommendations. The matching ROOT
file is still missing from this checkout, so the nominal-fit input check continues to block.
The source config's `UseMinos: all` is retained in adapted configs; this is specifically part
of the documented fix for the earlier covariance issue.

No extra 0.7% lineshape uncertainty is added: §0 says the revised SigModel already covers it.
The stability and counting differences are validation checks, not extra independent errors.
The quoted fiducial uncertainties are cross-checks, not inputs replacing the binned likelihood.

## Still needed

- The ROOT templates matching the updated config and normalization metadata.
- Signed acceptance responses by source (PDF, scale, etc.) and cross-channel correlations
  for the external refits. The last slide supplies nominal A, not these responses.
- Agreement with ee and tautau on the common truth-window definition and event orthogonality.
  No acceptance uncertainty or correlation has been inferred from the rounded slide numbers.
