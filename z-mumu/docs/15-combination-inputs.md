# What the combination gets from this channel

See `fitting/CONVENTIONS.md` (the contract) and `handoff.md` (the numbers).

- `fit/fitinputs/zmumu.root` -- every template with Sumw2, names `mumu_SR__<sample>[__<NP>Up|Down]`,
  `mumu_CRemu__...`; `meta_json` inside the file (and `zmumu.root.meta.json` next to it) with
  L, C, sigma_fid^pred, A_60_120, A_m50, the acceptance uncertainties and the counting numbers.
- `fit/zmumu.config` -- the TRExFitter configuration (POI `mu_Z`, NP names and categories).
- `fit/results/zmumu/RooStats/zmumu_combined_zmumu_model.root` -- the workspace the MultiFit
  combines (`fitting/combination_skeleton.config`).
- `fit/results/zmumu_fit_result.json` -- mu_Z with stat/syst/lumi, grouped impacts, the three
  cross sections, n_obs, post-fit background.

Open decisions for the three channels: the acceptance denominator (60 < m < 120 GeV,
recommended, vs m > 50 GeV) and NLO vs LO for A (3.2% apart). Dressed and Born fiducial
definitions differ by 1.8%; this channel quotes dressed leptons (dR < 0.1).
