# What the combination gets from this channel

See `fitting/CONVENTIONS.md` (the contract) and `handoff.md` (the numbers).

- `fit/fitinputs/zmumu.root` (committed) -- every template with Sumw2, names `mumu_SR__<sample>[__<NP>Up|Down]`,
  `mumu_CRemu__...`; `meta_json` inside the file (and `zmumu.root.meta.json` next to it) with
  L, C, sigma_fid^pred, A_60_120, A_m50, the reconstruction SF, the counting numbers and the `acceptance` block
  (`zmumu/acceptance.py`: 60-120 GeV at the top level, m > 50 GeV under `m50`), which the MultiFit turns
  into `Acc_PDF`, `Acc_AlphaS`, `Acc_QCDScale`, `Acc_PS_FSR`, `Acc_PTZ_mumu`, `Acc_Generator_mumu`,
  `Acc_QEDFSR_mumu`, `AccStat_mumu` (`combination/combLieke/config/channels.json`).
- `fit/zmumu.config` -- the TRExFitter configuration (POI `mu_Z`, NP names and categories).
- `fit/results/zmumu/RooStats/zmumu_combined_zmumu_model.root` (committed) -- the channel's workspace (the
  MultiFit in `combination/combLieke/` rebuilds its own from the config and inputs).
- `fit/results/zmumu_fit_result.json` -- mu_Z with stat/syst/lumi, grouped impacts, the three
  cross sections, n_obs, post-fit background.

Settled for the three channels (`fitting/CONVENTIONS.md` §6, `handoff.md` "Cross-channel checks"):
the acceptance denominator is 60 < m_LHE < 120 GeV (Born level, NLO); the fiducial volume of this
channel uses dressed leptons (dR < 0.1; Born differs by 1.8%). Exact reference numbers:
sigma_fid^pred = 799.566 pb, A_60_120 = 0.409209, sigma^pred(60-120) = 1953.93 pb, A_m50 = 0.394703.
Orthogonality: the tau tau channel vetoes muons; the ee selection shares only ZZ -> 4l events with
the mumu SR (61 events, 0.0006 % of the SR).
