# Handoff – Combination

## Team

- Lieke Gijsen
- Caspar van Eck

## What we do

We take the final numbers from the three channel groups (`z-ee/`, `z-mumu/`, `z-tautau/`) and combine
them into one Z cross section.

Each channel gives us `n_obs`, `n_bkg`, `acc_eff` and its systematic uncertainties. From those we
compute a cross section per channel, build the covariance matrix between the three, and average them
with BLUE.

## How to run

```
jupyter notebook combination.ipynb
```

## Results so far

- Nothing yet: waiting for the numbers from the three channels.
