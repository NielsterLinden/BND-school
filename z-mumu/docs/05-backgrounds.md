# Step 3 -- background estimation

`scripts/step3_backgrounds.py`

Every background here is estimated **from the data**. No simulation is used,
because none is available locally (`/dcache/atlas/sjankovy/BND/mc/` is empty).
That constraint turns out to be workable, because the two background classes
that matter can each be measured in a control region.

## The two classes, and why you need both

| Class | Charge | Example | Caught by |
|---|---|---|---|
| non-prompt / fake | uncorrelated | b -> mu inside a jet, punch-through | same-sign region |
| prompt flavour-symmetric | opposite | `ttbar`, `WW`, `tW`, `Z -> tau tau` | e-mu region |

These are disjoint and complementary. Using only one of them is the classic way
to underestimate the background in this channel: the same-sign method is blind
to `ttbar`, which is the larger of the two.

## 1. Non-prompt muons -- same-sign

Covered in detail in [03-fake-leptons.md](03-fake-leptons.md). In short:

```
    N_fake(OS)  =  R_OS/SS * N(SS),     R_OS/SS = 1.0 +/- 50%
```

## 2. Flavour-symmetric backgrounds -- e-mu

`ttbar`, `tW`, `WW` and `Z -> tau tau` all decay to lepton pairs **without
preferring a flavour**. Each W (or tau) decays to `e` or `mu` independently, so
the underlying rates are in fixed combinatorial proportion:

```
    N(ee) : N(mu mu) : N(e mu)   =   1 : 1 : 2
```

The factor of two on `e mu` is just combinatorics -- `e+mu-` and `mu+e-` are
both `e mu`, while `mu+mu-` is one way. Folding in the per-flavour selection
efficiencies:

```
    N(e mu)   =  2 * N0 * eff_e * eff_mu
    N(mu mu)  =      N0 * eff_mu^2

        =>    N(mu mu)  =  0.5 * (eff_mu / eff_e) * N(e mu)
```

So counting `e mu` events under an otherwise identical selection measures the
`mu mu` contamination, **including the `Z -> tau tau` component**, which is a
genuine irreducible background to `Z -> mu mu` and would be easy to forget.

Implemented in `_fill_emu` in step 1: exactly one selected muon above the
trigger plateau, exactly one selected electron (`cutBased >= 3`,
`pfRelIso03_all < 0.15`, `pT > 20`, `|eta| < 2.5`), opposite sign.

### The assumption, stated plainly

`k = eff_mu / eff_e` is taken as **1.0 +/- 50%** (`config.FS_EFF_RATIO`,
`FS_REL_UNC`). It cannot be measured in this dataset: fixing it properly needs
`Z -> ee`, and the SingleMuon primary dataset does not contain an unbiased
electron sample.

This assumption is wrong at the tens-of-percent level -- electron and muon
efficiencies genuinely differ. It is tolerable only because the background is
0.12% of the signal, so a 50% error on it is 0.06% on the cross section. **If
you tighten the selection or move off the Z peak, fix this first.** The clean
way is to take `k` from the `SingleElectron` dataset (recids 30529/30562, which
the z-ee subgroup is already using) via the ratio of Z yields.

Two further approximations, both folded into the same 50%:

- The e-mu region requires the *muon* to fire the trigger and pass `pT > 26`,
  while the mu-mu region lets *either* muon fire. The trigger acceptance is
  therefore not identical between the regions.
- The e-mu region is not itself fake-free (`W+jets` with a fake electron), so
  there is a small double-count with the same-sign estimate.

## What is not estimated

`WZ` and `ZZ` where both muons are prompt and same-flavour. These are neither
charge-uncorrelated nor flavour-symmetric, so neither control region sees them.
They are ~0.1% of the Z peak in this selection and are absorbed into the
background systematic rather than subtracted. With MC available they would be
taken from simulation, which is the normal treatment.

## Result

Measured over the full dataset, in `60 < m < 120` GeV:

| Component | Yield | Fraction of observed |
|---|---:|---:|
| non-prompt (same-sign) | see `RESULTS.md` | ~0.03% |
| flavour-symmetric (e-mu) | see `RESULTS.md` | ~0.12% |
| **total** | | **~0.15%** |

A 0.15% background is very small, and that is the expected outcome of a tight
isolated-muon selection on the Z peak. The value of doing this properly is not
that the number is large -- it is that the number is now *shown* to be small,
with a control region behind each piece.

## Outputs

- `output/plots/step3_control_regions.png` -- same-sign and e-mu mass spectra
- `output/plots/step3_data_vs_background.png` -- data with the stacked estimates
  and a ratio panel
- `output/plots/step3_subtracted.png` -- the Z peak before and after subtraction

Check the same-sign spectrum for a Z peak. If one is visible, charge
misassignment is leaking signal into the control region and `R_OS/SS` is no
longer the right transfer factor.
