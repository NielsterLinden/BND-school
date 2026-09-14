# Final-state radiation (FSR) photons

## The problem

A muon is a charged particle, and charged particles radiate. In roughly 15% of
`Z -> mu mu` decays at least one muon emits a photon hard enough to matter, and
in a few percent the photon carries away more than a GeV. The detector measures
the muon *after* it radiated, so the reconstructed pair carries less energy than
the Z did:

```
    Z  ->  mu+  mu-  gamma            m(mu mu)  <  m(Z)
```

This does two things, and it is worth keeping them separate because they matter
to different parts of the measurement:

1. **It distorts the lineshape.** Events migrate from the peak down into a
   low-mass tail. The peak is no longer a symmetric Breit-Wigner; it acquires a
   shoulder on the left. This is why the Voigt-profile fit in step 4 has a poor
   chi2 -- a Voigt has no tail, so the fitted background bends up to absorb it.

2. **It moves events across the window boundary.** Events whose true mass is
   inside 60-120 GeV can radiate down below 60 GeV and leave the fiducial
   volume. This is the part that affects the *cross section*, and it is much
   smaller than the lineshape distortion, because 60 GeV is a long way below
   the peak.

Measured here: FSR recovery changes the integrated 60-120 GeV yield by
**0.007%**, while visibly reshaping the peak. If you take one thing from this
page: for a *wide* mass window FSR is a lineshape issue, not a yield issue. It
would matter far more with a narrow window such as 80-100 GeV.

## What NanoAOD gives you

NanoAOD ships a dedicated `FsrPhoton` collection. The association work is
already done: `FsrPhoton_muonIdx` is the index of the muon the photon was
matched to during reconstruction. You do **not** have to do the dR matching
yourself, and you should not -- the reconstruction used tracker-level
information that is not in NanoAOD.

| Branch | Meaning |
|---|---|
| `FsrPhoton_pt`, `_eta`, `_phi` | photon kinematics (massless) |
| `FsrPhoton_muonIdx` | index into the `Muon` collection |
| `FsrPhoton_relIso03` | photon isolation, rejects photons from jets |
| `FsrPhoton_dROverEt2` | `dR(mu, gamma) / pT(gamma)^2`, the collinearity measure |

`dROverEt2` is the key discriminant. Genuine FSR is *collinear* with the muon
and *soft*, so real FSR photons sit at small `dR` and/or large `pT`, giving a
small ratio. Photons from pileup or from a nearby jet are uncorrelated with the
muon direction and land at large values.

## What this analysis does

In `zmumu/objects.py`:

```python
def fsr_photon_mask(events):
    return (
        (events.FsrPhoton_pt > 2.0)              # FSR_PT_MIN
        & (events.FsrPhoton_relIso03 < 1.8)      # FSR_REL_ISO_MAX
        & (events.FsrPhoton_dROverEt2 < 0.012)   # FSR_DR_OVER_ET2_MAX
    )
```

These are the standard CMS `H -> 4l` FSR recovery cuts, which is where the
algorithm and its tuning come from. Then, in `dimuon_mass(..., with_fsr=True)`,
every surviving photon whose `muonIdx` points at one of the two signal muons has
its four-vector added to the pair:

```
    m(mu mu gamma...)  =  | p(mu1) + p(mu2) + sum_i p(gamma_i) |
```

All cuts live in `config.py` (`FSR_*`), so they can be varied without touching
the analysis code.

## Results

- **4.6%** of selected events have at least one recovered FSR photon
  (`output/plots/step1_kinematics.png`, bottom-right panel).
- The recovered mass shift is shown in
  `output/plots/step1_fsr_mass_shift.png` -- a steeply falling spectrum, mostly
  below 2 GeV, with a tail out to ~15 GeV.
- The before/after comparison of the peak is
  `output/plots/step1_mass_fsr.png`. The peak sharpens and the low-mass
  shoulder shrinks.

The measurement uses FSR-recovered masses everywhere. The difference from the
bare-muon mass is carried as the "FSR recovery" systematic in step 5.

## What this does *not* do

Recovering the photon corrects back to the **dressed** muon -- the muon plus its
collinear radiation. It does **not** correct back to the **Born** level (the
muon before any radiation). Those differ by the non-collinear radiation the
algorithm cannot catch, and the difference can only be computed with a
generator.

This matters when comparing to a theory prediction: make sure the prediction is
for dressed leptons, or the comparison is inconsistent at the ~1% level. The
fiducial definition in `docs/06-cross-section.md` is stated at dressed level for
exactly this reason.

## If you want to go further

- Vary the three FSR cuts and watch the yield move; that is a more honest
  systematic than the bare-vs-recovered difference used now.
- Tighten the mass window to 80-100 GeV and repeat. The FSR systematic should
  grow by roughly an order of magnitude, which is a good check that the
  machinery is doing what this page claims.
