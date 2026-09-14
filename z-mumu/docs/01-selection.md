# Step 1 -- event selection

`scripts/step1_selection.py`

## What it does

Reads the full 31.5 GB skim **once** and fills three regions in the same pass:

| Region | Definition | Used for |
|---|---|---|
| **OS** | two opposite-sign selected muons | signal |
| **SS** | two same-sign selected muons | fake-muon background ([03](03-fake-leptons.md)) |
| **e-mu** | one selected muon + one selected electron, opposite sign | flavour-symmetric background ([05](05-backgrounds.md)) |

One pass matters: each pass over the skim costs a few minutes, and filling all
three together guarantees they share every cut except the pairing, which is what
makes the background subtraction valid.

## The cuts, in order

| Cut | Why |
|---|---|
| golden JSON | Only certified lumisections are in the `brilcalc` luminosity. The numerator and denominator of a cross section must cover the same data. Non-negotiable. |
| `HLT_IsoMu24 \|\| HLT_IsoTkMu24` | The analysis trigger. `IsoTkMu24` is the tracker-muon variant; in 2016 both are needed to avoid a few-percent hole. |
| MET filters + `PV_npvsGood >= 1` | Removes beam halo, HCAL noise and events without a good primary vertex. |
| exactly 2 selected muons | See below. |
| leading muon `pT > 26` GeV | The IsoMu24 plateau. Sitting at 24 would put the measurement on the trigger turn-on, where the efficiency is steep and badly known. |
| opposite sign | Z decays to opposite-sign pairs. |
| `60 < m(mu mu) < 120` GeV | The fiducial mass window. |

Muon object selection (`zmumu/objects.good_muon_mask`):

```
Muon_mediumId  and  pT > 20 GeV  and  |eta| < 2.4
and  pfRelIso04_all < 0.15
and  |dxy| < 0.2 cm  and  |dz| < 0.5 cm
```

## Two choices worth defending

**Exactly two muons, not at least two.** "At least two" leaves the fiducial
volume ambiguous -- which pair is the Z? -- and lets `WZ` and `ZZ` in. "Exactly
two" is unambiguous and suppresses them. The cost is a small inefficiency when a
real Z event picks up a third muon, which is rare enough to ignore at this
precision. If you change this, the fiducial definition in
[06-cross-section.md](06-cross-section.md) changes with it.

**Isolation at 0.15, not 0.25.** Tighter isolation costs ~5% efficiency (which
is measured, so it costs nothing in accuracy) and buys a factor of a few in fake
rejection. With the background at 0.15% either way, this is not a critical
choice -- but the tighter cut makes the same-sign control region cleaner.

## The trigger-menu wrinkle

The 2016 HLT menu changed during data-taking. Of the ~328 trigger branches in
the skim, **39 are absent from some files** -- 24 exist only in Run2016G, 15
only in Run2016H, and the boundaries move even within era H.

`zmumu/io.read_chunks` handles this by intersecting the requested branches with
what each file actually contains, and reporting anything missing as all-False
rather than raising. **`HLT_IsoMu24` and `HLT_IsoTkMu24` are present in all 152
files**, so the analysis trigger is unaffected; the reference triggers used in
step 2 were chosen from the universal set for the same reason.

If you add a trigger to `config.TRIGGERS`, check it exists everywhere first:

```python
from zmumu import io
import uproot, collections
counter = collections.Counter()
for f in io.find_files():
    counter.update(uproot.open(f"{f}:Events").keys())
# any path with count < 152 is not universal
```

## Outputs

- `output/data/step1_selection.pkl` -- histograms, cutflow, region yields, and
  the (pT, eta) maps of the signal muons that step 5 needs
- `output/plots/step1_cutflow.png` -- absolute and relative survival per cut
- `output/plots/step1_mass_regions.png` -- the three regions overlaid
- `output/plots/step1_mass_fsr.png` -- FSR recovery, before and after
- `output/plots/step1_mass_wide.png` -- 10-200 GeV, showing the peak in context
- `output/plots/step1_kinematics.png` -- muon pT/eta, Z pT/rapidity, pileup, FSR multiplicity
- `output/plots/step1_fsr_mass_shift.png` -- mass recovered per event

## Performance

12 workers over 152 files, reading ~30 branches. The branches this step needs
are a small fraction of the skim, so it is CPU-bound on the awkward operations
rather than I/O-bound.
