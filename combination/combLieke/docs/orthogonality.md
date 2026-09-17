# Are the three channels orthogonal?

**Yes, to 2 × 10⁻⁴ or better in every pair, measured event by event on the 2016 G+H data.**
The measurement is `checks/orthogonality.py`. Its output, including the (run, lumi, event) of every
overlapping event, is `checks/orthogonality.json`.

## Why the selections, not the primary datasets, have to be checked

CMS writes an event into every primary dataset whose trigger it fires. `SingleMuon` (μμ),
`SingleElectron` (ee) and `Tau` (ττ) are therefore not disjoint by construction. Two channels share an
event only if that event passes both offline selections, so each pair is tested on the data of one
channel with the other channel's offline selection applied on top.

| | μμ | ee | ττ |
|---|---|---|---|
| trigger | `HLT_IsoMu24 \|\| HLT_IsoTkMu24` | `HLT_Ele27_WPTight_Gsf` | `HLT_Double*IsoPFTau35/40*` |
| leptons | exactly 2 tight muons, p_T > 26/20 GeV, iso < 0.15 | exactly 2 cut-based medium electrons, p_T > 20 GeV | 2 τh, DeepTau Tight, opposite sign |
| e veto | no | – | yes: MVA noIso WP90, p_T > 10, \|d_xy\| < 0.045, \|d_z\| < 0.2, iso < 0.3, conversion veto, ≤ 1 lost hit |
| μ veto | – | no | yes: medium ID, p_T > 10, \|d_xy\| < 0.045, \|d_z\| < 0.2, iso < 0.3 |

The other channel's trigger is not required: the ee trigger bit is not in the μμ skim, and the μ/e
triggers are not in the ττ skim. Dropping a requirement can only add events, so every count below is an
upper bound.

## Result

| pair | sample tested | events tested | overlap | fraction of channel A | fraction of channel B |
|---|---|---:|---:|---:|---:|
| μμ ∩ ee | μμ signal region, SingleMuon skim (reproduces the 10 378 567 fitted events exactly) | 10 378 567 | ≤ **97** | 9 × 10⁻⁶ of μμ | 1.5 × 10⁻⁵ of ee |
| ττ ∩ μμ | Tau skim after the ττ lepton veto | 5 011 196 | **0** | 0 | 0 |
| ττ ∩ ee | the same, looked up in the ττ signal region (OS, both τh Tight) | 21 160 in the SR | ≤ **4** | 1.9 × 10⁻⁴ of ττ | 6 × 10⁻⁷ of ee |

* **μμ ∩ ee.** Neither channel vetoes the other flavour. An event with exactly two tight muons *and*
  exactly two opposite-sign medium electrons in 60–120 GeV is a four-lepton event (ZZ → 2e2μ). There are
  at most 97 such events. That induces a correlation between the two data-statistical uncertainties of
  about √(9·10⁻⁶ × 1.5·10⁻⁵) ≈ 10⁻⁵, and data statistics are 0.5 pb of the 33 pb total.
* **ττ ∩ μμ = ∅, exactly.** None of the 5.0 M events that pass the ττ lepton veto (a superset of the
  fitted ττ events, since the veto is applied when the ττ ntuples are built) has two muons passing the μμ
  signal-muon definition, even without any charge, mass or trigger requirement.
* **ττ ∩ ee is not empty by construction, but negligible.** The ττ electron veto is not a superset of the
  z-ee electron: a cut-based medium electron can fail MVA WP90, the impact-parameter cuts or the lost-hit
  requirement. 7 339 veto-passing events in the Tau skim pass the full z-ee offline selection. Only **4**
  of them are in the ττ signal region: (280242, 48, 85602248), (282037, 725, 1310902895),
  (283059, 315, 588096663) and (282035, 29, 65662226).

The fit treats data statistics as independent between channels, which is exact to this precision.

## Simulation

The signal templates of the three channels come from disjoint generator events: DYee, DYmumu and
DYtautau are split by LHE lepton flavour. Z → ττ appears as a background in the ee (6.7k events) and μμ
(11k) signal regions, while the ττ channel's signal uses the same simulated process. The MC-statistical
correlation this induces is below 10⁻³ of either measurement and is not modelled.

## Reproduce

```bash
source ../../setup.sh
python checks/orthogonality.py --workers 8     # ~3 min: 152 SingleMuon + 100 Tau skim files
```

Inputs: `$BND_SKIM_DIR/data_2016{G,H}` (z-mumu v2 skims) and
`$BND_TAUTAU_CACHE/skims_v1/data_2016{G,H}` plus `ntuples_v1/data_2016{G,H}.root` (z-tautau). The
selections are imported from the channels' own code (`z-mumu/zmumu/regions.py`,
`z-tautau/ztautau/objects.extra_lepton_veto`, `z-tautau/ztautau/analysis.regions`). The ee selection is
copied from `z-ee/z-ee.ipynb` (`process_sample`), which has no importable module.
