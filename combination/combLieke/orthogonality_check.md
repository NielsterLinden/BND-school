# Quick channel orthogonality check

Quick source-code and handoff review, not an event-level data overlap measurement. No zero-overlap claim is made.

| Pair | Assessment |
|---|---|
| mumu / ee | Overlap allowed by selections |
| mumu / tautau | Veto suppresses overlap; strict exclusivity not proven |
| ee / tautau | Veto suppresses overlap; strict exclusivity not proven |

## mumu / ee

Exactly two tight muons versus exactly two medium electrons; neither SR vetoes the other flavour. Both may accept a 2e2mu event.
mumu handoff reports 61 ZZ_4L events satisfying ee-like cuts in mumu SR (0.0006% of 10.38M SR events); not independently reproduced and not a full data-ID intersection.

## mumu / tautau

tau veto uses mediumId, pt>10, |eta|<2.4, |dxy|<0.045, |dz|<0.2, iso<0.3. mumu uses tightId, pt>20, |dxy|<0.2, |dz|<0.5, iso<0.15. Different ID branches and impact-parameter cuts prevent a simple subset proof.


## ee / tautau

ee uses cutBased>=3, pt>20, |eta|<2.5. tau veto instead requires mvaFall17V2noIso_WP90, pt>10, impact-parameter/isolation cuts, convVeto and lostHits<=1. An ee-selected electron is not guaranteed by these expressions to satisfy the tau veto.


The μμ eμ region is marked VALIDATION, so the SR is the relevant μμ fitted selection. The ττ BDT categories partition one selected sample. Shared control-sample and MC-statistical correlations are separate checks.

No selected-event ID intersection was performed. The ee notebook reads event IDs but returns mass/weight arrays without retaining those IDs; the histogram delivery cannot establish overlap counts. Raw event caches exist, but regenerating the selected ID lists is beyond this quick review.

Next step: Compare unique (run, luminosityBlock, event) IDs from the actual fitted selections, with identical certification/run coverage. Include tau category assignment and fitted-bin cuts. Review shared data-driven control samples separately.

Sources:
- z-mumu/zmumu/regions.py: muon_masks and dimuon_regions
- z-mumu/zmumu/config.py: MU_DXY_MAX and MU_DZ_MAX
- z-ee/z-ee.ipynb: process_sample
- z-tautau/ztautau/objects.py: extra_lepton_veto
- z-tautau/ztautau/config.py: VETO_MU and VETO_EL
- z-tautau/scripts/step2_ntuples.py: veto applied at ntuple construction
- z-mumu/handoff.md: Event orthogonality
