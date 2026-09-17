# Resume here (17 Sep 2026, node reboot)

Committed locally as `b860653` on `main`, **not pushed** (no GitHub credentials in the batch session:
`git push` asks for a username). Push from an interactive shell.

## Done
Repository restructured (four-channel measurement is the only one, canonical `fit/` + `output/`, job
`ztautau`), every code fix of `REVIEW_v4.md` applied, templates rebuilt
(`fit/fitinputs/ztautau.root`, 17 Sep 10:51), per-channel exports and **all TRExFitter configs written**.
The combination can already build workspaces from those.

## Not done — run this next, in this order
```bash
cd z-tautau && source ../setup.sh
python run_all.py --from 5 --to 6      # all fits (~4 h: the ranking of the combined fit dominates), report,
                                       # and scripts/update_docs.py injects the numbers into the docs
# then, to regenerate the tau_h tau_h base plots that were dropped with the v3 output/ directory
python scripts/step3_fakefactors.py && python scripts/step3b_bdt.py --no-train && python scripts/step4a_tautau_base.py
# and check that fit/fitinputs/tautau_base.root still gives the yields of the committed one before committing it
```
Then commit `output/`, `fit/results/*_fit_result.json`, the workspaces
(`git add -f fit/results/*/RooStats/*_combined_*_model.root`) and the updated docs.

## Two things to fix in the documentation
1. `docs/11-combination-inputs.md` §6 was written against the **BLUE / covariance** combination
   (`combination/comb/inputs.py`, `ChannelResult`, `combination/docs/01-inputs.md`). Commit `b6aeeb2`
   removed that combination: the live one is `combination/combLieke/`, a TRExFitter MultiFit. Route (a)
   of §6 is therefore the route; route (b) and the `for_combination` block of `output/results.json`
   should be kept (they cost nothing and the numbers are right) but presented as the secondary route,
   and the references to `combination/docs/` and to `load_tautau` updated to `combination/combLieke/`.
2. `slides/` is marked superseded (it was built against the v3 result schema) and has to be rebuilt.
