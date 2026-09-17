# State of z-tautau (17 September 2026)

Everything `REVIEW_v4.md` asked for is done and committed; `REVIEW_v4_RESPONSE.md` is the point-by-point
answer and `output/RESULTS.md` the numbers. Nothing is pending in the analysis chain.

**Not pushed.** The commits are local on `main` (this batch session has no GitHub credentials: `git push`
asks for a username). Push from an interactive shell.

## Two things left for other people

1. **The combination** (`combination/combLieke/`) must change three entries of `config/channels.json`
   before it can read our inputs: `signal` (the signal is 15 `DYtautau_tDM*` templates, listed in the
   metadata as `signal_samples`), `acceptance` (must become `null` -- A x epsilon is inside our fit now and
   adding `Acc_*` on top double counts), and `drop_empty_bins`. `docs/11-combination-inputs.md` section 6
   says exactly how, and section 7 what has to be understood before our number is used.
2. **The slide deck** in `slides/` was built against the v3 result schema and is marked superseded in
   `build_deck.py` and `make_figures.py`. It has to be rebuilt before it is shown again.

## How to re-run

```bash
source ../setup.sh
python run_tautau_base.py --from 3        # tau_h tau_h fake factors, BDT, fit/fitinputs/tautau_base.root
python run_all.py --from 4                # templates (~45 min), then the fits and the report
nohup setsid bash condor/orchestrator.sh > $BND_TAUTAU_CACHE/logs/condor_chain.log 2>&1 &   # step 5+6 on
                                          # HTCondor instead: 12 fits in parallel, the ranking one job per
                                          # parameter (26 min + 7 min, against ~4 h serially)
```
