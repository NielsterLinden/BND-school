#!/usr/bin/env bash
# Freeze every Z -> tautau number the presentation prints: runs the nine extractors in order.
#
#   cd <BND-school> && source setup.sh && bash presentation/data/extract_ztautau_all.sh [--check-only]
#
# LCG env only (uproot / numpy; never ROOT, never pip). Reads z-tautau/ read-only (committed outputs + the git-ignored
# output/data/*.json, fit/results/ztautau/Plots/*.yaml and the ntuples of /data/atlas/users/sjankovy/BND-school-cache/ztautau/ntuples_v1).
# Every extractor prints its assertion list against the channel anchors (handoff.md, output/RESULTS.md, docs/) and exits non-zero
# on any hard mismatch (numbers are never adjusted); soft checks (numbers quoted in docs from an earlier pass) are reported only.
# `--check-only` re-runs the assertions on the committed JSON files without writing anything.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"
if ! python -c "import uproot, awkward" 2>/dev/null; then
    echo "extract_ztautau_all.sh: run 'source setup.sh' first (LCG env with uproot/awkward)" >&2
    exit 2
fi
MODE="${1:-}"
SEED="${SEED:-20260916}"
t0=$(date +%s)
python "$HERE/extract_ztautau_fit.py" $MODE                                  # fit result + 3 x 14-bin plots (seconds)
python "$HERE/extract_ztautau_sr_stack.py" $MODE                             # fit-input histograms (seconds)
python "$HERE/extract_ztautau_bdt.py" $MODE                                  # BDT score, categories (seconds)
python "$HERE/extract_ztautau_bdt_inputs.py" $MODE                           # the 16 BDT input shapes, signal vs fakes (~1 min)
python "$HERE/extract_ztautau_fakes.py" $MODE                                # FF table + closure 'after' from the ntuples (~1-2 min)
python "$HERE/extract_ztautau_corrections.py" $MODE                          # TauPOG SFs + mean weights on the fiducial signal (~30 s)
python "$HERE/extract_ztautau_mass.py" $MODE                                 # m_vis / m_tt / m_col shapes of the signal (~30 s)
python "$HERE/extract_ztautau_events.py" --seed "$SEED" $MODE                # one SR event + posterior, ten SS pairs (~1 min)
python "$HERE/extract_ztautau_reference.py" $MODE                            # published CMS / ATLAS sigma, this work, prediction (seconds)
echo "[extract_ztautau_all] done in $(( $(date +%s) - t0 )) s: $(ls -1 "$HERE"/ztautau_*.json | xargs -n1 basename | tr '\n' ' ')"
