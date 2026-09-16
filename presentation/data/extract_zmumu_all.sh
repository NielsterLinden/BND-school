#!/usr/bin/env bash
# Freeze every Z -> mumu number the presentation prints: runs the six extractors in dependency order.
#
#   cd /project/atlas/Users/nterlind/BND-school && source setup.sh && bash presentation/data/extract_zmumu_all.sh [--check-only]
#
# LCG env only (uproot/awkward; never ROOT, never pip). Reads z-mumu/ read-only; part files of the MC pass go to
# presentation/work/extract/sr_stack (git-ignored, resumable: delete a part to redo it). Every extractor prints its
# assertion list against the channel anchors and exits non-zero on any mismatch (numbers are never adjusted).
# `--check-only` re-runs the assertions on the committed JSON files without writing anything.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"
if ! python -c "import uproot, awkward" 2>/dev/null; then
    echo "extract_zmumu_all.sh: run 'source setup.sh' first (LCG env with uproot/awkward)" >&2
    exit 2
fi
MODE="${1:-}"
WORKERS="${WORKERS:-8}"
t0=$(date +%s)
python "$HERE/extract_zmumu_fit.py" $MODE
python "$HERE/extract_zmumu_fakes.py" $MODE
if [ "$MODE" = "--check-only" ]; then
    python "$HERE/extract_zmumu_sr_stack.py" --check-only
else
    python "$HERE/extract_zmumu_sr_stack.py" --workers "$WORKERS"          # ~2-4 min on local skims
fi
python "$HERE/extract_zmumu_corrections.py" $MODE                          # needs zmumu_sr_stack.json
python "$HERE/extract_zmumu_events.py" --files 6 --seed 20260915 $MODE
python "$HERE/extract_zmumu_tnp.py" $MODE                                  # tag-and-probe cell, maps, two T&P events
echo "[extract_zmumu_all] done in $(( $(date +%s) - t0 )) s: $(ls -1 "$HERE"/zmumu_*.json | xargs -n1 basename | tr '\n' ' ')"
