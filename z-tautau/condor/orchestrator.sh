#!/bin/bash
# Full step 5 + 6 of the Z -> tautau measurement on HTCondor.
#
#   nohup setsid bash condor/orchestrator.sh > $BND_TAUTAU_CACHE/logs/condor_chain.log 2>&1 &
#
# Stage 1  the TRExFitter jobs of step 5, one Condor job each (batch ztt_fit): first the cross-checks, among them
#          ztautau_flatsf and ztautau_ptsplit, then the fits that read those two for the size of TauIDpT_tautau
#          (the measurement and ztautau_emutrig2x; scripts/step5_fit.py pt_model)
# Stage 2  the nuisance-parameter ranking of the combined fit, one Condor job per parameter (ztt_rank)
# Stage 3  merge the ranking, re-read the result, MultiFit of the four channel workspaces
# Stage 4  the report and the result block injected into the documents
# NB: no `set -u` before the LCG view is sourced -- its setup.sh reads unset variables (nikhef-condor skill)
CH=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)      # this checkout: the jobs run where the orchestrator is
D=$CH/condor
mkdir -p "$D"/{out,err,logs}
cd "$CH" || exit 1
source ../fitting/setup.sh
set -o pipefail

wait_for_batch() {           # $1 = JobBatchName, $2 = stage label
    local prefix="$1" stage="$2" n held
    sleep 10
    while :; do
        n=$(condor_q -af JobBatchName JobStatus 2>/dev/null | grep -c "^${prefix} " || true)
        [[ "$n" -eq 0 ]] && { echo "[$stage] batch $prefix drained"; return 0; }
        held=$(condor_q -af JobBatchName JobStatus 2>/dev/null | awk -v p="$prefix" '$1==p && $2==5' | wc -l)
        if [[ "$held" -gt 0 ]]; then
            echo "[$stage] $held HELD job(s) in $prefix -- stopping"
            condor_q -held -af JobBatchName HoldReason | grep "^${prefix}" | head
            return 2
        fi
        echo "[$stage] $n job(s) queued/running in $prefix ($(date +%H:%M:%S))"
        sleep 60
    done
}

echo "######## stage 1: the cross-check fits ($(date))"
for j in ztautau ztautau_flatsf ztautau_ptsplit ztautau_emutrig2x; do      # never read a result of an earlier chain
    rm -f "$CH/fit/results/${j}_fit_result.json"
done
condor_submit "CH=$CH" "params=$D/params_step5.txt" "$D/step5.sub" || exit 1
wait_for_batch ztt_fit stage1 || exit 2
for j in ztautau_flatsf ztautau_ptsplit; do
    [[ -f "$CH/fit/results/${j}_fit_result.json" ]] || { echo "MISSING result for $j -- see condor/err"; exit 3; }
done
echo "######## stage 1b: the measurement and ztautau_emutrig2x, with TauIDpT_tautau ($(date))"
condor_submit "CH=$CH" "params=$D/params_step5_measurement.txt" "$D/step5.sub" || exit 1
wait_for_batch ztt_fit stage1b || exit 2
for j in ztautau ztautau_flatsf ztautau_tautau ztautau_mutau ztautau_etau ztautau_emu ztautau_taulep ztautau_ptsplit ztautau_emutrig2x; do
    [[ -f "$CH/fit/results/${j}_fit_result.json" ]] || { echo "MISSING result for $j -- see condor/err"; exit 3; }
done
echo "stage 1 results:"
grep -h '"poi_value"' "$CH"/fit/results/*_fit_result.json | head -20

echo "######## stage 2: the ranking, one Condor job per parameter ($(date))"
python condor/make_rank_params.py || exit 4
# the ranking jobs are single-core: one CPU per parameter beats six CPUs per fit here
sed 's/^  NumCPU: .*/  NumCPU: 1/' "$CH/fit/ztautau.config" > "$CH/fit/ztautau_rank.config"
condor_submit "CH=$CH" "$D/rank.sub" || exit 5
wait_for_batch ztt_rank stage2 || exit 6
echo "ranking files: $(ls "$CH"/fit/results/ztautau/Fits/NPRanking_* 2>/dev/null | wc -l)"

echo "######## stage 3: merge, re-read, MultiFit ($(date))"
( cd "$CH/fit" && trex-fitter r ztautau.config Ranking=plot ) > "$D/out/rank_plot.log" 2>&1 || echo "WARNING: Ranking=plot failed, see $D/out/rank_plot.log"
python scripts/step5_fit.py --summarise-only || exit 7
python scripts/step5b_multifit.py || echo "WARNING: MultiFit failed"

echo "######## stage 4: report ($(date))"
python scripts/step6_report.py || exit 8
python scripts/update_docs.py || exit 9
echo "######## chain finished ($(date))"
