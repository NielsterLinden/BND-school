#!/bin/bash
# The ranking of ONE nuisance parameter of the combined fit (`trex-fitter r ... Ranking=<name>`).
# TRExFitter writes fit/results/ztautau/Fits/NPRanking_<name>_mu_Z.txt; the orchestrator merges them
# afterwards with Ranking=plot. $1 is the z-tautau checkout, $2 the parameter name, $3 the config to use.
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
cd "$1" || exit 1
shift
source ../fitting/setup.sh
set -e
set -o pipefail
cd fit
echo "==== ranking $1 @ $(date) on $(hostname) ===="
trex-fitter r "${2:-ztautau_rank.config}" "Ranking=$1"
echo "==== done @ $(date) ===="
