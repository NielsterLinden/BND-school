#!/bin/bash
# The ranking of ONE nuisance parameter of the combined fit (`trex-fitter r ... Ranking=<name>`).
# TRExFitter writes fit/results/ztautau/Fits/NPRanking_<name>_mu_Z.txt; the orchestrator merges them
# afterwards with Ranking=plot. $1 is the parameter name, $2 the config to use.
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
cd /project/atlas/users/sjankovy/BND/BND-school/z-tautau || exit 1
source ../setup.sh
set -e
set -o pipefail
cd fit
echo "==== ranking $1 @ $(date) on $(hostname) ===="
trex-fitter r "${2:-ztautau_rank.config}" "Ranking=$1"
echo "==== done @ $(date) ===="
