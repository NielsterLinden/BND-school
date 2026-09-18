#!/bin/bash
# One TRExFitter job of step 5 on a worker node: $1 is the z-tautau checkout, the rest goes straight to
# scripts/step5_fit.py.
# Workers are Debian: the submit file wraps this in the ATLAS almalinux9 image so that the LCG_110
# el9 view of ../fitting/setup.sh works (nikhef-condor skill).
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
cd "$1" || exit 1
shift
source ../fitting/setup.sh
set -e
set -o pipefail
echo "==== step5 $* @ $(date) on $(hostname) ===="
python scripts/step5_fit.py "$@"
echo "==== done @ $(date) ===="
