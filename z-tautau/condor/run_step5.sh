#!/bin/bash
# One TRExFitter job of step 5 on a worker node. Arguments are passed straight to scripts/step5_fit.py.
# Workers are Debian: the submit file wraps this in the ATLAS almalinux9 image so that the LCG_110
# el9 view of ../setup.sh works (nikhef-condor skill).
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
cd /project/atlas/users/sjankovy/BND/BND-school/z-tautau || exit 1
source ../setup.sh
set -e
set -o pipefail
echo "==== step5 $* @ $(date) on $(hostname) ===="
python scripts/step5_fit.py "$@"
echo "==== done @ $(date) ===="
