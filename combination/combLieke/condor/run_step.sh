#!/bin/bash
# One job of the combination DAG on a worker node:  run_step.sh <combLieke directory> <run.py arguments>
# The workers are Debian; the submit files (mf/condor.py) wrap this in the ATLAS almalinux9 image so that the
# el9 LCG_110 view of ../../fitting/setup.sh and the TRExFitter v1.8.0 build work (nikhef-condor skill).
# No `set -e` before the environment is up: the LCG setup returns non-zero during normal operation.
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
cd "$1" || exit 1
shift
source ../../fitting/setup.sh
set -e
set -o pipefail
echo "==== run.py $* @ $(date) on $(hostname) ===="
python run.py "$@"
echo "==== done @ $(date) ===="
