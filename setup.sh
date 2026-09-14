#!/bin/bash
# Shared environment for every channel of the BND-school Z analysis.
#
#     source setup.sh          # from anywhere; safe to source twice
#
# Provides: LCG_110 (ROOT 6.40 + python 3.13 + uproot/awkward/hist/mplhep/iminuit),
# the TRExFitter v1.8.0 build (trex-fitter on PATH once fitting/build_trexfitter.sh ran),
# and the environment variables the analysis code reads.
source /cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc14-opt/setup.sh

export BND_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TREXFITTER_HOME="${BND_ROOT}/TRExFitter-v1.8.0"
case ":${PATH}:" in *":${TREXFITTER_HOME}/build/bin:"*) ;; *) export PATH="${TREXFITTER_HOME}/build/bin:${PATH}" ;; esac
case ":${LD_LIBRARY_PATH}:" in *":${TREXFITTER_HOME}/build/lib:"*) ;; *) export LD_LIBRARY_PATH="${TREXFITTER_HOME}/build/lib:${LD_LIBRARY_PATH}" ;; esac

# Where the v2 skims live (bulk copy on /data; laptop bundle documented in z-mumu/handoff.md).
export BND_SKIM_DIR="${BND_SKIM_DIR:-/data/atlas/users/nterlind/BND-school-cache/skims_v2}"

if [ ! -x "${TREXFITTER_HOME}/build/bin/trex-fitter" ]; then
    echo "[setup.sh] trex-fitter is not built yet: run  bash ${BND_ROOT}/fitting/build_trexfitter.sh"
fi
