#!/bin/bash
# Get and build TRExFitter v1.8.0 (the mandatory version) against the LCG_110 ROOT.
#
#     bash fitting/build_trexfitter.sh [jobs]
#
# TRExFitter is not part of this repository: the source is cloned at tag v1.8.0 into
# $TREXFITTER_HOME (default: ../TRExFitter-v1.8.0, next to the repository) if it is not there yet.
#
# ROOT 6.40 in LCG_110 is compiled with C++20, so the build must use the same standard
# (TRExFitter's CMakeLists defaults to 17 through a CACHE variable, hence the override).
set -eo pipefail   # no -u: the LCG setup.sh reads unset variables
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/setup.sh" >/dev/null
SRC="${TREXFITTER_HOME}"
if [ ! -f "${SRC}/CMakeLists.txt" ]; then
    git clone --branch v1.8.0 --depth 1 --recurse-submodules --shallow-submodules \
        https://gitlab.cern.ch/TRExStats/TRExFitter.git "${SRC}"
fi
if [ ! -f "${SRC}/xroofit/CMakeLists.txt" ]; then
    git -C "${SRC}" submodule update --init --recursive
fi
# v1.8.0 was written for C++17: std::filesystem::path::u8string() returns std::u8string under
# C++20, which does not convert to std::string. Use string() instead (paths are ASCII here).
# The sed is idempotent, so a fresh checkout is patched again on rebuild.
sed -i 's/\.u8string()/.string()/g' "${SRC}/Root/Common.cc"
cmake -S "${SRC}" -B "${SRC}/build" -DCMAKE_CXX_STANDARD=20 -DUSE_LOCAL_XROOFIT=ON
cmake --build "${SRC}/build" -j "${1:-12}" 2>&1 | tee "${SRC}/build/build.log" | grep -E "error|Error|Built target|warning: .*deprecated" | grep -vE "^\s*$" | tail -40
ls -l "${SRC}/build/bin/trex-fitter"
