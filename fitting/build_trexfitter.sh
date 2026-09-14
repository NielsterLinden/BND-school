#!/bin/bash
# Build TRExFitter v1.8.0 (the git submodule TRExFitter-v1.8.0/) against the LCG_110 ROOT.
#
#     bash fitting/build_trexfitter.sh [jobs]
#
# ROOT 6.40 in LCG_110 is compiled with C++20, so the build must use the same standard
# (TRExFitter's CMakeLists defaults to 17 through a CACHE variable, hence the override).
set -eo pipefail   # no -u: the LCG setup.sh reads unset variables
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${HERE}/setup.sh" >/dev/null
SRC="${TREXFITTER_HOME}"
if [ ! -f "${SRC}/CMakeLists.txt" ] || [ ! -f "${SRC}/xroofit/CMakeLists.txt" ]; then
    echo "submodules missing: run  git submodule update --init --recursive"; exit 1
fi
# v1.8.0 was written for C++17: std::filesystem::path::u8string() returns std::u8string under
# C++20, which does not convert to std::string. Use string() instead (paths are ASCII here).
# The sed is idempotent, so a fresh `git submodule update` checkout is patched again on rebuild.
sed -i 's/\.u8string()/.string()/g' "${SRC}/Root/Common.cc"
cmake -S "${SRC}" -B "${SRC}/build" -DCMAKE_CXX_STANDARD=20 -DUSE_LOCAL_XROOFIT=ON
cmake --build "${SRC}/build" -j "${1:-12}" 2>&1 | tee "${SRC}/build/build.log" | grep -E "error|Error|Built target|warning: .*deprecated" | grep -vE "^\s*$" | tail -40
ls -l "${SRC}/build/bin/trex-fitter"
