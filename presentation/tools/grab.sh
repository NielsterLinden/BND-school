#!/usr/bin/env bash
# grab.sh <clip.mp4> <out.png> [time_s | end]   -- one frame, for inspection
set -euo pipefail
if [[ "${3:-end}" == "end" ]]; then
  ffmpeg -y -loglevel error -sseof -0.15 -i "$1" -frames:v 1 -update 1 "$2"
else
  ffmpeg -y -loglevel error -ss "$3" -i "$1" -frames:v 1 -update 1 "$2"
fi
echo "$2"
