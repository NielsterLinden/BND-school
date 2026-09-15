#!/bin/bash
# Build the summary deck (cvmfs TeX Live 2025; the user's TeX Live 2026 lacks translator.sty).
set -e
cd "$(dirname "$0")"
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2025/bin/x86_64-linux:$PATH
pdflatex -interaction=nonstopmode -halt-on-error zmumu_summary.tex > build.log 2>&1 || { grep -A3 "^!" build.log | head -40; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error zmumu_summary.tex >> build.log 2>&1
rm -f zmumu_summary.aux zmumu_summary.nav zmumu_summary.out zmumu_summary.snm zmumu_summary.toc zmumu_summary.log
echo "built zmumu_summary.pdf ($(du -h zmumu_summary.pdf | cut -f1))"
