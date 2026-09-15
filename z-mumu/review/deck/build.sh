#!/bin/bash
# Build the review deck. The user's TeX Live 2026 lacks translator.sty (beamer); the cvmfs TeX Live 2025 is complete.
set -e
cd "$(dirname "$0")"
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2025/bin/x86_64-linux:$PATH
pdflatex -interaction=nonstopmode -halt-on-error zmumu_review.tex > build.log 2>&1 || { grep -A3 "^!" build.log | head -40; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error zmumu_review.tex >> build.log 2>&1
rm -f zmumu_review.aux zmumu_review.nav zmumu_review.out zmumu_review.snm zmumu_review.toc zmumu_review.log
echo "built zmumu_review.pdf ($(du -h zmumu_review.pdf | cut -f1))"
