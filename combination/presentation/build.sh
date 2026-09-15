#!/bin/bash
# Build the combination deck.
#
#     bash build.sh            # figures (if missing) + pdflatex twice + checks
#     bash build.sh --figures  # regenerate figures/*.pdf first
#
# The system TeX Live 2020 in /bin has no metropolis/pgfopts/siunitx; the cvmfs TeX Live 2025 is
# complete. Same route the z-mumu review deck uses.
set -e
cd "$(dirname "$0")"
JOB=bnd_z_combination

if [ "$1" = "--figures" ] || [ ! -f figures/precision_budget.pdf ]; then
    echo "== deck figures"
    ( source ../../setup.sh >/dev/null 2>&1; python figures.py )
fi
if [ ! -f ../output/plots/forest.pdf ]; then
    echo "!! ../output/plots is empty -- run  cd .. && python run_combination.py"
    exit 1
fi

export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2025/bin/x86_64-linux:$PATH
command -v pdflatex >/dev/null || { echo "!! no pdflatex: is /cvmfs/sft.cern.ch mounted?"; exit 1; }

echo "== pdflatex (pass 1/2)"
pdflatex -interaction=nonstopmode -halt-on-error $JOB.tex > build.log 2>&1 \
    || { echo "!! LaTeX error:"; grep -A4 "^!" build.log | head -40; exit 1; }
echo "== pdflatex (pass 2/2)"
pdflatex -interaction=nonstopmode -halt-on-error $JOB.tex >> build.log 2>&1

# report anything that would show as a hole or a '??' on a slide
if grep -qE "Warning: (Reference|Citation).*undefined" build.log; then
    echo "!! undefined references:"; grep -E "Warning: (Reference|Citation).*undefined" build.log
fi
if grep -q "File .* not found" build.log; then
    echo "!! missing figures:"; grep "File .* not found" build.log
fi
OVERFULL=$(grep -c "Overfull \\\\hbox" build.log || true)
[ "$OVERFULL" -gt 0 ] && echo "   ($OVERFULL overfull hboxes -- check the slides they name in build.log)"

rm -f $JOB.aux $JOB.nav $JOB.out $JOB.snm $JOB.toc $JOB.vrb
echo "== built $JOB.pdf ($(du -h $JOB.pdf | cut -f1), $(pdfinfo $JOB.pdf 2>/dev/null | awk '/^Pages/{print $2}') pages)"
