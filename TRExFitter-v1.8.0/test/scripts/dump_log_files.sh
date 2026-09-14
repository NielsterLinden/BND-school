#!/bin/bash

echo ""
echo "#########################################################################"
echo "TRExFitter validation logfiles"
echo "#########################################################################"
echo ""
echo "You are about to change the logfiles used online to validate the results."
echo "Note that ANY difference with the current logfiles must be justified in the"
echo "merge request form to make the bookmark of changes easier."
echo ""
echo "Are you sure you want to continue ? [y/n]"
read -n 1 ok_to_continue
echo ""
if [[ "${ok_to_continue}" != "y" ]]; then
  echo "Stopping the execution of the script !"
  return 0
fi
echo ""

##
## Checks if ROOT is setup up by any other way as the setup.sh of the folder
##
which_root=`which root`
which_root=""
if [[ "$which_root" != "" ]]; then
  echo "Looks like ROOT has been setup :( Stopping the execution of the script now !"
  return 0;
fi

##
## Sources ROOT with THE good way
##
echo ""
source setup.sh #setting up ROOT as defined in the package
which_root=`which root`
which_root="toto"#FIXME
if [[ "$which_root" == "" ]]; then
  echo "For some reasons, the setup failed (access to /cvmfs/ ?) Stopping the "
  echo "execution of the script now !"
  return 0;
fi

##
## Remove the output files
##

echo "Removing the output files, 'rm -rf FitExample*'"
rm -rf FitExample*

##
## Compiling the code
##
echo ""
echo "Compiling the code. Cleaning it first, and then recompile."
rm -rf build/
mkdir build && cd build/
cmake -DTREXFITTER_AVX_OPTIMISATION=FALSE ../
make -j4
cd ..
if [[ ! -f build/bin/trex-fitter ]]; then
  echo "!!ERROR!! The binary file is not found. Need to investigate !!"
  return 0
fi
echo ""

##
## Now, actually runs the thing
##
echo ""
echo "Things seem to be in order now. I am about to produce the logfiles. Note that "
echo "you still can abort the process at any moment."
echo ""
for step in h w f l s r d p ; do
  echo "==> $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExample.config >& LOG_$step
  cat LOG_$step | grep -v "TRExFitter" >& test/logs/FitExample/LOG_$step
  rm -f LOG_$step
done

echo "Copying files..."
cp FitExample/PruningText.txt test/reference/FitExample/
cp FitExample/Tables/*.txt test/reference/FitExample/Tables/
cp FitExample/Fits/FitExample.txt test/reference/FitExample/Fits/
cp FitExample/Fits/NPRanking_SigXsecOverSM.txt test/reference/FitExample/Fits/
cp FitExample/Fits/FitExample_errDecomp_SigXsecOverSM.txt test/reference/FitExample/Fits/
cp FitExample/Fits/FitExample_group_errDecomp_*.txt test/reference/FitExample/Fits/

for step in l s ; do
  echo "==> $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleToys.config >& LOG_${step}_toys
  cat LOG_${step}_toys | grep -v "TRExFitter" >& test/logs/FitExample/LOG_${step}_toys
  rm -f LOG_${step}_toys
done

for step in h w f ; do
  echo "==> stat only $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleStatOnly.config  >& LOG_STATONLY_$step
  cat LOG_STATONLY_$step | grep -v "TRExFitter" >& test/logs/FitExampleStatOnly/LOG_STATONLY_$step
  rm -f LOG_STATONLY_$step
done

echo "Copying files..."
cp FitExampleStatOnly/Fits/FitExampleStatOnly.txt test/reference/FitExampleStatOnly/Fits/
cp FitExampleStatOnly/Fits/FitExampleStatOnly_errDecomp_SigXsecOverSM.txt test/reference/FitExampleStatOnly/Fits/
cp FitExampleStatOnly/Fits/FitExampleStatOnly_group_errDecomp_*.txt test/reference/FitExampleStatOnly/Fits/

for step in h w f d p ; do
  echo "==> Morphing $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleMorphing.config >& LOG_MORPH_$step
  cat LOG_MORPH_$step | grep -v "TRExFitter" >& test/logs/FitExampleMorphing/LOG_MORPH_$step
  rm -f LOG_MORPH_$step
done

echo "Copying files..."
cp FitExampleMorphing/PruningText.txt test/reference/FitExampleMorphing/
cp FitExampleMorphing/Tables/*.txt test/reference/FitExampleMorphing/Tables/
cp FitExampleMorphing/Fits/FitExampleMorphing.txt test/reference/FitExampleMorphing/Fits/
cp FitExampleMorphing/Fits/FitExampleMorphing_errDecomp_topWidth.txt test/reference/FitExampleMorphing/Fits/
cp FitExampleMorphing/Fits/FitExampleMorphing_group_errDecomp_*.txt test/reference/FitExampleMorphing/Fits/

for step in h w f d p ; do
  echo "==> Morphing dilepton $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleMorphingDilep.config >& LOG_MORPH_DILEP_$step
  cat LOG_MORPH_DILEP_$step | grep -v "TRExFitter" >& test/logs/FitExampleMorphingDilep/LOG_MORPH_DILEP_$step
  rm -f LOG_MORPH_DILEP_$step
done

echo "Copying files..."
cp FitExampleMorphingDilep/PruningText.txt test/reference/FitExampleMorphingDilep/
cp FitExampleMorphingDilep/Tables/*.txt test/reference/FitExampleMorphingDilep/Tables/
cp FitExampleMorphingDilep/Fits/FitExampleMorphingDilep.txt test/reference/FitExampleMorphingDilep/Fits/
cp FitExampleMorphingDilep/Fits/FitExampleMorphingDilep_errDecomp_topWidth.txt test/reference/FitExampleMorphingDilep/Fits/
cp FitExampleMorphingDilep/Fits/FitExampleMorphingDilep_group_errDecomp_*.txt test/reference/FitExampleMorphingDilep/Fits/

for step in h t w f d p ; do
  echo "==> Template_1D $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleTemplate_1D.config >& LOG_TEMPLATE_1D_$step
  cat LOG_TEMPLATE_1D_$step | grep -v "TRExFitter" >& test/logs/FitExampleTemplate_1D/LOG_TEMPLATE_1D_$step
  rm -f LOG_TEMPLATE_1D_$step
done

echo "Copying files..."
cp FitExampleTemplate_1D/Tables/*.txt test/reference/FitExampleTemplate_1D/Tables/
cp FitExampleTemplate_1D/Fits/FitExampleTemplate_1D.txt test/reference/FitExampleTemplate_1D/Fits/
cp FitExampleTemplate_1D/Fits/FitExampleTemplate_1D_errDecomp_topMass.txt test/reference/FitExampleTemplate_1D/Fits/
cp FitExampleTemplate_1D/Templates/Parametrization.txt test/reference/FitExampleTemplate_1D/Templates/
cp FitExampleTemplate_1D/Fits/FitExampleTemplate_1D_group_errDecomp_*.txt test/reference/FitExampleTemplate_1D/Fits/

for step in h t w f d p ; do
  echo "==> Template_2D $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleTemplate_2D.config >& LOG_TEMPLATE_2D_$step
  cat LOG_TEMPLATE_2D_$step | grep -v "TRExFitter" >& test/logs/FitExampleTemplate_2D/LOG_TEMPLATE_2D_$step
  rm -f LOG_TEMPLATE_2D_$step
done

echo "Copying files..."
cp FitExampleTemplate_2D/Tables/*.txt test/reference/FitExampleTemplate_2D/Tables/
cp FitExampleTemplate_2D/Fits/FitExampleTemplate_2D.txt test/reference/FitExampleTemplate_2D/Fits/
cp FitExampleTemplate_2D/Fits/FitExampleTemplate_2D_errDecomp_topWidth.txt test/reference/FitExampleTemplate_2D/Fits/
cp FitExampleTemplate_2D/Templates/Parametrization.txt test/reference/FitExampleTemplate_2D/Templates/
cp FitExampleTemplate_2D/Fits/FitExampleTemplate_2D_group_errDecomp_*.txt test/reference/FitExampleTemplate_2D/Fits/

for step in h t w f ; do
  echo "==> Template_2D_AD $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleTemplate_2D_AD.config >& LOG_TEMPLATE_2D_AD_$step
  cat LOG_TEMPLATE_2D_AD_$step | grep -v "TRExFitter" >& test/logs/FitExampleTemplate_2D_AD/LOG_TEMPLATE_2D_AD_$step
  rm -f LOG_TEMPLATE_2D_AD_$step
done

echo "Copying files..."
cp FitExampleTemplate_2D_AD/Fits/FitExampleTemplate_2D_AD.txt test/reference/FitExampleTemplate_2D_AD/Fits/
cp FitExampleTemplate_2D_AD/Fits/FitExampleTemplate_2D_AD_errDecomp_topWidth.txt test/reference/FitExampleTemplate_2D_AD/Fits/
cp FitExampleTemplate_2D_AD/Templates/Parametrization.txt test/reference/FitExampleTemplate_2D_AD/Templates/
cp FitExampleTemplate_2D_AD/Fits/FitExampleTemplate_2D_AD_group_errDecomp_*.txt test/reference/FitExampleTemplate_2D_AD/Fits/

for step in n b w f d p i a ; do
  echo "==> Ntuple $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleNtuple.config >& LOG_NTUPLE_$step
  cat LOG_NTUPLE_$step | grep -v "TRExFitter" >& test/logs/FitExampleNtuple/LOG_NTUPLE_$step
  rm -f LOG_NTUPLE_$step
done

echo "Copying files..."
cp FitExampleNtuple/PruningText.txt test/reference/FitExampleNtuple/
cp FitExampleNtuple/Tables/*.txt test/reference/FitExampleNtuple/Tables/
cp FitExampleNtuple/Fits/FitExampleNtuple.txt test/reference/FitExampleNtuple/Fits/
cp FitExampleNtuple/Fits/GroupedImpact_mu_XS_ttH.txt test/reference/FitExampleNtuple/Fits/
cp FitExampleNtuple/Fits/FitExampleNtuple_errDecomp_mu_XS_ttH.txt test/reference/FitExampleNtuple/Fits/
cp FitExampleNtuple/Fits/FitExampleNtuple_group_errDecomp_*.txt test/reference/FitExampleNtuple/Fits/

echo "==> Ntuple mixed f step ongoing"
./build/bin/trex-fitter f test/configs/FitExampleNtupleMixed.config >& LOG_NTUPLE_MIXED_f
cat LOG_NTUPLE_MIXED_f | grep -v "TRExFitter" >& test/logs/FitExampleNtuple/LOG_NTUPLE_MIXED_f
rm -f LOG_NTUPLE_MIXED_f

echo "==> Ntuple combined nwfdp step ongoing"
./build/bin/trex-fitter nwfdp test/configs/FitExampleNtupleAllSteps.config >& LOG_NTUPLE_nwfdp
cat LOG_NTUPLE_nwfdp | grep -v "TRExFitter" >& test/logs/FitExampleNtupleAllSteps/LOG_NTUPLE_nwfdp
rm -f LOG_NTUPLE_nwfdp

for step in w f; do
  echo "==> Morphing multifit $step step ongoing"
  ./build/bin/trex-fitter m$step test/configs/FitExampleMorphMultifit.config >& LOG_MORPH_MULTI_$step
  cat LOG_MORPH_MULTI_$step | grep -v "TRExFitter" >& test/logs/FitExampleMorphMultifit/LOG_MORPH_MULTI_$step
  rm -f LOG_MORPH_MULTI_$step
done

echo "Copying files..."
cp FitExampleMorphMultifit/Fits/FitExampleMorphMultifit.txt test/reference/FitExampleMorphMultifit/Fits/
cp FitExampleMorphMultifit/Fits/FitExampleMorphMultifit_errDecomp_topWidth.txt test/reference/FitExampleMorphMultifit/Fits/
cp FitExampleMorphMultifit/Fits/FitExampleMorphMultifit_group_errDecomp_*.txt test/reference/FitExampleMorphMultifit/Fits/

for step in u h w f; do
  echo "==> Unfolding $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleUnfolding.config >& LOG_UNFOLDING_$step
  cat LOG_UNFOLDING_$step | grep -v "TRExFitter" >& test/logs/FitExampleUnfolding/LOG_UNFOLDING_$step
  rm -f LOG_UNFOLDING_$step
done

echo "Copying files..."
cp FitExampleUnfolding/PruningText.txt test/reference/FitExampleUnfolding/
cp FitExampleUnfolding/Fits/FitExampleUnfolding.txt test/reference/FitExampleUnfolding/Fits/
cp FitExampleUnfolding/Fits/Unfolding_UnfoldedResults.txt test/reference/FitExampleUnfolding/Fits/
cp FitExampleUnfolding/Fits/FitExampleUnfolding_errDecomp_Unfolding_Bin_*_mu.txt test/reference/FitExampleUnfolding/Fits/
cp FitExampleUnfolding/Fits/FitExampleUnfolding_group_errDecomp_*.txt test/reference/FitExampleUnfolding/Fits/

for step in u h w f d p; do
  echo "==> Normalized unfolding $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleUnfoldingNormalized.config >& LOG_UNFOLDING_NORM_$step
  cat LOG_UNFOLDING_NORM_$step | grep -v "TRExFitter" >& test/logs/FitExampleUnfoldingNormalized/LOG_UNFOLDING_NORM_$step
  rm -f LOG_UNFOLDING_NORM_$step
done

echo "Copying files..."
cp FitExampleUnfoldingNormalized/PruningText.txt test/reference/FitExampleUnfoldingNormalized/
cp FitExampleUnfoldingNormalized/Fits/FitExampleUnfoldingNormalized.txt test/reference/FitExampleUnfoldingNormalized/Fits/
cp FitExampleUnfoldingNormalized/Fits/Unfolding_UnfoldedResults.txt test/reference/FitExampleUnfoldingNormalized/Fits/
cp FitExampleUnfoldingNormalized/Fits/FitExampleUnfoldingNormalized_errDecomp_Unfolding_Bin_*_mu.txt test/reference/FitExampleUnfoldingNormalized/Fits/
cp FitExampleUnfoldingNormalized/Tables/*.txt test/reference/FitExampleUnfoldingNormalized/Tables/
cp FitExampleUnfoldingNormalized/Fits/FitExampleUnfoldingNormalized_group_errDecomp_*.txt test/reference/FitExampleUnfoldingNormalized/Fits/

for step in n e w f d p; do
  echo "==> EFT $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleEFT.config >& LOG_EFT_$step
  cat LOG_EFT_$step | grep -v "TRExFitter" >& test/logs/FitExampleEFT/LOG_EFT_$step
  rm -f LOG_EFT_$step
done

echo "Copying files..."
cp FitExampleEFT/PruningText.txt test/reference/FitExampleEFT/
cp FitExampleEFT/EFT/*.txt test/reference/FitExampleEFT/EFT/
cp FitExampleEFT/Fits/FitExampleEFT.txt test/reference/FitExampleEFT/Fits/
cp FitExampleEFT/Fits/FitExampleEFT_errDecomp_citW.txt test/reference/FitExampleEFT/Fits/
cp FitExampleEFT/Fits/FitExampleEFT_errDecomp_ctW.txt test/reference/FitExampleEFT/Fits/
cp FitExampleEFT/Fits/FitExampleEFT_errDecomp_ttbar_NORM.txt test/reference/FitExampleEFT/Fits/
cp FitExampleEFT/Tables/*.txt test/reference/FitExampleEFT/Tables/
cp FitExampleEFT/Fits/FitExampleEFT_group_errDecomp_*.txt test/reference/FitExampleEFT/Fits/

for step in h e w f; do
  echo "==> EFT shape factor $step step ongoing"
  ./build/bin/trex-fitter $step test/configs/FitExampleEFTShapeFactor.config >& LOG_EFT_SHAPE_FACTOR_$step
  cat LOG_EFT_SHAPE_FACTOR_$step | grep -v "TRExFitter" >& test/logs/FitExampleEFTShapeFactor/LOG_EFT_SHAPE_FACTOR_$step
  rm -f LOG_EFT_SHAPE_FACTOR_$step
done

echo "Copying files..."
cp FitExampleEFTShapeFactor/PruningText.txt test/reference/FitExampleEFTShapeFactor/
cp FitExampleEFTShapeFactor/EFT/*.txt test/reference/FitExampleEFTShapeFactor/EFT/
cp FitExampleEFTShapeFactor/Fits/*.txt test/reference/FitExampleEFTShapeFactor/Fits/

##
## Making a git status and asks if the files have to be added
##
echo ""
echo "Logfiles have been produced. So far, nothing has been added to git. This is "
echo "your responsability to do so. A few advices: "
echo "  - use git diff <path to the file> to see what has changed."
echo "  - if you understand ALL changes, do git add <path to the file>"
echo ""
echo "When doing git commit: "
echo "  - please mention which test files have been updated and why in the commit"
echo "    message"
echo ""
echo "Push the change to your branch on the git server, and check the result of the"
echo "build, here: https://gitlab.cern.ch/TRExStats/TRExFitter/pipelines"
echo ""
