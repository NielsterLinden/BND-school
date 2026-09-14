#!/bin/bash
diff -w -I "Developed by" -I "Version:" -I "Creation of NLL object" -I "Evaluated function and gradient" -I "MnHesse Done after" -I "VariableMetricBuilder Stop iterating after" FitExampleUnfolding/Fits/FitExampleUnfolding.txt test/reference/FitExampleUnfolding/Fits/FitExampleUnfolding.txt && diff -w FitExampleUnfolding/Fits/Unfolding_UnfoldedResults.txt test/reference/FitExampleUnfolding/Fits/Unfolding_UnfoldedResults.txt && for file in `ls FitExampleUnfolding/Fits/*.txt`; do diff -w $file test/reference/$file; done

