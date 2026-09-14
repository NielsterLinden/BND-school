#!/bin/bash
diff -w -I "Developed by" -I "Real time" -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "library compiled with" -I "Version:" -I "NumericIntegration" -I "using CPU computation library" -I "png file FitExampleEFTShapeFactor" -I "Creation of NLL object" -I "Evaluated function and gradient" -I "MnHesse Done after" -I "VariableMetricBuilder Stop iterating after" LOG_EFT_SHAPE_FACTOR_f test/logs/FitExampleEFTShapeFactor/LOG_EFT_SHAPE_FACTOR_f && for file in `ls FitExampleEFTShapeFactor/Fits/*.txt`; do diff -w $file test/reference/$file; done

