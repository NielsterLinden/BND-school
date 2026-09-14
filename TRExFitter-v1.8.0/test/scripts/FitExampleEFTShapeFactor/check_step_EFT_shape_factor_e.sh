#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleEFTShapeFactor" -I "library compiled with" -I "Version:" LOG_EFT_SHAPE_FACTOR_e test/logs/FitExampleEFTShapeFactor/LOG_EFT_SHAPE_FACTOR_e && for file in `ls FitExampleEFTShapeFactor/EFT/*.txt`; do diff -w $file test/reference/$file; done


