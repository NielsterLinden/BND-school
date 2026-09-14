#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleUnfoldingNormalized" -I "library compiled with" -I "Version:" LOG_UNFOLDING_NORM_p test/logs/FitExampleUnfoldingNormalized/LOG_UNFOLDING_NORM_p && for file in `ls FitExampleUnfoldingNormalized/Tables/*postFit.txt`; do diff -w $file test/reference/$file; done
