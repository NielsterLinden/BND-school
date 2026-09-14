#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleTemplate_1D" -I "library compiled with" -I "Version:" -I "Inferring initial errors" LOG_TEMPLATE_1D_p test/logs/FitExampleTemplate_1D/LOG_TEMPLATE_1D_p && for file in `ls FitExampleTemplate_1D/Tables/*postFit.txt`; do diff -w $file test/reference/$file; done
