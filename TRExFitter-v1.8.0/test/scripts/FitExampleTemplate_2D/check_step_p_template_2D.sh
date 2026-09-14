#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleTemplate_2D" -I "library compiled with" -I "Version:" -I "Inferring initial errors" LOG_TEMPLATE_2D_p test/logs/FitExampleTemplate_2D/LOG_TEMPLATE_2D_p && for file in `ls FitExampleTemplate_2D/Tables/*postFit.txt`; do diff -w $file test/reference/$file; done
