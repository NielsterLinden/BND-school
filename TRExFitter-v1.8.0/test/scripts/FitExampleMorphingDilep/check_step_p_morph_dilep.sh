#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphingDilep" -I "library compiled with" -I "Version:" LOG_MORPH_DILEP_p test/logs/FitExampleMorphingDilep/LOG_MORPH_DILEP_p && for file in `ls FitExampleMorphingDilep/Tables/*postFit.txt`; do diff -w $file test/reference/$file; done
