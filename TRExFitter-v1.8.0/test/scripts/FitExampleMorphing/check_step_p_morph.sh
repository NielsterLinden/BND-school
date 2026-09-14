#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphing" -I "library compiled with" -I "Version:" LOG_MORPH_p test/logs/FitExampleMorphing/LOG_MORPH_p && for file in `ls FitExampleMorphing/Tables/*postFit.txt`; do diff -w $file test/reference/$file; done
