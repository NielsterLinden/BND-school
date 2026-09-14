#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleEFT" -I "library compiled with" -I "Version:" LOG_EFT_p test/logs/FitExampleEFT/LOG_EFT_p && for file in `ls FitExampleEFT/Tables/*postFit.txt`; do diff -w $file test/reference/$file; done
