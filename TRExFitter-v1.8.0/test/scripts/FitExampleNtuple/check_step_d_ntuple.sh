#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" -I "library compiled with" -I "Version:" LOG_NTUPLE_d test/logs/FitExampleNtuple/LOG_NTUPLE_d && for file in `ls FitExampleNtuple/Tables/*.txt`; do diff -w $file test/reference/$file; done
