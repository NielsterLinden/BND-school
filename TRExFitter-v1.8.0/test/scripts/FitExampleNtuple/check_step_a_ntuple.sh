#!/bin/bash
diff -w -I "Developed by" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" -I "library compiled with" -I "Version:" LOG_NTUPLE_a test/logs/FitExampleNtuple/LOG_NTUPLE_a
