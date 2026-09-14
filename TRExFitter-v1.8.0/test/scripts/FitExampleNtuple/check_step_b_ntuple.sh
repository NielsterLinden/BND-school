#!/bin/bash
diff -w -I "Developed by" -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" -I "Version:" LOG_NTUPLE_b test/logs/FitExampleNtuple/LOG_NTUPLE_b
