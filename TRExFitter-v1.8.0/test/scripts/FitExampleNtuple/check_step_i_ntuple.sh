#!/bin/bash
diff -w -I "Developed by" -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" -I "Version:" LOG_NTUPLE_i test/logs/FitExampleNtuple/LOG_NTUPLE_i && diff -w FitExampleNtuple/Fits/GroupedImpact_mu_XS_ttH.txt test/reference/FitExampleNtuple/Fits/GroupedImpact_mu_XS_ttH.txt
