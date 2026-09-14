#!/bin/bash
diff -w -I "Developed by" -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" -I "using CPU computation library" -I "Version:" -I "Creation of NLL object" -I "Evaluated function and gradient" -I "MnHesse Done after" -I "VariableMetricBuilder Stop iterating after" LOG_NTUPLE_MIXED_f test/logs/FitExampleNtuple/LOG_NTUPLE_MIXED_f
