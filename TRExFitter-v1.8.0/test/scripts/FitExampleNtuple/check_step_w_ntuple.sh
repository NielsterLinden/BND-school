 #!/bin/bash
diff -w -I "Developed by" -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple/Pruning.png" -I "using CPU computation library" -I "Version:" LOG_NTUPLE_w test/logs/FitExampleNtuple/LOG_NTUPLE_w && diff FitExampleNtuple/PruningText.txt test/reference/FitExampleNtuple/PruningText.txt
