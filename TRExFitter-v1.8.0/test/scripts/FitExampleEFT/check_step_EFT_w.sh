 #!/bin/bash
diff -w -I "Developed by" -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExampleEFT/Pruning.png" -I "png file FitExampleEFT/NormalisationPlot.png" -I "using CPU computation library" -I "Version:" LOG_EFT_w test/logs/FitExampleEFT/LOG_EFT_w && diff FitExampleEFT/PruningText.txt test/reference/FitExampleEFT/PruningText.txt
