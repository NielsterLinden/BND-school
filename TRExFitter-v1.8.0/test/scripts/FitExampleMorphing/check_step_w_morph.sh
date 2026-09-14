 #!/bin/bash
diff -w -I "Developed by" -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphing/Pruning.png" -I "using CPU computation library" -I "Version:" LOG_MORPH_w test/logs/FitExampleMorphing/LOG_MORPH_w && diff FitExampleMorphing/PruningText.txt test/reference/FitExampleMorphing/PruningText.txt
