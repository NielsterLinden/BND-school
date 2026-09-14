 #!/bin/bash
diff -w -I "Developed by" -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphingDilep/Pruning.png" -I "using CPU computation library" -I "Version:" LOG_MORPH_DILEP_w test/logs/FitExampleMorphingDilep/LOG_MORPH_DILEP_w && diff FitExampleMorphingDilep/PruningText.txt test/reference/FitExampleMorphingDilep/PruningText.txt
