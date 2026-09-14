#include "TRandom3.h"

double BootstrapDefault(int seed){
    gRandom->SetSeed(seed);     // set event specific seed which is different for every bootstrapid/DSID/eventnumber
    return gRandom->Poisson(1);
}
