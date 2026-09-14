#ifndef __TREXLABELS_H
#define __TREXLABELS_H

#include "Rtypes.h"

void TRExLabel(Double_t x,Double_t y, const char* experiment, const char* text, Color_t color = 1);

void TRExLabelNew(Double_t x,Double_t y, const char* experiment, const char* text, Color_t color, float text_size=0, float delx = 0.16);

#endif // __TREXLABELS_H
