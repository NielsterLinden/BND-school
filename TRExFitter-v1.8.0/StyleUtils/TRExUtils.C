#include "TLatex.h"

void myText(Double_t x,Double_t y,Color_t color, const char *text) {

  TLatex l; 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
