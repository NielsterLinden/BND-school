#include "TRExLabels.h"

#include "TLatex.h"
#include "TPad.h"

void TRExLabel(Double_t x,Double_t y, const char* experiment, const char* text, Color_t color) {
  TLatex l; 
  l.SetNDC();
  l.SetTextFont(73);
  l.SetTextColor(color);
  double delx = 0.11*550*gPad->GetWh()/(472*gPad->GetWw());
  l.DrawLatex(x,y,experiment);
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(43);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
  }
}

void TRExLabelNew(Double_t x,Double_t y, const char* experiment, const char* text, Color_t color, float text_size, float delx) {

  TLatex l;
  l.SetNDC();
  l.SetTextFont(73);
  l.SetTextColor(color);
  l.SetTextSize(text_size);

  l.DrawLatex(x,y,experiment);
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(43);
    p.SetTextColor(color);
    p.SetTextSize(text_size);
    p.DrawLatex(x+delx,y,text);
  }
  return;

}
