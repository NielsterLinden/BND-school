#include <iostream>

#include "TRExStyle.h"

#include "TROOT.h"

void SetTRExStyle ()
{
  static TStyle* trexStyle = 0;
  if ( trexStyle==0 ) trexStyle = TRExStyle();
  gROOT->SetStyle("TREx");
  gROOT->ForceStyle();
  gStyle->SetErrorX(0.5);
}

TStyle* TRExStyle() 
{
  TStyle *trexStyle = new TStyle("TREx","TRExFitter style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  trexStyle->SetFrameBorderMode(icol);
  trexStyle->SetFrameFillColor(icol);
  trexStyle->SetCanvasBorderMode(icol);
  trexStyle->SetCanvasColor(icol);
  trexStyle->SetPadBorderMode(icol);
  trexStyle->SetPadColor(icol);
  trexStyle->SetStatColor(icol);

  // set the paper & margin sizes
  trexStyle->SetPaperSize(20,26);

  // set margin sizes
  trexStyle->SetPadTopMargin(0.05);
  trexStyle->SetPadRightMargin(0.05);
  trexStyle->SetPadBottomMargin(0.16);
  trexStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  trexStyle->SetTitleXOffset(1.4);
  trexStyle->SetTitleYOffset(1.4);

  // use large fonts
  // --- Michele
  Int_t font=43;
  Double_t tsize=21;
  // ---
  trexStyle->SetTextFont(font);

  trexStyle->SetTextSize(tsize);
  trexStyle->SetLabelFont(font,"x");
  trexStyle->SetTitleFont(font,"x");
  trexStyle->SetLabelFont(font,"y");
  trexStyle->SetTitleFont(font,"y");
  trexStyle->SetLabelFont(font,"z");
  trexStyle->SetTitleFont(font,"z");

  trexStyle->SetLabelSize(tsize,"x");
  trexStyle->SetTitleSize(tsize,"x");
  trexStyle->SetLabelSize(tsize,"y");
  trexStyle->SetTitleSize(tsize,"y");
  trexStyle->SetLabelSize(tsize,"z");
  trexStyle->SetTitleSize(tsize,"z");

  // --- Michele
  trexStyle->SetLegendFont(font);
  // ---
  
  // use bold lines and markers
  trexStyle->SetMarkerStyle(20);
  trexStyle->SetMarkerSize(1.2);
  trexStyle->SetHistLineWidth(2.);
  trexStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  // get rid of error bar caps
  trexStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  trexStyle->SetOptTitle(0);
  trexStyle->SetOptStat(0);
  trexStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  trexStyle->SetPadTickX(1);
  trexStyle->SetPadTickY(1);

  return trexStyle;

}

