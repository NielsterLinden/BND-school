#include "../StyleUtils/TRExStyle.C"
#include "../StyleUtils/TRExLabels.C"

TString ParseSplitColon(TString line)
{
    TObjArray* splitString = line.Tokenize(":");
    line = ((TObjString*) splitString->At(1))->String();
    line = line.Strip(TString::kBoth, '\t');
    line = line.Strip(TString::kBoth, ' ');
    line = line.Strip(TString::kBoth, '\"');
    return line;
}

TString Clean(std::string input)
{
  TString line = TString(input);
  line = line.Strip(TString::kBoth, '\t');
  line = line.Strip(TString::kBoth, ' ');
  return line;
}

void SmallText(Double_t x,Double_t y,Color_t color, const char *text) {
  TLatex l;
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}


TH1F* MakeContour(TH2 *h1, double* cl_value, Color_t col, Style_t ls, Width_t lw){
  TH1F *hist = (TH1F*) h1->Clone();
  hist->SetContour(1, cl_value);
  hist->SetLineColor(col);
  hist->SetLineStyle(ls);
  hist->SetLineWidth(lw);  
  return hist;
}



void MergeLikelihoodScan2D(TString configFileName, float zMaximum=-1.0, TString mode="XY")
{
  if (!(mode == "XY" or mode == "X" or mode == "Y") ) {
    std::cout << "Error: optional parameter mode must be X, Y, or XY" << std::endl;
    return;
  }


  // Parse the config file to identify directory and 2D variables
  TString baseDir = "";
  TString varXname = "";
  TString varYname = "";
  int nScanSteps = 30;
  int nScanStepsY = 30;
  
  struct NormFactor{
    TString name;
    TString title;
  };
  std::vector<NormFactor> NormFactorStore;
  


  std::string templine;
  ifstream config(configFileName);
  if (!config.is_open()) {
    std::cout << "ERROR! Could not open file" << std::endl;
    exit;
  }
  
  while ( getline (config, templine) ) {
    TString line = Clean(templine);
    if (line.BeginsWith("%")) continue;
    
    if (line.BeginsWith("Job:")){
      baseDir = ParseSplitColon(line);
    }

    if (line.BeginsWith("OutputDir:")){
	    baseDir = ParseSplitColon(line) + "/" + baseDir;
    }

    if (line.BeginsWith("Fit:")) {
      while ( getline(config, templine) ) {
	line = Clean(templine);
	if (line.BeginsWith("%")) continue;
	
	if (line.BeginsWith("do2DLHscan:")) {
	  TString varList = ParseSplitColon(line);
	  TObjArray* varListArray = varList.Tokenize(",");
	  varXname = ((TObjString*)varListArray->At(0))->String().Strip(TString::kBoth, ' ');
	  varYname = ((TObjString*)varListArray->At(1))->String().Strip(TString::kBoth, ' ');
	  delete varListArray;
	}

	if (line.BeginsWith("LHscanSteps:")) {
	  nScanSteps = ParseSplitColon(line).Atoi();
	}

	if (line.BeginsWith("LHscanStepsY:")) {
	  nScanStepsY = ParseSplitColon(line).Atoi();
	}

	if (line == "") break;
      }
    }
    
    if (line.BeginsWith("NormFactor:")) {
      NormFactor tempNF;
      TString normFactorName = ParseSplitColon(line);
      tempNF.name = normFactorName;
      tempNF.title = "";
      
      while ( getline(config, templine) ) {
	line = Clean(templine);
	if (line.BeginsWith("%")) continue;
	
	if (line.BeginsWith("Title:")) {
	  TString title = ParseSplitColon(line);
	  tempNF.title = title;
	}
	
	if (line == "") {
	  NormFactorStore.push_back(tempNF);
	  break;
	}
      }
    }
  }
  config.close();
  



  // Manually set TREx style so the 
  std::cout << "\nApplying TREx style settings...\n" << std::endl ;
  static TStyle* trexStyle  = TRExStyle();
  trexStyle->SetPadRightMargin(0.16);
  gROOT->SetStyle("TREx");
  gROOT->ForceStyle();

  TH2D* h1 = NULL;
  int minx = (mode == "Y") ? -1 : 0;
  int maxx = (mode == "Y") ?  0 : nScanSteps;
  int miny = (mode == "X") ? -1 : 0;
  int maxy = (mode == "X") ?  0 : nScanStepsY;
  
  for (auto i=minx; i<maxx; i++) {
    for (auto j=miny; j<maxy; j++) {
      TString filename = "NLLscan_"+varXname+"_"+varYname;
      if (minx != -1) {
        filename += Form("_StepX%d", i);
      }
      if (miny != -1) {
        filename += Form("_StepY%d", j);
      }
      filename += "_histo.root";

      if(gSystem->AccessPathName(baseDir + "/LHoodPlots/" + filename)){
        std::cout << "ERROR! Attempting to load file that does not exist:" << baseDir + "/LHoodPlots/" + filename << std::endl;
        std::cout << "Are you running in the correct mode?" << std::endl;
        return;
      }
      
      TFile *inputfile = new TFile(baseDir + "/LHoodPlots/" + filename);
      TH2D *inputhisto = (TH2D*)inputfile->Get("NLL");
      inputhisto->SetDirectory(NULL);

      if (h1 == NULL) {
	h1 = inputhisto;
	delete inputfile;
      } else {
	h1->Add(inputhisto);
	delete inputhisto;
	delete inputfile;
      }
    }
  }
  

  gStyle->SetPalette(kBird);
  TCanvas *c1 = new TCanvas("c1","",200,10,800,600);
  h1->SetTitle("");

  // Assign x axis title
  for (auto nf : NormFactorStore) {
    if (nf.name == varXname) {
      if (nf.title == "") {
	h1->GetXaxis()->SetTitle(nf.name);
      } else {
	h1->GetXaxis()->SetTitle(nf.title);
      }
      break;
    }
  }

  // Assign x axis title
  for (auto nf : NormFactorStore) {
    if (nf.name == varYname) {
      if (nf.title == "") {
	h1->GetYaxis()->SetTitle(nf.name);
      } else {
	h1->GetYaxis()->SetTitle(nf.title);
      }
      break;
    }
  }

  h1->GetZaxis()->SetTitle("-#DeltalogL");
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->GetZaxis()->CenterTitle();

  //// define CL values
  double d_cl68[1]; d_cl68[0] = 1.15; // 68% CL
  double d_cl9545[1]; d_cl9545[0] = 3.09; // 95.45% CL
  double d_cl99[1]; d_cl99[0] = 4.605; // 99% CL
  double d_cl90[1]; d_cl90[0] = 2.305; // 90% CL
  double d_cl95[1]; d_cl95[0] = 2.995; // 95% CL
  double d_cl9973[1]; d_cl9973[0] = 5.915; // 99.73% CL

  Int_t binxx, binyy, bin; 
  Int_t xbin = h1->GetXaxis()->GetNbins();
  Int_t ybin = h1->GetYaxis()->GetNbins();
  double min = h1->GetMinimum();
  std::cout << "hist minimum value: " << min << std::endl;
  for(int i=1; i<xbin+1; i++) {
    for(int j=1; j<ybin+1; j++) {
      double_t binContent = h1->GetBinContent(i, j);
      h1->SetBinContent(i, j, binContent - min);
      double_t newbinContent = h1->GetBinContent(i, j);
    }
  }  

  TH1F *h_d_cl68   = MakeContour(h1, d_cl68,   kRed+1, 2, 3);
  TH1F *h_d_cl9545 = MakeContour(h1, d_cl9545, kRed+1, 3, 3);
  TH1F *h_d_cl9973 = MakeContour(h1, d_cl9973, kRed+1, 4, 3);
  
  if (zMaximum < 0) zMaximum = h1->GetMaximum();
  h1        ->GetZaxis()->SetRangeUser(0, zMaximum);
  h_d_cl68  ->GetZaxis()->SetRangeUser(0, zMaximum);
  h_d_cl9545->GetZaxis()->SetRangeUser(0, zMaximum);
  h_d_cl9973->GetZaxis()->SetRangeUser(0, zMaximum);

  //h1->Smooth(1,"R");
  h1->DrawCopy("col z");
  h1->SetContour(1);
  h1->Draw("cont3 same");

  /// Draw contours ///
  h_d_cl68->Draw("CONT3 same");
  h_d_cl9545->Draw("CONT3 same");
  h_d_cl9973->Draw("CONT3 same");

  SmallText(0.18, 0.86, kBlack,"#sqrt{s}=13TeV, 139 fb^{-1}");  
  TRExLabel(0.18,0.9,"ATLAS","Internal");

  auto legend = new TLegend(0.7,0.83,0.85,0.93); // for 3 best contours

  legend->AddEntry(h_d_cl68,   "1#sigma", "l");
  legend->AddEntry(h_d_cl9545, "2#sigma", "l");
  legend->AddEntry(h_d_cl9973, "3#sigma", "l");    
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();

  TString outputFilePath = baseDir + "/LHoodPlots/";
  c1->SaveAs(outputFilePath + "/NLLscan_" + varXname + "_" + varYname + "_" + "Merged_histo_2D.png");

  TFile *output = new TFile(outputFilePath + "/NLLscan_" + varXname + "_" + varYname + "_" + "Merged_histo_2D.root", "RECREATE");
  output->cd();
  h1->Write();
  output->Close();
  delete output;
} 
