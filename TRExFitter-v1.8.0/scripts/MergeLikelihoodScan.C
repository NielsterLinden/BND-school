#include "../StyleUtils/TRExStyle.C"

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


void MergeLikelihoodScan(TString configFileName)
{
 // Parse the config file to identify directory and 2D variables
  TString baseDir = "";
  TString varXname = "";
  int nScanSteps = 30;
  
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

	if (line.BeginsWith("doLHscan:")) {
	  varXname = ParseSplitColon(line);
	}

	if (line.BeginsWith("LHscanSteps:")) {
	  nScanSteps = ParseSplitColon(line).Atoi();
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

  SetTRExStyle();
  std::vector<double> x_values;
  std::vector<double> y_values;

  for (auto i=0; i<nScanSteps; i++) {
    TString filename = "NLLscan_" + varXname + Form("Step%d", i) + "_curve.root";
    TFile *inputfile = new TFile(baseDir + "/LHoodPlots/" + filename);
    TGraph *inputgraph = (TGraph*)inputfile->Get("LHscan");
    x_values.push_back(inputgraph->GetPointX(0));
    y_values.push_back(inputgraph->GetPointY(0));
    delete inputgraph;
    delete inputfile;
  }
  
  // Compute minimum
  double mnll = 9999999;
  for (auto &y : y_values) {
    if (y < mnll) mnll = y;
  }
  for (auto &y : y_values) {
    y = y - mnll;
  }
  
  
  // Plot the graph
  // Code lifted from TRexFit.C
  // Minor edits to variable names, that's all
  TCanvas can("NLLscan");
    
  TGraph graph(x_values.size(), &x_values[0], &y_values[0]);
  graph.Draw("ALP");
  graph.GetXaxis()->SetRangeUser(x_values[0],x_values.back());

  // Assign x axis title
  for (auto nf : NormFactorStore) {
    if (nf.name == varXname) {
      if (nf.title == "") {
	graph.GetXaxis()->SetTitle(nf.name);
      } else {
	graph.GetXaxis()->SetTitle(nf.title);
      }
      break;
    }
  }

  graph.GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");

  can.SetTitle("NLLscan_merged");
  can.SetName("NLLscan_merged");
  can.cd();

  TLatex tex{};
  tex.SetTextColor(kGray+2);

  TLine l1s(x_values[0],0.5,x_values.back(),0.5);
  l1s.SetLineStyle(kDashed);
  l1s.SetLineColor(kGray);
  l1s.SetLineWidth(2);
  if(graph.GetMaximum()>2){
    l1s.Draw();
    tex.DrawLatex(x_values.back(),0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
  }

  if(graph.GetMaximum()>2){
    TLine l2s(x_values[0],2,x_values.back(),2);
    l2s.SetLineStyle(kDashed);
    l2s.SetLineColor(kGray);
    l2s.SetLineWidth(2);
    l2s.Draw();
    tex.DrawLatex(x_values.back(),2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
  }
   
  if(graph.GetMaximum()>4.5){
    TLine l3s(x_values[0],4.5,x_values.back(),4.5);
    l3s.SetLineStyle(kDashed);
    l3s.SetLineColor(kGray);
    l3s.SetLineWidth(2);
    l3s.Draw();
    tex.DrawLatex(x_values.back(),4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
  }
  
  TLine lv0(0,graph.GetMinimum(),0,graph.GetMaximum());
  lv0.Draw();
  
  TLine lh0(x_values[0],0,x_values.back(),0);
  lh0.Draw();

  can.RedrawAxis();

  TString outputFilePath = baseDir + "/LHoodPlots/";
  can.SaveAs(outputFilePath + "/NLLscan_" + varXname + "_Merged_curve.png");
  TFile *outputFile = new TFile(outputFilePath + "/NLLscan_"+varXname+"_Merged_curve.root", "RECREATE");
  outputFile->cd();
  graph.Write("LHscan");
  outputFile->Close();
} 
