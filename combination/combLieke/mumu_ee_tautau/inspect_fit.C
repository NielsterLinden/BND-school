#include <TFile.h>
#include <TKey.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooStats/ModelConfig.h>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

// Run from this directory: root -l -b -q inspect_fit.C
void inspect_fit(const char* job="combLieke", const char* workspace_path="results/combLieke/ws_combined.root", const char* output_path="fit_diagnostics.json") {
  TFile file((std::string("results/")+job+"/Fits/"+job+".root").c_str());
  RooFitResult *fit = nullptr;
  for (auto key : *file.GetListOfKeys()) {
    fit = dynamic_cast<RooFitResult*>(file.Get(key->GetName()));
    if (fit) break;
  }
  if (!fit) throw std::runtime_error("Missing RooFitResult");
  TFile wsfile(workspace_path);
  auto *ws = dynamic_cast<RooWorkspace*>(wsfile.Get("combWS"));
  if (!ws) {
    for (auto key : *wsfile.GetListOfKeys()) {
      ws = dynamic_cast<RooWorkspace*>(wsfile.Get(key->GetName()));
      if (ws) break;
    }
  }
  if (!ws) throw std::runtime_error("Missing combined workspace");
  auto *mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
  auto *sim = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
  auto *data = ws->data("obsData");
  if (!sim || !data) throw std::runtime_error("Missing simultaneous PDF or observed data");
  std::ofstream out(output_path); out << std::setprecision(16);
  out << "{\n  \"status\": " << fit->status() << ",\n  \"covariance_quality\": " << fit->covQual()
      << ",\n  \"edm\": " << fit->edm() << ",\n  \"min_nll\": " << fit->minNll() << ",\n  \"status_history\": [";
  for (unsigned i=0; i<fit->numStatusHistory(); ++i) {
    if (i) out << ',';
    out << "{\"label\": \"" << fit->statusLabelHistory(i) << "\", \"code\": " << fit->statusCodeHistory(i) << '}';
  }
  out << "],\n  \"parameters\": [";
  bool first = true;
  for (auto arg : fit->floatParsFinal()) {
    auto *p = dynamic_cast<RooRealVar*>(arg); if (!p) continue;
    if (!first) out << ','; first = false;
    bool boundary = (p->hasMin() && std::abs(p->getVal()-p->getMin()) < 1.e-5)
                 || (p->hasMax() && std::abs(p->getVal()-p->getMax()) < 1.e-5);
    out << "\n    {\"name\": \"" << p->GetName() << "\", \"value\": " << p->getVal()
        << ", \"hesse_error\": " << p->getError() << ", \"minos_up\": " << p->getAsymErrorHi()
        << ", \"minos_down\": " << p->getAsymErrorLo() << ", \"at_boundary\": " << (boundary?"true":"false") << '}';
  }
  out << "\n  ],\n  \"workspace_regions\": ["; first = true;
  for (auto const& state : sim->indexCat()) {
    if (!first) out << ','; first = false;
    const std::string cut = std::string(sim->indexCat().GetName()) + "==" + sim->indexCat().GetName() + "::" + state.first;
    std::unique_ptr<RooAbsData> channelData(data->reduce(cut.c_str()));
    out << "{\"name\": \"" << state.first << "\", \"entries\": " << channelData->numEntries()
        << ", \"sum_weights\": " << channelData->sumEntries() << '}';
  }
  out << "]\n}\n";
}
