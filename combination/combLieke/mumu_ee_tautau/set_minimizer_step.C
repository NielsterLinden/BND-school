#include <TFile.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <stdexcept>
void set_minimizer_step() {
 TFile f("results/combLieke/ws_combined.root", "UPDATE");
 auto *w=dynamic_cast<RooWorkspace*>(f.Get("combWS"));
 if(!w) throw std::runtime_error("No combined workspace");
 auto *g=w->var("gamma_stat_tautau_SR2_bin_0");
 if(!g) throw std::runtime_error("No target gamma");
 std::cout << "Initial minimizer step: " << g->getError() << " -> 1; value, range, constant flag, and likelihood unchanged" << std::endl;
 g->setError(1.0);
 w->Write("combWS",TObject::kOverwrite);
}
