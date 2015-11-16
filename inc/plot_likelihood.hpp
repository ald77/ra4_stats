#ifndef H_PLOT_LIKELIHOOD
#define H_PLOT_LIKELIHOOD

#include <vector>
#include <string>

#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooRealVar.h"

void PlotVars(RooWorkspace &w, bool bkg_only);
void GetOptions(int argc, char *argv[]);
RooRealVar * SetVariables(RooWorkspace &w,
                          const RooFitResult &f);
std::vector<std::string> GetVarNames(const RooWorkspace &w);

#endif
