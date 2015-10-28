#ifndef H_GET_DILEPTON_UNCERTAINTIES
#define H_GET_DILEPTON_UNCERTAINTIES

#include <vector>
#include <string>

#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooRealVar.h"

std::vector<std::string> GetBinNames(const RooWorkspace &w);

double GetBkgPred(const RooWorkspace &w,
                  const std::string &bin_name);

double GetBkgPredErr(const RooWorkspace &w,
                     const RooFitResult &f,
                     const std::string &bin_name);

RooRealVar * SetVariables(RooWorkspace &w,
                          const RooFitResult &f);

#endif
