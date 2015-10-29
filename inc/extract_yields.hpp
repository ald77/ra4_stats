#ifndef H_EXTRACT_YIELDS
#define H_EXTRACT_YIELDS

#include <string>
#include <vector>
#include <limits>

#include "TH1D.h"
#include "TGraphErrors.h"

#include "RooWorkspace.h"
#include "RooFitResult.h"

std::string GetSignalName(const RooWorkspace &w);

std::string TexFriendly(const std::string &s);

void PrintTable(RooWorkspace &w,
                const RooFitResult &f,
                const std::string &file_name);

double GetMCYield(const RooWorkspace &w,
                  const std::string &bin_name,
                  const std::string &prc_name);

double GetMCTotal(const RooWorkspace &w,
                  const std::string &bin_name);
double GetMCTotalErr(const RooWorkspace &w,
                     const RooFitResult &f,
                     const std::string &bin_name);

double GetBkgPred(const RooWorkspace &w,
                  const std::string &bin_name);
double GetBkgPredErr(const RooWorkspace &w,
                     const RooFitResult &f,
                     const std::string &bin_name);

double GetSigPred(const RooWorkspace &w,
                  const std::string &bin_name);
double GetSigPredErr(const RooWorkspace &w,
                     const RooFitResult &f,
                     const std::string &bin_name);

double GetTotPred(const RooWorkspace &w,
                  const std::string &bin_name);
double GetTotPredErr(const RooWorkspace &w,
                     const RooFitResult &f,
                     const std::string &bin_name);

double GetObserved(const RooWorkspace &w,
                   const std::string &bin_name);

double GetLambda(const RooWorkspace &w,
                 const std::string &bin_name);
double GetLambdaErr(const RooWorkspace &w,
                    const RooFitResult &f,
                    const std::string &bin_name);

RooRealVar * SetVariables(RooWorkspace &w,
                          const RooFitResult &f);

void MakeYieldPlot(RooWorkspace &w,
                   const RooFitResult &f,
                   const std::string &file_name);

std::vector<std::string> GetBinNames(const RooWorkspace &w);
std::vector<std::string> GetPlainBinNames(const RooWorkspace &w);
std::vector<std::string> GetProcessNames(const RooWorkspace &w);

std::vector<std::vector<double> > GetComponentYields(const RooWorkspace &w,
                                                     const std::vector<std::string> &bin_names,
                                                     const std::vector<std::string> &process_names);

std::vector<TH1D> MakeBackgroundHistos(const std::vector<std::vector<double> > &component_yields,
                                       const std::vector<std::string> &bin_names,
                                       const std::vector<std::string> &prc_names);

TH1D MakeTotalHisto(const RooWorkspace &w,
                    const RooFitResult &f,
                    const std::vector<std::string> &bin_names);

TH1D MakeObserved(const RooWorkspace &w,
                  const std::vector<std::string> &bin_names);

void SetBounds(TH1D &a,
               TH1D &b,
               std::vector<TH1D> &cs);

double GetMaximum(const TH1D &a,
                  const TH1D &b,
                  const std::vector<TH1D> &cs);

double GetMinimum(const TH1D &a,
                  const TH1D &b,
                  const std::vector<TH1D> &cs);

double GetMaximum(const TH1D &h, double y = std::numeric_limits<double>::max());
double GetMinimum(const TH1D &h, double y = -std::numeric_limits<double>::max());

TGraphErrors MakeErrorBand(const TH1D &h);
TGraphErrors MakeRatio(const TH1D &num, const TH1D &den);

void MakeCorrectionPlot(RooWorkspace &w,
			const RooFitResult &f,
			const std::string &file_name);

#endif
