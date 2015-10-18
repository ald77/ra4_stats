#ifndef H_EXTRACT_YIELDS
#define H_EXTRACT_YIELDS

#include <string>
#include <vector>
#include <limits>

#include "TH1D.h"
#include "TGraphErrors.h"

#include "RooWorkspace.h"
#include "RooFitResult.h"

void MakeYieldPlot(RooWorkspace &w,
		   const RooFitResult &f,
		   const std::string &file_name);

std::vector<std::string> GetBinNames(const RooWorkspace &w);
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

#endif
