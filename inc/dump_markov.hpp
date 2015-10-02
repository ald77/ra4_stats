#ifndef H_DUMP_MARKOV
#define H_DUMP_MARKOV

#include <map>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"

#include "RooArgSet.h"
#include "RooDataSet.h"

RooDataSet * GetData(TFile &file);

void GetValues(RooWorkspace &w,
	       RooDataSet &data,
	       std::map<std::string, std::vector<double> > &values,
	       std::vector<double> &weights);

void FillValues(std::map<std::string, std::vector<double> > &values,
                const RooWorkspace &w);

void FillValues(std::map<std::string, std::vector<double> > &values,
                const RooArgSet &args);

void MakePlots(const std::map<std::string, std::vector<double> > &values,
	       const std::vector<double> &weights);

void MakePlot(const std::string &name, const std::vector<double> &values,
	      const std::vector<double> &weights);

#endif
