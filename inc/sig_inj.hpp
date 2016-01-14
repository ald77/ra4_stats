#ifndef H_SIG_INJ
#define H_SIG_INJ

#include <vector>
#include <string>
#include <utility>

#include "RooAbsReal.h"
#include "RooFitResult.h"

void InjectSignal(const std::string id_string, double inject, std::size_t index);
std::pair<double,double> ExtractSignal(const std::string id_string, std::size_t index, std::size_t toy, bool is_nc);
double GetAsymError(const RooRealVar &var, double &ehi, double &elo);
double GetError(const RooAbsReal &var, const RooFitResult &f);
void GetOptions(int argc, char *argv[]);
void GetStats(const std::vector<double> &vals, double &center, double &up, double &down);
double GetValue(std::vector<double> vals, double fraction);
double GetMode(const std::vector<double> &v, double frac);
double GetMedian(std::vector<double> v);
std::vector<double>  GetSmallestRange(const std::vector<double> &v, double frac);
void MakePlot(const std::vector<double> &injections, const std::vector<std::vector<double> > &yvals,
	      bool is_nc, bool is_pull);
void Plot1D(const std::string &base, double inj, const std::vector<double> &vals,
	    double center, double up, double down);
void RemoveBadResults(std::vector<double> &vals, std::vector<double> &pulls);

#endif
