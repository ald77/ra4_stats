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
void GetStats(const std::vector<double> &vals,
              double &mean, double &median,
              double &up, double &down,
              double &up2, double &down2);
double GetValue(std::vector<double> vals, double fraction);
double GetMode(const std::vector<double> &v, double frac);
std::vector<double>  GetSmallestRange(const std::vector<double> &v, double frac);
void MakePlot(const std::vector<double> &injections, const std::vector<std::vector<double> > &yvals,
	      bool is_nc, bool is_pull);
void Plot1D(double inj, const std::vector<double> &vals,
            double mean, double median,
            double up, double down,
            double up2, double down2,
            bool is_nc, bool is_pull);
void RemoveBadResults(std::vector<double> &vals, std::vector<double> &pulls);
void MergeWithText(std::vector<double> &inj,
                   std::vector<std::vector<double> > &yvals,
                   std::vector<std::vector<double> > &pulls,
                   bool is_nc);
std::size_t GetIndex(const std::vector<double> &v, double x);
void SortByInjectionStrength(std::vector<double> &inj,
                             std::vector<std::vector<double> > &yvals,
                             std::vector<std::vector<double> > &pulls);

#endif
