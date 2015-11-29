#ifndef H_SIG_INJ
#define H_SIG_INJ

#include <vector>
#include <utility>

void InjectSignal(double inject, std::size_t index);
std::pair<double,double> ExtractSignal(std::size_t index, std::size_t toy, bool is_nc);
void GetOptions(int argc, char *argv[]);
void GetStats(std::vector<double> vals, double &center, double &up, double &down);
double GetMode(const std::vector<double> &v, double frac);
double GetMedian(std::vector<double> v);
std::vector<double>  GetSmallestRange(const std::vector<double> &v, double frac);
void MakePlot(const std::vector<double> &injections, const std::vector<std::vector<double> > &yvals,
	      bool is_nc, bool is_pull);
void RemoveBadResults(std::vector<double> &vals, std::vector<double> &pulls);

#endif
