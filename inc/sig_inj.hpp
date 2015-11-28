#ifndef H_SIG_INJ
#define H_SIG_INJ

#include <vector>

void InjectSignal(double inject, std::size_t index);
double ExtractSignal(std::size_t index, std::size_t toy, bool is_nc);
void GetOptions(int argc, char *argv[]);
void GetStats(std::vector<double> vals, double &center, double &up, double &down);
double GetMode(const std::vector<double> &v, double frac);
std::vector<double>  GetSmallestRange(const std::vector<double> &v, double frac);
void MakePlot(const std::vector<double> &injections, const std::vector<std::vector<double> > &yvals, bool is_nc);

#endif
