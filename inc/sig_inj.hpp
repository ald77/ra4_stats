#ifndef H_SIG_INJ
#define H_SIG_INJ

#include <vector>

void GetOptions(int argc, char *argv[]);
void GetStats(std::vector<double> vals, double &center, double &up, double &down);
double GetMode(const std::vector<double> &v, double frac);
std::vector<double>  GetSmallestRange(const std::vector<double> &v, double frac);

#endif
