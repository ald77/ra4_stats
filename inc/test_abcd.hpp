#ifndef H_TEST_ABCD
#define H_TEST_ABCD

#include <ios>
#include <string>
#include <vector>
#include <random>

#include "TH1D.h"

std::string GetParamString(int argc, char *argv[]);

double GetTestStatistic(double a, double b, double c, double d);

std::vector<double> SampleTestStatistic(std::mt19937 &prng, size_t n, double a, double b, double c, double d);

double QToP(double q);
double QToZ(double q);
double ZToQ(double z);
double ZToP(double z);
double PToQ(double p);
double PToZ(double p);

std::vector<double> GetSignificances(const std::vector<double> &b_qs, const std::vector<double> &sb_qs);

std::vector<double> GetSignificances(const std::vector<double> &sb_qs);

double FindFraction(const std::vector<double> &qs, double q);

void GetPredictions(double a_in, double b_in, double c_in, double d_in,
                    double &a_out, double &b_out, double &c_out, double &d_out);

double BinLogLikelihood(double obs, double pred);

void PrintRates(std::ostream &out, const std::string &name,
                double a, double b, double c, double d);

void InitializePRNG(std::mt19937 &prng);

void PlotQDistributions(const std::string &params,
                        const std::vector<double> &b_qs,
                        const std::vector<double> &sb_qs);

void PlotZDistributions(const std::string &params,
                        const std::vector<double> &b_qs,
                        const std::vector<double> &sb_qs,
                        double z = -1.);

void NormalizeWithOverflow(TH1D &h);

double AddInQuadrature(double x, double y);

double RealMax(const TH1D &h);

#endif
