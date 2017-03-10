#ifndef H_LIMIT_SCAN
#define H_LIMIT_SCAN

#include <string>
#include <vector>

#include "TGraph2D.h"
#include "TGraph.h"
#include "TH2D.h"

class TLegend;

void ReadPoints(std::vector<double> &vmx,
		std::vector<double> &vmy,
		std::vector<double> &vxsec,
		std::vector<double> &vobs,
		std::vector<double> &vobsup,
		std::vector<double> &vobsdown,
		std::vector<double> &vexp,
		std::vector<double> &vup,
		std::vector<double> &vdown,
		std::vector<double> &vsigobs,
		std::vector<double> &vsigexp);

TH2D MakeObservedSignificancePlot(std::vector<double> vmx,
                                  std::vector<double> vmy,
                                  std::vector<double> vobs);

TH2D MakeExpectedSignificancePlot(std::vector<double> vmx,
                                  std::vector<double> vmy,
                                  std::vector<double> vobs);

void MakeLimitPlot(std::vector<double> vmx,
                   std::vector<double> vmy,
                   std::vector<double> vlim,
                   std::vector<double> vobs,
                   std::vector<double> vobsup,
                   std::vector<double> vobsdown,
                   std::vector<double> vexp,
                   std::vector<double> vup,
                   std::vector<double> vdown,
                   const TH2D &hsigobs,
                   const TH2D &hsigexp);

int GetNumBins(const std::vector<double> &pts, double width);
void GetParticleNames(std::string &xparticle, std::string &yparticle);

void Style(TGraph *c, int color, int style);
		
TGraph DrawContours(TGraph2D &g2, int color, int style, double width,
		    int n_smooth, double val = 1.);

void FixGraph(TGraph &graph);
void ReverseGraph(TGraph &graph);

void SetupColors();
void SetupSignedColors();

void GetOptions(int argc, char *argv[]);

#endif
