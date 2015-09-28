#include <cmath>

#include <limits>

#include "TH2D.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;

int main(){
  const int bands = 999;
  int colors[bands];
  const unsigned num = 6;
  double red[num] =   {1.,0.,0.,0.,1.,1.};
  double green[num] = {0.,0.,1.,1.,1.,0.};
  double blue[num] =  {1.,1.,1.,0.,0.,0.};
  double stops[num] = {0.,0.2,0.4,0.6,0.8,1.};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    colors[i] = fi+i;
  }
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);

  double inf = numeric_limits<double>::infinity();
  TH2D h("h", "Combined Significance;Significance A;Significance B", 55, 0., 5., 55, 0., 5.);
  h.SetMinimum(0.);
  h.SetMaximum(7.);
  h.SetStats(false);
  for(int ix = 1; ix <= h.GetNbinsX(); ++ix){
    double zx = h.GetXaxis()->GetBinCenter(ix);
    double px = erfc(zx/sqrt(2.));
    double qx = -2.*log(px);
    for(int iy = 1; iy <= h.GetNbinsY(); ++iy){
      double zy = h.GetYaxis()->GetBinCenter(iy);
      double py = erfc(zy/sqrt(2.));
      double qy = -2.*log(py);

      double q = qx + qy;
      double p = TMath::Prob(q, 4);
      double cp = 1.-0.5*p;
      double z = cp <= 0. ? -inf
	: cp >= 1. ? inf
	: TMath::NormQuantile(cp);
      h.SetBinContent(ix, iy, z);
    }
  }
  TCanvas c;
  h.Draw("colz");
  c.Print("combined_significance.pdf");
}
