#include "test_abcd.hpp"

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <string>
#include <random>
#include <functional>
#include <array>
#include <algorithm>

#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"

using namespace std;

namespace{
  size_t num_bkg_toys = 1e7;
  size_t num_sig_toys = 1e5;

  double inf = numeric_limits<double>::infinity();
}

int main(int argc, char *argv[]){
  if(argc != 9){
    cerr << "Must supply 8 rates as arguments: bkg A, bkg B, bkg C, bkg D, sig A, sig B, sig C, sig D" << endl;
    return EXIT_FAILURE;
  }

  //Get a string to record input args
  string params = GetParamString(argc, argv);

  //Read background parameters from command line
  double bkg_a = atof(argv[1]);
  double bkg_b = atof(argv[2]);
  double bkg_c = atof(argv[3]);
  double bkg_d = atof(argv[4]);

  //Fix background to close under ABCD
  GetPredictions(bkg_a, bkg_b, bkg_c, bkg_d,
                 bkg_a, bkg_b, bkg_c, bkg_d);

  //Read background parameters from command line
  double sig_a = atof(argv[5]);
  double sig_b = atof(argv[6]);
  double sig_c = atof(argv[7]);
  double sig_d = atof(argv[8]);

  //Compute S+B yields for printing
  double tot_a = bkg_a + sig_a;
  double tot_b = bkg_b + sig_b;
  double tot_c = bkg_c + sig_c;
  double tot_d = bkg_d + sig_d;

  //Get ABCD predictions for printing
  double pred_a, pred_b, pred_c, pred_d;
  GetPredictions(tot_a, tot_b, tot_c, tot_d,
                 pred_a, pred_b, pred_c, pred_d);

  //Print rates
  PrintRates(cout, "Background", bkg_a, bkg_b, bkg_c, bkg_d);
  PrintRates(cout, "    Signal", sig_a, sig_b, sig_c, sig_d);
  PrintRates(cout, "     Total", tot_a, tot_b, tot_c, tot_d);
  PrintRates(cout, "      ABCD", pred_a, pred_b, pred_c, pred_d);

  //Compute -2 ln(l_null/l_alt) and significance for Asimov data set
  double q = GetTestStatistic(tot_a, tot_b, tot_c, tot_d);
  double z = QToZ(q);
  cout << "Z-score: " << setprecision(6) << z << endl;

  //PRNG for toys
  mt19937 prng;
  InitializePRNG(prng);

  //Generate toys
  vector<double> b_qs = SampleTestStatistic(prng, num_bkg_toys, bkg_a, bkg_b, bkg_c, bkg_d);
  vector<double> sb_qs = SampleTestStatistic(prng, num_sig_toys, tot_a, tot_b, tot_c, tot_d);

  //Make plots
  PlotQDistributions(params, b_qs, sb_qs);
  PlotZDistributions(params, b_qs, sb_qs, z);
}

string GetParamString(int argc, char *argv[]){
  ostringstream oss;
  for(int i = 1; i < argc - 1; ++i){
    oss << argv[i] << "_";
  }
  oss << argv[argc-1];

  return oss.str();
}

double GetTestStatistic(double a, double b, double c, double d){
  //Get ABCD predictions with signal added
  double pred_a, pred_b, pred_c, pred_d;
  GetPredictions(a, b, c, d,
                 pred_a, pred_b, pred_c, pred_d);

  //Compute log-likelihoods
  double null_ll = BinLogLikelihood(a, pred_a)
    + BinLogLikelihood(b, pred_b)
    + BinLogLikelihood(c, pred_c)
    + BinLogLikelihood(d, pred_d);
  double alt_ll = BinLogLikelihood(a, a)
    + BinLogLikelihood(b, b)
    + BinLogLikelihood(c, c)
    + BinLogLikelihood(d, d);

  //Compute the likelihood ratio test statistic
  double q = 2.*(alt_ll-null_ll);
  return q >= 0. ? q : 0.;
}

double QToP(double q){
  return erfc(sqrt(0.5*q));
}

double QToZ(double q){
  return sqrt(q);
}

double ZToQ(double z){
  return z*z;
}

double ZToP(double z){
  return erfc(z/sqrt(2.));
}

double PToZ(double p){
  double cp = 1.-0.5*p;
  return cp <= 0. ? -inf
    : cp >= 1. ? inf
    : TMath::NormQuantile(cp);
}

double PToQ(double p){
  return ZToQ(PToZ(p));
}

vector<double> GetSignificances(const vector<double> &b_qs, const vector<double> &sb_qs){
  vector<double> out(sb_qs.size());
  for(size_t i = 0; i < out.size(); ++i){
    out.at(i) = PToZ(FindFraction(b_qs, sb_qs.at(i)));
  }
  return out;
}

vector<double> GetSignificances(const vector<double> &sb_qs){
  vector<double> out(sb_qs.size());
  for(size_t i = 0; i < out.size(); ++i){
    out.at(i) = QToZ(sb_qs.at(i));
  }
  return out;
}

double FindFraction(const vector<double> &qs, double q){
  return static_cast<double>(distance(lower_bound(qs.begin(), qs.end(), q), qs.end()))/qs.size();
}

vector<double> SampleTestStatistic(mt19937 & prng, size_t n, double a, double b, double c, double d){
  vector<double> out(n);
  poisson_distribution<unsigned> pa(a), pb(b), pc(c), pd(d);
  for(size_t i = 0; i < n; ++i){
    out.at(i) = GetTestStatistic(a > 0. ? pa(prng) : 0.,
                                 b > 0. ? pb(prng) : 0.,
                                 c > 0. ? pc(prng) : 0.,
                                 d > 0. ? pd(prng) : 0.);
  }
  sort(out.begin(), out.end());
  return out;
}

void GetPredictions(double a_in, double b_in, double c_in, double d_in,
                    double &a_out, double &b_out, double &c_out, double &d_out){
  double total = a_in + b_in + c_in + d_in;
  a_out = total ? (a_in+b_in)*(a_in+c_in)/total : 0.;
  b_out = total ? (a_in+b_in)*(b_in+d_in)/total : 0.;
  c_out = total ? (a_in+c_in)*(c_in+d_in)/total : 0.;
  d_out = total ? (b_in+d_in)*(c_in+d_in)/total : 0.;
}

double BinLogLikelihood(double obs, double pred){
  return (pred > 0.) ? (obs*log(pred) - pred) : (obs > 0. ? -inf : 0.);
}

void PrintRates(ostream &out, const string &name,
                double a, double b, double c, double d){
  out << fixed << showpoint << setprecision(2);
  out
    << name << ": "
    << ' ' << setw(8) << a
    << ' ' << setw(8) << b
    << ' ' << setw(8) << c
    << ' ' << setw(8) << d
    << endl;
}

void InitializePRNG(mt19937 &prng){
  array<int, mt19937::state_size> vals;
  random_device rd;
  generate(vals.begin(), vals.end(), ref(rd));
  seed_seq ss(vals.begin(), vals.end());
  prng.seed(ss);
}

void PlotQDistributions(const string &params,
                        const vector<double> &b_qs,
                        const vector<double> &sb_qs){
  TH1D hb("hb", (params+";-2*ln(L_{null}/L_{alt});PDF").c_str(), 100, 0., 30.);
  hb.SetLineColor(kGreen+2);
  hb.SetLineStyle(1);
  hb.SetLineWidth(4);
  hb.SetStats(false);
  TH1D hsb("hsb", (params+";-2*ln(L_{null}/L_{alt});PDF").c_str(), 100, 0., 30.);
  hsb.SetLineColor(kRed+2);
  hsb.SetLineStyle(1);
  hsb.SetLineWidth(4);
  hsb.SetStats(false);

  for(const double &q: b_qs) hb.Fill(q);
  for(const double &q: sb_qs) hsb.Fill(q);

  NormalizeWithOverflow(hb);
  NormalizeWithOverflow(hsb);

  TF1 pdf("pdf", "ROOT::Math::chisquared_pdf(x,1.)", 0., 30.);
  pdf.SetTitle((params+";-2*ln(L_{null}/L_{alt});PDF").c_str());
  pdf.SetNpx(1000);
  pdf.SetLineColor(kBlack);
  pdf.SetLineStyle(1);
  pdf.SetLineWidth(4);
  pdf.SetMaximum(1.);
  pdf.SetMinimum(1e-8);

  TLegend l(0.1, 0.1, 0.4, 0.4);
  l.AddEntry(&hsb, "S+B", "le");
  l.AddEntry(&hb, "B-only", "le");
  l.AddEntry(&pdf, "#chi^{2}(1)", "l");

  TCanvas c;
  c.SetLogy();
  pdf.Draw();
  hb.Draw("e0 same");
  hsb.Draw("e0 same");
  l.Draw("same");
  c.Print(("q_"+params+".pdf").c_str());
}

void PlotZDistributions(const string &params,
                        const vector<double> &b_qs,
                        const vector<double> &sb_qs,
                        double z){
  TH1D htoy("htoy", (params+";Significance;PDF").c_str(), 100, 0., 10.);
  htoy.SetLineColor(kRed+2);
  htoy.SetLineStyle(1);
  htoy.SetLineWidth(4);
  htoy.SetStats(false);
  TH1D hasympt("hasympt", (params+";Significance;PDF").c_str(), 100, 0., 10.);
  hasympt.SetLineColor(kBlack);
  hasympt.SetLineStyle(2);
  hasympt.SetLineWidth(4);
  hasympt.SetStats(false);
  TH1D hnull("hnull", (params+";Significance;PDF").c_str(), 100, 0., 10.);
  hnull.SetLineColor(kGreen+2);
  hnull.SetLineStyle(1);
  hnull.SetLineWidth(4);
  hnull.SetStats(false);
  TF1 pdf("pdf", "2*TMath::Gaus(x,0,1,1)", 0., 10.);
  pdf.SetNpx(1000);
  pdf.SetTitle((params+";Significance;PDF").c_str());
  pdf.SetLineColor(kBlack);
  pdf.SetLineStyle(1);
  pdf.SetLineWidth(4);

  vector<double> toy_zs = GetSignificances(b_qs, sb_qs);

  for(const double &zi: toy_zs) htoy.Fill(zi);
  for(const double &q: sb_qs) hasympt.Fill(QToZ(q));
  for(const double &q: b_qs) hnull.Fill(QToZ(q));

  NormalizeWithOverflow(htoy);
  NormalizeWithOverflow(hasympt);
  NormalizeWithOverflow(hnull);

  pdf.SetMinimum(0.);
  pdf.SetMaximum(1.);

  TLine line(z, 0., z, 1.);
  line.SetLineColor(kRed+2);
  line.SetLineStyle(3);
  line.SetLineWidth(2);

  TLegend leg(0.6, 0.6, 0.9, 0.9);
  leg.AddEntry(&htoy, "S+B (toys)", "l");
  leg.AddEntry(&hnull, "B-only (toys)", "l");
  leg.AddEntry(&hasympt, "S+B (asympt)", "l");
  leg.AddEntry(&pdf, "B-only (asympt)", "l");

  TCanvas c;
  pdf.Draw();
  hasympt.Draw("hist same");
  hnull.Draw("hist same");
  htoy.Draw("hist same");
  leg.Draw("same");
  if(z >= 0.) line.Draw("same");
  c.Print(("z_"+params+".pdf").c_str());
}

void NormalizeWithOverflow(TH1D &h){
  int nbins = h.GetNbinsX();
  double entries = h.GetEntries();
  h.Sumw2();
  h.SetBinContent(1, h.GetBinContent(0)+h.GetBinContent(1));
  h.SetBinError(1, hypot(h.GetBinError(0), h.GetBinError(1)));
  h.SetBinContent(nbins, h.GetBinContent(nbins)+h.GetBinContent(nbins+1));
  h.SetBinError(nbins, hypot(h.GetBinError(nbins), h.GetBinError(nbins+1)));

  h.SetBinContent(0, 0.);
  h.SetBinError(0, 0.);
  h.SetBinContent(nbins+1, 0.);
  h.SetBinError(nbins+1, 0.);

  h.Scale(1./h.Integral("width"));
  h.SetEntries(entries);
}

double RealMax(const TH1D &h){
  return h.GetBinContent(h.GetMaximumBin());
}
