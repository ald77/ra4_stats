#include "sig_inj.hpp"

#include <cmath>

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <getopt.h>

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFitResult.h"

#include "utilities.hpp"
#include "styles.hpp"

using namespace std;

namespace{
  int ntoys = 5;
  vector<double> injections;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  styles style("RA4");
  style.setDefaultStyle();

  vector<vector<double> > yvals(injections.size());
  for(size_t i = 0; i < injections.size(); ++i){
    cout << "Injecting signal with strength " << injections.at(i) << endl;
    double x = injections.at(i);
    ostringstream oss;
    oss << "./run/make_workspace.exe --method m1bk --lumi 1.264 --use_r4 --toys " << ntoys-1
	<< " --sig_strength " << x << flush;
    execute(oss.str());
    string file_name = execute("ls -rt m1bk*.root | tail -n 1");
    for(int itoy = 0; itoy < ntoys; ++itoy){
      cout << "Evaluating toy " << (itoy+1) << " of " << ntoys << endl;
      oss.str("");
      oss << "export blah=$(pwd); cd ~/cmssw/CMSSW_7_1_5/src; eval `scramv1 runtime -sh`; cd $blah; combine -M MaxLikelihoodFit --saveWorkspace --saveWithUncertainties --minos=all -w w";
      if(itoy != 0) oss << '_' << itoy;
      oss << " " << file_name << flush;
      execute(oss.str());
      TFile r_file("mlfit.root","read");
      if(!r_file.IsOpen()) continue;
      oss.str("");
      RooFitResult *f = static_cast<RooFitResult*>(r_file.Get("fit_s"));
      if(f == nullptr) continue;
      RooArgList pars = f->floatParsFinal();
      for(int ipar = 0; ipar < pars.getSize(); ++ipar){
	RooRealVar *var = static_cast<RooRealVar*>(pars.at(ipar));
	if(var == nullptr) continue;
	if(var->GetName() != string("r")) continue;
	yvals.at(i).push_back(var->getVal());
      }
    }
  }

  vector<double> centers(injections.size()), ups(injections.size()), downs(injections.size()), zeros(injections.size(), 0.);
  vector<double> bots(centers.size()), tops(centers.size());
  for(size_t i = 0; i < injections.size(); ++i){
    GetStats(yvals.at(i), centers.at(i), ups.at(i), downs.at(i));
    tops.at(i) = centers.at(i) + ups.at(i);
    bots.at(i) = centers.at(i) - downs.at(i);
  }

  TGraphAsymmErrors g(injections.size(), &injections.at(0), &centers.at(0),
		      &zeros.at(0), &zeros.at(0),
		      &downs.at(0), &ups.at(0));
  g.SetMarkerStyle(20);
  g.SetMarkerSize(2);
  g.SetMarkerColor(1);
  g.SetLineStyle(1);
  g.SetLineWidth(5);
  g.SetLineColor(1);
  double xmin = *min_element(injections.cbegin(), injections.cend());
  double xmax = *max_element(injections.cbegin(), injections.cend());
  double xdelta = xmax-xmin;
  double ymin = *min_element(bots.cbegin(), bots.cend());
  double ymax = *max_element(tops.cbegin(), tops.cend());
  double ydelta = ymax-ymin;
  double margin = 0.05;
  double themin = min(xmin-margin*xdelta, ymin-margin*ydelta);
  double themax = max(xmax+margin*xdelta, ymax+margin*ydelta);
  TH1D h("", ";Injected Signal Strength;Extracted Signal Strength", 1, themin, themax);
  h.SetMinimum(themin);
  h.SetMaximum(themax);
  h.Fill(0.5, 1.0);
  TF1 f("", "x", themin, themax);
  f.SetLineColor(2);
  f.SetLineWidth(4);
  f.SetLineStyle(2);
  TCanvas c;
  h.Draw("axis");
  f.Draw("same");
  g.Draw("p same");
  c.Print("siginj.pdf");

  cout << "Extracted signal strengths:" << endl;
  for(size_t i = 0; i < injections.size(); ++i){
    cout << injections.at(i) << ": " << centers.at(i) << " + " << ups.at(i) << " - " << downs.at(i) << ": ";
    for(const auto &y: yvals.at(i)){
      cout << y << " ";
    }
    cout << endl;
  }
}

void GetStats(vector<double> vals, double &center, double &up, double &down){
  sort(vals.begin(), vals.end());
  double se = erf(1./sqrt(2.));
  vector<double> vband = GetSmallestRange(vals, se);
  center = GetMode(vband, se);
  up = vband.back()-center;
  down = center-vband.front();
}

double GetMode(const vector<double> &v, double frac){
  if(v.size() == 0){
    return 0.;
  }else if( v.size() == 1){
    return v.front();
  }else{
    vector<double> temp = GetSmallestRange(v, frac);
    if(temp.size() == v.size()){
      if(v.size() == 2){
	return v.front()+0.5*(v.back()-v.front());
      }else{
	double fdiff = v.at(1)-v.at(0);
	double bdiff = v.at(2)-v.at(1);
	if(fdiff>bdiff){
	  return GetMode(vector<double>(v.cbegin()+1, v.cend()), frac);
	}else{
	  return GetMode(vector<double>(v.cbegin(), v.cend()-1), frac);
	}
      }
    }else{
      return GetMode(temp, frac);
    }
  }
}

vector<double> GetSmallestRange(const vector<double> &v, double frac){
  if(v.size() <= 1){
    return v;
  }else{
    if(frac < 0.) frac = 0.;
    if(frac > 1.) frac = 1.;
    size_t dist = ceil(frac*(v.size()-1));
    if(dist == v.size()-1) return v;
    size_t best_start = 0;
    double best_sep = v.back() - v.front();
    for(size_t i = 0; i < v.size()-dist; ++i){
      double sep = v.at(i+dist)-v.at(i);
      if(sep < best_sep){
	best_start = i;
	best_sep = sep;
      }
    }
    return vector<double>(v.cbegin()+best_start, v.cbegin()+best_start+dist);
  }
}

void GetOptions(int argc, char *argv[]){
  set<double> injection_set;
  while(true){
    static struct option long_options[] = {
      {"toys", required_argument, 0, 't'},
      {"inject", required_argument, 0, 'i'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "t:i:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 't':
      ntoys = atoi(optarg);
      break;
    case 'i':
      injection_set.insert(atof(optarg));
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
  if(injection_set.size() == 0){
    injections.resize(2);
    injections.at(0) = 0.;
    injections.at(1) = 1.;
  }else{
    injections = vector<double>(injection_set.cbegin(), injection_set.cend());
  }
}
