#include "sig_inj.hpp"

#include <cmath>

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <future>
#include <utility>
#include <mutex>

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

  mutex global_mutex;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  styles style("RA4");
  style.setDefaultStyle();

  cout << "Creating workspaces..." << endl;
  vector<future<void> > injected(injections.size());
  for(size_t i = 0; i < injections.size(); ++i){
    injected.at(i) =  async(launch::async, InjectSignal, injections.at(i), i);
  }
  for(size_t i = 0; i < injected.size(); ++i){
    injected.at(i).get();
  }

  cout << "Extracting signal strength from toys..." << endl;
  vector<vector<future<pair<double,double> > > > toyed(ntoys);
  for(size_t i = 0; i < toyed.size(); ++i){
    vector<future<pair<double,double> > > (injections.size()).swap(toyed.at(i));
  }
  vector<vector<double> > yvals_nc(injections.size(), vector<double>(ntoys));
  vector<vector<double> > yvals_c(injections.size(), vector<double>(ntoys));
  vector<vector<double> > pulls_nc(injections.size(), vector<double>(ntoys));
  vector<vector<double> > pulls_c(injections.size(), vector<double>(ntoys));
  for(int toy = 0; toy < ntoys; ++toy){
    //Right now, running one job in parallel per injected signal strength. GCC seems to default to always using deferred on lxplus when launch::async is dropped. If this is ever improved, could just do all the asyncs up front and have a completely separate loop for the gets.
    for(size_t i = 0; i < injections.size(); ++i){
      toyed.at(toy).at(i) = async(launch::async, ExtractSignal, i, toy, true);
    }
    for(size_t i = 0; i < injections.size(); ++i){
      auto res = toyed.at(toy).at(i).get();
      yvals_nc.at(i).at(toy) = res.first;
      pulls_nc.at(i).at(toy) = res.second;
    }
    for(size_t i = 0; i < injections.size(); ++i){
      toyed.at(toy).at(i) = async(launch::async, ExtractSignal, i, toy, false);
    }
    for(size_t i = 0; i < injections.size(); ++i){
      auto res = toyed.at(toy).at(i).get();
      yvals_c.at(i).at(toy) = res.first;
      pulls_c.at(i).at(toy) = res.second;
    }
  }
  execute("rm -f *_sig_inj_*.root");

  cout << "Generating plots..." << endl;
  MakePlot(injections, yvals_nc, true, false);
  MakePlot(injections, yvals_c, false, false);
  MakePlot(injections, pulls_nc, true, true);
  MakePlot(injections, pulls_c, false, true);
}

void InjectSignal(double inject, size_t index){
  {
    lock_guard<mutex> lock(global_mutex);
    cout << "Starting to inject signal with strength " << inject << endl;
  }
  ostringstream oss;
  oss << "./run/make_workspace.exe --method m1bk --lumi 2.094 --use_r4 --toys " << ntoys
      << " --sig_strength " << inject << " --identifier sig_inj_" << index << " &> /dev/null" << flush;
  {
    lock_guard<mutex> lock(global_mutex);
    cout << "Executing " << oss.str() << endl;
  }
  execute(oss.str());
  {
    lock_guard<mutex> lock(global_mutex);
    cout << "Done injecting signal with strength " << inject << endl;
  }
}

pair<double, double> ExtractSignal(size_t index, size_t toy, bool is_nc){
  {
    lock_guard<mutex> lock(global_mutex);
    cout << "Extracting signal given strength " << injections.at(index) << " in toy " << toy << endl;
  }
  ostringstream oss;
  oss << "ls -rt m1bk*" << (is_nc ? "_nc_" : "_c_") << "*sig_inj_" << index << ".root | tail -n 1" << flush;
  string file_name = execute(oss.str());
  while(file_name.back() == '\n' || file_name.back() == ' '){
    file_name.pop_back();
  }
  oss.str("");
  string workdir = MakeDir("sig_inj_");
  oss << "export origdir=$(pwd); "
      << "cd ~/cmssw/CMSSW_7_1_5/src; "
      << "eval `scramv1 runtime -sh` &> /dev/null; "
      << "cd $origdir; "
      << "cp " << file_name << ' ' << workdir << "; "
      << "cd " << workdir << "; "
      << "combine -M MaxLikelihoodFit --skipBOnlyFit --dataset data_obs_" << toy << " " << file_name << " &> /dev/null; "
      << "cd $origdir; "
      << flush;
  {
    lock_guard<mutex> lock(global_mutex);
    cout << "Executing " << oss.str() << endl;
  }
  string output = execute(oss.str());
  {
    lock_guard<mutex> lock(global_mutex);
    TFile r_file((workdir+"/mlfit.root").c_str(),"read");
    if(!r_file.IsOpen()) return pair<double, double>(-1., -1.);
    RooFitResult *f = static_cast<RooFitResult*>(r_file.Get("fit_s"));
    RooArgList pars = f->floatParsFinal();
    for(int ipar = 0; ipar < pars.getSize(); ++ipar){
      RooRealVar *var = static_cast<RooRealVar*>(pars.at(ipar));
      if(var == nullptr) continue;
      if(var->GetName() != string("r")) continue;
      double val = var->getVal();
      double delta = val - injections.at(index);
      double ehi = fabs(var->getErrorHi());
      double elo = fabs(var->getErrorLo());
      double pull = delta > 0. ? delta/elo : delta/ehi;
      cout << "Extracted strength " << val << " + " << ehi << " - " << elo << " given strength " << injections.at(index) << " in toy " << toy << ". pull=" << pull << endl;
      r_file.Close();
      execute("rm -rf "+workdir);
      return pair<double, double>(val, pull);
    }
    execute("rm -rf "+workdir);
    r_file.Close();
  }
  return pair<double, double>(-1., -1.);
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

void MakePlot(const vector<double> &injections,
	      const vector<vector<double> > &yvals,
	      bool is_nc, bool is_pull){
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
  ostringstream oss;
  oss << ";Injected " << (is_nc ? "NC" : "C") << " Signal Strength;Extracted " << (is_nc ? "NC" : "C") << "Signal Strength" << flush;
  TH1D h("", oss.str().c_str(), 1, is_pull ? xmin-margin*xdelta : themin, is_pull ? xmax+margin*xdelta : themax);
  if(is_pull){
    h.SetMinimum(-5.);
    h.SetMaximum(5.);
  }else{
    h.SetMinimum(themin);
    h.SetMaximum(themax);
  }
  h.Fill(0.5, 1.0);
  TF1 f("", "x", themin, themax);
  f.SetLineColor(2);
  f.SetLineWidth(4);
  f.SetLineStyle(2);
  TCanvas c;
  h.Draw("axis");
  if(!is_pull) f.Draw("same");
  g.Draw("p same");
  oss.str("");
  if(!is_pull){
    oss << "siginj_" << (is_nc ? "nc" : "c") << ".pdf" << flush;
  }else{
    oss << "pull_" << (is_nc ? "nc" : "c") << ".pdf" << flush;
  }
  c.Print(oss.str().c_str());

  if(!is_pull){
    cout << "Extracted signal strengths for " << (is_nc ? "NC" : "C") << ':' << endl;
  }else{
    cout << "Pulls for " << (is_nc ? "NC" : "C") << ':' << endl;
  }
  for(size_t i = 0; i < injections.size(); ++i){
    cout << injections.at(i) << ": " << centers.at(i) << " + " << ups.at(i) << " - " << downs.at(i) << ": ";
    for(const auto &y: yvals.at(i)){
      cout << y << " ";
    }
    cout << endl;
  }
}
