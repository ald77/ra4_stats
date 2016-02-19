#include "sig_inj.hpp"

#include <cmath>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <future>
#include <utility>
#include <mutex>
#include <limits>
#include <numeric>
#include <functional>

#include <stdlib.h>
#include <getopt.h>

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFitResult.h"

#include "utilities.hpp"
#include "styles.hpp"
#include "thread_pool.hpp"

using namespace std;

namespace{
  int ntoys = 5;
  vector<double> injections;
  double lumi = 2.246;
  bool do_asymmetric_error = false;
  bool do_systematics = false;
  bool draw_only = false;

  mutex global_mutex;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  styles style("RA4");
  style.setDefaultStyle();

  ostringstream oss;
  oss
    << "toys_" << ntoys
    << "_lumi_" << lumi
    << (do_asymmetric_error ? "_asym_error" : "")
    << (do_systematics ? "_with_syst" : "")
    << flush;
  string id_string = oss.str();

  ThreadPool thread_pool;
  if(!draw_only){
    cout << "Creating workspaces..." << endl;
    vector<future<void> > injected(injections.size());
    for(size_t i = 0; i < injections.size(); ++i){
      injected.at(i) =  thread_pool.Push(InjectSignal, id_string, injections.at(i), i);
    }
    for(size_t i = 0; i < injected.size(); ++i){
      injected.at(i).get();
    }
  }

  vector<vector<double> > yvals_nc(injections.size(), vector<double>(ntoys, -9876543210.));
  vector<vector<double> > yvals_c(injections.size(), vector<double>(ntoys, -9876543210.));
  vector<vector<double> > pulls_nc(injections.size(), vector<double>(ntoys, -9876543210.));
  vector<vector<double> > pulls_c(injections.size(), vector<double>(ntoys, -9876543210.));

  if(!draw_only){
    cout << "Extracting signal strength from toys..." << endl;
    vector<vector<future<pair<double,double> > > > toyed_nc(ntoys);
    vector<vector<future<pair<double,double> > > > toyed_c(ntoys);
    for(int toy = 0; toy < ntoys; ++toy){
      vector<future<pair<double,double> > > (injections.size()).swap(toyed_nc.at(toy));
      vector<future<pair<double,double> > > (injections.size()).swap(toyed_c.at(toy));
      for(size_t i = 0; i < injections.size(); ++i){
        toyed_nc.at(toy).at(i) = thread_pool.Push(ExtractSignal, id_string, i, toy, true);
        toyed_c.at(toy).at(i) = thread_pool.Push(ExtractSignal, id_string, i, toy, false);
      }
    }
    for(int toy = 0; toy < ntoys; ++toy){
      for(size_t i = 0; i < injections.size(); ++i){
        auto res_nc = toyed_nc.at(toy).at(i).get();
        yvals_nc.at(i).at(toy) = res_nc.first;
        pulls_nc.at(i).at(toy) = res_nc.second;
        auto res_c = toyed_c.at(toy).at(i).get();
        yvals_c.at(i).at(toy) = res_c.first;
        pulls_c.at(i).at(toy) = res_c.second;
      }
    }
    execute("rm -f *_sig_inj_*.root");
  }

  cout << "Generating plots..." << endl;
  for(size_t i = 0; i < injections.size(); ++i){
    RemoveBadResults(yvals_nc.at(i), pulls_nc.at(i));
    RemoveBadResults(yvals_c.at(i), pulls_c.at(i));
  }

  vector<double> inj_nc = injections;
  vector<double> inj_c = injections;
  
  MergeWithText(inj_nc, yvals_nc, pulls_nc, true);
  MergeWithText(inj_c, yvals_c, pulls_c, false);

  MakePlot(inj_nc, yvals_nc, true, false);
  MakePlot(inj_c, yvals_c, false, false);
  MakePlot(inj_nc, pulls_nc, true, true);
  MakePlot(inj_c, pulls_c, false, true);
}

void InjectSignal(const string id_string, double inject, size_t index){
  {
    lock_guard<mutex> lock(global_mutex);
    cout << "Starting to inject signal with strength " << inject << endl;
  }
  ostringstream oss;
  oss << "./run/make_workspace.exe --method m1bk"
      << (do_systematics ? "" : " --no_syst") << " --lumi " << lumi << " --use_r4 --toys " << ntoys
      << " --sig_strength " << inject << " --identifier sig_inj_" << id_string << "_" << index
      << " < /dev/null &> /dev/null" << flush;
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

pair<double, double> ExtractSignal(const string id_string, size_t index, size_t toy, bool is_nc){
  {
    lock_guard<mutex> lock(global_mutex);
    cout << "Extracting signal given strength " << injections.at(index) << " in toy " << toy << endl;
  }
  ostringstream oss;
  oss << "ls -rt m1bk*" << (is_nc ? "_nc_" : "_c_")
      << "*sig_inj_" << id_string << "_" << index << ".root | tail -n 1" << flush;
  string file_name = execute(oss.str());
  while(file_name.back() == '\n' || file_name.back() == ' '){
    file_name.pop_back();
  }
  oss.str("");
  string workdir = MakeDir("sig_inj_");
  {
    lock_guard<mutex> lock(global_mutex);
    oss << "export origdir=$(pwd); "
	<< "cd ~/cmssw/CMSSW_7_1_5/src; "
	<< "eval `scramv1 runtime -sh` &> /dev/null; "
	<< "cd $origdir; "
	<< "cp " << file_name << ' ' << workdir << "; "
	<< "cd " << workdir << "; "
	<< "combine -M MaxLikelihoodFit --skipBOnlyFit --dataset data_obs_" << toy << " --preFitValue " << max(injections.at(index),0.01) << ' ' << file_name << " &> /dev/null; "
	<< "cd $origdir; "
	<< flush;
    cout << "Executing " << oss.str() << endl;
  }
  string output = execute(oss.str());
  {
    lock_guard<mutex> lock(global_mutex);
    TFile r_file((workdir+"/mlfit.root").c_str(),"read");
    if(!r_file.IsOpen()) return pair<double, double>(-1., -1.);
    RooFitResult *f = static_cast<RooFitResult*>(r_file.Get("fit_s"));
    if(f == nullptr) return pair<double, double>(-1., -1.);
    RooArgList pars = f->floatParsFinal();
    for(int ipar = 0; ipar < pars.getSize(); ++ipar){
      RooRealVar *var = static_cast<RooRealVar*>(pars.at(ipar));
      if(var == nullptr) continue;
      if(var->GetName() != string("r")) continue;
      double val = var->getVal();
      double delta = val - injections.at(index);
      double ehi, elo;
      if(do_asymmetric_error){
	GetAsymError(*var, ehi, elo);
      }else{
	ehi = GetError(*var, *f);
	elo = ehi;
      }
      double pull = 0.;
      if(delta < 0.){
	pull = -fabs(delta)/fabs(ehi);
      }else{
	pull = fabs(delta)/fabs(elo);
      }
      cout << "Extracted strength " << val << " + " << ehi << " - " << fabs(elo) << " given strength " << injections.at(index) << " in toy " << toy << ". pull=" << pull << endl;
      r_file.Close();
      execute("rm -rf "+workdir);
      return pair<double, double>(val, pull);
    }
    execute("rm -rf "+workdir);
    r_file.Close();
  }
  return pair<double, double>(-1., -1.);
}

double GetAsymError(const RooRealVar& var, double &ehi, double &elo){
  if(var.hasAsymError()){
    ehi = var.getAsymErrorHi();
    elo = var.getAsymErrorLo();
    double aehi = fabs(ehi);
    double aelo = fabs(elo);
    double delta = aehi - aelo;
    return aelo + 0.5*delta;
  }else if(var.hasError()){
    ehi = var.getError();
    elo = ehi;
    return ehi;
  }else{
    ehi = 0.;
    elo = 0.;
    return 0.;
  }
}

double GetError(const RooAbsReal &var,
                const RooFitResult &f){
  // Clone self for internal use
  RooAbsReal* cloneFunc = static_cast<RooAbsReal*>(var.cloneTree());
  RooArgSet* errorParams = cloneFunc->getObservables(f.floatParsFinal());
  RooArgSet* nset = cloneFunc->getParameters(*errorParams);

  // Make list of parameter instances of cloneFunc in order of error matrix
  RooArgList paramList;
  const RooArgList& fpf = f.floatParsFinal();
  vector<int> fpf_idx;
  for (int i=0; i<fpf.getSize(); i++) {
    RooAbsArg* par = errorParams->find(fpf[i].GetName());
    if (par) {
      paramList.add(*par);
      fpf_idx.push_back(i);
    }
  }

  vector<double> errors(paramList.getSize());
  for (Int_t ivar=0; ivar<paramList.getSize(); ivar++) {
    RooRealVar& rrv = static_cast<RooRealVar&>(fpf[fpf_idx[ivar]]);

    double cenVal = rrv.getVal();
    double errVal = rrv.getError();

    // Make Plus variation
    static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal+0.5*errVal);
    double up = cloneFunc->getVal(nset);

    // Make Minus variation
    static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal-0.5*errVal);
    double down = cloneFunc->getVal(nset);

    errors.at(ivar) = (up-down);

    static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal);
  }

  vector<double> right(errors.size());
  for(size_t i = 0; i < right.size(); ++i){
    right.at(i) = 0.;
    for(size_t j = 0; j < errors.size(); ++j){
      right.at(i) += f.correlation(paramList.at(i)->GetName(),paramList.at(j)->GetName())*errors.at(j);
    }
  }
  double sum = 0.;
  for(size_t i = 0; i < right.size(); ++i){
    sum += errors.at(i)*right.at(i);
  }

  if(cloneFunc != nullptr){
    delete cloneFunc;
    cloneFunc = nullptr;
  }
  if(errorParams != nullptr){
    delete errorParams;
    errorParams = nullptr;
  }
  if(nset != nullptr){
    delete nset;
    nset = nullptr;
  }

  return sqrt(sum);
}

void GetStats(const vector<double> &vals, double &mean, double &median,
              double &up, double &down,
              double &up2, double &down2){
  double tail = erfc(1./sqrt(2.))*0.5; //0.159 = (1-0.683)/2
  double tail2 = erfc(sqrt(2.))*0.5; //0.023 = (1-0.954)/2
  median = GetValue(vals, 0.5);
  mean = vals.size() > 0 ? accumulate(vals.cbegin(), vals.cend(), 0.)/vals.size() : 0.;
  double low = GetValue(vals, tail);
  double high = GetValue(vals, 1.-tail);
  double low2 = GetValue(vals, tail2);
  double high2 = GetValue(vals, 1.-tail2);
  up = max(high-median,0.);
  down = max(median-low,0.);
  up2 = max(high2-median,0.);
  down2 = max(median-low2,0.);
}

double GetValue(vector<double> vals, double fraction){
  if(vals.size() == 0) return 0.;
  if(vals.size() == 1) return vals.front();
  sort(vals.begin(), vals.end());
  double frac_index = fraction*vals.size()-0.5;
  long lo_index = static_cast<long>(floor(frac_index));
  long hi_index = static_cast<long>(ceil(frac_index));
  if(lo_index < 0){
    lo_index = 0;
    hi_index = 1;
  }
  if(static_cast<size_t>(hi_index) >= vals.size()){
    lo_index = vals.size() - 2;
    hi_index = vals.size() - 1;
  }
  double lo_value = vals.at(lo_index);
  double hi_value = vals.at(hi_index);
  return lo_value+(hi_value-lo_value)*(frac_index - lo_index);
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
      {"lumi", required_argument, 0, 'l'},
      {"asym", no_argument, 0, 'a'},
      {"syst", no_argument, 0, 's'},
      {"draw", no_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "t:i:l:asd", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 't':
      ntoys = atoi(optarg);
      break;
    case 'i':
      injection_set.insert(atof(optarg));
      break;
    case 'l':
      lumi = atof(optarg);
      break;
    case 'a':
      do_asymmetric_error = true;
      break;
    case 's':
      do_systematics = true;
      break;
    case 'd':
      draw_only = true;
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

void MakePlot(const vector<double> &injections_list,
	      const vector<vector<double> > &yvals,
	      bool is_nc, bool is_pull){
  vector<double> zeros(injections_list.size(), 0.);
  vector<double> means(injections_list.size()), medians(injections_list.size());
  vector<double> ups(injections_list.size()), downs(injections_list.size());
  vector<double> ups2(injections_list.size()), downs2(injections_list.size());
  vector<double> tops(injections_list.size()), bots(injections_list.size());
  vector<double> tops2(injections_list.size()), bots2(injections_list.size());
  for(size_t i = 0; i < injections_list.size(); ++i){
    GetStats(yvals.at(i), means.at(i), medians.at(i),
             ups.at(i), downs.at(i),
             ups2.at(i), downs2.at(i));
    tops.at(i) = medians.at(i) + ups.at(i);
    bots.at(i) = medians.at(i) - downs.at(i);
    tops2.at(i) = medians.at(i) + ups2.at(i);
    bots2.at(i) = medians.at(i) - downs2.at(i);
  }

  TGraphAsymmErrors g(injections_list.size(), &injections_list.at(0), &medians.at(0),
		      &zeros.at(0), &zeros.at(0),
		      &downs.at(0), &ups.at(0));
  TGraphAsymmErrors g2(injections_list.size(), &injections_list.at(0), &medians.at(0),
                       &zeros.at(0), &zeros.at(0),
                       &downs2.at(0), &ups2.at(0));
  TGraph gm(injections_list.size(), &injections_list.at(0), &means.at(0));
  g.SetMarkerStyle(20);
  g.SetMarkerSize(2);
  g.SetMarkerColor(kRed+2);
  g.SetLineStyle(1);
  g.SetLineWidth(5);
  g.SetLineColor(kRed+2);
  g2.SetMarkerStyle(0);
  g2.SetMarkerSize(0);
  g2.SetMarkerColor(0);
  g2.SetLineStyle(2);
  g2.SetLineWidth(3);
  g2.SetLineColor(kRed+2);
  gm.SetMarkerStyle(20);
  gm.SetMarkerSize(2);
  gm.SetMarkerColor(kBlue);
  double xmax = *max_element(injections_list.cbegin(), injections_list.cend());
  double xright = 1.0625*xmax;
  double ymax = *max_element(tops.cbegin(), tops.cend());
  ostringstream oss;
  oss << ";Injected " << (is_nc ? "NC" : "C") << " Signal Strength;";
  if(is_pull){
    oss << "Pull";
  }else{
    oss << "Extracted " << (is_nc ? "NC" : "C") << " Signal Strength";
  }
  oss << flush;
  TH1D h("", oss.str().c_str(), 1, 0., xright);
  if(is_pull){
    h.SetMinimum(-3.);
    h.SetMaximum(3.);
  }else{
    h.SetMinimum(0.);
    h.SetMaximum(ymax);
  }
  h.Fill(0.5, 1.0); 
  TF1 f("", "x", 0., xright);
  f.SetLineColor(kBlack);
  f.SetLineWidth(1);
  f.SetLineStyle(3);
  TCanvas c;
  h.Draw("axis");
  if(!is_pull){
    f.Draw("same");
  }
  g.Draw("p 0 same");
  g2.Draw("p 0 same");
  gm.Draw("p same");
  TLine hline;
  hline.SetLineColor(kBlack);
  hline.SetLineWidth(1);
  hline.SetLineStyle(3);
  if(is_pull){
    hline.DrawLine(0, -2., xright, -2.);
    hline.DrawLine(0, -1., xright, -1.);
    hline.DrawLine(0, 0., xright, 0.);
    hline.DrawLine(0, 1., xright, 1.);
    hline.DrawLine(0, 2., xright, 2.);
  }
  oss.str("");
  oss
    << "siginj"
    << (is_nc ? "_nc" : "_c")
    << "_toys_" << ntoys
    << "_lumi_" << lumi
    << (do_asymmetric_error ? "_asym_error" : "")
    << (do_systematics ? "_with_syst" : "")
    << ".pdf"
    << flush;
  string plot_name = oss.str();
  if(is_pull){
    ReplaceAll(plot_name, "siginj", "pull");
  }
  c.Print(plot_name.c_str());

  for(size_t i = 0; i < injections_list.size(); ++i){
    Plot1D(injections_list.at(i), yvals.at(i),
	   means.at(i), medians.at(i),
           ups.at(i), downs.at(i),
           ups2.at(i), downs2.at(i),
           is_nc, is_pull);
  }
}

void Plot1D(double inj, const vector<double> &vals,
	    double mean, double median,
            double up, double down,
            double up2, double down2,
            bool is_nc, bool is_pull){
  if(vals.size() == 0) return;
  ostringstream oss;
  oss
    << (is_pull ? "pull" : "siginj")
    << "_" << (is_nc ? "nc" : "c")
    << "_lumi_" << lumi
    << "_inj_" << inj
    << ".pdf"
    << flush;
  string name = oss.str();
  oss.str("");
  oss << "Inj. " << (is_nc ? "NC" : "C") << " Sig.=" << inj << ", N_{toys}=" << vals.size() << ";";
  if(is_pull){
    oss << "Pull";
  }else{
    if(is_nc){
      oss << "Extracted NC Signal Strength";
    }else{
      oss << "Extracted C Signal Strength";
    }
  }
  oss << ";# of Toys" << flush;
    
  TH1D h("", oss.str().c_str(), floor(sqrt(vals.size())),
	 *min_element(vals.cbegin(), vals.cend())-0.001,
	 *max_element(vals.cbegin(), vals.cend())+0.001);
  h.SetLineWidth(5);
  h.SetLineColor(kBlack);
  h.Sumw2();
  for(const auto &val: vals){
    h.Fill(val);
  }
  double hmax = 0.;
  for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
    double content = h.GetBinContent(bin)+h.GetBinError(bin);
    if(content > hmax) hmax = content;
  }
  h.SetMinimum(0.);
  h.SetMaximum(hmax);
  TLine ldown2(median-down2, 0., median-down2, hmax);
  ldown2.SetLineColor(kRed+2);
  ldown2.SetLineStyle(3);
  ldown2.SetLineWidth(2);
  TLine ldown(median-down, 0., median-down, hmax);
  ldown.SetLineColor(kRed+2);
  ldown.SetLineStyle(2);
  ldown.SetLineWidth(4);
  TLine lmedian(median, 0., median, hmax);
  lmedian.SetLineColor(kRed+2);
  lmedian.SetLineWidth(5);
  TLine lmean(mean, 0., mean, hmax);
  lmean.SetLineColor(kBlue);
  lmean.SetLineWidth(5);
  TLine lup(median+up, 0., median+up, hmax);
  lup.SetLineColor(kRed+2);
  lup.SetLineStyle(2);
  lup.SetLineWidth(4);
  TLine lup2(median+up2, 0., median+up2, hmax);
  lup2.SetLineColor(kRed+2);
  lup2.SetLineStyle(3);
  lup2.SetLineWidth(2);
  TCanvas c;
  h.Draw("e1p");
  ldown2.Draw("same");
  ldown.Draw("same");
  lmedian.Draw("same");
  lmean.Draw("same");
  lup.Draw("same");
  lup2.Draw("same");
  c.Print(name.c_str());
}

void RemoveBadResults(vector<double> &vals, vector<double> &pulls){
  if(vals.size() != pulls.size()) ERROR("Vals and pull must have same length");
  vector<size_t> bad_indices;
  for(size_t i = 0; i < vals.size(); ++i){
    if(vals.at(i) < 0.
       || fabs(pulls.at(i)) >= 1000.){
      bad_indices.push_back(i);
    }
  }
  sort(bad_indices.begin(), bad_indices.end(), greater<size_t>());
  for(size_t ii = 0; ii < bad_indices.size(); ++ii){
    size_t i = bad_indices.at(ii);
    vals.erase(vals.begin()+i);
    pulls.erase(pulls.begin()+i);
  }
}

void MergeWithText(vector<double> &injs,
                   vector<vector<double> > &yvals,
                   vector<vector<double> > &pulls,
                   bool is_nc){
  ostringstream oss;
  oss << "txt/siginj/lumi_" << lumi << "_" << (is_nc ? "nc" : "c") << ".txt" << flush;
  string path = oss.str();

  string line;
  ifstream infile(path);
  while(getline(infile, line)){
    istringstream iss(line);
    double inj, y, pull;
    iss >> inj >> y >> pull;
    size_t i = GetIndex(injs, inj);
    if(i == static_cast<size_t>(-1)){
      injs.push_back(inj);
      yvals.push_back(vector<double>());
      pulls.push_back(vector<double>());
      i = injs.size() - 1;
    }
    yvals.at(i).push_back(y);
    pulls.at(i).push_back(pull);
  }
  SortByInjectionStrength(injs, yvals, pulls);

  ofstream outfile(path);
  outfile.precision(std::numeric_limits<double>::max_digits10);
  for(size_t iinj = 0; iinj < injs.size(); ++iinj){
    for(size_t itoy = 0; itoy < yvals.at(iinj).size() || itoy < pulls.at(iinj).size(); ++itoy){
      outfile
        << setw(32) << injs.at(iinj)
        << ' ' << setw(32) << (itoy < yvals.at(iinj).size() ? yvals.at(iinj).at(itoy) : -9876543210.)
        << ' ' << setw(32) << (itoy < pulls.at(iinj).size() ? pulls.at(iinj).at(itoy) : -9876543210.)
        << endl;
    }
  }
  outfile.close();
}

size_t GetIndex(const vector<double> &v, double x){
  for(size_t i = 0; i < v.size(); ++i){
    if(v.at(i) == x) return i;
  }
  return -1;
}

void SortByInjectionStrength(vector<double> &inj,
                             vector<vector<double> > &yvals,
                             vector<vector<double> > &pulls){
  vector<size_t> vi(inj.size());
  for(size_t i = 0; i < vi.size(); ++i) vi.at(i) = i;
  sort(vi.begin(), vi.end(), [&](size_t a, size_t b){return inj.at(a)<inj.at(b);});
  auto inj_temp = inj;
  auto yvals_temp = yvals;
  auto pulls_temp = pulls;
  for(size_t i = 0; i < vi.size(); ++i){
    inj.at(i) = inj_temp.at(vi.at(i));
    yvals.at(i) = yvals_temp.at(vi.at(i));
    pulls.at(i) = pulls_temp.at(vi.at(i));
    sort(yvals.at(i).begin(), yvals.at(i).end());
    sort(pulls.at(i).begin(), pulls.at(i).end());
  }
}
