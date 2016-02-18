#include "plot_likelihood.hpp"

#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <limits>

#include <stdlib.h>
#include <getopt.h>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"

#include "utilities.hpp"
#include "styles.hpp"

using namespace std;

namespace{
  string file_name = "";
  string workspace_name = "w";
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  if(file_name == ""){
    ERROR("Specify name of file containing workspace with -f.");
  }
  execute("export blah=$(pwd); cd ~/cmssw/CMSSW_7_1_5/src; eval `scramv1 runtime -sh`; cd $blah; combine -M MaxLikelihoodFit --saveWorkspace --saveWithUncertainties --minos=all -w "+workspace_name+" "+file_name);

  styles style("RA4");
  style.setDefaultStyle();

  string fit_workspace_file_name("higgsCombineTest.MaxLikelihoodFit.mH120.root");
  TFile fit_workspace_file(fit_workspace_file_name.c_str(),"read");
  if(!fit_workspace_file.IsOpen()){
    ERROR("File "+fit_workspace_file_name+" was not produced.");
  }
  RooWorkspace *w = static_cast<RooWorkspace*>(fit_workspace_file.Get(workspace_name.c_str()));
  if(w == nullptr) {
    ERROR("Workspace "+workspace_name+" not found in file "+fit_workspace_file_name);
  }
    
  string fit_file_name("mlfit.root");
  TFile fit_file(fit_file_name.c_str(), "read");
  if(!fit_file.IsOpen()){
    ERROR("File "+fit_file_name+" was not produced.");
  }

  string out_file_name = ChangeExtension(file_name, "_profiles.root");
  TFile out_file(out_file_name.c_str(), "recreate");
  if(!out_file.IsOpen()){
    ERROR("Could not open output file "+out_file_name);
  }
  out_file.cd();

  RooFitResult *fit_b = static_cast<RooFitResult*>(fit_file.Get("fit_b"));
  if(fit_b != nullptr){
    SetVariables(*w, *fit_b);
    PlotVars(*w, true);
  }
  RooFitResult *fit_s = static_cast<RooFitResult*>(fit_file.Get("fit_s"));
  if(fit_s != nullptr){
    SetVariables(*w, *fit_s);
    PlotVars(*w, false);
  }

  out_file.Close();
  fit_file.Close();
  fit_workspace_file.Close();
}

void PlotVars(RooWorkspace &w, bool bkg_only){
  RooAbsPdf *model = w.pdf(bkg_only ? "model_b" : "model_s");
  if(model == nullptr) return;
  RooDataSet *data = static_cast<RooDataSet*>(w.data("data_obs"));
  if(data ==nullptr) return;
  RooAbsReal *nll = model->createNLL(*data);
  if(nll == nullptr) return;
  RooMinuit m(*nll);
  m.migrad();
  TIter iter(w.allVars().createIterator());
  int size = w.allVars().getSize();
  RooRealVar *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooRealVar*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    if(arg->isConstant()) continue;
    string name = arg->GetName();
    double vmin = arg->getMin();
    double vmax = arg->getMax();
    double vval = arg->getVal();
    double vehi = arg->getErrorHi();
    double velo = arg->getErrorLo();
    double verr = arg->getError();
    cout << name << ": " << vmin << " " << (vval+velo) << " " << vval << " (" << verr << ") " << (vval+vehi) << " " << vmax << endl;
    double low = max(vmin, vval+5*velo);
    double high = min(vmax, vval+5*vehi);
    TH1D h("", (";"+name+";-2 log(L)").c_str(), 11, low, high);
    double minval=numeric_limits<double>::max();
    for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
      arg->setVal(h.GetBinCenter(bin));
      arg->setConstant(true);
      m.migrad();
      RooFitResult *f = m.save();
      if(f == nullptr) continue;
      double val = f->minNll();
      h.SetBinContent(bin, val);
      if(val<minval) minval = val;
    }
    for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
      h.SetBinContent(bin, 2.*(h.GetBinContent(bin)-minval));
    }
    arg->setVal(vval);
    arg->setConstant(false);
    TCanvas c;
    h.Draw();
    c.Print((name+".pdf").c_str());
  }
  iter.Reset();
}

vector<string> GetVarNames(const RooWorkspace &w){
  vector<string> names;
  TIter iter(w.allVars().createIterator());
  int size = w.allVars().getSize();
  TObject *obj = nullptr;
  int i = 0;
  while((obj = iter()) && i < size){
    ++i;
    if(obj == nullptr) continue;
    string name = obj->GetName();
    Append(names, name);
  }
  iter.Reset();
  sort(names.begin(), names.end());
  return names;
}

RooRealVar * SetVariables(RooWorkspace &w,
                          const RooFitResult &f){
  bool set_r = false;
  RooArgList pars = f.floatParsFinal();
  for(int ipar = 0; ipar < pars.getSize(); ++ipar){
    RooRealVar *fit_var = static_cast<RooRealVar*>(pars.at(ipar));
    if(fit_var == nullptr) continue;
    RooRealVar *w_var = w.var(fit_var->GetName());
    if(w_var == nullptr) continue;
    w_var->setMin(fit_var->getMin());
    w_var->setMax(fit_var->getMax());
    w_var->setVal(fit_var->getVal());
    w_var->setError(fit_var->getError());
    if(fit_var->GetName() == string("r")) set_r = true;
  }

  RooRealVar *r_var = static_cast<RooRealVar*>(w.var("r"));
  if(r_var != nullptr){
    if(!set_r){
      r_var->setVal(0);
      r_var->setConstant(true);
    }else{
      r_var->setConstant(false);
    }
  }
  return r_var;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"file", required_argument, 0, 'f'},
      {"workspace", required_argument, 0, 'w'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:w:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'w':
      workspace_name = optarg;
      break;
    case 'f':
      file_name = optarg;
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
