#include "get_dilepton_uncertainties.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "TFile.h"

#include "RooArgList.h"
#include "RooRealVar.h"

#include "utilities.hpp"

using namespace std;

int main(int argc, char *argv[]){
  for(int argi = 1; argi < argc; ++argi){
    execute("export blah=$(pwd); cd ~/cmssw/CMSSW_7_1_5/src; eval `scramv1 runtime -sh`; cd $blah; combine -M MaxLikelihoodFit --saveWorkspace "+string(argv[argi]));

    TFile w_file("higgsCombineTest.MaxLikelihoodFit.mH120.root","read");
    if(!w_file.IsOpen()) continue;
    RooWorkspace *w = static_cast<RooWorkspace*>(w_file.Get("w"));
    if(w == nullptr) continue;
    
    TFile fit_file("mlfit.root","read");
    if(!fit_file.IsOpen()) continue;
    RooFitResult *fit_b = static_cast<RooFitResult*>(fit_file.Get("fit_b"));
    if(fit_b == nullptr) continue;

    string outname = ChangeExtension(argv[argi], "_dilep.txt");
    ostringstream out;
    out << "SYSTEMATIC dilep_closure\n";
    out << "  PROCESS ttbar,other\n";
    vector<string> bin_names = GetBinNames(*w);
    for(const auto &bin_name: bin_names){
      double yield = GetBkgPred(*w, bin_name);
      double uncert = GetBkgPredErr(*w, *fit_b, bin_name);

      double frac = yield > 0. ? uncert/yield : 2.;
      out << "    " << bin_name << ' ' << frac << '\n';
    }
    out << flush;
    ofstream outfile(outname);
    outfile << out.str();
    cout << "Dilepton closure written to " << outname << ":\n";
    cout << out.str() << endl;
  }
}

vector<string> GetBinNames(const RooWorkspace &w){
  vector<string> names;
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  TObject *obj;
  int i = 0;
  while((obj = iter()) && i < size){
    ++i;
    RooAbsArg *arg = static_cast<RooAbsArg*>(obj);
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nexp_BLK_") continue;
    auto bpos = name.find("_BIN_");
    auto ppos = name.find("_PRC_");
    if(bpos == string::npos) continue;
    string bin_name = name.substr(bpos+5, ppos-bpos-5);
    Append(names, bin_name);
  }
  iter.Reset();
  reverse(names.begin(), names.end());
  return names;
}

double GetBkgPred(const RooWorkspace &w,
                  const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  TObject *obj;
  int i = 0;
  while((obj = iter()) && i < size){
    ++i;
    RooAbsArg *arg = static_cast<RooAbsArg*>(obj);
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nbkg_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetBkgPredErr(const RooWorkspace &w,
                     const RooFitResult &f,
                     const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  TObject *obj;
  int i = 0;
  while((obj = iter()) && i < size){
    ++i;
    RooAbsArg *arg = static_cast<RooAbsArg*>(obj);
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nbkg_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getPropagatedError(f);
  }
  iter.Reset();
  return -1.;
}
