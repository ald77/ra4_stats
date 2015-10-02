#include "dump_markov.hpp"

#include <cmath>

#include <iostream>
#include <string>
#include <map>
#include <utility>
#include <algorithm>

#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"
#include "TList.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "RooStats/MarkovChain.h"

#include "timer.hpp"

using namespace std;

int main(int argc, char *argv[]){
  for(int iarg = 1; iarg < argc; ++iarg){
    TFile file(argv[iarg], "read");
    if(!file.IsOpen()) continue;
    RooWorkspace *w = static_cast<RooWorkspace*>(file.Get("w"));
    if(w == nullptr) continue;
    RooDataSet *data = GetData(file);
    if(data == nullptr) continue;

    map<string, vector<double> > values;
    vector<double> weights;
    GetValues(*w, *data, values, weights);
    MakePlots(values, weights);
  }
}

RooDataSet * GetData(TFile &file){
  if(!file.IsOpen()) return nullptr;
  TDirectory *dir = static_cast<TDirectory*>(file.Get("toys"));
  if(dir == nullptr) return nullptr;
  TList *keys = dir->GetListOfKeys();
  if(keys == nullptr) return nullptr;
  TIter next_key(keys);
  TKey *key = nullptr;
  RooDataSet *data = nullptr;
  while((key = static_cast<TKey*>(next_key()))){
    if(key->GetClassName() != string("RooStats::MarkovChain")) continue;
    RooStats::MarkovChain *mcmc = static_cast<RooStats::MarkovChain*>(key->ReadObj());
    if(mcmc == nullptr) continue;
    RooDataSet *this_data = static_cast<RooDataSet*>(mcmc->GetAsDataSet());
    if(this_data == nullptr) continue;
    if(data == nullptr){
      data = this_data;
    }else{
      data->append(*this_data);
    }
  }
  return data;
}

void GetValues(RooWorkspace &w,
	       RooDataSet &data,
	       map<string, vector<double> > &values,
	       vector<double> &weights){
  values.clear();
  weights.clear();

  int num_entries = data.numEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(int entry = 0; entry < num_entries; ++entry){
    timer.Iterate();
    const RooArgSet *args = static_cast<const RooArgSet*>(data.get(entry));
    TIterator *iter = args->createIterator();
    for(; iter != nullptr && *(*iter) != nullptr; iter->Next()){
      RooRealVar *var = static_cast<RooRealVar*>(*(*iter));
      if(var == nullptr) continue;
      RooRealVar *wvar = w.var(var->GetName());
      if(wvar == nullptr) continue;
      wvar->setVal(var->getVal());

      FillValues(values, w);
      weights.push_back(data.weight());
    }
  }
}

void FillValues(map<string, vector<double> > &values,
                const RooWorkspace &w){
  FillValues(values, w.allVars());
  FillValues(values, w.allFunctions());
  FillValues(values, w.allPdfs());
}

void FillValues(map<string, vector<double> > &values,
                const RooArgSet &args){
  TIterator *iter = args.createIterator();
  for(; iter != nullptr && *(*iter) != nullptr; iter->Next()){
    RooRealVar *var = static_cast<RooRealVar*>(*(*iter));
    if(var == nullptr) continue;
    values[var->GetName()].push_back(var->getVal());
  }
  if(iter != nullptr){
    delete iter;
    iter = 0;
  }
}

void MakePlots(const map<string, vector<double> > &values,
	       const vector<double> &weights){
  for(auto var = values.cbegin();
      var != values.cend();
      ++var){
    MakePlot(var->first, var->second, weights);
  }
}

void MakePlot(const string &name, const vector<double> &values,
	      const vector<double> &weights){
  int num_bins = TMath::Nint(sqrt(static_cast<double>(values.size())));
  if(num_bins < 1) num_bins = 1;
  if(num_bins > 100) num_bins = 100;
  double min_val = *min_element(values.cbegin(), values.cend());
  double max_val = *max_element(values.cbegin(), values.cend());

  TH1D h(("hist_"+name).c_str(), (name+";"+name+";").c_str(),
         num_bins, min_val, max_val);
  h.Sumw2();
  size_t num_entries = min(values.size(), weights.size());
  for(size_t entry = 0; entry < num_entries; ++entry){
    h.Fill(values.at(entry), weights.at(entry));
  }

  TCanvas canvas;
  h.Draw("e1p");
  canvas.Print(("hist_"+name+"_lin.pdf").c_str());
  canvas.SetLogy();
  h.Draw("e1p");
  canvas.Print(("hist_"+name+"_log.pdf").c_str());
}
