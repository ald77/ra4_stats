#include "dump_markov.hpp"

#include <cmath>

#include <iostream>
#include <string>
#include <map>

#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"
#include "TList.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

#include "RooStats/MarkovChain.h"

using namespace std;

int main(int argc, char *argv[]){
  for(int iarg = 1; iarg < argc; ++iarg){
    TFile file(argv[iarg], "read");
    if(!file.IsOpen()) continue;
    TDirectory *dir = static_cast<TDirectory*>(file.Get("toys"));
    if(dir == nullptr) continue;
    TList *keys = dir->GetListOfKeys();
    if(keys == nullptr) continue;
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

      // const RooArgSet *args = static_cast<const RooArgSet*>(mcmc->Get());
      // if(args == nullptr) continue;
      // TIterator *iter = args->createIterator();
      // for(; iter != nullptr && *(*iter) != nullptr; iter->Next()){
      // 	RooAbsArg *arg = static_cast<RooAbsArg*>(*(*iter));
      // 	if(arg == nullptr) continue;
      // 	arg->Print();
      // }
    }
    if(data == nullptr) continue;
    int num_entries = data->numEntries();
    bool setup_done = false;
    map<unsigned long,TH1D> hmap;
    for(int i = 0; i< num_entries; ++i){
      const RooArgSet *args  = static_cast<const RooArgSet*>(data->get(i));
      if(args == nullptr) continue;
      if(!setup_done){
	SetupHistos(hmap, *args, *data);
	setup_done = true;
      }
      TIterator *iter = args->createIterator();
      for(; iter != nullptr && *(*iter) != nullptr; iter->Next()){
	RooRealVar *var = static_cast<RooRealVar*>(*(*iter));
	if(var == nullptr) continue;
	hmap.at(var->Hash()).Fill(var->getVal(), data->weight());
      }
    }

    TCanvas canvas;
    for(auto histo = hmap.begin();
	histo != hmap.end();
	++histo){
      canvas.SetLogy(false);
      histo->second.Draw("e1p");
      canvas.Print((histo->second.GetName()+string("_lin.pdf")).c_str());
      canvas.SetLogy(true);
      histo->second.Draw("e1p");
      canvas.Print((histo->second.GetName()+string("_log.pdf")).c_str());
    }
  }
}

void SetupHistos(map<unsigned long,TH1D> &hmap,
		 const RooArgSet &args,
		 const RooDataSet &data){
  hmap.clear();
  TIterator *iter = args.createIterator();
  int nbins = sqrt(data.numEntries());
  if(nbins <= 0) nbins = 1;
  if(nbins > 100) nbins = 100;
  for(; iter != nullptr && *(*iter) != nullptr; iter->Next()){
    RooRealVar *var = static_cast<RooRealVar*>(*(*iter));
    if(var == nullptr) continue;
    if(hmap.find(var->Hash()) != hmap.end()) continue;
    double low, high;
    data.getRange(*var, low, high);
    hmap[var->Hash()] = TH1D((string("hist_")+var->GetName()).c_str(),
		      var->GetName(),
		      nbins, low, high);
    hmap[var->Hash()].Sumw2();
  }
}
