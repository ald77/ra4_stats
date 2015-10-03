#include "workspace_generator.hpp"

#include <iostream>
#include <sstream>
#include <array>

using namespace std;

map<YieldKey, GammaParams> WorkspaceGenerator::yields_ = map<YieldKey, GammaParams>();

WorkspaceGenerator::WorkspaceGenerator(const Cut &baseline,
				       const set<Block> &blocks,
				       const set<Process> &backgrounds,
				       const Process &signal,
				       const Process &data):
  signal_(signal),
  data_(data),
  blocks_(blocks),
  backgrounds_(backgrounds),
  baseline_(baseline){
  }

void WorkspaceGenerator::WriteToFile(const string &file_name){
  GetYields();
  cout << "Wrote workspace to file " << file_name << endl;
}

void WorkspaceGenerator::GetYields(){
  for(const auto &block: blocks_){
    for(const auto &vbin: block.Bins()){
      for(const auto &bin: vbin){
	StoreYield(bin, data_);
	StoreYield(bin, signal_);
	for(const auto &bkg: backgrounds_){
	  StoreYield(bin, bkg);
	}
      }
    }
  }
}

void WorkspaceGenerator::StoreYield(const Bin &bin, const Process &process){
  cout << "Getting yields for bin " << bin.Name()
       << ", process " << process.Name() << endl;

  GammaParams gps;
  YieldKey key(bin, process, baseline_);
  if(yields_.find(key) != yields_.end()){
    cout << "Recycling already computed yield." << endl;
    gps = yields_.at(key);
  }else if(process.GetEntries() == 0){
    cout << "No entries found." << endl;
    gps.SetNEffectiveAndWeight(0., 0.);
  }else{
    ostringstream oss;
    oss << 3.0 << flush;//QQQ Set lumi
    Cut lumi_weight = Cut(oss.str()+"*weight");

    array<Cut, 6> cuts;
    cuts.at(0) = lumi_weight*(baseline_ && bin.Cut() && process.Cut());
    cuts.at(1) = lumi_weight*(baseline_ && process.Cut());
    cuts.at(2) = lumi_weight*(process.Cut());
    cuts.at(3) = lumi_weight;
    cuts.at(4) = Cut(oss.str());
    cuts.at(5) = Cut();

    for(size_t icut = 0; icut < cuts.size() && gps.Weight()<=0.; ++icut){
      if(icut > 0 && !process.CountZeros()){
	gps.SetNEffectiveAndWeight(0., 0.);
	break;
      }
      Cut &cut = cuts.at(icut);
      cout << "Trying cut " << cut << endl;
      GammaParams temp_gps = process.GetYield(cut);
      if(icut == 0) gps = temp_gps;
      else gps.SetNEffectiveAndWeight(0., temp_gps.Weight());
    }
  }
  cout
    << "Found yield=" << gps.Yield()
    << ", uncertainty=" << gps.CorrectedUncertainty()
    << ", raw sqrt(n) uncertainty=" << gps.Uncertainty()
    << ", N_eff=" << gps.NEffective()
    << ", weight=" << gps.Weight()
    << "\n" << endl;
  yields_[key] =  gps;
}
