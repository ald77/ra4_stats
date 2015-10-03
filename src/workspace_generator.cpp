#include "workspace_generator.hpp"

#include <iostream>
#include <sstream>
#include <array>
#include <algorithm>
#include <numeric>

#include "yield_key.hpp"

using namespace std;

bool blinded = true;
double lumi = 3.;

map<YieldKey, GammaParams> WorkspaceGenerator::yields_ = map<YieldKey, GammaParams>();

WorkspaceGenerator::WorkspaceGenerator(const Cut &baseline,
				       const set<Block> &blocks,
				       const set<Process> &backgrounds,
				       const Process &signal,
				       const Process &data):
  baseline_(baseline),
  backgrounds_(backgrounds),
  signal_(signal),
  data_(data),
  blocks_(blocks),
  w_("w"),
  poi_(),
  observables_(),
  nuisances_(),
  systematics_(){
  w_.cd();
  }

void WorkspaceGenerator::WriteToFile(const string &file_name){
  GetYields();
  AddPOI();

  for(const auto block: blocks_){
    AddData(block);
    AddBackgroundFractions(block);
    BlockYields by = AddABCDParameters(block);
    AddBackgroundPredictions(block, by);
  }

  cout << "Wrote workspace to file " << file_name << endl;
}

void WorkspaceGenerator::GetYields() const{
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

void WorkspaceGenerator::StoreYield(const Bin &bin, const Process &process) const{
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
    oss << lumi << flush;
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

void WorkspaceGenerator::AddPOI(){
  w_.factory("r[1.,0.,20.]");
  poi_.insert(poi_.end(), "r");
}

void WorkspaceGenerator::AddData(const Block &block){
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      GammaParams gps(0., 0.);
      if(!blinded){
	YieldKey key(bin, data_, baseline_);
	gps = yields_.at(key);
      }else{
	for(const auto &bkg: backgrounds_){
	  YieldKey key(bin, bkg, baseline_);
	  gps += yields_.at(key);
	}
      }
      ostringstream oss;
      oss << "nobs_BLK_" << block.Name()
          << "_BIN_" << bin.Name() << flush;
      observables_.insert(observables_.end(), oss.str());
      oss << "[" << gps.Yield() << "]" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddBackgroundFractions(const Block &block){
  ostringstream oss;
  if(backgrounds_.size()>1){
    map<Process, double> bkg_fracs = GetBackgroundFractions(block);
    set<string> names;
    auto bkg = backgrounds_.cbegin();
    for(++bkg; bkg != backgrounds_.cend(); ++bkg){
      oss.str("");
      oss << "frac_BLK_" << block.Name()
          << "_PRC_" << bkg->Name()
          << flush;
      nuisances_.insert(nuisances_.end(), oss.str());
      names.insert(names.end(), oss.str());
      oss << "[" << bkg_fracs.at(*bkg) << ",0.,1.]" << flush;
      w_.factory(oss.str().c_str());
    }
    oss.str("");
    oss << "expr::frac_BLK_" << block.Name() << "_PRC_"
	<< backgrounds_.cend()->Name()
	<< "('";
    for(const auto &name: names) oss << "-" << name;
    oss << ")" << flush;
    w_.factory(oss.str().c_str());
  }else{
    oss << "frac_BLK_" << block.Name() << "_PRC_"
        << backgrounds_.begin()->Name()
        << "[1]" << flush;
    w_.factory(oss.str().c_str());
  }
}

map<Process, double> WorkspaceGenerator::GetBackgroundFractions(const Block &block) const{
  map<Process, double> output;
  for(const auto &bkg: backgrounds_){
    for(const auto &vbin: block.Bins()){
      for(const auto &bin: vbin){
	YieldKey key(bin, bkg, baseline_);
	output[bkg] += yields_.at(key).Yield();
      }
    }
  }
  double scale = accumulate(output.cbegin(), output.cend(), 0.,
			    [](double result, map<Process, double>::value_type x){return result+x.second;});
  scale = 1./scale;
  for_each(output.begin(), output.end(), [scale](map<Process, double>::value_type x){x.second*=scale;});

  return output;
}

BlockYields WorkspaceGenerator::AddABCDParameters(const Block &block){
  BlockYields by(block, backgrounds_, baseline_, yields_);

  ostringstream rxss, ryss;
  rxss << "sum::rxnorm_BLK_" << block.Name() << "(1.,";
  ryss << "sum::rynorm_BLK_" << block.Name() << "(1.,";
  ostringstream oss;
  oss << "norm_BLK_" << block.Name() << flush;
  nuisances_.insert(nuisances_.end(), oss.str());
  oss << "[" << max(1., by.Total().Yield()) << ",0.,"
      << max(5.*by.Total().Yield(), 20.) << "]" << flush;
  w_.factory(oss.str().c_str());

  for(size_t irow = 0; irow < by.RowSums().size(); ++irow){
    if(irow == by.MaxRow()) continue;
    oss.str("");
    oss << "ry" << (irow+1) << (by.MaxRow()+1) << "_BLK_" << block.Name() << flush;
    ryss << "," << oss.str();
    nuisances_.insert(nuisances_.end(), oss.str());
    oss << "[" << by.RowSums().at(irow).Yield()/by.RowSums().at(by.MaxRow()).Yield()
        << ",0.,10.]" << flush;
    w_.factory(oss.str().c_str());
  }
  ryss << ")" << flush;
  w_.factory(ryss.str().c_str());
  for(size_t icol = 0; icol < by.ColSums().size(); ++icol){
    if(icol == by.MaxCol()) continue;
    oss.str("");
    oss << "rx" << (icol+1) << (by.MaxCol()+1) << "_BLK_" << block.Name() << flush;
    rxss << "," << oss.str();
    nuisances_.insert(nuisances_.end(), oss.str());
    oss << "[" << by.ColSums().at(icol).Yield()/by.ColSums().at(by.MaxCol()).Yield()
        << ",0.,10.]" << flush;
    w_.factory(oss.str().c_str());
  }
  rxss << ")" << flush;
  w_.factory(rxss.str().c_str());
  oss.str("");
  oss << "prod::rnorm_BLK_" << block.Name()
      << "(rxnorm_BLK_" << block.Name()
      << ",rynorm_BLK_" << block.Name() << ")" << flush;
  w_.factory(oss.str().c_str());
  oss.str("");
  oss << "expr::rscale_BLK_" << block.Name()
      << "('norm_BLK_" << block.Name() << "/rnorm_BLK_" << block.Name()
      << "',norm_BLK_" << block.Name()
      << ",rnorm_BLK_" << block.Name() << ")" << flush;
  w_.factory(oss.str().c_str());

  return by;
}

void AddBackgroundPredictions(const Block &block, const BlockYields &by){
  cout << block.Bins().size() << by.MaxRow();
}
