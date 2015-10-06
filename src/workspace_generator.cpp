#include "workspace_generator.hpp"

#include <iostream>
#include <sstream>
#include <array>
#include <algorithm>
#include <numeric>

#include "RooPoisson.h"
#include "RooDataSet.h"

#include "RooStats/ModelConfig.h"

#include "yield_key.hpp"

using namespace std;

bool blinded = true;
bool do_syst = false;
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
  AddSystematicsGenerators();

  for(const auto &block: blocks_){
    AddData(block);
    AddBackgroundFractions(block);
    AddABCDParameters(block);
    AddRawBackgroundPredictions(block);
    AddFullBackgroundPredictions(block);
    AddSignalPredictions(block);
    AddPdfs(block);
  }

  AddFullPdf();
  AddParameterSets();
  AddModels();

  cout << *this << endl;
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

void WorkspaceGenerator::AddSystematicsGenerators(){
  for(const auto &block: blocks_){
    for(const auto &vbin: block.Bins()){
      for(const auto &bin: vbin){
        for(const auto &syst: bin.Systematics()){
          AddSystematicGenerator(syst.Name());
          string full_name = syst.Name()+"_BLK_"+block.Name()+"_BIN_"+bin.Name();
          ostringstream oss;
          oss << "strength_" << full_name << "[" << syst.Strength() << "]" << flush;
          w_.factory(oss.str().c_str());
          oss << "expr::" << full_name
              << "('exp(strength_" << full_name << "*" << syst.Name() << ")',"
              << "strength_" << full_name << "," << syst.Name() << ")" << flush;
          w_.factory(oss.str().c_str());
        }
      }
    }
  }
}

void WorkspaceGenerator::AddSystematicGenerator(const string &name){
  if(systematics_.find(name) != systematics_.end()) return;
  w_.factory(("RooGaussian::CONSTRAINT_"+name+"("+name+"[0.,-10.,10.],0.,1.)").c_str());
  nuisances_.insert(nuisances_.end(), name);
  systematics_.insert(systematics_.end(),name);
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
  for_each(output.begin(), output.end(),
           [scale](map<Process, double>::value_type x){x.second*=scale;});

  return output;
}

void WorkspaceGenerator::AddABCDParameters(const Block &block){
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
}

void WorkspaceGenerator::AddRawBackgroundPredictions(const Block &block){
  BlockYields by(block, backgrounds_, baseline_, yields_);
  size_t max_row = by.MaxRow();
  size_t max_col = by.MaxCol();
  for(size_t irow = 0; irow < block.Bins().size(); ++irow){
    for(size_t icol = 0; icol < block.Bins().at(irow).size(); ++icol){
      const Bin &bin = block.Bins().at(irow).at(icol);
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      vector<string> prod_list;
      for(const auto &bkg: backgrounds_){
        string prod_name = "rate_"+bb_name+"_PRC_"+bkg.Name();
        prod_list.push_back(prod_name);
        string factory_string = "prod::"+prod_name+"(rscale_BLK_"+block.Name();
        if(icol != max_col){
          ostringstream oss;
          oss << ",rx" << (icol+1) << (max_col+1) << "_BLK_" << block.Name() << flush;
          factory_string += oss.str();
        }
        if(irow != max_row){
          ostringstream oss;
          oss << ",ry" << (irow+1) << (max_row+1) << "_BLK_" << block.Name() << flush;
          factory_string += oss.str();
        }
        factory_string += (",frac_BLK_"+block.Name()+"_PRC_"+bkg.Name()+")");
        w_.factory(factory_string.c_str());
      }
      string factory_string="sum::nbkg_raw_"+bb_name+"(";
      for(auto prod = prod_list.cbegin(); prod != prod_list.cend(); ++prod){
        if(prod != prod_list.cbegin()) factory_string += ",";
        factory_string += *prod;
      }
      factory_string += ")";
      w_.factory(factory_string.c_str());
    }
  }
}

void WorkspaceGenerator::AddFullBackgroundPredictions(const Block &block){
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      ostringstream oss;
      oss << "prod::nbkg_" << bb_name << "("
          << "nbkg_raw_" << bb_name;
      for(const auto &syst: bin.Systematics()){
        oss << "," << syst.Name() << "_" << bb_name;
      }
      oss << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddSignalPredictions(const Block &block){
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      YieldKey yk(bin, signal_, baseline_);
      double yield = yields_.at(yk).Yield();
      ostringstream oss;
      oss << "rate_BLK_" << block.Name()
          << "_BIN_" << bin.Name()
          << "_PRC_" << signal_.Name()
          << "[" << yield << "]" << flush;
      w_.factory(oss.str().c_str());
      oss.str("");
      oss << "prod::nsig_BLK_" << block.Name()
          << "_BIN_" << bin.Name()
          << "(r,"
          << "rate_BLK_" << block.Name()
          << "_BIN_" << bin.Name()
          << "_PRC_" << signal_.Name()
          << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddPdfs(const Block &block){
  string null_list = "", alt_list = "";
  bool is_first = true;
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      if(!is_first){
        null_list += ",";
        alt_list += ",";
      }
      string bb_name = "_BLK_"+block.Name() +"_BIN_"+bin.Name();
      string null_name = "pdf_null"+bb_name;
      string alt_name = "pdf_alt"+bb_name;
      null_list += null_name;
      alt_list += alt_name;
      w_.factory(("sum::nexp"+bb_name+"(nbkg"+bb_name+",nsig"+bb_name+")").c_str());
      w_.factory(("RooPoisson::pdf_null"+bb_name+"(nobs"+bb_name+",nbkg"+bb_name+")").c_str());
      (static_cast<RooPoisson*>(w_.pdf(null_name.c_str())))->setNoRounding();
      w_.factory(("RooPoisson::pdf_alt"+bb_name+"(nobs"+bb_name+",nexp"+bb_name+")").c_str());
      (static_cast<RooPoisson*>(w_.pdf(alt_name.c_str())))->setNoRounding();
      is_first = false;
    }
  }
  w_.factory(("PROD:pdf_null_BLK_"+block.Name()+"("+null_list+")").c_str());
  w_.factory(("PROD:pdf_alt_BLK_"+block.Name()+"("+alt_list+")").c_str());
}

void WorkspaceGenerator::AddFullPdf(){
  if(blocks_.size() == 0){
    w_.factory("RooPoisson::model_b(0,0)");
    w_.factory("RooPoisson::model_s(0,0)");
  }else{
    string null_list = "pdf_null_BLK_"+blocks_.cbegin()->Name();
    string alt_list = "pdf_alt_BLK_"+blocks_.cbegin()->Name();
    auto block = blocks_.cbegin();
    for(++block; block != blocks_.cend(); ++block){
      null_list += (",pdf_null_BLK_"+block->Name());
      alt_list += (",pdf_alt_BLK_"+block->Name());
    }
    if(do_syst){
      for(const auto &syst: systematics_){
        null_list += (",CONSTRAINT_"+syst);
        alt_list += (",CONSTRAINT_"+syst);
      }
    }
    w_.factory(("PROD::model_b("+null_list+")").c_str());
    w_.factory(("PROD::model_s("+alt_list+")").c_str());
  }
}

void WorkspaceGenerator::AddParameterSets(){
  DefineParameterSet("nuisances", nuisances_);
  DefineParameterSet("observables", observables_);
  RooDataSet data_obs{"data_obs", "data_obs", *w_.set("observables")};
  data_obs.add(*w_.set("observables"));
}

void WorkspaceGenerator::DefineParameterSet(const string &set_name,
                                            const set<string> &var_names){
  if(var_names.size()==0){
    w_.defineSet(set_name.c_str(), "");
  }else{
    string cat_names = *var_names.cbegin();
    auto name = var_names.cbegin();
    for(++name; name != var_names.cend(); ++name){
      cat_names += ("," + *name);
    }
    w_.defineSet(set_name.c_str(), cat_names.c_str());
  }
}

void WorkspaceGenerator::AddModels(){
  RooStats::ModelConfig model_config("ModelConfig", &w_);
  model_config.SetPdf(*w_.pdf("model_s"));
  model_config.SetParametersOfInterest(*w_.set("POI"));
  model_config.SetObservables(*w_.set("observables"));
  model_config.SetNuisanceParameters(*w_.set("nuisances"));

  RooStats::ModelConfig model_config_bonly("ModelConfig_bonly", &w_);
  model_config_bonly.SetPdf(*w_.pdf("model_b"));
  model_config_bonly.SetParametersOfInterest(*w_.set("POI"));
  model_config_bonly.SetObservables(*w_.set("observables"));
  model_config_bonly.SetNuisanceParameters(*w_.set("nuisances"));

  w_.import(model_config);
  w_.import(model_config_bonly);
}

ostream & operator<<(ostream& stream, const WorkspaceGenerator &wg){
  stream << wg.signal_.Name() << endl;
  return stream;
}
