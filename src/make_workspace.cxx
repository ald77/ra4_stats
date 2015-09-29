#include "make_workspace.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <initializer_list>
#include <map>
#include <array>
#include <vector>
#include <string>

#include <unistd.h>
#include <getopt.h>

#include "TChain.h"
#include "TH1D.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPoisson.h"

#include "RooStats/ModelConfig.h"

#include "gamma_params.hpp"

using namespace std;
using namespace RooStats;

namespace{
  double lumi = 10.;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  //Define processes. Try to minimize splitting
  Process ttbar{"ttbar", {
      {"archive/2015_09_28_ana/skim/*TTJets*Lept*.root/tree"}
    }};
  Process other{"other", {
      {"archive/2015_09_28_ana/skim/*DYJetsToLL*.root/tree"},
        {"archive/2015_09_28_ana/skim/*QCD_Pt*.root/tree"},
	  {"archive/2015_09_28_ana/skim/*_ST_*.root/tree"},
	    {"archive/2015_09_28_ana/skim/*WJetsToLNu*.root/tree"},
	      {"archive/2015_09_28_ana/skim/*_WWTo*.root/tree"},
		{"archive/2015_09_28_ana/skim/*ggZH_HToBB*.root/tree"},
		  {"archive/2015_09_28_ana/skim/*ttHJetTobb*.root/tree"}
    }};
  Process signal{"signal", {
      {"archive/2015_09_28_ana/skim/*T1tttt*1500*100*.root/tree"}
    }};
  Process data{"data", {}};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  vector<reference_wrapper<Process> > backgrounds{ref(ttbar), ref(other)};

  //Baseline selection applied to all bins and processes
  string baseline{"ht>500&&met>200&njets>=7&&nbm>=2&&(nels+nmus)==1"};

  //Declare bins
  //Method 3
  Bin r1_lowmet_lownb{"r1_lowmet_lownb", "mt<=140&&mj<=400&&met<=400&&nbm<=2"};
  Bin r1_lowmet_highnb{"r1_lowmet_highnb", "mt<=140&&mj<=400&&met<=400&&nbm>2"};
  Bin r1_highmet{"r1_highmet", "mt<=140&&mj<=400&&met>400"};

  Bin r2_lowmet_lownj_lownb{"r2_lowmet_lownj_lownb", "mt<=140&&mj>400&&met<=400&&njets<=8&&nbm<=2"};
  Bin r2_lowmet_lownj_highnb{"r2_lowmet_lownj_highnb", "mt<=140&&mj>400&&met<=400&&njets<=8&&nbm>2"};
  Bin r2_lowmet_highnj_lownb{"r2_lowmet_highnj_lownb", "mt<=140&&mj>400&&met<=400&&njets>8&&nbm<=2"};
  Bin r2_lowmet_highnj_highnb{"r2_lowmet_highnj_highnb", "mt<=140&&mj>400&&met<=400&&njets>8&&nbm>2"};
  Bin r2_highmet_lownj{"r2_highmet_lownj", "mt<=140&&mj>400&&met>400&&njets<=8"};
  Bin r2_highmet_highnj{"r2_highmet_highnj", "mt<=140&&mj>400&&met>400&&njets>8"};

  Bin r3_lowmet_lownb{"r3_lowmet_lownb", "mt>140&&mj<=400&&met<=400&&nbm<=2"};
  Bin r3_lowmet_highnb{"r3_lowmet_highnb", "mt>140&&mj<=400&&met<=400&&nbm>2"};
  Bin r3_highmet{"r3_highmet", "mt>140&&mj<=400&&met>400"};

  Bin r4_lowmet_lownj_lownb{"r4_lowmet_lownj_lownb", "mt>140&&mj>400&&met<=400&&njets<=8&&nbm<=2"};
  Bin r4_lowmet_lownj_highnb{"r4_lowmet_lownj_highnb", "mt>140&&mj>400&&met<=400&&njets<=8&&nbm>2"};
  Bin r4_lowmet_highnj_lownb{"r4_lowmet_highnj_lownb", "mt>140&&mj>400&&met<=400&&njets>8&&nbm<=2"};
  Bin r4_lowmet_highnj_highnb{"r4_lowmet_highnj_highnb", "mt>140&&mj>400&&met<=400&&njets>8&&nbm>2"};
  Bin r4_highmet_lownj{"r4_highmet_lownj", "mt>140&&mj>400&&met>400&&njets<=8"};
  Bin r4_highmet_highnj{"r4_highmet_highnj", "mt>140&&mj>400&&met>400&&njets>8"};

  //Method 1
  Bin m1_r1_lowmet_lownj{"m1_r1_lowmet_lownj", "mt<=140&&mj<=600&&met<=400&&njets<=8"};
  Bin m1_r1_lowmet_highnj{"m1_r1_lowmet_highnj", "mt<=140&&mj<=600&&met<=400&&njets>8"};
  Bin m1_r1_highmet_lownj{"m1_r1_highmet_lownj", "mt<=140&&mj<=600&&met>400&&njets<=8"};
  Bin m1_r1_highmet_highnj{"m1_r1_highmet_highnj", "mt<=140&&mj<=600&&met>400&&njets>8"};

  Bin m1_r2_lowmet_lownj{"m1_r2_lowmet_lownj", "mt<=140&&mj>600&&met<=400&&njets<=8"};
  Bin m1_r2_lowmet_highnj{"m1_r2_lowmet_highnj", "mt<=140&&mj>600&&met<=400&&njets>8"};
  Bin m1_r2_highmet_lownj{"m1_r2_highmet_lownj", "mt<=140&&mj>600&&met>400&&njets<=8"};
  Bin m1_r2_highmet_highnj{"m1_r2_highmet_highnj", "mt<=140&&mj>600&&met>400&&njets>8"};

  Bin m1_r3_lowmet_lownj{"m1_r3_lowmet_lownj", "mt>140&&mj<=600&&met<=400&&njets<=8"};
  Bin m1_r3_lowmet_highnj{"m1_r3_lowmet_highnj", "mt>140&&mj<=600&&met<=400&&njets>8"};
  Bin m1_r3_highmet_lownj{"m1_r3_highmet_lownj", "mt>140&&mj<=600&&met>400&&njets<=8"};
  Bin m1_r3_highmet_highnj{"m1_r3_highmet_highnj", "mt>140&&mj<=600&&met>400&&njets>8"};

  Bin m1_r4_lowmet_lownj{"m1_r4_lowmet_lownj", "mt>140&&mj>600&&met<=400&&njets<=8"};
  Bin m1_r4_lowmet_highnj{"m1_r4_lowmet_highnj", "mt>140&&mj>600&&met<=400&&njets>8"};
  Bin m1_r4_highmet_lownj{"m1_r4_highmet_lownj", "mt>140&&mj>600&&met>400&&njets<=8"};
  Bin m1_r4_highmet_highnj{"m1_r4_highmet_highnj", "mt>140&&mj>600&&met>400&&njets>8"};

  //Specify ABCD constraints
  vector<Block> blocks_m3{
    {"lowmet_lownb", {{r1_lowmet_lownb, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb},
          {r3_lowmet_lownb, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb}}},
      {"lowmet_highnb", {{r1_lowmet_highnb, r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb},
            {r3_lowmet_highnb, r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb}}},
        {"highmet", {{r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
              {r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  vector<Block> blocks_m1{
    {"lowmet_lownj", {{m1_r1_lowmet_lownj, m1_r2_lowmet_lownj},
          {m1_r3_lowmet_lownj, m1_r4_lowmet_lownj}}},
      {"lowmet_highnj", {{m1_r1_lowmet_highnj, m1_r2_lowmet_highnj},
            {m1_r3_lowmet_highnj, m1_r4_lowmet_highnj}}},
        {"highmet_lownj", {{m1_r1_highmet_lownj, m1_r2_highmet_lownj},
              {m1_r3_highmet_lownj, m1_r4_highmet_lownj}}},
          {"highmet_highnj", {{m1_r1_highmet_highnj, m1_r2_highmet_highnj},
                {m1_r3_highmet_highnj, m1_r4_highmet_highnj}}}
  };

  vector<Block> blocks_test1{
    {"lowmet_lownj", {{m1_r1_lowmet_lownj, m1_r2_lowmet_lownj},
          {m1_r3_lowmet_lownj, m1_r4_lowmet_lownj}}}
  };
  vector<Block> blocks_test2{
    {"lowmet_highnj", {{m1_r1_lowmet_highnj, m1_r2_lowmet_highnj},
          {m1_r3_lowmet_highnj, m1_r4_lowmet_highnj}}},
      };
  vector<Block> blocks_test3{
    {"highmet_lownj", {{m1_r1_highmet_lownj, m1_r2_highmet_lownj},
          {m1_r3_highmet_lownj, m1_r4_highmet_lownj}}},
      };
  vector<Block> blocks_test4{
    {"highmet_highnj", {{m1_r1_highmet_highnj, m1_r2_highmet_highnj},
          {m1_r3_highmet_highnj, m1_r4_highmet_highnj}}}
  };

  ostringstream oss;
  oss << (10.*lumi) << flush;
  string lumi_string = oss.str();
  auto decimal = lumi_string.find('.');
  while( decimal != string::npos ) {
    lumi_string.erase(decimal,1);
    decimal = lumi_string.find('.');
  }

  map<BinProc, GammaParams> yields;
  MakeWorkspace("methodtest1_"+lumi_string+".root", baseline, blocks_test1, data, signal, backgrounds, yields);
  MakeWorkspace("methodtest2_"+lumi_string+".root", baseline, blocks_test2, data, signal, backgrounds, yields);
  MakeWorkspace("methodtest3_"+lumi_string+".root", baseline, blocks_test3, data, signal, backgrounds, yields);
  MakeWorkspace("methodtest4_"+lumi_string+".root", baseline, blocks_test4, data, signal, backgrounds, yields);
  MakeWorkspace("method3_"+lumi_string+".root", baseline, blocks_m3, data, signal, backgrounds, yields);
  MakeWorkspace("method1_"+lumi_string+".root", baseline, blocks_m1, data, signal, backgrounds, yields);
}

Bin::Bin(const string &name, const string &cut):
  name_(name),
  cut_(cut){
  }

bool Bin::operator<(const Bin &b) const{
  return name_ < b.name_
    || (name_ == b.name_
        && cut_ < b.cut_);
}

Process::Process(const string &name,
                 const vector<string> &file_names,
                 const string &cut,
                 bool count_zeros):
  chain_("tree", "tree"),
  name_(name),
  cut_(cut),
  count_zeros_(count_zeros){
  for(auto file_name = file_names.cbegin();
      file_name != file_names.cend();
      ++file_name){
    chain_.Add(file_name->c_str());
  }
  }

Process::Process(const string &name,
                 initializer_list<string> file_names,
                 const string &cut,
                 bool count_zeros):
  chain_("tree","tree"),
  name_(name),
  cut_(cut),
  count_zeros_(count_zeros){
  for(auto file_name = file_names.begin();
      file_name != file_names.end();
      ++file_name){
    chain_.Add(file_name->c_str());
  }
  }

bool Process::operator<(const Process &p) const{
  return name_ < p.name_
    || (name_ == p.name_
        && (cut_ < p.cut_
            || (cut_ == p.cut_
                && (count_zeros_ < p.count_zeros_
                    || (count_zeros_ == p.count_zeros_
                        && chain_.Hash() < p.chain_.Hash())))));
}

Block::Block(const string &name, const vector<vector<Bin> > &bins):
  bins_(bins),
  name_(name){
  }

Block::Block(const string &name, initializer_list<vector<Bin> > bins):
  bins_(bins),
  name_(name){
  }

BinProc::BinProc(const Bin &bin, Process &process):
  process_(process),
  bin_(bin){
  }

bool BinProc::operator<(const BinProc &bp) const{
  return bin_ < bp.bin_
    || (!(bp.bin_ < bin_)
        && process_ < bp.process_);
}

void GetYields(const vector<Block> &blocks,
               const string &baseline,
               Process &data,
               Process &signal,
               vector<reference_wrapper<Process> > &backgrounds,
               map<BinProc, GammaParams> &yields){
  for(auto block = blocks.cbegin();
      block != blocks.cend();
      ++block){
    for(auto vbin = block->bins_.cbegin();
        vbin != block->bins_.cend();
        ++vbin){
      for(auto bin = vbin->cbegin();
          bin != vbin->cend();
          ++bin){
        BinProc bp_data{*bin, data};
        StoreYield(bp_data, baseline, yields);
        BinProc bp_sig{*bin, signal};
        StoreYield(bp_sig, baseline, yields);
        for(auto bkg = backgrounds.cbegin();
            bkg != backgrounds.cend();
            ++bkg){
          BinProc bp{*bin, *bkg};
          StoreYield(bp, baseline, yields);
        }
      }
    }
  }
}

void StoreYield(const BinProc &bp,
                const string &baseline,
                map<BinProc, GammaParams> &yields){
  cout << "Getting yields for bin " << bp.bin_.name_
       << ", process " << bp.process_.name_ << endl;

  GammaParams gps;

  if(yields.find(bp) != yields.end()){
    cout << "Recycling already computed yield." << endl;
    gps = yields.at(bp);
  }else if(bp.process_.chain_.GetEntries() == 0){
    cout << "No entries found." << endl;
    gps.SetNEffectiveAndWeight(0., 0.);
  }else{
    ostringstream oss;
    oss << lumi << flush;
    string lumi_string = oss.str();
    array<string, 6> cuts;
    cuts.at(0) = "("+lumi_string+"*weight)*(("+baseline+")&&("+bp.bin_.cut_+")&&("+bp.process_.cut_+"))";
    cuts.at(1) = "("+lumi_string+"*weight)*(("+baseline+")&&("+bp.process_.cut_+"))";
    cuts.at(2) = "("+lumi_string+"*weight)*(&&("+bp.process_.cut_+"))";
    cuts.at(3) = "("+lumi_string+"*weight)";
    cuts.at(4) = lumi_string;
    cuts.at(5) = "1";

    for(size_t icut = 0;
        icut < cuts.size() && gps.Weight()<=0.;
        ++icut){
      if(icut > 0 && !bp.process_.count_zeros_){
        gps.SetNEffectiveAndWeight(0., 0.);
        break;
      }
      const string &cut = cuts.at(icut);
      cout << "Trying cut " << cut << endl;
      double count, uncertainty;
      GetCountAndUncertainty(bp.process_.chain_, cut, count, uncertainty);
      GammaParams temp_gps;
      temp_gps.SetYieldAndUncertainty(count, uncertainty);
      if(icut == 0){
        gps = temp_gps;
      }else{
        gps.SetNEffectiveAndWeight(0., temp_gps.Weight());
      }
    }
  }

  cout
    << "Found yield=" << gps.Yield()
    << ", uncertainty=" << gps.CorrectedUncertainty()
    << ", raw sqrt(n) uncertainty=" << gps.Uncertainty()
    << ", N_eff=" << gps.NEffective()
    << ", weight=" << gps.Weight()
    << "\n" << endl;
  yields[bp] =  gps;
}

void GetCountAndUncertainty(TTree &tree,
                            const string &cut,
                            double &count,
                            double &uncertainty){
  const string hist_name{"temp"};
  TH1D temp{hist_name.c_str(), "", 1, -1.0, 1.0};
  temp.Sumw2();
  tree.Project(hist_name.c_str(), "0.", cut.c_str());
  count=temp.IntegralAndError(0,2,uncertainty);
}

void MakeWorkspace(const string &file_name,
                   const string &baseline,
                   const vector<Block> &blocks,
                   Process &data,
                   Process &signal,
                   vector<reference_wrapper<Process> > &backgrounds,
                   map<BinProc, GammaParams> &yields){
  GetYields(blocks, baseline, data,
            signal, backgrounds, yields);

  RooWorkspace w{"w"};
  w.cd();
  w.factory("r[1.,0.,10.]");
  w.defineSet("POI","r");

  vector<string> obs_names, nuis_names;
  for(auto block = blocks.cbegin();
      block != blocks.cend();
      ++block){
    AddData(w, *block, data, yields, obs_names);
    AddBackgroundFractions(w, *block, backgrounds, yields, nuis_names);
    size_t max_col, max_row;
    AddABCDParams(w, *block, backgrounds, yields, nuis_names, max_col, max_row);
    AddBackgroundPreds(w, *block, backgrounds, max_col, max_row);
    AddSignalPreds(w, *block, signal, yields);
    AddBinPdfs(w, *block);
  }

  AddModels(w, blocks);

  DefineSet(w, "nuisances", nuis_names);
  DefineSet(w, "observables", obs_names);
  RooDataSet data_obs{"data_obs", "data_obs", *w.set("observables")};
  data_obs.add(*w.set("observables"));
  w.import(data_obs);

  ModelConfig model_config("ModelConfig", &w);
  model_config.SetPdf(*w.pdf("model_s"));
  model_config.SetParametersOfInterest(*w.set("POI"));
  model_config.SetObservables(*w.set("observables"));
  model_config.SetNuisanceParameters(*w.set("nuisances"));

  ModelConfig model_config_bonly("ModelConfig_bonly", &w);
  model_config_bonly.SetPdf(*w.pdf("model_b"));
  model_config_bonly.SetParametersOfInterest(*w.set("POI"));
  model_config_bonly.SetObservables(*w.set("observables"));
  model_config_bonly.SetNuisanceParameters(*w.set("nuisances"));

  w.import(model_config);
  w.import(model_config_bonly);

  w.writeToFile(file_name.c_str());
  PrintDiagnostics(w, blocks, data, signal, backgrounds, yields);
}

vector<double> GetBackgroundFractions(const Block &block,
                                      vector<reference_wrapper<Process> > &backgrounds,
                                      const map<BinProc, GammaParams> &yields){
  vector<double> output(backgrounds.size(), 0.);

  for(size_t ibkg = 0; ibkg < backgrounds.size(); ++ibkg){
    Process & bkg = backgrounds.at(ibkg);
    for(auto vbin = block.bins_.cbegin();
        vbin != block.bins_.cend();
        ++vbin){
      for(auto bin = vbin->cbegin();
          bin != vbin->cend();
          ++bin){
        BinProc bp{*bin, bkg};
        output.at(ibkg) += yields.at(bp).Yield();
      }
    }
  }

  double scale = 1./accumulate(output.cbegin(), output.cend(), 0.);
  transform(output.begin(), output.end(), output.begin(),
            bind1st(multiplies<double>(), scale));

  return output;
}

void AddBackgroundFractions(RooWorkspace &w,
                            const Block &block,
                            vector<reference_wrapper<Process> > &backgrounds,
                            const map<BinProc, GammaParams> &yields,
                            vector<string> &nuis_names){
  ostringstream oss;
  if(backgrounds.size()>1){
    vector<double> bkg_fracs = GetBackgroundFractions(block,
                                                      backgrounds,
                                                      yields);
    vector<string> list_of_names(bkg_fracs.size()-1);
    for(size_t ibkg = 0; ibkg < bkg_fracs.size()-1; ++ibkg){
      oss.str("");
      oss << "frac_BLK_" << block.name_
          << "_PRC_" << static_cast<Process&>(backgrounds.at(ibkg)).name_
          << flush;
      list_of_names.at(ibkg) = oss.str();
      nuis_names.push_back(oss.str());
      oss << "[" << bkg_fracs.at(ibkg) << ",0.,1.]" << flush;
      w.factory(oss.str().c_str());
    }
    oss.str("");
    oss << "expr::frac_BLK_" << block.name_ << "_PRC_"
        << static_cast<Process&>(backgrounds.at(backgrounds.size()-1)).name_
        << "('1";
    for(auto name = list_of_names.cbegin();
        name != list_of_names.cend();
        ++name){
      oss << "-" << (*name);
    }
    oss << "'";
    for(auto name = list_of_names.cbegin();
        name != list_of_names.cend();
        ++name){
      oss << "," << (*name);
    }
    oss << ")" << flush;
    w.factory(oss.str().c_str());
  }else{
    oss << "frac_BLK_" << block.name_ << "_PRC_"
        << static_cast<Process&>(backgrounds.back()).name_
        << "[1]" << flush;
    w.factory(oss.str().c_str());
  }
}

void AddABCDParams(RooWorkspace &w,
                   const Block &block,
                   vector<reference_wrapper<Process> > &backgrounds,
                   const map<BinProc, GammaParams> &yields,
                   vector<string> &nuis_names,
                   size_t &max_col, size_t &max_row){
  vector<double> row_sums(block.bins_.size());
  vector<double> col_sums(block.bins_.size() ? block.bins_.at(0).size() : 0);

  for(size_t irow = 0; irow < block.bins_.size(); ++irow){
    for(size_t icol = 0; icol < block.bins_.at(0).size(); ++icol){
      for(auto bkg = backgrounds.cbegin();
          bkg != backgrounds.cend();
          ++bkg){
        BinProc bp{block.bins_.at(irow).at(icol), *bkg};
        double yield = yields.at(bp).Yield();
        row_sums.at(irow) += yield;
        col_sums.at(icol) += yield;
      }
    }
  }

  max_col = MaxIndex(col_sums);
  max_row = MaxIndex(row_sums);

  double total = accumulate(row_sums.cbegin(), row_sums.cend(), 0.);

  ostringstream rxss, ryss;
  rxss << "sum::rxnorm_BLK_" << block.name_ << "(1.,";
  ryss << "sum::rynorm_BLK_" << block.name_ << "(1.,";
  ostringstream oss;
  oss << "norm_BLK_" << block.name_ << flush;
  nuis_names.push_back(oss.str());
  oss << "[" << max(1., total) << ",0.," << max(5.*total, 20.) << "]" << flush;
  w.factory(oss.str().c_str());

  for(size_t irow = 0; irow < row_sums.size(); ++irow){
    if(irow == max_row) continue;
    oss.str("");
    oss << "ry" << (irow+1) << (max_row+1) << "_BLK_" << block.name_ << flush;
    rxss << "," << oss.str();
    nuis_names.push_back(oss.str());
    oss << "[" << row_sums.at(irow)/row_sums.at(max_row)
        << ",0.,10.]" << flush;
    w.factory(oss.str().c_str());
  }
  rxss << ")" << flush;
  w.factory(rxss.str().c_str());
  for(size_t icol = 0; icol < col_sums.size(); ++icol){
    if(icol == max_col) continue;
    oss.str("");
    oss << "rx" << (icol+1) << (max_col+1) << "_BLK_" << block.name_ << flush;
    ryss << "," << oss.str();
    nuis_names.push_back(oss.str());
    oss << "[" << col_sums.at(icol)/col_sums.at(max_col)
        << ",0.,10.]" << flush;
    w.factory(oss.str().c_str());
  }
  ryss << ")" << flush;
  w.factory(ryss.str().c_str());
  oss.str("");
  oss << "prod::rnorm_BLK_" << block.name_
      << "(rxnorm_BLK_" << block.name_
      << ",rynorm_BLK_" << block.name_ << ")" << flush;
  w.factory(oss.str().c_str());
  oss.str("");
  oss << "expr::rscale_BLK_" << block.name_
      << "('norm_BLK_" << block.name_ << "/rnorm_BLK_" << block.name_
      << "',norm_BLK_" << block.name_
      << ",rnorm_BLK_" << block.name_ << ")" << flush;
  w.factory(oss.str().c_str());
}

void AddBackgroundPreds(RooWorkspace &w,
                        const Block &block,
                        const vector<reference_wrapper<Process> > &backgrounds,
                        size_t max_col, size_t max_row){
  for(size_t irow = 0; irow < block.bins_.size(); ++irow){
    bool no_ry = (irow == max_row);
    for(size_t icol = 0; icol < block.bins_.at(0).size(); ++icol){
      bool no_rx = (icol == max_col );
      const Bin & bin = block.bins_.at(irow).at(icol);
      vector<string> prod_list;
      for(size_t iprocess = 0; iprocess < backgrounds.size(); ++iprocess){
        const Process & bkg = backgrounds.at(iprocess);
        string prod_name = "rate_BLK_"+block.name_+"_BIN_"+bin.name_+"_PRC_"+bkg.name_;
        prod_list.push_back(prod_name);
        string fact_str = "prod::"+prod_name+"(rscale_BLK_"+block.name_;
        if(!no_rx){
          ostringstream oss;
          oss << ",rx" << (icol+1) << (max_col+1) << "_BLK_" << block.name_ << flush;
          fact_str += oss.str();
        }
        if(!no_ry){
          ostringstream oss;
          oss << ",ry" << (irow+1) << (max_row+1) << "_BLK_" << block.name_ << flush;
          fact_str += oss.str();
        }
        fact_str += (",frac_BLK_"+block.name_+"_PRC_"+bkg.name_+")");
        w.factory(fact_str.c_str());
      }
      string fact_str="sum::nbkg_BLK_"+block.name_+"_BIN_"+bin.name_+"(";
      for(size_t iprod = 0; iprod < prod_list.size(); ++iprod){
        const string &prod_name = prod_list.at(iprod);
        if(iprod != 0) fact_str += ",";
        fact_str += prod_name;
      }
      fact_str += ")";
      w.factory(fact_str.c_str());
    }
  }
}

void AddSignalPreds(RooWorkspace &w,
                    const Block &block,
                    Process &signal,
                    const map<BinProc, GammaParams> &yields){
  for(auto vbin = block.bins_.cbegin();
      vbin != block.bins_.cend();
      ++vbin){
    for(auto bin = vbin->cbegin();
        bin != vbin->cend();
        ++bin){
      BinProc bp{*bin, signal};
      double yield = yields.at(bp).Yield();
      ostringstream oss;
      oss << "rate_BLK_" << block.name_
          << "_BIN_" << bin->name_
          << "_PRC_" << signal.name_
          << "[" << yield << "]" << flush;
      w.factory(oss.str().c_str());
      oss.str("");
      oss << "prod::nsig_BLK_" << block.name_
          << "_BIN_" << bin->name_
          << "(r,"
          << "rate_BLK_" << block.name_
          << "_BIN_" << bin->name_
          << "_PRC_" << signal.name_
          << ")" << flush;
      w.factory(oss.str().c_str());
    }
  }
}

void AddBinPdfs(RooWorkspace &w,
                const Block &block){
  string null_list{""}, alt_list{""};
  for(auto vbin = block.bins_.cbegin();
      vbin != block.bins_.cend();
      ++vbin){
    for(auto bin = vbin->cbegin();
        bin != vbin->cend();
        ++bin){
      string bb_name = "_BLK_"+block.name_ +"_BIN_"+bin->name_;
      string null_name = "pdf_null"+bb_name;
      string alt_name = "pdf_alt"+bb_name;
      if(vbin != block.bins_.cbegin() || bin != vbin->cbegin()){
        null_list += ",";
        alt_list += ",";
      }
      null_list += null_name;
      alt_list += alt_name;
      w.factory(("sum::nexp"+bb_name+"(nbkg"+bb_name+",nsig"+bb_name+")").c_str());
      w.factory(("RooPoisson::pdf_null"+bb_name+"(nobs"+bb_name+",nbkg"+bb_name+")").c_str());
      (static_cast<RooPoisson*>(w.pdf(("pdf_null"+bb_name).c_str())))->setNoRounding();
      w.factory(("RooPoisson::pdf_alt"+bb_name+"(nobs"+bb_name+",nexp"+bb_name+")").c_str());
      (static_cast<RooPoisson*>(w.pdf(("pdf_alt"+bb_name).c_str())))->setNoRounding();
    }
  }
  w.factory(("PROD:pdf_null_BLK_"+block.name_+"("+null_list+")").c_str());
  w.factory(("PROD:pdf_alt_BLK_"+block.name_+"("+alt_list+")").c_str());
}

void AddData(RooWorkspace &w,
             const Block &block,
             Process &data,
             const map<BinProc, GammaParams> &yields,
             vector<string> &obs_names){
  for(auto vbin = block.bins_.cbegin();
      vbin != block.bins_.cend();
      ++vbin){
    for(auto bin = vbin->cbegin();
        bin != vbin->cend();
        ++bin){
      BinProc bp{*bin, data};
      double yield = yields.at(bp).Yield();
      ostringstream oss;
      oss << "nobs_BLK_" << block.name_
          << "_BIN_" << bin->name_ << flush;
      obs_names.push_back(oss.str());
      oss << "[" << yield << "]" << flush;
      w.factory(oss.str().c_str());
    }
  }
}

void DefineSet(RooWorkspace &w,
               const string &set_name,
               const vector<string> &var_names){
  if(var_names.size()==0){
    w.defineSet(set_name.c_str(), "");
  }else{
    string cat_names = var_names.at(0);
    for(size_t ivar = 1; ivar < var_names.size(); ++ivar){
      cat_names += ("," + var_names.at(ivar));
    }
    w.defineSet(set_name.c_str(), cat_names.c_str());
  }
}

void AddModels(RooWorkspace &w,
               const vector<Block> & blocks){
  if(blocks.size() == 0){
    w.factory("RooPoisson::model_b(0,0)");
    w.factory("RooPoisson::model_s(0,0)");
  }else{
    string null_list = "pdf_null_BLK_"+blocks.at(0).name_;
    string alt_list = "pdf_alt_BLK_"+blocks.at(0).name_;
    for(size_t iblock = 1; iblock < blocks.size(); ++iblock){
      null_list += (",pdf_null_BLK_"+blocks.at(iblock).name_);
      alt_list += (",pdf_alt_BLK_"+blocks.at(iblock).name_);
    }
    w.factory(("PROD::model_b("+null_list+")").c_str());
    w.factory(("PROD::model_s("+alt_list+")").c_str());
  }
}

size_t MaxIndex(const vector<double> &v){
  if(v.size() == 0) return -1;
  size_t imax = 0;
  for(size_t i = 1; i < v.size(); ++i){
    if(v.at(i) > v.at(imax)){
      imax = i;
    }
  }
  return imax;
}

void PrintDiagnostics(const RooWorkspace &w,
                      const vector<Block> &blocks,
                      Process &data,
                      Process &signal,
                      vector<reference_wrapper<Process> > &backgrounds,
                      const map<BinProc, GammaParams> &yields){
  for(auto block = blocks.cbegin();
      block != blocks.cend();
      ++block){
    for(auto vbin = block->bins_.cbegin();
        vbin != block->bins_.cend();
        ++vbin){
      for(auto bin = vbin->cbegin();
          bin != vbin->cend();
          ++bin){

        BinProc bp_data{*bin, data};
        PrintComparison(w, *block, bp_data, yields, true);
        BinProc bp_sig{*bin, signal};
        PrintComparison(w, *block, bp_sig, yields);
        for(auto bkg = backgrounds.cbegin();
            bkg != backgrounds.cend();
            ++bkg){
          BinProc bp{*bin, *bkg};
          PrintComparison(w, *block, bp, yields);
        }
      }
    }
  }
}

void PrintComparison(const RooWorkspace &w,
                     const Block &block,
                     BinProc &bp,
                     const std::map<BinProc, GammaParams> &yields,
                     bool is_data){
  GammaParams gp = yields.at(bp);
  ostringstream name;
  name << (is_data ? "nobs" : "rate")
       << "_BLK_" << block.name_
       << "_BIN_" << bp.bin_.name_;
  if(!is_data){
    name << "_PRC_" << bp.process_.name_;
  }
  name << flush;

  cout << fixed << setprecision(2);
  cout << setw(64) << name.str() << ": "
       << setw(8) << gp.Yield()
       << " +- " << setw(8) << gp.CorrectedUncertainty()
       << " => ";
  RooAbsReal *fp = w.function(name.str().c_str());
  if(fp){
    cout << setw(8) << fp->getVal() << endl;
  }else{
    cout << "Not found" << endl;
  }
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:d:l:s:tv", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
