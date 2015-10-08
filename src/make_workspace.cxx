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

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPoisson.h"

#include "RooStats/ModelConfig.h"

#include "gamma_params.hpp"
#include "bin.hpp"
#include "process.hpp"
#include "bin_proc.hpp"
#include "utilities.hpp"
#include "systematic.hpp"

using namespace std;
using namespace RooStats;

namespace{
  double lumi = 3.;
  bool blinded = true;
  bool do_syst = true;
}

int main(int argc, char *argv[]){
  cout << fixed << setprecision(2);
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
  Process signal_nc{"signal_nc", {
      {"archive/2015_09_28_ana/skim/*T1tttt*1500*100*.root/tree"}
    }};
  Process signal_c{"signal_c", {
      {"archive/2015_09_28_ana/skim/*T1tttt*1200*800*.root/tree"}
    }};
  Process data{"data", {}};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  vector<reference_wrapper<Process> > backgrounds{ref(ttbar), ref(other)};

  //Baseline selection applied to all bins and processes
  string baseline{"ht>500&&met>200&njets>=7&&nbm>=2&&(nels+nmus)==1"};
  string baseline_david{"ht>500&&met>250&njets>=9&&nbm>=2&&(nels+nmus)==1"};

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

  //Method David
  Bin r1{"r1", "mt<=140&&mj<=400"};
  Bin r2{"r2", "mt<=140&&mj>400"};
  Bin r3{"r3", "mt>140&&mj<=400"};
  Bin r4{"r4", "mt>140&&mj>400"};

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

  vector<Block> blocks_david{
    {"all", {{r1, r2}, {r3, r4}}}
  };

  ostringstream oss;
  oss << (10.*lumi) << flush;
  string lumi_string = oss.str();
  auto decimal = lumi_string.find('.');
  while( decimal != string::npos ) {
    lumi_string.erase(decimal,1);
    decimal = lumi_string.find('.');
  }
  string no_syst = do_syst ? "" : "_nosyst";

  map<BinProc, GammaParams> yields, dilep_yields;
  MakeWorkspace("method3nc_"+lumi_string+no_syst+".root",
		baseline, blocks_m3, data, signal_nc, backgrounds, yields);
  MakeWorkspace("method3c_"+lumi_string+no_syst+".root",
		baseline, blocks_m3, data, signal_c, backgrounds, yields);
  MakeWorkspace("method1nc_"+lumi_string+no_syst+".root",
		baseline, blocks_m1, data, signal_nc, backgrounds, yields);
  MakeWorkspace("method1c_"+lumi_string+no_syst+".root",
  baseline, blocks_m1, data, signal_c, backgrounds, yields);
  MakeWorkspace("methoddavidnc_"+lumi_string+no_syst+".root",
		baseline_david, blocks_david, data, signal_nc, backgrounds, yields);
  MakeWorkspace("methoddavidc_"+lumi_string+no_syst+".root",
		baseline_david, blocks_david, data, signal_c, backgrounds, yields);
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
	Bin dilep_bin = *bin;
	string dilep_baseline = baseline;
	bool do_dilep = NeedsDileptonBin(*bin, baseline);
	if(do_dilep) MakeDileptonBin(*bin, baseline,
				     dilep_bin, dilep_baseline);
        BinProc bp_data{*bin, data};
        if(!blinded) StoreYield(bp_data, baseline, yields);
        BinProc bp_sig{*bin, signal};
        StoreYield(bp_sig, baseline, yields);
        for(auto bkg = backgrounds.cbegin();
            bkg != backgrounds.cend();
            ++bkg){
          BinProc bp{*bin, *bkg};
          StoreYield(bp, baseline, yields);
	  if(do_dilep) StoreYield(BinProc{dilep_bin, *bkg}, dilep_baseline, yields);
        }
      }
    }
  }
}

bool NeedsDileptonBin(const Bin &bin, const string &baseline){
  return Contains(bin.cut_, "mt>")
    && (Contains(bin.cut_, "(nels+nmus)==1")
	|| Contains(bin.cut_, "(nmus+nels)==1")
	|| Contains(bin.cut_, "nels+nmus==1")
	|| Contains(bin.cut_, "nmus+nels==1")
	|| Contains(bin.cut_, "nleps==1")
	|| Contains(baseline, "(nels+nmus)==1")
	|| Contains(baseline, "(nmus+nels)==1")
	|| Contains(baseline, "nels+nmus==1")
	|| Contains(baseline, "nmus+nels==1")
	|| Contains(baseline, "nleps==1"));
}

void MakeDileptonBin(const Bin &bin, const string &baseline,
		     Bin &dilep_bin, string &dilep_baseline){
  dilep_bin = bin;
  dilep_bin.name_ = "DILEPTON_"+dilep_bin.name_;
  dilep_baseline = baseline;
  ReplaceAll(dilep_bin.cut_, "(nels+nmus)==1", "(nels+nmus)==2");
  ReplaceAll(dilep_bin.cut_, "(nmus+nels)==1", "(nmus+nels)==2");
  ReplaceAll(dilep_bin.cut_, "nels+nmus==1", "nels+nmus==2");
  ReplaceAll(dilep_bin.cut_, "nmus+nels==1", "nmus+nels==2");
  ReplaceAll(dilep_bin.cut_, "nleps==1", "nleps==2");
  RmCutOn(dilep_bin.cut_, "nbm", "nbm>=1&&nbm<=2");
  RmCutOn(dilep_bin.cut_, "met", "met>200&&met<=400");
  ReplaceAll(dilep_baseline, "(nels+nmus)==1", "(nels+nmus)==2");
  ReplaceAll(dilep_baseline, "(nmus+nels)==1", "(nmus+nels)==2");
  ReplaceAll(dilep_baseline, "nels+nmus==1", "nels+nmus==2");
  ReplaceAll(dilep_baseline, "nmus+nels==1", "nmus+nels==2");
  ReplaceAll(dilep_baseline, "nleps==1", "nleps==2");
  RmCutOn(dilep_baseline, "nbm", "nbm>=1&&nbm<=2");
  RmCutOn(dilep_baseline, "met", "met>200&&met<=400");
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

void MakeWorkspace(const string &file_name,
                   const string &baseline,
                   vector<Block> blocks,
                   Process &data,
                   Process &signal,
                   vector<reference_wrapper<Process> > &backgrounds,
                   map<BinProc, GammaParams> &yields){
  string baseline_fix = baseline;
  ReplaceAll(baseline_fix, " ", "");
  GetYields(blocks, baseline_fix, data,
            signal, backgrounds, yields);

  RooWorkspace w{"w"};
  w.cd();
  w.factory("r[1.,0.,10.]");
  w.defineSet("POI","r");

  vector<string> obs_names, nuis_names;
  set<string> syst_generators;
  for(auto block = blocks.begin();
      block != blocks.end();
      ++block){
    if(!blinded){
      AddData(w, *block, data, yields, obs_names);
    }else{
      AddMockData(w, *block, backgrounds, yields, obs_names);
    }
    AddBackgroundFractions(w, *block, backgrounds, yields, nuis_names);
    size_t max_col, max_row;
    AddABCDParams(w, *block, backgrounds, yields, nuis_names, max_col, max_row);
    if(do_syst) AddDileptonSystematics(*block, baseline, backgrounds, yields);
    AddBackgroundPreds(w, *block, backgrounds, max_col, max_row, syst_generators, nuis_names);
    AddSignalPreds(w, *block, signal, yields);
    AddBinPdfs(w, *block);
  }

  AddModels(w, blocks, syst_generators);

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
  cout << "Wrote workspace to " << file_name << "!\n" << endl;
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

void AddDileptonSystematics(Block &block,
			    const string &baseline,
			    vector<reference_wrapper<Process> > &backgrounds,
			    const map<BinProc, GammaParams> &yields){
  for(auto vbin = block.bins_.begin();
      vbin != block.bins_.end();
      ++vbin){
    for(auto bin = vbin->begin();
	bin != vbin->end();
	++bin){
      Bin dilep_bin = *bin;
      string dilep_baseline = baseline;
      MakeDileptonBin(*bin, baseline,
		      dilep_bin, dilep_baseline);
      GammaParams gp{0., 0.};
      GammaParams dilep_gp{0., 0.};
      bool found_one = false;
      for(auto bkg = backgrounds.cbegin();
	  bkg != backgrounds.cend();
	  ++bkg){
	BinProc bp{*bin, *bkg};
	BinProc dilep_bp{dilep_bin, *bkg};
	if(yields.find(dilep_bp) == yields.end()
	   || yields.find(bp) == yields.end()) continue;
	found_one = true;
	gp += yields.at(bp);
	dilep_gp += yields.at(dilep_bp);
      }
      if(!found_one) continue;
      double syst = 1.;
      string name = "DILEP_CLOSE_"+bin->name_;
      if(dilep_gp.Yield()>1.){
	syst = 1./sqrt(dilep_gp.Yield());
      }
      bin->systematics_.push_back(Systematic{name, syst});
    }
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

void AddSystGenerator(RooWorkspace &w,
		      set<string> &syst_generators,
		      vector<string> &nuis_names,
		      const string &name){
  if(syst_generators.find(name) != syst_generators.end()) return;
  w.factory(("RooGaussian::CONSTRAINT_"+name+"("+name+"[0.,-10.,10.],0.,1.)").c_str());
  nuis_names.push_back(name);
  syst_generators.insert(name);
}

void AddBackgroundPreds(RooWorkspace &w,
                        const Block &block,
                        const vector<reference_wrapper<Process> > &backgrounds,
                        size_t max_col, size_t max_row,
			set<string> &syst_generators,
			vector<string> &nuis_names){
  for(size_t irow = 0; irow < block.bins_.size(); ++irow){
    bool no_ry = (irow == max_row);
    for(size_t icol = 0; icol < block.bins_.at(0).size(); ++icol){
      bool no_rx = (icol == max_col );
      const Bin &bin = block.bins_.at(irow).at(icol);
      string bb_name = "BLK_"+block.name_+"_BIN_"+bin.name_;
      vector<string> prod_list;
      for(size_t iprocess = 0; iprocess < backgrounds.size(); ++iprocess){
        const Process & bkg = backgrounds.at(iprocess);
        string prod_name = "rate_"+bb_name+"_PRC_"+bkg.name_;
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
      string fact_str="sum::nbkg_raw_"+bb_name+"(";
      for(size_t iprod = 0; iprod < prod_list.size(); ++iprod){
        const string &prod_name = prod_list.at(iprod);
        if(iprod != 0) fact_str += ",";
        fact_str += prod_name;
      }
      fact_str += ")";
      w.factory(fact_str.c_str());

      fact_str="prod::nbkg_"+bb_name+"("
	+"nbkg_raw_"+bb_name;
      for(auto syst = bin.systematics_.cbegin();
	  syst != bin.systematics_.cend();
	  ++syst){
	AddSystGenerator(w, syst_generators, nuis_names, syst->base_name_);
	string full_name = syst->base_name_ + "_" + bb_name;
	string sname = "strength_"+full_name;
	ostringstream oss;
	oss << sname << "[" << syst->multiplier_ << "]" << flush;
	w.factory(oss.str().c_str());
	oss.str("");
	oss << "expr::" << full_name
	    << "('exp(" << sname << "*" << syst->base_name_ << ")',"
	    << sname << "," << syst->base_name_ << ")" << flush;
	w.factory(oss.str().c_str());
	fact_str+=","+full_name;
      }
      fact_str+=")";
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

void AddMockData(RooWorkspace &w,
		 const Block &block,
		 vector<reference_wrapper<Process> > &backgrounds,
		 const map<BinProc, GammaParams> &yields,
		 vector<string> &obs_names){
  for(auto vbin = block.bins_.cbegin();
      vbin != block.bins_.cend();
      ++vbin){
    for(auto bin = vbin->cbegin();
        bin != vbin->cend();
        ++bin){
      GammaParams gp;
      for(auto bkg = backgrounds.cbegin();
	  bkg != backgrounds.cend();
	  ++bkg){
	BinProc bp{*bin, *bkg};
	gp += yields.at(bp);
      }
      ostringstream oss;
      oss << "nobs_BLK_" << block.name_
          << "_BIN_" << bin->name_ << flush;
      obs_names.push_back(oss.str());
      oss << "[" << gp.Yield() << "]" << flush;
      w.factory(oss.str().c_str());
    }
  }
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

void AddModels(RooWorkspace &w,
               const vector<Block> & blocks,
	       const set<string> &syst_generators){
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
    if(do_syst){
      for(auto syst = syst_generators.cbegin();
	  syst != syst_generators.cend();
	  ++syst){
	null_list += (",CONSTRAINT_"+(*syst));
	alt_list += (",CONSTRAINT_"+(*syst));
      }
    }
    w.factory(("PROD::model_b("+null_list+")").c_str());
    w.factory(("PROD::model_s("+alt_list+")").c_str());
  }
}

void PrintDiagnostics(const RooWorkspace &w,
                      const vector<Block> &blocks,
                      Process &data,
                      Process &signal,
                      vector<reference_wrapper<Process> > &backgrounds,
                      const map<BinProc, GammaParams> &yields){
  w.Print();
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
                     const map<BinProc, GammaParams> &yields,
                     bool is_data){
  GammaParams gp{0., 0.};
  if(yields.find(bp) != yields.end()) gp = yields.at(bp);

  ostringstream name;
  name << (is_data ? "nobs" : "rate")
       << "_BLK_" << block.name_
       << "_BIN_" << bp.bin_.name_;
  if(!is_data){
    name << "_PRC_" << bp.process_.name_;
  }
  name << flush;

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
      {"unblind", no_argument, 0, 'u'},
      {"no_syst", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:u", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 'u':
      blinded = false;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "no_syst"){
	do_syst = false;
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
