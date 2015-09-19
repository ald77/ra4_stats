#include "make_workspace.hpp"

#include <iostream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <initializer_list>
#include <map>
#include <array>
#include <vector>
#include <string>

#include "TChain.h"
#include "TH1D.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"

#include "RooStats/ModelConfig.h"

#include "gamma_params.hpp"

using namespace std;
using namespace RooStats;

int main(){
  //Define processes. Try to minimize splitting
  Process ttbar{"ttbar", {
      {"archive/2015_08_13/*TTJets*.root/tree"}
    }};
  Process other{"other", {
      {"archive/2015_08_13/*_ST_*.root/tree"},
	{"archive/2015_08_13/*WJetsToLNu*.root/tree"}
	  //{"archive/2015_08_13/*_WWTo*.root/tree"},
	  //{"archive/2015_08_13/*ttHJetTobb*.root/tree"},
	  //{"archive/2015_08_13/*DYJetsToLL*.root/tree"},
	  //{"archive/2015_08_13/*QCD_Pt*.root/tree"}
    }};
  Process signal{"signal", {
      {"archive/2015_07_22/small_quick_SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_MINIAODSIM_UCSB2377_v78.root/tree"}
    }};
  Process data{"data", {}};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  vector<reference_wrapper<Process> > backgrounds{ref(ttbar), ref(other)};

  //Baseline selection applied to all bins and processes
  string baseline{"ht>500&&met>200&njets>=7&&nbm>=2"};

  //Declare bins
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

  //Specify ABCD constraints
  vector<Block> blocks{
    {"lowmet_lownb", {{r1_lowmet_lownb, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb},
	  {r3_lowmet_lownb, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb}}},
      {"lowmet_highnb", {{r1_lowmet_highnb, r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb},
	    {r3_lowmet_highnb, r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb}}},
	{"highmet", {{r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
	      {r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
  };
  vector<Block> blocks2{
    {"lowmet_lownb", {{r1_lowmet_lownb, r2_lowmet_lownj_lownb},
	  {r3_lowmet_lownb, r4_lowmet_lownj_lownb}}},
      {"lowmet_highnb", {{r1_lowmet_highnb, r2_lowmet_lownj_highnb},
	    {r3_lowmet_highnb, r4_lowmet_lownj_highnb}}},
	{"highmet", {{r1_highmet, r2_highmet_lownj},
	      {r3_highmet, r4_highmet_lownj}}}
  };

  MakeWorkspace("method3.root", baseline, blocks, data, signal, backgrounds);
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

map<BinProc, GammaParams> GetYields(const vector<Block> &blocks,
				    const string &baseline,
				    Process &data,
				    Process &signal,
				    vector<reference_wrapper<Process> > &backgrounds){
  map<BinProc, GammaParams> yields;
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
	yields[bp_data] = GetYield(bp_data, baseline);
	BinProc bp_sig{*bin, signal};
	yields[bp_sig] = GetYield(bp_sig, baseline);
	for(auto bkg = backgrounds.cbegin();
	    bkg != backgrounds.cend();
	    ++bkg){
	  BinProc bp{*bin, *bkg};
	  yields[bp] = GetYield(bp, baseline);
	}
      }
    }
  }
  return yields;
}

GammaParams GetYield(const BinProc &bp,
		     const string &baseline){
  cout << "Getting yields for bin " << bp.bin_.name_
       << ", process " << bp.process_.name_ << endl;
  if(bp.process_.chain_.GetEntries() == 0){
    cout << "No entries found.\n" << endl;
    return {0., 0.};
  }
  array<string, 4> cuts;
  cuts.at(0) = "3*weight*(("+baseline+")&&("+bp.bin_.cut_+")&&("+bp.process_.cut_+"))";
  cuts.at(1) = "3*weight*(("+baseline+")&&("+bp.process_.cut_+"))";
  cuts.at(2) = "3*weight*(&&("+bp.process_.cut_+"))";
  cuts.at(3) = "3*weight";

  GammaParams gps;

  for(size_t icut = 0;
      icut < cuts.size() && gps.NEffective()<=0. && gps.Weight()<=0.;
      ++icut){
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

  cout
    << "Found yield=" << gps.Yield()
    << ", uncertainty=" << gps.Uncertainty()
    << ", N_eff=" << gps.NEffective()
    << ", weight=" << gps.Weight()
    << "\n" << endl;
  return gps;
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
		   vector<reference_wrapper<Process> > &backgrounds){
  map<BinProc, GammaParams> yields = GetYields(blocks, baseline, data,
					       signal, backgrounds);

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
    AddABCDParams(w, *block, backgrounds, yields, nuis_names); 
    AddBackgroundPreds(w, *block, backgrounds);
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
  w.Print();
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
		   vector<string> &nuis_names){
  vector<double> row_sums(block.bins_.size());
  vector<double> col_sums(block.bins_.size() ? block.bins_.at(0).size() : 0);

  for(size_t irow = 0; irow < block.bins_.size(); ++irow){
    for(size_t icol = 0; icol < block.bins_.size(); ++icol){
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

  double total = accumulate(row_sums.cbegin(), row_sums.cend(), 0.);

  ostringstream oss;
  oss << "norm_BLK_" << block.name_ << flush;
  nuis_names.push_back(oss.str());
  oss << "[" << max(1., total) << ",0.," << max(5.*total, 20.) << "]" << flush;
  w.factory(oss.str().c_str());

  for(size_t irow = 0; irow < row_sums.size() - 1; ++irow){
    oss.str("");
    oss << "ry" << (irow+1) << "_BLK_" << block.name_ << flush;
    nuis_names.push_back(oss.str());
    oss << "[" << max(1.,row_sums.at(irow)/row_sums.at(row_sums.size()-1))
	<< ",0.,100.]" << flush;
    w.factory(oss.str().c_str());
  }
  for(size_t icol = 0; icol < col_sums.size() - 1; ++icol){
    oss.str("");
    oss << "rx" << (icol+1) << "_BLK_" << block.name_ << flush;
    nuis_names.push_back(oss.str());
    oss << "[" << max(1.,col_sums.at(icol)/col_sums.at(col_sums.size()-1))
	<< ",0.,100.]" << flush;
    w.factory(oss.str().c_str());
  }
}

void AddBackgroundPreds(RooWorkspace &w,
			const Block &block,
			const vector<reference_wrapper<Process> > &backgrounds){
  for(size_t irow = 0; irow < block.bins_.size(); ++irow){
    bool no_ry = (irow == (block.bins_.size()-1) );
    for(size_t icol = 0; icol < block.bins_.at(0).size(); ++icol){
      bool no_rx = (icol == (block.bins_.at(0).size()-1) );
      const Bin & bin = block.bins_.at(irow).at(icol);
      vector<string> prod_list;
      for(size_t iprocess = 0; iprocess < backgrounds.size(); ++iprocess){
	const Process & bkg = backgrounds.at(iprocess);
	string prod_name = "rate_BLK_"+block.name_+"_BIN_"+bin.name_+"_PRC_"+bkg.name_;
	prod_list.push_back(prod_name);
	string fact_str = "prod::"+prod_name+"(norm_BLK_"+block.name_;
	if(!no_rx){
	  ostringstream oss;
	  oss << ",rx" << (icol+1) << "_BLK_" << block.name_ << flush;
	  fact_str += oss.str();
	}
	if(!no_ry){
	  ostringstream oss;
	  oss << ",ry" << (irow+1) << "_BLK_" << block.name_ << flush;
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
      w.factory(("RooPoisson::pdf_alt"+bb_name+"(nobs"+bb_name+",nexp"+bb_name+")").c_str());
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
