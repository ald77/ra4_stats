#include "workspace_generator.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>
#include <algorithm>
#include <numeric>

#include "RooPoisson.h"
#include "RooDataSet.h"

#include "RooStats/ModelConfig.h"

#include "utilities.hpp"

using namespace std;

YieldManager WorkspaceGenerator::yields_ = YieldManager(4.);

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
  systematics_(),
  free_systematics_(),
  luminosity_(4.),
  print_level_(PrintLevel::normal),
  blind_level_(BlindLevel::blinded),
  do_systematics_(true),
  do_dilepton_(true),
  do_mc_kappa_correction_(true),
  w_is_valid_(false){
  w_.cd();
}

void WorkspaceGenerator::WriteToFile(const string &file_name){
  if(print_level_ >= PrintLevel::everything){
    cout << "WriteToFile(" << file_name << ")" << endl;
  }
  if(!w_is_valid_) UpdateWorkspace();

  w_.writeToFile(file_name.c_str());
  if(print_level_ >= PrintLevel::everything){
    w_.Print();
  }
  if(print_level_ >= PrintLevel::normal){
    cout << *this << endl;
  }
  if(print_level_ >= PrintLevel::important){
    cout << "Wrote workspace to file " << file_name << endl;
  }
}

double WorkspaceGenerator::GetLuminosity() const{
  return luminosity_;
}

WorkspaceGenerator & WorkspaceGenerator::SetLuminosity(double luminosity){
  if(luminosity != luminosity_){
    luminosity_ = luminosity;
    w_is_valid_ = false;
  }
  return *this;
}

WorkspaceGenerator::BlindLevel WorkspaceGenerator::GetBlindLevel() const{
  return blind_level_;
}

WorkspaceGenerator & WorkspaceGenerator::SetBlindLevel(BlindLevel blind_level){
  if(blind_level != blind_level_){
    blind_level_ = blind_level;
    w_is_valid_ = false;
  }
  return *this;
}

bool WorkspaceGenerator::GetDoSystematics() const{
  return do_systematics_;
}

WorkspaceGenerator & WorkspaceGenerator::SetDoSystematics(bool do_systematics){
  if(do_systematics != do_systematics_){
    do_systematics_ = do_systematics;
    w_is_valid_ = false;
  }
  return *this;
}

bool WorkspaceGenerator::GetDoDilepton() const{
  return do_dilepton_;
}

WorkspaceGenerator & WorkspaceGenerator::SetDoDilepton(bool do_dilepton){
  if(do_dilepton != do_dilepton_){
    do_dilepton_ = do_dilepton;
    w_is_valid_ = false;
  }
  return *this;
}

WorkspaceGenerator::PrintLevel WorkspaceGenerator::GetPrintLevel() const{
  return print_level_;
}

WorkspaceGenerator & WorkspaceGenerator::SetPrintLevel(WorkspaceGenerator::PrintLevel print_level){
  print_level_ = print_level;
  return *this;
}

bool WorkspaceGenerator::GetKappaCorrected() const{
  return do_mc_kappa_correction_;
}

WorkspaceGenerator & WorkspaceGenerator::SetKappaCorrected(bool do_kappa_correction){
  if(do_mc_kappa_correction_ != do_kappa_correction){
    do_mc_kappa_correction_ = do_kappa_correction;
    w_is_valid_ = false;
  }
  return *this;
}

GammaParams WorkspaceGenerator::GetYield(const YieldKey &key) const{
  yields_.Luminosity() = luminosity_;
  return yields_.GetYield(key);
}

GammaParams WorkspaceGenerator::GetYield(const Bin &bin,
					 const Process &process,
					 const Cut &cut) const{
  return GetYield(YieldKey(bin, process, cut));
}

GammaParams WorkspaceGenerator::GetYield(const Bin &bin,
					 const Process &process) const{
  return GetYield(bin, process, baseline_);
}

void WorkspaceGenerator::UpdateWorkspace(){
  if(print_level_ >= PrintLevel::everything){
    cout << "UpdateWorkspace()" << endl;
  }
  w_ = RooWorkspace("w");
  w_.cd();

  if(do_dilepton_){
    AddDileptonSystematic();
  }
  AddPOI();
  AddSystematicsGenerators();

  for(const auto &block: blocks_){
    AddData(block);
    AddMCYields(block);
    AddMCPdfs(block);
    AddMCProcessSums(block);
    AddBackgroundFractions(block);
    AddABCDParameters(block);
    AddRawBackgroundPredictions(block);
    if(do_mc_kappa_correction_) AddKappas(block);
    AddFullBackgroundPredictions(block);
    AddSignalPredictions(block);
    AddPdfs(block);
  }

  AddFullPdf();
  AddParameterSets();
  AddModels();

  w_is_valid_ = true;
}

void WorkspaceGenerator::AddPOI(){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddPOI()" << endl;
  }
  w_.factory("r[1.,0.,20.]");
  Append(poi_, "r");
}

void WorkspaceGenerator::AddDileptonSystematic(){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddDileptonSystematic()" << endl;
  }

  set<Block> new_blocks;
  for(const auto &block: blocks_){
    Block new_block = block;
    for(auto &vbin: new_block.Bins()){
      for(auto &bin: vbin){
	if(!NeedsDileptonBin(bin)) continue;
	Bin dilep_bin = bin;
	Cut dilep_baseline = baseline_;
	MakeDileptonBin(bin, dilep_bin, dilep_baseline);
	GammaParams dilep_gp(0., 0.);
	if(blind_level_ != BlindLevel::unblinded){
	  for(const auto &bkg: backgrounds_){
	    dilep_gp += GetYield(dilep_bin, bkg, dilep_baseline);
	  }
	}else{
	  dilep_gp = GetYield(dilep_bin, data_, dilep_baseline);
	}
	double strength = 1.;
	string name = "dilep_"+bin.Name();
	if(dilep_gp.Yield()>1.){
	  strength = 1./sqrt(dilep_gp.Yield());
	}
	Systematic syst(name, strength);
	bin.AddSystematic(syst);
      }
    }
    Append(new_blocks, new_block);
  }
  blocks_ = new_blocks;
}

bool WorkspaceGenerator::NeedsDileptonBin(const Bin &bin) const{
  if(print_level_ >= PrintLevel::everything){
    cout << "NeedsDileptonBin(" << bin << ")" << endl;
  }
  return Contains(static_cast<string>(bin.Cut()), "mt>")
    && (Contains(static_cast<string>(bin.Cut()), "(nels+nmus)==1")
	|| Contains(static_cast<string>(bin.Cut()), "(nmus+nels)==1")
	|| Contains(static_cast<string>(bin.Cut()), "nels+nmus==1")
	|| Contains(static_cast<string>(bin.Cut()), "nmus+nels==1")
	|| Contains(static_cast<string>(bin.Cut()), "nleps==1")
	|| Contains(static_cast<string>(baseline_), "(nels+nmus)==1")
	|| Contains(static_cast<string>(baseline_), "(nmus+nels)==1")
	|| Contains(static_cast<string>(baseline_), "nels+nmus==1")
	|| Contains(static_cast<string>(baseline_), "nmus+nels==1")
	|| Contains(static_cast<string>(baseline_), "nleps==1"));
}

void WorkspaceGenerator::MakeDileptonBin(const Bin &bin, Bin &dilep_bin, Cut &dilep_cut) const{
  if(print_level_ >= PrintLevel::everything){
    cout << "MakeDileptonBin(" << bin << ", " << dilep_bin << ", " << dilep_cut << ")" << endl;
  }
  dilep_bin = bin;
  dilep_bin.Name("dilep_"+dilep_bin.Name());
  dilep_cut = baseline_;
  dilep_bin.Cut().Replace("(nels+nmus)==1", "(nels+nmus)==2");
  dilep_bin.Cut().Replace("(nmus+nels)==1", "(nmus+nels)==2");
  dilep_bin.Cut().Replace("nels+nmus==1", "nels+nmus==2");
  dilep_bin.Cut().Replace("nmus+nels==1", "nmus+nels==2");
  dilep_bin.Cut().Replace("nleps==1", "nleps==2");
  dilep_bin.Cut().RmCutOn("nbm", "nbm>=1&&nbm<=2");
  dilep_bin.Cut().RmCutOn("met", "met>200&&met<=400");
  dilep_cut.Replace("(nels+nmus)==1", "(nels+nmus)==2");
  dilep_cut.Replace("(nmus+nels)==1", "(nmus+nels)==2");
  dilep_cut.Replace("nels+nmus==1", "nels+nmus==2");
  dilep_cut.Replace("nmus+nels==1", "nmus+nels==2");
  dilep_cut.Replace("nleps==1", "nleps==2");
  dilep_cut.RmCutOn("nbm", "nbm>=1&&nbm<=2");
  dilep_cut.RmCutOn( "met", "met>200&&met<=400");
}

void WorkspaceGenerator::AddSystematicsGenerators(){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddSystematicsGenerators()" << endl;
  }
  for(const auto &block: blocks_){
    for(const auto &vbin: block.Bins()){
      for(const auto &bin: vbin){
        for(const auto &syst: bin.Systematics()){
          AddSystematicGenerator(syst.Name());
          string full_name = syst.Name()+"_BLK_"+block.Name()+"_BIN_"+bin.Name();
          ostringstream oss;
          oss << "strength_" << full_name << "[" << syst.Strength() << "]" << flush;
          w_.factory(oss.str().c_str());
	  oss.str("");
          oss << "expr::" << full_name
              << "('exp(strength_" << full_name << "*" << syst.Name() << ")',"
              << "strength_" << full_name << "," << syst.Name() << ")" << flush;
          w_.factory(oss.str().c_str());
        }
      }
    }
  }

  auto all_prcs = backgrounds_;
  Append(all_prcs, signal_);
  for(const auto &bkg: all_prcs){
    for(const auto &syst: bkg.Systematics()){
      AddSystematicGenerator(syst.Name());
      string full_name = syst.Name()+"_PRC_"+bkg.Name();
      ostringstream oss;
      oss << "strength_" << full_name << "[" << syst.Strength() << "]" << flush;
      w_.factory(oss.str().c_str());
      oss.str("");
      oss << "expr::" << full_name
	  << "('exp(strength_" << full_name << "*" << syst.Name() << ")',"
	  << "strength_" << full_name << "," << syst.Name() << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }

  for(const auto &syst: free_systematics_){
    AddSystematicGenerator(syst.Name());
    for(const auto &block: blocks_){
      for(const auto &vbin: block.Bins()){
	for(const auto &bin: vbin){
	  for(const auto &prc: all_prcs){
	    if(!syst.HasEntry(bin, prc)) continue;
	    string full_name = syst.Name()+"_BIN_"+bin.Name()+"_PRC_"+prc.Name();
	    ostringstream oss;
	    oss << "strength_" << full_name << "[" << syst.Strength(bin, prc) << "]" << flush;
	    w_.factory(oss.str().c_str());
	    oss.str("");
	    oss << "expr::" << full_name
		<< "('exp(strength_" << full_name << "*" << syst.Name() << ")" << flush;
	    w_.factory(oss.str().c_str());
	  }
	}
      }
    }
  }
}

void WorkspaceGenerator::AddSystematicGenerator(const string &name){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddSystematicGenerator(" << name << ")" << endl;
  }
  if(systematics_.find(name) != systematics_.end()) return;
  w_.factory(("RooGaussian::constraint_"+name+"("+name+"[0.,-10.,10.],0.,1.)").c_str());
  Append(nuisances_, name);
  Append(systematics_, name);
}

void WorkspaceGenerator::AddData(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddData(" << block << ")" << endl;
  }
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      GammaParams gps(0., 0.);
      if(blind_level_ == BlindLevel::unblinded){
        gps = GetYield(bin, data_);
      }else{
        for(const auto &bkg: backgrounds_){
	  gps += GetYield(bin, bkg);
        }
      }
      ostringstream oss;
      oss << "nobs_BLK_" << block.Name()
          << "_BIN_" << bin.Name() << flush;
      Append(observables_, oss.str());
      oss << "[" << gps.Yield() << "]" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddBackgroundFractions(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddBackgroundFractions(" << block << ")" << endl;
  }
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      for(const auto &bkg: backgrounds_){
	string var = "expr::frac_BIN_"+bin.Name()+"_PRC_"+bkg.Name()
	  +"('ymc_BLK_"+block.Name()+"_BIN_"+bin.Name()+"_PRC_"+bkg.Name()
	  +"/ymc_BLK_"+block.Name()+"_BIN_"+bin.Name()
	  +"',ymc_BLK_"+block.Name()+"_BIN_"+bin.Name()+"_PRC_"+bkg.Name()
	  +",ymc_BLK_"+block.Name()+"_BIN_"+bin.Name()+")";
	w_.factory(var.c_str());	
      }
    }
  }
}

map<Process, double> WorkspaceGenerator::GetBackgroundFractions(const Block &block) const{
  if(print_level_ >= PrintLevel::everything){
    cout << "GetBackgroundFractions(" << block << ")" << endl;
  }
  map<Process, double> output;
  for(const auto &bkg: backgrounds_){
    for(const auto &vbin: block.Bins()){
      for(const auto &bin: vbin){
        output[bkg] += GetYield(bin, bkg).Yield();
      }
    }
  }
  double scale = 0.;
  for(const auto &process_value: output) scale += process_value.second;
  scale = 1./scale;
  for(auto &process_value: output) process_value.second *= scale;

  return output;
}

void WorkspaceGenerator::AddABCDParameters(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddABCDParameters(" << block << ")" << endl;
  }
  BlockYields by(block, backgrounds_, baseline_, yields_);

  ostringstream rxss, ryss;
  rxss << "sum::rxnorm_BLK_" << block.Name() << "(1.,";
  ryss << "sum::rynorm_BLK_" << block.Name() << "(1.,";
  ostringstream oss;
  oss << "norm_BLK_" << block.Name() << flush;
  Append(nuisances_, oss.str());
  oss << "[" << max(1., by.Total().Yield()) << ",0.,"
      << max(5.*by.Total().Yield(), 20.) << "]" << flush;
  w_.factory(oss.str().c_str());
  for(size_t irow = 0; irow < by.RowSums().size(); ++irow){
    if(irow == by.MaxRow()) continue;
    oss.str("");
    oss << "ry" << (irow+1) << (by.MaxRow()+1) << "_BLK_" << block.Name() << flush;
    ryss << "," << oss.str();
    Append(nuisances_, oss.str());
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
    Append(nuisances_, oss.str());
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
  if(print_level_ >= PrintLevel::everything){
    cout << "AddRawBackgroundPredictions(" << block << ")" << endl;
  }
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
	Append(prod_list, prod_name);
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
        factory_string += (",frac_BIN_"+bin.Name()+"_PRC_"+bkg.Name());
	if(do_systematics_){
	  for(const auto &syst: bkg.Systematics()){
	    factory_string += (","+syst.Name()+"_PRC_"+bkg.Name());
	  }
	  for(const auto &syst: free_systematics_){
	    if(syst.HasEntry(bin, bkg)){
	      factory_string += (","+syst.Name()+"_BIN_"+bin.Name()+"_PRC_"+bkg.Name());
	    }
	  }
	}
	factory_string += ")";
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

void WorkspaceGenerator::AddKappas(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddKappas(" << block << ")" << endl;
  }
  AddMCRowSums(block);
  AddMCColSums(block);
  AddMCTotal(block);
  AddMCPrediction(block);
  AddMCKappa(block);
}

void WorkspaceGenerator::AddMCYields(const Block & block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddMCYields(" << block << ")" << endl;
  }
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      ostringstream oss;
      for(const auto &bkg: backgrounds_){
	GammaParams gp = GetYield(bin, bkg);
	string bbp_name = bb_name + "_PRC_"+bkg.Name();
	oss.str("");
	oss << "nobsmc_" << bbp_name << flush;
	Append(observables_, oss.str());
	oss << "[" << gp.NEffective() << "]" << flush;
	w_.factory(oss.str().c_str());
	oss.str("");
	oss << "nmc_" << bbp_name << flush;
	Append(nuisances_, oss.str());
	oss << "[" << gp.NEffective()
	    << ",0.," << max(5.*gp.NEffective(), 20.) << "]" << flush;
	w_.factory(oss.str().c_str());
	oss.str("");
	oss << "wmc_" << bbp_name << "[" << gp.Weight() << "]" << flush;
	w_.factory(oss.str().c_str());
	oss.str("");
	oss << "prod::ymc_" << bbp_name
	    << "(nmc_" << bbp_name
	    << ",wmc_" << bbp_name << ")" << flush;
	w_.factory(oss.str().c_str());
      }
      oss.str("");
    }
  }
}

void WorkspaceGenerator::AddMCPdfs(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddKappaPdfs(" << block << ")" << endl;
  }
  bool first = true;
  string factory_string = "PROD::pdf_mc_"+block.Name()+"(";
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      for(const auto &bkg: backgrounds_){
	string bbp_name = "BLK_"+block.Name()+"_BIN_"+bin.Name()+"_PRC_"+bkg.Name();
	w_.factory(("RooPoisson::pdf_mc_"+bbp_name
		    +"(nobsmc_"+bbp_name
		    +",nmc_"+bbp_name+")").c_str());
	(static_cast<RooPoisson*>(w_.pdf(("pdf_mc_"+bbp_name).c_str())))->setNoRounding();
	if(first) first = false;
	else factory_string += ",";
	factory_string += "pdf_mc_"+bbp_name;
      }
    }
  }
  factory_string += ")";
  w_.factory(factory_string.c_str());
}

void WorkspaceGenerator::AddMCProcessSums(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddMCProcessSums(" << block << ")" << endl;
  }
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      ostringstream oss;
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      oss << "sum::ymc_" << bb_name << "(";
      if(backgrounds_.size() > 0){
	auto bkg = backgrounds_.cbegin();
	string name = "ymc_"+bb_name+"_PRC_"+bkg->Name();
	oss << name;
	for(++bkg; bkg != backgrounds_.cend(); ++bkg){
	  name = "ymc_"+bb_name+"_PRC_"+bkg->Name();
	  oss << "," << name;
	}
      }
      oss << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddMCRowSums(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddMCRowSums(" << block << ")" << endl;
  }
  for(size_t irow = 0; irow < block.Bins().size(); ++irow){
    ostringstream oss;
    oss << "sum::rowmc" << (irow+1) << "_BLK_" << block.Name() << "(";
    auto bin = block.Bins().at(irow).cbegin();
    oss << "ymc_BLK_" << block.Name() << "_BIN_" << bin->Name();
    for(++bin; bin != block.Bins().at(irow).cend(); ++bin){
      oss << ",ymc_BLK_" << block.Name() << "_BIN_" << bin->Name();
    }
    oss << ")" << flush;
    w_.factory(oss.str().c_str());
  }
}

void WorkspaceGenerator::AddMCColSums(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddMCColSums(" << block << ")" << endl;
  }
  if(block.Bins().size() > 0 && block.Bins().at(0).size() > 0){
    for(size_t icol = 0; icol < block.Bins().at(0).size(); ++icol){
      ostringstream oss;
      oss << "sum::colmc" << (icol+1) << "_BLK_" << block.Name() << "(";
      size_t irow = 0;
      oss << "ymc_BLK_" << block.Name() << "_BIN_" << block.Bins().at(irow).at(icol).Name();
      for(irow = 1; irow < block.Bins().size(); ++irow){
	if(icol < block.Bins().at(irow).size()){
	  oss << ",ymc_BLK_" << block.Name() << "_BIN_" << block.Bins().at(irow).at(icol).Name();
	}
      }
      oss << ")" << flush;;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddMCTotal(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddMCTotalSums(" << block << ")" << endl;
  }
  ostringstream oss;
  oss << "sum::totmc_BLK_" << block.Name() << "(";
  if(block.Bins().size() > 0){
    oss << "rowmc1_BLK_" << block.Name();
    for(size_t irow = 1; irow < block.Bins().size(); ++irow){
      oss << ",rowmc" << (irow+1) << "_BLK_" << block.Name();
    }
  }
  oss << ")" << flush;
  w_.factory(oss.str().c_str());
}

void WorkspaceGenerator::AddMCPrediction(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddMCPrediction(" << block << ")" << endl;
  }
  for(size_t irow = 0; irow < block.Bins().size(); ++irow){
    for(size_t icol = 0; icol < block.Bins().at(irow).size(); ++icol){
      const Bin &bin = block.Bins().at(irow).at(icol);
      ostringstream oss;
      oss << "expr::predmc_BLK_" << block.Name() << "_BIN_" << bin.Name() << "('"
	  << "(rowmc" << (irow+1) << "_BLK_" << block.Name()
	  << "*colmc" << (icol+1) << "_BLK_" << block.Name()
	  << ")/totmc_BLK_" << block.Name()
	  << "',rowmc" << (irow+1) << "_BLK_" << block.Name()
	  << ",colmc" << (icol+1) << "_BLK_" << block.Name()
	  << ",totmc_BLK_" << block.Name() << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddMCKappa(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddMCPrediction(" << block << ")" << endl;
  }
  for(size_t irow = 0; irow < block.Bins().size(); ++irow){
    for(size_t icol = 0; icol < block.Bins().at(irow).size(); ++icol){
      const Bin &bin = block.Bins().at(irow).at(icol);
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      ostringstream oss;
      oss << "expr::kappamc_" << bb_name << "('"
	  << "ymc_" << bb_name
	  << "/predmc_" << bb_name
	  << "',ymc_" << bb_name
	  << ",predmc_" << bb_name
	  << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddFullBackgroundPredictions(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddFullBackgroundPredictions(" << block << ")" << endl;
  }
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      ostringstream oss;
      oss << "prod::nbkg_" << bb_name << "("
          << "nbkg_raw_" << bb_name;
      if(do_systematics_ || do_dilepton_){
	for(const auto &syst: bin.Systematics()){
	  oss << "," << syst.Name() << "_" << bb_name;
	}
      }
      if(do_mc_kappa_correction_){
	oss << ",kappamc_" << bb_name;
      }
      oss << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddSignalPredictions(const Block &block){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddSignalPredictions(" << block << ")" << endl;
  }
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      double yield = GetYield(bin, signal_).Yield();
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
  if(print_level_ >= PrintLevel::everything){
    cout << "AddSignalPredictions(" << block << ")" << endl;
  }
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
  if(print_level_ >= PrintLevel::everything){
    cout << "AddFullPdf()" << endl;
  }
  if(blocks_.size() == 0){
    w_.factory("RooPoisson::model_b(0,0)");
    w_.factory("RooPoisson::model_s(0,0)");
  }else{
    string null_list = "pdf_null_BLK_"+blocks_.cbegin()->Name();
    string alt_list = "pdf_alt_BLK_"+blocks_.cbegin()->Name();
    auto blockp = blocks_.cbegin();
    for(++blockp; blockp != blocks_.cend(); ++blockp){
      null_list += (",pdf_null_BLK_"+blockp->Name());
      alt_list += (",pdf_alt_BLK_"+blockp->Name());
    }
    if(do_systematics_ || do_dilepton_){
      for(const auto &syst: systematics_){
        null_list += (",constraint_"+syst);
        alt_list += (",constraint_"+syst);
      }
    }
    if(do_mc_kappa_correction_){
      for(const auto & block: blocks_){
	null_list += (",pdf_mc_"+block.Name());
	alt_list += (",pdf_mc_"+block.Name());
      }
    }
    w_.factory(("PROD::model_b("+null_list+")").c_str());
    w_.factory(("PROD::model_s("+alt_list+")").c_str());
  }
}

void WorkspaceGenerator::AddParameterSets(){
  if(print_level_ >= PrintLevel::everything){
    cout << "AddParameterSets()" << endl;
  }
  DefineParameterSet("POI", poi_);
  DefineParameterSet("nuisances", nuisances_);
  DefineParameterSet("observables", observables_);
  RooDataSet data_obs{"data_obs", "data_obs", *w_.set("observables")};
  data_obs.add(*w_.set("observables"));
  w_.import(data_obs);
}

void WorkspaceGenerator::DefineParameterSet(const string &set_name,
                                            const set<string> &var_names){
  if(print_level_ >= PrintLevel::everything){
    cout << "DefineParameterSet(" << set_name << ",[var_names])" << endl;
  }
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
  if(print_level_ >= PrintLevel::everything){
    cout << "AddModels()" << endl;
  }
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
  for(const auto &block: wg.blocks_){
    for(const auto &vbin: block.Bins()){
      for(const auto &bin: vbin){
        wg.PrintComparison(stream, bin, wg.data_, block, true);
        wg.PrintComparison(stream, bin, wg.signal_, block, false);
        for(const auto &bkg: wg.backgrounds_){
          wg.PrintComparison(stream, bin, bkg, block, false);
        }
      }
    }
  }
  return stream;
}

void WorkspaceGenerator::PrintComparison(ostream &stream, const Bin &bin, const Process &process,
					 const Block &block, bool is_data) const{
  if(print_level_ >= PrintLevel::everything){
    cout << "PrintComparison([stream], " << bin << ", " << process
	 << ", " << block << ", " << is_data << ")" << endl;
  }
  GammaParams gp(0., 0.);
  if(!is_data || blind_level_ == BlindLevel::unblinded){
    gp = GetYield(bin, process);
  }

  ostringstream name;
  name << (is_data ? "nobs" : "rate")
       << "_BLK_" << block.Name()
       << "_BIN_" << bin.Name();
  if(!is_data){
    name << "_PRC_" << process.Name();
  }
  name << flush;

  stream << setw(64) << name.str() << ": "
         << setw(8) << gp.Yield()
         << " +- " << setw(8) << gp.CorrectedUncertainty()
         << " => ";
  RooAbsReal *fp = w_.function(name.str().c_str());
  if(fp){
    stream << setw(8) << fp->getVal() << endl;
  }else{
    stream << "Not found" << endl;
  }
}
