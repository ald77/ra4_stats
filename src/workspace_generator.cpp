#include "workspace_generator.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <array>
#include <algorithm>
#include <numeric>

#include "TDirectory.h"

#include "RooPoisson.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"

#include "utilities.hpp"

using namespace std;

YieldManager WorkspaceGenerator::yields_ = YieldManager(4.);
mt19937_64 WorkspaceGenerator::prng_ = WorkspaceGenerator::InitializePRNG();
poisson_distribution<> WorkspaceGenerator::dist_(1.);

WorkspaceGenerator::WorkspaceGenerator(const Cut &baseline,
                                       const set<Block> &blocks,
                                       const set<Process> &backgrounds,
                                       const Process &signal,
                                       const Process &data,
                                       const string &systematics_file,
                                       const bool use_r4, const double sig_strength, const double sig_xsec_f):
  baseline_(baseline),
  backgrounds_(backgrounds),
  signal_(signal),
  data_(data),
  blocks_(blocks),
  obs_vals_(),
  obs_gens_(),
  systematics_file_(systematics_file),
  use_r4_(use_r4),
  sig_strength_(sig_strength),
  sig_xsec_f_(sig_xsec_f),
  rmax_(20.),
  w_("w"),
  poi_(),
  observables_(),
  nuisances_(),
  systematics_(),
  free_systematics_(),
  luminosity_(4.),
  print_level_(PrintLevel::important),
  do_systematics_(true),
  do_dilepton_(false),
  do_mc_kappa_correction_(true),
  num_toys_(0),
  w_is_valid_(false){
  w_.cd();
}

void WorkspaceGenerator::WriteToFile(const string &file_name){
  if(print_level_ >= PrintLevel::everything) DBG(file_name);
  if(!w_is_valid_) UpdateWorkspace();
  w_.writeToFile(file_name.c_str());
  if(print_level_ >= PrintLevel::everything){
    DBG("");
    w_.Print();
  }
  if(print_level_ >= PrintLevel::normal){
    cout << *this << endl;
  }
  if(print_level_ >= PrintLevel::important){
    cout << endl << "Wrote workspace to file " << file_name << endl<< endl;
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

double WorkspaceGenerator::GetRMax() const{
  return rmax_;
}

WorkspaceGenerator & WorkspaceGenerator::SetRMax(double rmax){
  if(rmax != rmax_){
    rmax_ = rmax;
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

size_t WorkspaceGenerator::AddToys(size_t num_toys){
  if(print_level_ >= PrintLevel::everything) DBG(num_toys);
  if(!w_is_valid_) UpdateWorkspace();
  if(num_toys == 0) return num_toys_;
  const RooArgSet *obs_orig = w_.set("observables");
  if(obs_orig == nullptr) ERROR("Could not get observables list for toy generation");
  RooArgSet obs(*obs_orig);
  SetupToys(obs);
  for(size_t itoy = num_toys_; itoy < num_toys_+num_toys; ++itoy){
    GenerateToys(obs);
    RooDataSet data_new(("data_obs_"+to_string(itoy)).c_str(), ("data_obs_"+to_string(itoy)).c_str(), obs);
    data_new.add(obs);
    w_.import(data_new);
  }
  ResetToys(obs);
  num_toys_ += num_toys;
  return num_toys_;
}

void WorkspaceGenerator::SetupToys(const RooArgSet &obs){
  if(print_level_ >= PrintLevel::everything) DBG("");
  obs_vals_.clear();
  obs_gens_.clear();
  TIterator *iter_ptr = obs.createIterator();
  if(iter_ptr == nullptr) ERROR("Could not generator iterator to set up toys");
  for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
    RooAbsReal *arg = static_cast<RooAbsReal*>(*(*iter_ptr));
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(Contains(name,"nobsmc")) continue;
    double value = arg->getVal();
    obs_vals_[name] = value;
    obs_gens_[name] = poisson_distribution<>(value);
  }
  iter_ptr->Reset();
}

void WorkspaceGenerator::GenerateToys(RooArgSet &obs){
  if(print_level_ >= PrintLevel::everything) DBG("");
  TIterator *iter_ptr = obs.createIterator();
  if(iter_ptr == nullptr) ERROR("Could not generator iterator to set up toys");
  for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
    RooRealVar *arg = static_cast<RooRealVar*>(*(*iter_ptr));
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(Contains(name,"nobsmc")) continue;
    arg->setVal(obs_gens_[name](prng_));
  }
  iter_ptr->Reset();
}

void WorkspaceGenerator::ResetToys(RooArgSet &obs){
  if(print_level_ >= PrintLevel::everything) DBG("");
  TIterator *iter_ptr = obs.createIterator();
  if(iter_ptr == nullptr) ERROR("Could not generator iterator to set up toys");
  for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
    RooRealVar *arg = static_cast<RooRealVar*>(*(*iter_ptr));
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(Contains(name,"nobsmc")) continue;
    double value = obs_vals_.at(name);
    arg->setVal(value);
  }
  iter_ptr->Reset();
}

mt19937_64 WorkspaceGenerator::InitializePRNG(){
  array<int, 128> sd;
  random_device r;
  generate_n(sd.begin(), sd.size(), ref(r));
  seed_seq ss(begin(sd), end(sd));
  return mt19937_64(ss);
}

int WorkspaceGenerator::GetPoisson(double rate){
  using param_t = decltype(dist_)::param_type;
  dist_.param(param_t{rate});
  return dist_(prng_);
}

void WorkspaceGenerator::UpdateWorkspace(){
  if(print_level_ >= PrintLevel::everything) DBG("");
  string old_name = w_.GetName();
  w_.Delete();
  gDirectory->Delete(old_name.c_str());
  w_.Clear();
  w_ = RooWorkspace("w");
  w_.SetName("w");
  w_.cd();

  if(do_dilepton_){
    AddDileptonSystematic();
  }
  if(do_systematics_){
    ReadSystematicsFile();
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

  AddDummyNuisance();
  AddFullPdf();
  AddParameterSets();
  AddModels();

  w_is_valid_ = true;
}

void WorkspaceGenerator::AddPOI(){
  if(print_level_ >= PrintLevel::everything) DBG("");
  w_.factory(("r[1.,0.,"+to_string(rmax_)+"]").c_str());
  Append(poi_, "r");
}

void WorkspaceGenerator::ReadSystematicsFile(){
  if(systematics_file_ == "") return;
  ifstream file(systematics_file_);
  vector<string> lines;
  string one_line;
  while(getline(file, one_line)){
    CleanLine(one_line);
    if(one_line != ""){
      lines.push_back(one_line);
    }
  }

  vector<vector<string> > words;
  for(const auto &line: lines){
    vector<string> these_words = Tokenize(line);
    words.push_back(these_words);
  }

  auto all_prc = backgrounds_;
  Append(all_prc, signal_);
  string syst_name;
  auto process_list = all_prc;
  process_list.clear();

  free_systematics_.clear();
  FreeSystematic this_systematic("BADBADBADBADBADBADBADBADBAD");
  bool ready = false;
  for(const auto &line: words){
    if(line.size() < 2){
      string out;
      for(const auto &word: line) out += word;
      ERROR("Bad systematics line: "+out);
    }
    if(line.at(0) == "SYSTEMATIC"){
      if(ready){
        Append(free_systematics_, this_systematic);
      }
      this_systematic = FreeSystematic(line.at(1));
      ready = true;
    }else if(line[0] == "PROCESSES"){
      process_list.clear();
      for(size_t iword = 1; iword < line.size(); ++iword){
        bool found = false;
        vector<string> names = Tokenize(line.at(iword), ", ");
        for(const auto &name: names){
          for(const auto &prc: all_prc){
            if(name == prc.Name()){
              Append(process_list, prc);
              found = true;
            }
          }
          if(!found){
            ERROR("Systematic "+this_systematic.Name()
                                +" could not be applied to process "+line.at(iword));
          }
        }
      }
    }else{
      bool found = false;
      string clean_line(line.at(0));
      ReplaceAll(clean_line, " ", "");
      ReplaceAll(clean_line, "\t", "");
      for(const auto &block: blocks_){
        for(const auto &vbin: block.Bins()){
          for(const auto &bin: vbin){
            if(bin.Name() != clean_line) continue;
            for(const auto &prc: process_list){
	      float syst(atof(line.at(1).c_str()));
              if(isnan(syst)){
		DBG("Systematic " << this_systematic.Name() << " is NaN for bin "
		    << bin.Name() << ", process " << prc.Name());
		syst = 0.;
	      }
              if(isinf(syst)){
		DBG("Systematic " << this_systematic.Name() << " is infinite for bin "
		    << bin.Name() << ", process " << prc.Name());
		if(syst>0.){
		  syst = 1.;
		}else{
		  syst = -1.;
		}
	      }
              if(syst>=0) this_systematic.Strength(bin, prc) = log(1+syst);
	      else this_systematic.Strength(bin, prc) = -log(1+fabs(syst));
              found = true;
            }

          }
        }
      }
      if(!found){
        ERROR("Systematic "+this_systematic.Name()
                            +" could not be applied to bin "+line.at(0));
      }
    }
  }
  if(ready){
    Append(free_systematics_, this_systematic);
  }
}

void WorkspaceGenerator::CleanLine(string &line){
  ReplaceAll(line, "=", "");
  string old = line;
  do{
    old = line;
    ReplaceAll(line, "  ", " ");
  } while (line != old);
  while(line.size() > 0 && line[0] == ' '){
    line = line.substr(1);
  }
  if(line.size() > 0 && line[0] == '#'){
    line = "";
  }
}

void WorkspaceGenerator::AddDileptonSystematic(){
  if(print_level_ >= PrintLevel::everything) DBG("");

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
        if(dilep_bin.Blind()){
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
  if(print_level_ >= PrintLevel::everything) DBG("");
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
    DBG(bin << ", " << dilep_bin << ", " << dilep_cut);
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
  dilep_bin.Cut().RmCutOn("mt");
  dilep_cut.Replace("(nels+nmus)==1", "(nels+nmus)==2");
  dilep_cut.Replace("(nmus+nels)==1", "(nmus+nels)==2");
  dilep_cut.Replace("nels+nmus==1", "nels+nmus==2");
  dilep_cut.Replace("nmus+nels==1", "nmus+nels==2");
  dilep_cut.Replace("nleps==1", "nleps==2");
  dilep_cut.RmCutOn("nbm", "nbm>=1&&nbm<=2");
  dilep_cut.RmCutOn("met", "met>200&&met<=400");
  dilep_cut.RmCutOn("mt");
}

void WorkspaceGenerator::AddSystematicsGenerators(){
  if(print_level_ >= PrintLevel::everything) DBG("");
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
                << "('exp(strength_" << full_name << "*" << syst.Name() << ")',"
                << "strength_" << full_name << "," << syst.Name() << ")" << flush;
            w_.factory(oss.str().c_str());
          }
        }
      }
    }
  }
}

void WorkspaceGenerator::AddSystematicGenerator(const string &name){
  if(print_level_ >= PrintLevel::everything) DBG(name);
  if(systematics_.find(name) != systematics_.end()) return;
  w_.factory(("RooGaussian::constraint_"+name+"("+name+"[0.,-10.,10.],0.,1.)").c_str());
  Append(nuisances_, name);
  Append(systematics_, name);
}

void WorkspaceGenerator::AddData(const Block &block){
  if(print_level_ >= PrintLevel::everything) DBG(block);
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      GammaParams gps(0., 0.);
      if(bin.Blind()){
        for(const auto &bkg: backgrounds_){
          gps += GetYield(bin, bkg);
        }
        // Injecting signal
        gps += sig_strength_*GetYield(bin, signal_);
      }else{
        gps = GetYield(bin, data_);
      }

      ostringstream oss;
      oss << "nobs_BLK_" << block.Name()
          << "_BIN_" << bin.Name() << flush;
      if(use_r4_ || !Contains(bin.Name(), "4")){
        Append(observables_, oss.str());
      }
      oss << "[" << gps.Yield() << "]" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddBackgroundFractions(const Block &block){
  if(print_level_ >= PrintLevel::everything) DBG(block);
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

void WorkspaceGenerator::AddABCDParameters(const Block &block){
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
  AddMCRowSums(block);
  AddMCColSums(block);
  AddMCTotal(block);
  AddMCPrediction(block);
  AddMCKappa(block);
}

void WorkspaceGenerator::AddMCYields(const Block & block){
  if(print_level_ >= PrintLevel::everything) DBG(block);
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      ostringstream oss;
      auto all_prcs = backgrounds_;
      Append(all_prcs, signal_);
      for(const auto &bkg: all_prcs){
        GammaParams gp = GetYield(bin, bkg);
	if(Contains(bkg.Name(), "sig")) gp *= sig_xsec_f_;
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
  bool first = true;
  string factory_string = "PROD::pdf_mc_"+block.Name()+"(";
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      auto all_prcs = backgrounds_;
      Append(all_prcs, signal_);
      for(const auto &bkg: all_prcs){
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      string bb_name = "BLK_"+block.Name()+"_BIN_"+bin.Name();
      ostringstream oss;
      oss << "prod::nbkg_" << bb_name << "("
          << "nbkg_raw_" << bb_name;
      for(const auto &syst: bin.Systematics()){
        if(do_systematics_ && syst.Name().substr(0,6) != "dilep_"){
          oss << "," << syst.Name() << "_" << bb_name;
        }
      }
      if(do_systematics_){
        for(const auto &prc: backgrounds_){
          for(const auto &syst: prc.Systematics()){
            oss << "," << syst.Name() << "_BIN_" << bin.Name() << "_PRC_" << prc.Name();
          }
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
  if(print_level_ >= PrintLevel::everything) DBG(block);
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      ostringstream oss;
      oss.str("");
      oss << "prod::nsig_BLK_" << block.Name()
          << "_BIN_" << bin.Name()
          << "(r,"
          << "ymc_BLK_" << block.Name()
          << "_BIN_" << bin.Name()
          << "_PRC_" << signal_.Name();
      if(do_systematics_){
        for(const auto &syst: signal_.Systematics()){
          oss << "," << syst.Name() << "_BIN_" << bin.Name() << "_PRC_" << signal_.Name();
        }
        for(const auto &syst: free_systematics_){
          if(!syst.HasEntry(bin, signal_)) continue;
          oss << "," << syst.Name() << "_BIN_" << bin.Name() << "_PRC_" << signal_.Name();
        }
      }
      oss << ")" << flush;
      w_.factory(oss.str().c_str());
    }
  }
}

void WorkspaceGenerator::AddPdfs(const Block &block){
  if(print_level_ >= PrintLevel::everything) DBG(block);
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
      w_.factory(("sum::nexp"+bb_name+"(nbkg"+bb_name+",nsig"+bb_name+")").c_str());
      if(use_r4_ || !Contains(bb_name, "4")){
        null_list += null_name;
        alt_list += alt_name;
        w_.factory(("RooPoisson::pdf_null"+bb_name+"(nobs"+bb_name+",nbkg"+bb_name+")").c_str());
        (static_cast<RooPoisson*>(w_.pdf(null_name.c_str())))->setNoRounding();
        w_.factory(("RooPoisson::pdf_alt"+bb_name+"(nobs"+bb_name+",nexp"+bb_name+")").c_str());
        (static_cast<RooPoisson*>(w_.pdf(alt_name.c_str())))->setNoRounding();
        is_first = false;
      }
    }
  }
  w_.factory(("PROD:pdf_null_BLK_"+block.Name()+"("+null_list+")").c_str());
  w_.factory(("PROD:pdf_alt_BLK_"+block.Name()+"("+alt_list+")").c_str());
}

void WorkspaceGenerator::AddDummyNuisance(){
  w_.factory("RooGaussian::pdf_dummy_nuisance(dummy_nuisance[0.,-10.,10.],0.,1.)");
  Append(nuisances_, "dummy_nuisance");
}

void WorkspaceGenerator::AddFullPdf(){
  if(print_level_ >= PrintLevel::everything) DBG("");
  string null_list = "pdf_dummy_nuisance";
  string alt_list = "pdf_dummy_nuisance";
  for(const auto &block: blocks_){
    null_list += (",pdf_null_BLK_"+block.Name());
    alt_list += (",pdf_alt_BLK_"+block.Name());
  }
  if(do_systematics_ || do_dilepton_){
    for(const auto &syst: systematics_){
      null_list += (",constraint_"+syst);
      alt_list += (",constraint_"+syst);
    }
  }
  for(const auto & block: blocks_){
    null_list += (",pdf_mc_"+block.Name());
    alt_list += (",pdf_mc_"+block.Name());
  }
  w_.factory(("PROD::model_b("+null_list+")").c_str());
  w_.factory(("PROD::model_s("+alt_list+")").c_str());
}

void WorkspaceGenerator::AddParameterSets(){
  if(print_level_ >= PrintLevel::everything) DBG("");
  DefineParameterSet("POI", poi_);
  DefineParameterSet("nuisances", nuisances_);
  DefineParameterSet("observables", observables_);
  DefineParameterSet("globalObservables", set<string>());
  auto xxx = w_.set("observables");
  RooDataSet data_obs{"data_obs", "data_obs", *xxx};
  data_obs.add(*w_.set("observables"));
  w_.import(data_obs);
}

void WorkspaceGenerator::DefineParameterSet(const string &set_name,
                                            const set<string> &var_names){
  if(print_level_ >= PrintLevel::everything) DBG(set_name);
  if(var_names.size()==0){
    w_.defineSet(set_name.c_str(), "");
  }else{
    w_.defineSet(set_name.c_str(), var_names.cbegin()->c_str());
    for(auto name = ++var_names.cbegin();
	name != var_names.cend();
	++name){
      w_.extendSet(set_name.c_str(), name->c_str());
    }
  }
}

void WorkspaceGenerator::AddModels(){
  if(print_level_ >= PrintLevel::everything) DBG("");
  RooStats::ModelConfig model_config("ModelConfig", &w_);
  model_config.SetPdf(*w_.pdf("model_s"));
  model_config.SetParametersOfInterest(*w_.set("POI"));
  model_config.SetObservables(*w_.set("observables"));
  model_config.SetNuisanceParameters(*w_.set("nuisances"));
  model_config.SetGlobalObservables(*w_.set("globalObservables"));

  RooStats::ModelConfig model_config_bonly("ModelConfig_bonly", &w_);
  model_config_bonly.SetPdf(*w_.pdf("model_b"));
  model_config_bonly.SetParametersOfInterest(*w_.set("POI"));
  model_config_bonly.SetObservables(*w_.set("observables"));
  model_config_bonly.SetNuisanceParameters(*w_.set("nuisances"));
  model_config_bonly.SetGlobalObservables(*w_.set("globalObservables"));

  w_.import(model_config);
  w_.import(model_config_bonly);
}

ostream & operator<<(ostream& stream, const WorkspaceGenerator &wg){
  for(const auto &block: wg.blocks_){
    for(const auto &vbin: block.Bins()){
      for(const auto &bin: vbin){
        wg.PrintComparison(stream, bin, wg.data_, block);
        wg.PrintComparison(stream, bin, wg.signal_, block);
        for(const auto &bkg: wg.backgrounds_){
          wg.PrintComparison(stream, bin, bkg, block);
        }
      }
    }
  }
  return stream;
}

void WorkspaceGenerator::PrintComparison(ostream &stream, const Bin &bin,
                                         const Process &process, const Block &block) const{
  if(print_level_ >= PrintLevel::everything){
    DBG(bin << ", " << process << ", " << block);
  }

  GammaParams gp(0., 0.);
  if(!(process.IsData() && bin.Blind())){
    gp = GetYield(bin, process);
  }

  ostringstream name;
  name << (process.IsData() ? "nobs" : process.IsSignal() ? "ymc" : "rate")
       << "_BLK_" << block.Name()
       << "_BIN_" << bin.Name();
  if(!process.IsData()){
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
