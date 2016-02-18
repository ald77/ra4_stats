#include "process.hpp"

#include <string>
#include <set>
#include <initializer_list>
#include <algorithm>

#include "TChain.h"

#include "utilities.hpp"

using namespace std;

Process::Process(const string &name,
                 const set<string> &file_names,
                 const class Cut &cut,
                 bool is_data,
                 bool is_signal,
                 bool count_zeros,
                 const SystCollection &systematics):
  file_names_(file_names),
  chain_(make_shared<TChain>("tree", "tree")),
  cut_(cut),
  name_(name),
  is_data_(is_data),
  is_signal_(is_signal),
  count_zeros_(count_zeros),
  systematics_(systematics){
  CleanName();
  AddFiles();
  }

Process::Process(const string &name,
                 initializer_list<string> file_names,
                 const class Cut &cut,
                 bool is_data,
                 bool is_signal,
                 bool count_zeros,
                 const SystCollection &systematics):
  file_names_(file_names),
  chain_(make_shared<TChain>("tree", "tree")),
  cut_(cut),
  name_(name),
  is_data_(is_data),
  is_signal_(is_signal),
  count_zeros_(count_zeros),
  systematics_(systematics){
  CleanName();
  AddFiles();
  }

const string & Process::Name() const{
  return name_;
}

Process & Process::Name(const string &name){
  name_ = name;
  CleanName();
  return *this;
}

const class Cut & Process::Cut() const{
  return cut_;
}

class Cut & Process::Cut(){
  return cut_;
}

const set<string> & Process::FileNames() const{
  return file_names_;
}

long Process::GetEntries() const{
  return chain_->GetEntries();
}

GammaParams Process::GetYield(const class Cut &cut) const{
  double count, uncertainty;
  ::GetCountAndUncertainty(*chain_, cut*cut_, count, uncertainty);
  GammaParams gps;
  gps.SetYieldAndUncertainty(count, uncertainty);
  return gps;
}

const bool & Process::IsData() const{
  return is_data_;
}

bool & Process::IsData(){
  return is_data_;
}

const bool & Process::IsSignal() const{
  return is_signal_;
}

bool & Process::IsSignal(){
  return is_signal_;
}

const bool & Process::CountZeros() const{
  return count_zeros_;
}

bool & Process::CountZeros(){
  return count_zeros_;
}

const Process::SystCollection & Process::Systematics() const{
  return systematics_;
}

Process & Process::Systematics(const SystCollection &systematics){
  systematics_ = systematics;
  return *this;
}

Process & Process::AddSystematic(const Systematic &systematic){
  if(!HasSystematic(systematic)){
    systematics_.insert(systematics_.end(), systematic);
  }
  return *this;
}

Process & Process::AddSystematics(const SystCollection &systematics){
  for(const auto& systematic: systematics){
    AddSystematic(systematic);
  }
  return *this;
}

bool Process::HasSystematic(const Systematic &systematic) const{
  return find(systematics_.cbegin(), systematics_.cend(), systematic) != systematics_.cend();
}

Process & Process::RemoveSystematic(const Systematic &systematic){
  try{
    systematics_.erase(find(systematics_.begin(), systematics_.end(), systematic));
  }catch(const out_of_range &e){
    ERROR(string(e.what())+": bin "+name_+" does not contain systematic "+systematic.Name()+".");
  }
  return *this;
}

Process & Process::RemoveSystematics(){
  systematics_.clear();
  return *this;
}

Process & Process::SetSystematicStrength(const std::string &name, double strength){
  bool found_it = false;
  for(auto systematic = systematics_.cbegin(); systematic != systematics_.cend(); ++systematic){
    if(systematic->Name() == name){
      Systematic new_syst = *systematic;
      new_syst.Strength() = strength;
      found_it = true;
      systematics_.erase(systematic);
      systematics_.insert(systematics_.end(), new_syst);
    }
  }
  if(!found_it){
    ERROR("Process "+name_+" does not contain systematic "+name);
  }
  return *this;
}

bool Process::operator<(const Process &p) const{
  return tie(cut_, file_names_, count_zeros_, systematics_)
    < tie(p.cut_, p.file_names_, p.count_zeros_, p.systematics_);
}

bool Process::operator==(const Process &p) const{
  return tie(cut_, file_names_, count_zeros_, systematics_)
    == tie(p.cut_, p.file_names_, p.count_zeros_, p.systematics_);
}

void Process::CleanName(){
  ReplaceAll(name_, " ", "");
}

void Process::AddFiles(){
  for(const auto &file_name: file_names_){
    chain_->Add(file_name.c_str());
  }
}

ostream & operator<<(ostream &stream, const Process &proc){
  stream << "Process::" << proc.Name()
         << "(cut=" << proc.Cut()
         << ",count_zeros=" << proc.CountZeros()
         << ",is_data=" << proc.IsData()
         << ",is_signal=" << proc.IsSignal()
         << ")";
  return stream;
}
