#include "process.hpp"

#include <string>
#include <set>
#include <initializer_list>

#include "TChain.h"

#include "utilities.hpp"

using namespace std;

Process::Process(const string &name,
                 const set<string> &file_names,
                 const class Cut &cut,
                 bool count_zeros):
  file_names_(file_names),
  chain_(make_shared<TChain>("tree", "tree")),
  cut_(cut),
  name_(name),
  count_zeros_(count_zeros){
  CleanName();
  AddFiles();
  }

Process::Process(const string &name,
                 initializer_list<string> file_names,
                 const class Cut &cut,
                 bool count_zeros):
  file_names_(file_names),
  chain_(make_shared<TChain>("tree", "tree")),
  cut_(cut),
  name_(name),
  count_zeros_(count_zeros){
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

const bool & Process::CountZeros() const{
  return count_zeros_;
}

bool & Process::CountZeros(){
  return count_zeros_;
}

bool Process::operator<(const Process &p) const{
  return ComparisonTuple() < p.ComparisonTuple();
}

void Process::CleanName(){
  ReplaceAll(name_, " ", "");
}

void Process::AddFiles(){
  for(const auto &file_name: file_names_){
    chain_->Add(file_name.c_str());
  }
}
