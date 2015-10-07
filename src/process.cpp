#include "process.hpp"

#include <string>
#include <vector>
#include <initializer_list>

#include "TChain.h"

#include "utilities.hpp"

using namespace std;

Process::Process(const string &name,
                 const vector<string> &file_names,
                 const class Cut &cut,
                 bool count_zeros):
  chain_(make_shared<TChain>("tree", "tree")),
  cut_(cut),
  name_(name),
  count_zeros_(count_zeros){
  CleanName();
  for(auto file_name = file_names.cbegin();
      file_name != file_names.cend();
      ++file_name){
    chain_->Add(file_name->c_str());
  }
  }

Process::Process(const string &name,
                 initializer_list<string> file_names,
                 const class Cut &cut,
                 bool count_zeros):
  chain_(make_shared<TChain>("tree","tree")),
  cut_(cut),
  name_(name),
  count_zeros_(count_zeros){
  CleanName();
  for(auto file_name = file_names.begin();
      file_name != file_names.end();
      ++file_name){
    chain_->Add(file_name->c_str());
  }
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
