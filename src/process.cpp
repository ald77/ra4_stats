#include "process.hpp"

#include <string>
#include <vector>
#include <initializer_list>

#include "TChain.h"

#include "utilities.hpp"

using namespace std;

Process::Process(const string &name,
                 const vector<string> &file_names,
                 const string &cut,
                 bool count_zeros):
  chain_("tree", "tree"),
  name_(name),
  cut_(cut),
  count_zeros_(count_zeros){
  CleanName();
  CleanCut();
  ReplaceAll(name_, " ", "");
  ReplaceAll(cut_, " ", "");
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
  CleanName();
  CleanCut();
  ReplaceAll(name_, " ", "");
  ReplaceAll(cut_, " ", "");
  for(auto file_name = file_names.begin();
      file_name != file_names.end();
      ++file_name){
    chain_.Add(file_name->c_str());
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

const string & Process::Cut() const{
  return cut_;
}

Process & Process::Cut(const string &cut){
  cut_ = cut;
  CleanCut();
  return *this;
}

long Process::GetEntries() const{
  return chain_.GetEntries();
}

void Process::GetCountAndUncertainty(double &count, double &uncertainty,
				     const std::string &cut) const{
  return ::GetCountAndUncertainty(chain_, "("+cut+")*("+cut_+")",
				  count, uncertainty);
}

bool Process::CountZeros() const{
  return count_zeros_;
}

Process & Process::CountZeros(bool count_zeros){
  count_zeros_ = count_zeros;
  return *this;
}

bool Process::operator<(const Process &p) const{
  return ComparisonTuple() < p.ComparisonTuple();
}

void Process::CleanName(){
  ReplaceAll(name_, " ", "");
}

void Process::CleanCut(){
  ReplaceAll(cut_, " ", "");
}
