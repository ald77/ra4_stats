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
  ReplaceAll(name_, " ", "");
  ReplaceAll(cut_, " ", "");
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
