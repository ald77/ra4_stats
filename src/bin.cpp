#include "bin.hpp"

#include <string>
#include <vector>

#include "systematic.hpp"
#include "utilities.hpp"

using namespace std;

Bin::Bin(const string &name, const string &cut,
	 const vector<Systematic> &systematics):
  name_(name),
  cut_(cut),
  systematics_(systematics){
  ReplaceAll(name_, " ", "");
  ReplaceAll(cut_, " ", "");
  }

bool Bin::operator<(const Bin &b) const{
  return name_ < b.name_
    || (name_ == b.name_
        && cut_ < b.cut_)
    || (name_ == b.name_
	&& cut_ == b.cut_
	&& systematics_ < b.systematics_);
}
