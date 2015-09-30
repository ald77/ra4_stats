#include "systematic.hpp"

#include <string>

using namespace std;

Systematic::Systematic(const string &base_name,
		       double multiplier):
  base_name_(base_name),
  multiplier_(multiplier){
}

bool Systematic::operator<(const Systematic &s) const{
  return multiplier_ < s.multiplier_
    || (multiplier_ == s.multiplier_
	&& base_name_ == s.base_name_);
}
