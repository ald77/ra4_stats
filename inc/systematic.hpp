#ifndef H_SYSTEMATIC
#define H_SYSTEMATIC

#include <string>

struct Systematic{
  Systematic(const std::string &base_name,
	     double multiplier);

  std::string base_name_;
  double multiplier_;

  bool operator<(const Systematic &s) const;
};

#endif
