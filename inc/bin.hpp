#ifndef H_BIN
#define H_BIN

#include <string>
#include <vector>

#include "systematic.hpp"

struct Bin{
  Bin(const std::string &name, const std::string &cut,
      const std::vector<Systematic> &systematics = std::vector<Systematic>());

  std::string name_, cut_;
  std::vector<Systematic> systematics_;

  bool operator<(const Bin &b) const;
};

#endif
