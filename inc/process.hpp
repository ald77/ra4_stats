#ifndef H_PROCESS
#define H_PROCESS

#include <string>
#include <vector>
#include <initializer_list>
#include <tuple>

#include "TChain.h"

#include "cut.hpp"

class Process{
public:
  Process(const std::string &name,
          const std::vector<std::string> &file_names,
          const Cut &cut = ::Cut(),
          bool count_zeros = true);
  Process(const std::string &name,
          std::initializer_list<std::string> file_names,
          const Cut &cut = ::Cut(),
          bool count_zeros = true);

  const std::string & Name() const;
  Process & Name(const std::string &name);

  const Cut & Cut() const;
  Process & Cut(const class Cut &cut);

  bool CountZeros() const;
  Process & CountZeros(bool count_zeros);

  long GetEntries() const;
  void GetCountAndUncertainty(double &count, double &uncertainty,
			      const class Cut &cut = ::Cut("1")) const;

  bool operator<(const Process &p) const;

private:
  mutable TChain chain_;
  class Cut cut_;
  std::string name_;
  bool count_zeros_;

  Process(Process&& p) = delete;
  Process(const Process &p) = delete;
  Process& operator=(const Process &p) = delete;

  void CleanName();

  auto ComparisonTuple() const{
    return std::make_tuple(cut_, chain_.Hash(), count_zeros_);
  }
};

#endif
