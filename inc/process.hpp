#ifndef H_PROCESS
#define H_PROCESS

#include <string>
#include <set>
#include <initializer_list>
#include <tuple>
#include <memory>

#include "TChain.h"

#include "cut.hpp"
#include "gamma_params.hpp"

class Process{
public:
  Process(const std::string &name,
          const std::set<std::string> &file_names,
          const Cut &cut = ::Cut(),
          bool count_zeros = true);
  Process(const std::string &name,
          std::initializer_list<std::string> file_names,
          const Cut &cut = ::Cut(),
          bool count_zeros = true);

  const std::string & Name() const;
  Process & Name(const std::string &name);

  const class Cut & Cut() const;
  class Cut & Cut();
  
  const bool & CountZeros() const;
  bool & CountZeros();

  const std::set<std::string> & FileNames() const;

  long GetEntries() const;
  GammaParams GetYield(const class Cut &cut = ::Cut("1")) const;

  bool operator<(const Process &p) const;

private:
  std::set<std::string> file_names_;
  mutable std::shared_ptr<TChain> chain_;
  class Cut cut_;
  std::string name_;
  bool count_zeros_;

  void CleanName();
  void AddFiles();

  using CompType = std::tuple<const class Cut&, const std::set<std::string> &, const bool &>;
  CompType ComparisonTuple() const{
    return CompType(cut_, file_names_, count_zeros_);
  }
};

#endif
