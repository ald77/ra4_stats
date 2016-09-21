#ifndef H_PROCESS
#define H_PROCESS

#include <string>
#include <set>
#include <initializer_list>
#include <tuple>
#include <memory>
#include <ostream>

#include "TChain.h"

#include "cut.hpp"
#include "gamma_params.hpp"
#include "systematic.hpp"

class Process{
  typedef std::set<Systematic> SystCollection;
public:
  Process() = default;
  Process(const std::string &name,
          const std::set<std::string> &file_names,
          const Cut &cut = ::Cut(),
          bool is_data = false,
          bool is_signal = false,
          bool count_zeros = true,
          const SystCollection &systematics = SystCollection());
  Process(const std::string &name,
          std::initializer_list<std::string> file_names,
          const Cut &cut = ::Cut(),
          bool is_data = false,
          bool is_signal = false,
          bool count_zeros = true,
          const SystCollection &systematics = SystCollection());

  const std::string & Name() const;
  Process & Name(const std::string &name);

  const class Cut & Cut() const;
  class Cut & Cut();

  const bool & IsData() const;
  bool & IsData();
  
  const bool & IsSignal() const;
  bool & IsSignal();
  
  const bool & CountZeros() const;
  bool & CountZeros();

  const std::set<std::string> & FileNames() const;

  long GetEntries() const;
  GammaParams GetYield(const class Cut &cut = ::Cut("1")) const;

  const SystCollection & Systematics() const;
  Process & Systematics(const SystCollection &systematics);
  Process & AddSystematic(const Systematic &systematic);
  Process & AddSystematics(const SystCollection &systematic);
  bool HasSystematic(const Systematic &systematic) const;
  Process & RemoveSystematic(const Systematic &systematic);
  Process & RemoveSystematics();
  Process & SetSystematicStrength(const std::string &name, double strength);

  bool operator<(const Process &p) const;
  bool operator==(const Process &p) const;

private:
  std::set<std::string> file_names_;
  mutable std::shared_ptr<TChain> chain_;
  class Cut cut_;
  std::string name_;
  bool is_data_;
  bool is_signal_;
  bool count_zeros_;
  SystCollection systematics_;

  void CleanName();
  void AddFiles();
};

std::ostream & operator<<(std::ostream &stream, const Process &proc);

#endif
