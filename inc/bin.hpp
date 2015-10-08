#ifndef H_BIN
#define H_BIN

#include <ostream>
#include <string>
#include <set>
#include <tuple>

#include "systematic.hpp"
#include "cut.hpp"

class Bin{
  typedef std::set<Systematic> SystCollection;
public:
  Bin(const std::string &name, const class Cut &cut,
      const SystCollection &systematics = SystCollection());

  const std::string Name() const;
  Bin & Name(const std::string &name);

  const class Cut & Cut() const;
  class Cut & Cut();

  const SystCollection & Systematics() const;
  Bin & Systematics(const SystCollection &systematics);
  Bin & AddSystematic(const Systematic &systematic);
  Bin & AddSystematics(const SystCollection &systematic);
  bool HasSystematic(const Systematic &systematic) const;
  Bin & RemoveSystematic(const Systematic &systematic);
  Bin & RemoveSystematics();
  Bin & SetSystematicStrength(const std::string &name, double strength);

  bool operator<(const Bin &b) const;

private:
  class Cut cut_;
  std::string name_;
  SystCollection systematics_;

  using CompType = std::tuple<const class Cut&, const SystCollection&>;
  CompType ComparisonTuple() const{
    return CompType(cut_, systematics_);
  }
};

std::ostream & operator<<(std::ostream &stream, const Bin &bin);

#endif
