#ifndef H_FREE_SYSTEMATIC
#define H_FREE_SYSTEMATIC

#include <string>
#include <map>
#include <utility>
#include <ostream>

#include "bin.hpp"
#include "process.hpp"

class FreeSystematic{
public:
  FreeSystematic(const std::string &name);

  const std::string & Name() const;
  std::string & Name();

  bool HasEntry(const Bin &bin, const Process &process) const;

  double Strength(const Bin &bin, const Process &process) const;
  double & Strength(const Bin &bin, const Process &process);

  bool operator<(const FreeSystematic &syst) const;
  bool operator==(const FreeSystematic &syst) const;

private:
  std::string name_;
  std::map<std::pair<Bin, Process>, double> strengths_;
};

std::ostream & operator<<(std::ostream &stream, const FreeSystematic &syst);

#endif
