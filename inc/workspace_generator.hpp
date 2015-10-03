#ifndef H_WORKSPACE_GENERATOR
#define H_WORKSPACE_GENERATOR

#include <set>
#include <map>
#include <string>

#include "cut.hpp"
#include "block.hpp"
#include "process.hpp"
#include "yield_key.hpp"
#include "gamma_params.hpp"

class WorkspaceGenerator{
public:
  WorkspaceGenerator(const Cut &baseline,
		     const std::set<Block> &blocks,
		     const std::set<Process> &backgrounds,
		     const Process &signal,
		     const Process &data);

  void WriteToFile(const std::string &file_name);

private:
  Process signal_, data_;
  std::set<Block> blocks_;
  std::set<Process> backgrounds_;
  Cut baseline_;

  static std::map<YieldKey, GammaParams> yields_;

  void GetYields();
  void StoreYield(const Bin &bin, const Process &process);
};

#endif
