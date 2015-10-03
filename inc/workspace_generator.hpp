#ifndef H_WORKSPACE_GENERATOR
#define H_WORKSPACE_GENERATOR

#include <set>
#include <map>
#include <string>
#include <utility>

#include "RooWorkspace.h"

#include "cut.hpp"
#include "block.hpp"
#include "process.hpp"
#include "yield_key.hpp"
#include "gamma_params.hpp"
#include "block_yields.hpp"

class WorkspaceGenerator{
public:
  WorkspaceGenerator(const Cut &baseline,
		     const std::set<Block> &blocks,
		     const std::set<Process> &backgrounds,
		     const Process &signal,
		     const Process &data);

  void WriteToFile(const std::string &file_name);

private:
  Cut baseline_;
  std::set<Process> backgrounds_;
  Process signal_, data_;
  std::set<Block> blocks_;
  RooWorkspace w_;
  std::set<std::string>  poi_, observables_, nuisances_, systematics_;

  static std::map<YieldKey, GammaParams> yields_;

  void GetYields() const;
  void StoreYield(const Bin &bin, const Process &process) const;
  void AddPOI();
  void AddData(const Block &block);
  void AddBackgroundFractions(const Block &block);
  std::map<Process, double> GetBackgroundFractions(const Block &block) const;
  BlockYields AddABCDParameters(const Block &block);
  void AddBackgroundPredictions(const Block &block, const BlockYields &by);
};

#endif
