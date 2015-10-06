#ifndef H_WORKSPACE_GENERATOR
#define H_WORKSPACE_GENERATOR

#include <ostream>
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

  friend std::ostream & operator<<(std::ostream& stream, const WorkspaceGenerator &wg);
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
  void AddSystematicsGenerators();
  void AddSystematicGenerator(const std::string &name);
  void AddData(const Block &block);
  void AddBackgroundFractions(const Block &block);
  std::map<Process, double> GetBackgroundFractions(const Block &block) const;
  void AddABCDParameters(const Block &block);
  void AddRawBackgroundPredictions(const Block &block);
  void AddFullBackgroundPredictions(const Block &block);
  void AddSignalPredictions(const Block &block);
  void AddPdfs(const Block &block);
  void AddFullPdf();
  void AddParameterSets();
  void DefineParameterSet(const std::string &cat_name,
                          const std::set<std::string> &var_names);
  void AddModels();
};

#endif
