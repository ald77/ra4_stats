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

  enum class PrintLevel{silent, important, normal, everything};
  enum class BlindLevel{unblinded, r4_blinded, blinded};

  void WriteToFile(const std::string &file_name);

  double GetLuminosity() const;
  WorkspaceGenerator & SetLuminosity(double luminosity);

  BlindLevel GetBlindLevel() const;
  WorkspaceGenerator & SetBlindLevel(BlindLevel blind_level);

  bool GetDoSystematics() const;
  WorkspaceGenerator & SetDoSystematics(bool do_systematics);

  PrintLevel GetPrintLevel() const;
  WorkspaceGenerator & SetPrintLevel(PrintLevel print_level);

  static bool HaveYield(const YieldKey &key);
  GammaParams GetYield(const YieldKey &key) const;

  friend std::ostream & operator<<(std::ostream& stream, const WorkspaceGenerator &wg);

private:
  Cut baseline_;
  std::set<Process> backgrounds_;
  Process signal_, data_;
  std::set<Block> blocks_;
  RooWorkspace w_;
  std::set<std::string>  poi_, observables_, nuisances_, systematics_;
  double luminosity_;
  PrintLevel print_level_;
  BlindLevel blind_level_;
  bool do_systematics_;
  mutable bool w_is_valid_;

  static std::map<YieldKey, GammaParams> yields_;
  static const double yield_lumi_;

  void UpdateWorkspace();
  void ComputeYield(const Bin &bin, const Process &process) const;
  void ComputeYield(const YieldKey &key) const;
  void AddPOI();
  void AddDileptonSystematic();
  bool NeedsDileptonBin(const Bin &bin) const;
  void MakeDileptonBin(const Bin &bin, Bin &dilep_bin, Cut &dilep_cut) const;
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
  void PrintComparison(std::ostream &stream, const YieldKey &yk,
		       const Block &block, bool is_data) const;
};

#endif
