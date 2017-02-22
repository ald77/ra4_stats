#ifndef H_WORKSPACE_GENERATOR
#define H_WORKSPACE_GENERATOR

#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <random>

#include "RooWorkspace.h"

#include "cut.hpp"
#include "block.hpp"
#include "process.hpp"
#include "yield_key.hpp"
#include "gamma_params.hpp"
#include "block_yields.hpp"
#include "yield_manager.hpp"
#include "free_systematic.hpp"

class WorkspaceGenerator{
public:
  WorkspaceGenerator(const Cut &baseline,
                     const std::set<Block> &blocks,
                     const std::set<Process> &backgrounds,
                     const Process &signal,
                     const Process &data,
                     const std::string &systematics_file = "",
                     const bool use_r4 = true,
                     const double sig_strength = 0.,
                     const double sig_xsec_f = 1.);

  enum class PrintLevel{silent, important, normal, everything};

  void WriteToFile(const std::string &file_name);

  double GetLuminosity() const;
  WorkspaceGenerator & SetLuminosity(double luminosity);

  bool GetDoSystematics() const;
  WorkspaceGenerator & SetDoSystematics(bool do_systematics);

  bool GetDoDilepton() const;
  WorkspaceGenerator & SetDoDilepton(bool do_systematics);

  PrintLevel GetPrintLevel() const;
  WorkspaceGenerator & SetPrintLevel(PrintLevel print_level);

  bool GetKappaCorrected() const;
  WorkspaceGenerator & SetKappaCorrected(bool do_kappa_correction);

  double GetRMax() const;
  WorkspaceGenerator & SetRMax(double rmax);

  bool UseGausApprox() const;
  WorkspaceGenerator & UseGausApprox(bool use_gaus_approx);

  GammaParams GetYield(const YieldKey &key) const;
  GammaParams GetYield(const Bin &bin,
                       const Process &process,
                       const Cut &cut) const;
  GammaParams GetYield(const Bin &bin,
                       const Process &process) const;

  size_t AddToys(size_t num_toys = 0);

  const Process & GetInjectionModel() const;
  WorkspaceGenerator & SetInjectionModel(const Process &injection);
  bool GetDefaultInjectionModel() const;
  WorkspaceGenerator & SetDefaultInjectionModel();

  friend std::ostream & operator<<(std::ostream& stream, const WorkspaceGenerator &wg);

private:
  Cut baseline_;
  std::set<Process> backgrounds_;
  Process signal_, data_;
  Process injection_;
  bool inject_other_signal_;
  std::set<Block> blocks_;
  std::map<std::string, double> obs_vals_;
  std::map<std::string, std::poisson_distribution<> > obs_gens_;
  std::string systematics_file_;
  bool use_r4_;
  double sig_strength_, sig_xsec_f_;
  double rmax_;
  RooWorkspace w_;
  std::set<std::string>  poi_, observables_, nuisances_, systematics_;
  std::set<FreeSystematic> free_systematics_;
  double luminosity_;
  PrintLevel print_level_;
  bool do_systematics_;
  bool do_dilepton_;
  bool do_mc_kappa_correction_;
  size_t num_toys_;
  bool gaus_approx_;
  mutable bool w_is_valid_;

  static YieldManager yields_;
  static std::mt19937_64 prng_;
  static std::poisson_distribution<> dist_;

  static std::mt19937_64 InitializePRNG();
  static int GetPoisson(double rate);

  void SetupToys(const RooArgSet &obs);
  void GenerateToys(RooArgSet &obs);
  void ResetToys(RooArgSet &obs);
  void UpdateWorkspace();
  void AddPOI();
  void ReadSystematicsFile();
  static void CleanLine(std::string &line);
  void AddDileptonSystematic();
  bool NeedsDileptonBin(const Bin &bin) const;
  void MakeDileptonBin(const Bin &bin, Bin &dilep_bin, Cut &dilep_cut) const;
  void AddSystematicsGenerators();
  void AddSystematicGenerator(const std::string &name);
  void AddData(const Block &block);
  void AddBackgroundFractions(const Block &block);
  void AddABCDParameters(const Block &block);
  void AddRawBackgroundPredictions(const Block &block);
  void AddKappas(const Block &block);
  void AddMCYields(const Block &block);
  void AddMCPdfs(const Block &block);
  void AddMCProcessSums(const Block &block);
  void AddMCRowSums(const Block &block);
  void AddMCColSums(const Block &block);
  void AddMCTotal(const Block &block);
  void AddMCPrediction(const Block &block);
  void AddMCKappa(const Block &block);
  void AddFullBackgroundPredictions(const Block &block);
  void AddSignalPredictions(const Block &block);
  void AddPdfs(const Block &block);
  void AddDebug(const Block &block);
  void AddDummyNuisance();
  void AddFullPdf();
  void AddParameterSets();
  void DefineParameterSet(const std::string &cat_name,
                          const std::set<std::string> &var_names);
  void AddModels();
  void AddPoisson(const std::string &pdf_name,
		  const std::string &n_name,
		  const std::string &mu_name,
		  bool allow_approx);
  void PrintComparison(std::ostream &stream, const Bin &bin,
                       const Process &process, const Block &block) const;
};

#endif
