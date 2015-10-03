#ifndef H_MAKE_WORKSPACE
#define H_MAKE_WORKSPACE

#include <string>
#include <vector>
#include <functional>
#include <initializer_list>
#include <map>
#include <set>

#include "RooWorkspace.h"

#include "gamma_params.hpp"
#include "block.hpp"
#include "bin.hpp"
#include "process.hpp"
#include "bin_proc.hpp"

void GetYields(const std::vector<Block> &blocks,
               const std::string &baseline,
               const Process &data,
               const Process &signal,
               const std::vector<std::reference_wrapper<Process> > &backgrounds,
               std::map<BinProc, GammaParams> &yields);

bool NeedsDileptonBin(const Bin &bin, const std::string &baseline);

void MakeDileptonBin(const Bin &bin, const std::string &baseline,
		     Bin &dilep_bin, std::string &dilep_baseline);

void StoreYield(const BinProc &bp,
                const std::string &baseline,
                std::map<BinProc, GammaParams> &yields);

void MakeWorkspace(const std::string &file_name,
                   const std::string &baseline,
                   std::vector<Block> blocks,
                   const Process &data,
                   const Process &signal,
                   const std::vector<std::reference_wrapper<Process> > &backgrounds,
                   std::map<BinProc, GammaParams> &yields);

std::vector<double> GetBackgroundFractions(const Block &block,
                                           const std::vector<std::reference_wrapper<Process> > &backgrounds,
                                           const std::map<BinProc, GammaParams> &yields);

void AddBackgroundFractions(RooWorkspace &w,
                            const Block &block,
                            const std::vector<std::reference_wrapper<Process> > &backgrounds,
                            const std::map<BinProc, GammaParams> &yields,
                            std::vector<std::string> &nuis_names);

void AddDileptonSystematics(Block &block,
			    const std::string &baseline,
			    const std::vector<std::reference_wrapper<Process> > &backgrounds,
			    const std::map<BinProc, GammaParams> &yields);

void AddABCDParams(RooWorkspace &w,
                   const Block &block,
                   const std::vector<std::reference_wrapper<Process> > &backgrounds,
                   const std::map<BinProc, GammaParams> &yields,
                   std::vector<std::string> &nuis_names,
                   size_t &max_col, size_t &max_row);

void AddBackgroundPreds(RooWorkspace &w,
                        const Block &block,
                        const std::vector<std::reference_wrapper<Process> > &backgrounds,
                        size_t max_col, size_t max_row,
			std::set<std::string> &syst_generators,
			std::vector<std::string> &nuis_names);

void AddSignalPreds(RooWorkspace &w,
                    const Block &block,
                    const Process &signal,
                    const std::map<BinProc, GammaParams> &yields);

void AddBinPdfs(RooWorkspace &w,
                const Block &block);

void AddMockData(RooWorkspace &w,
		 const Block &block,
		 const std::vector<std::reference_wrapper<Process> > &backgrounds,
		 const std::map<BinProc, GammaParams> &yields,
		 std::vector<std::string> &obs_names);

void AddData(RooWorkspace &w,
             const Block &block,
             const Process &data,
             const std::map<BinProc, GammaParams> &yields,
             std::vector<std::string> &obs_names);

void AddModels(RooWorkspace &w,
               const std::vector<Block> &blocks,
	       const std::set<std::string> &syst_generators);

void PrintDiagnostics(const RooWorkspace &w,
                      const std::vector<Block> &blocks,
                      const Process &data,
                      const Process &signal,
                      const std::vector<std::reference_wrapper<Process> > &backgrounds,
                      const std::map<BinProc, GammaParams> &yields);

void PrintComparison(const RooWorkspace &w,
                     const Block &block,
                     const BinProc &bp,
                     const std::map<BinProc, GammaParams> &yields,
                     bool is_data = false);

void GetOptions(int argc, char *argv[]);

#endif
