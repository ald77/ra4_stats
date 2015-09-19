#ifndef H_MAKE_WORKSPACE
#define H_MAKE_WORKSPACE

#include <string>
#include <vector>
#include <functional>
#include <initializer_list>
#include <map>

#include "TChain.h"

#include "RooWorkspace.h"

#include "gamma_params.hpp"

struct Bin{
  Bin(const std::string &name, const std::string &cut);

  std::string name_, cut_;

  bool operator<(const Bin &b) const;
};

struct Process{
  Process(const std::string &name,
	  const std::vector<std::string> &file_names,
	  const std::string &cut = "1",
	  bool count_zeros = true);
  Process(const std::string &name,
	  std::initializer_list<std::string> file_names,
	  const std::string &cut = "1",
	  bool count_zeros = true);

  Process(Process&& p) = delete;
  Process(const Process &p) = delete;
  Process& operator=(const Process &p) = delete;

  bool operator<(const Process &p) const;

  TChain chain_;
  std::string name_, cut_;
  bool count_zeros_;
};

struct Block{
  Block(const std::string &name, const std::vector<std::vector<Bin> > &bins);
  Block(const std::string &name, std::initializer_list<std::vector<Bin> > bins);

  std::vector<std::vector<Bin> > bins_;
  std::string name_;
  };

struct BinProc{
  BinProc(const Bin &bin, Process &process);

  bool operator<(const BinProc &bp) const;

  Process &process_;
  Bin bin_;
};

std::map<BinProc, GammaParams> GetYields(const std::vector<Block> &blocks,
					 const std::string &baseline,
					 Process &data,
					 Process &signal,
					 std::vector<std::reference_wrapper<Process> > &backgrounds);

GammaParams GetYield(const BinProc &bp,
		     const std::string &baseline);

void GetCountAndUncertainty(TTree &tree,
			       const std::string &cut,
			       double &count,
			       double &uncertainty);

void MakeWorkspace(const std::string &file_name,
		   const std::string &baseline,
		   const std::vector<Block> &blocks,
		   Process &data,
		   Process &signal,
		   std::vector<std::reference_wrapper<Process> > &backgrounds);

std::vector<double> GetBackgroundFractions(const Block &block,
					   std::vector<std::reference_wrapper<Process> > &backgrounds,
					   const std::map<BinProc, GammaParams> &yields);

void AddBackgroundFractions(RooWorkspace &w,
			    const Block &block,
			    std::vector<std::reference_wrapper<Process> > &backgrounds,
			    const std::map<BinProc, GammaParams> &yields,
			    std::vector<std::string> &nuis_names);

void AddABCDParams(RooWorkspace &w,
		   const Block &block,
		   std::vector<std::reference_wrapper<Process> > &backgrounds,
		   const std::map<BinProc, GammaParams> &yields,
		   std::vector<std::string> &nuis_names);

void AddBackgroundPreds(RooWorkspace &w,
			const Block &block,
			const std::vector<std::reference_wrapper<Process> > &backgrounds);

void AddSignalPreds(RooWorkspace &w,
		    const Block &block,
		    Process &signal,
		    const std::map<BinProc, GammaParams> &yields);

void AddBinPdfs(RooWorkspace &w,
		const Block &block);

void AddData(RooWorkspace &w,
	     const Block &block,
	     Process &data,
	     const std::map<BinProc, GammaParams> &yields,
	     std::vector<std::string> &obs_names);

void DefineSet(RooWorkspace &w,
	       const std::string &set_name,
	       const std::vector<std::string> &var_names);

void AddModels(RooWorkspace &w,
	       const std::vector<Block> &blocks);

#endif
