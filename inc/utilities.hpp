#ifndef H_UTILITIES
#define H_UTILITIES

#include <string>
#include <vector>

#include "TTree.h"

#include "RooWorkspace.h"

#include "cut.hpp"

bool Contains(const std::string &str, const std::string &pat);

void ReplaceAll(std::string &str, const std::string &orig, const std::string &rep);

void RmCutOn(std::string &cut, const std::string &orig, const std::string &rep="1");

size_t MaxIndex(const std::vector<double> &v);

void DefineSet(RooWorkspace &w,
               const std::string &set_name,
               const std::vector<std::string> &var_names);

void GetCountAndUncertainty(TTree &tree,
                            const Cut &cut,
                            double &count,
                            double &uncertainty);

#endif
