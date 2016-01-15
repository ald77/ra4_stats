#ifndef H_UTILITIES
#define H_UTILITIES

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>

#include "TTree.h"

#include "RooWorkspace.h"

#include "cut.hpp"

#define ERROR(x) do{throw std::runtime_error(string("Error in file ")+__FILE__+" at line "+to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

bool Contains(const std::string &str, const std::string &pat);
bool StartsWith(const std::string &str, const std::string &pat);

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

std::string execute(const std::string &cmd);

std::vector<std::string> Tokenize(const std::string& input,
                                  const std::string& tokens=" ");

std::string ChangeExtension(std::string path, const std::string &new_ext);

std::string MakeDir(std::string prefix);

void parseMasses(const std::string &str, int &mglu, int &mlsp);

template<typename T>
void Append(T &collection, const typename T::value_type &value){
  collection.insert(collection.end(), value);
}

#endif
