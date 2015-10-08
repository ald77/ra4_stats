#ifndef H_SENSITIVITY
#define H_SENSITIVITY

#include <string>

#include "RooWorkspace.h"

double GetSignificance(const std::string &file, double lumi);
double GetLimit(const std::string &file_name, double lumi);
void ModifyLumi(const std::string &file_name, double lumi);
double ExtractNumber(const std::string &results, const std::string &key);
void FixSignificance(double &x);
void FixLimit(double &x);

#endif
