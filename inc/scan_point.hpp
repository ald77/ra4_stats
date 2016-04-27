#ifndef H_SCAN_POINT
#define H_SCAN_POINT

#include <string>

double GetSignif(const std::string &filename);
std::string GetBaseName(const std::string &path);
double ExtractNumber(const std::string &results, const std::string &key);
void GetOptions(int argc, char *argv[]);

#endif
