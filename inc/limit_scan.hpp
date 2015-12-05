#ifndef H_LIMIT_SCAN
#define H_LIMIT_SCAN

#include <string>

#include "TGraph2D.h"

class TLegend;

void DrawContours(TGraph2D &g2, int color, int style,
                  TLegend *l = nullptr, const std::string &name = "");
void GetOptions(int argc, char *argv[]);
void SetupColors();

#endif