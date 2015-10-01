#ifndef H_DUMP_MARKOV
#define H_DUMP_MARKOV

#include <map>

#include "TH1D.h"

#include "RooArgSet.h"
#include "RooDataSet.h"

void SetupHistos(std::map<unsigned long,TH1D> &hmap,
		 const RooArgSet &args,
		 const RooDataSet &data);

#endif
