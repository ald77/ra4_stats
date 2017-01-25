#include <iostream>
#include <string>

#include <unistd.h>
#include <getopt.h>

#include "TFile.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooFitResult.h"

using namespace std;

namespace{
  string file_path = "wspace_nosyst_nokappa_nor4_T1tttt_mGluino-1700_mLSP-100_xsecNom.root";
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      file_path = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  TFile in_file(file_path.c_str(), "read");
  RooWorkspace *w = static_cast<RooWorkspace*>(in_file.Get("w"));
  RooDataSet *data_obs = static_cast<RooDataSet*>(w->data("data_obs"));
  RooProdPdf *model_b = static_cast<RooProdPdf*>(w->pdf("model_b"));
  
  RooNLLVar nll("nll", "nll", *model_b, *data_obs);
  nll.Print();
  RooMinuit minuit(nll);
 
  minuit.setPrintLevel(99999);
  minuit.setVerbose(true);
  minuit.setStrategy(2);
  minuit.optimizeConst(true);
  minuit.hesse();
  minuit.migrad();

  nll.Print();

  minuit.save()->Print("v");
}
