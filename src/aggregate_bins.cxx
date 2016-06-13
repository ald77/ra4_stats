#include "aggregate_bins.hpp"

#include <string>
#include <set>
#include <limits>
#include <sstream>
#include <iomanip>

#include <unistd.h>
#include <getopt.h>

#include "bin.hpp"
#include "process.hpp"
#include "utilities.hpp"
#include "systematic.hpp"
#include "cut.hpp"
#include "cross_sections.hpp"

#include "workspace_generator.hpp"

using namespace std;

namespace{
  string out_dir = "/net/cms2/cms2r0/babymaker/wspaces/agg_bins";
  double lumi = 20.;
  bool do_track_veto = true;
  double met_low = 200.;
  double met_high = -1.;
  double njets_low = 5.5;
  double njets_high = -1.;
  double nbm_low = 0.5;
  double nbm_high = -1.;
  int mglu = 1500.;
  int mlsp = 100.;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  string mc = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_baseline/";
  Process ttbar{"ttbar", {
      {mc+"*TTJets*Lept*.root/tree"},
        {mc+"*TTJets*HT*.root/tree"}
    }, "stitch&&pass"};
  Process other{"other", {
      {mc+"*_WJetsToLNu*.root/tree"},
        {mc+"*_ST_*.root/tree"},
          {mc+"*_TTWJets*.root/tree"},
            {mc+"*_TTZTo*.root/tree"},
              {mc+"*DYJetsToLL*.root/tree"},
                {mc+"*_QCD_HT*.root/tree"},
                  {mc+"*_ZJET*.root/tree"},
                    {mc+"*_WWTo*.root/tree"},
                      {mc+"*ggH_HToBB*.root/tree"},
                        {mc+"*ttHJetTobb*.root/tree"},
                          {mc+"*_TTGJets*.root/tree"},
                            {mc+"*_TTTT_*.root/tree"},
                              {mc+"*_WH_HToBB*.root/tree"},
                                {mc+"*_WZTo*.root/tree"},
                                  {mc+"*_ZH_HToBB*.root/tree"},
                                    {mc+"*_ZZ_*.root/tree"},
    }, "stitch&&pass"};
  Process data{"data", {{"/net/cms27/cms27r0/babymaker/2016_04_29/data/merged_baseline/*.root/tree"}}, "(trig[4]||trig[8])&&pass", true};
  
  Process signal{"signal", {
      {"/net/cms27/cms27r0/babymaker/2016_04_29/mc/T1tttt/skim_baseline/*SMS-T1tttt_mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp)+"*.root/tree"}
    }, "stitch&&pass", false, true};

  Cut baseline("mj14>250.");

  string met = "&&met>"+to_string(met_low);
  if(met_high > met_low) met += "&&met<=" + to_string(met_high);
  string njets_nbm = "&&njets>"+to_string(njets_low)+"&&nbm>"+to_string(nbm_low);
  if(njets_high > njets_low) njets_nbm += "&&njets<=" + to_string(njets_high);
  if(nbm_high > nbm_low) njets_nbm += "&&nbm<=" + to_string(nbm_high);
  string tkveto = "&&nveto==0";

  Bin r1("r1", "mt<=140.&&mj14<=400."+met, true);
  Bin r2("r2", "mt<=140.&&mj14>400."+met+njets_nbm, true);
  Bin r3("r3", "mt>140.&&mj14<=400."+met+tkveto, true);
  Bin r4("r4", "mt>140.&&mj14>400."+met+njets_nbm+tkveto, true);

  set<Block> blocks = {{"all", {{r1, r2}, {r3, r4}}}};

  WorkspaceGenerator wg(baseline, blocks, {ttbar, other}, signal, data,
                        "txt/systematics/simple.txt", true, 0., 1.);
  wg.SetRMax(10.);
  wg.SetKappaCorrected(true);
  wg.SetDoSystematics(true);
  wg.SetLuminosity(lumi);

  ostringstream oss;
  oss << setprecision(numeric_limits<double>::digits10) << out_dir << "/wspace_aggbin_"
      << "t1tttt_" << mglu << '_' << mlsp << "_lumi_" << lumi
      << "_met_" << met_low << '_';
  if(met_high > met_low) oss << met_high;
  else oss << "inf";
  oss << "_njets_" << ceil(njets_low) << '_';
  if(njets_high > njets_low) oss << floor(njets_high);
  else oss << "inf";
  oss << "_nbm_" << ceil(nbm_low) << '_';
  if(nbm_high > nbm_low) oss << floor(nbm_high);
  else oss << "inf";
  oss << "_tkveto_" << (do_track_veto ? "true" : "false")
      << ".root" << flush;
  wg.WriteToFile(oss.str());
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"out_dir", required_argument, 0, 0},
      {"lumi", required_argument, 0, 0},
      {"mglu", required_argument, 0, 0},
      {"mlsp", required_argument, 0, 0},
      {"do_track_veto", required_argument, 0, 0},
      {"met_low", required_argument, 0, 0},
      {"met_high", required_argument, 0, 0},
      {"njets_low", required_argument, 0, 0},
      {"njets_high", required_argument, 0, 0},
      {"nbm_low", required_argument, 0, 0},
      {"nbm_high", required_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 0:
      optname = long_options[option_index].name;
      if(optname == "out_dir"){
        out_dir = optarg;
      }else if(optname == "lumi"){
        lumi = atof(optarg);
      }else if(optname == "do_track_veto"){
        do_track_veto = atoi(optarg);
      }else if(optname == "met_low"){
        met_low = atof(optarg);
      }else if(optname == "met_high"){
        met_high = atof(optarg);
      }else if(optname == "njets_low"){
        njets_low = atof(optarg);
      }else if(optname == "njets_high"){
        njets_high = atof(optarg);
      }else if(optname == "nbm_low"){
        nbm_low = atof(optarg);
      }else if(optname == "nbm_high"){
        nbm_high = atof(optarg);
      }else if(optname == "mglu"){
        mglu = atoi(optarg);
      }else if(optname == "mlsp"){
        mlsp = atoi(optarg);
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
