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
  string out_dir = "/net/top/homes/ald77/ra4_stats/wspaces";
  string syst_file = "txt/systematics/simple.txt";
  double lumi = 36.2;
  bool do_track_veto = true;
  double met_low = 200.;
  double met_high = -1.;
  double njets_low = 5.5;
  double njets_high = -1.;
  double nbm_low = 0.5;
  double nbm_high = -1.;
  int mglu = 1700.;
  int mlsp = 100.;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  string hostname = execute("echo $HOSTNAME");
  string basefolder("/net/cms2/cms2r0/babymaker/");
  if(Contains(hostname, "lxplus")) basefolder = "/afs/cern.ch/user/m/manuelf/work/";
  string foldermc(basefolder+"babies/2016_08_10/mc/merged_mcbase_abcd/"); 
  string foldersig(basefolder+"babies/2016_08_10/T1tttt/skim_abcd/"); 
  string folderdata(basefolder+"babies/2016_11_08/data/merged_database_abcd/");

  //Define processes. Try to minimize splitting
  string stitch_cuts("stitch&&pass");

  Process ttbar{"ttbar", {
      {foldermc+"/*_TTJets*SingleLept*.root/tree",
	  foldermc+"/*_TTJets*DiLept*.root/tree",
	  foldermc+"/*_TTJets_HT*.root/tree"}
    },stitch_cuts};
  Process other{"other", {
      {foldermc+"/*_WJetsToLNu*.root/tree",
	  foldermc+"/*_ST_*.root/tree",
	  foldermc+"/*_TTW*.root/tree",
	  foldermc+"/*_TTZ*.root/tree",
	  foldermc+"/*DYJetsToLL*.root/tree",
	  foldermc+"/*_ZJet*.root/tree",
	  foldermc+"/*_ttHJetTobb*.root/tree",
	  foldermc+"/*_TTGJets*.root/tree",
	  foldermc+"/*_TTTT*.root/tree",
	  foldermc+"/*_WH_HToBB*.root/tree",
	  foldermc+"/*_ZH_HToBB*.root/tree",
	  foldermc+"/*_WWTo*.root/tree",
	  foldermc+"/*_WZ*.root/tree",
	  foldermc+"/*_ZZ_*.root/tree",
	  foldermc+"/*QCD_HT*0_Tune*.root/tree",
	  foldermc+"/*QCD_HT*Inf_Tune*.root/tree"}
    },stitch_cuts};

  string data_cuts("trig_ra4&&pass");

  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};
  
  Process signal{"signal", {
      {foldersig+"/*SMS-T1tttt_mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp)+"_*.root/tree"}
    }, "stitch", false, true};

  Cut baseline("met/met_calo<5.&&pass_ra2_badmu&&st>500&&met>200&&nleps==1&&nbm>=1&&njets>=6&&mj14>250.");

  string met = "&&met>"+to_string(met_low);
  if(met_high > met_low) met += "&&met<=" + to_string(met_high);
  string njets_nbm = "&&njets>"+to_string(njets_low)+"&&nbm>"+to_string(nbm_low);
  if(njets_high > njets_low) njets_nbm += "&&njets<=" + to_string(njets_high);
  if(nbm_high > nbm_low) njets_nbm += "&&nbm<=" + to_string(nbm_high);
  string tkveto = "&&nveto==0";

  Bin r1("r1", "mt<=140.&&mj14<=400."+met+tkveto, false);
  Bin r2("r2", "mt<=140.&&mj14>400."+met+njets_nbm+tkveto, false);
  Bin r3("r3", "mt>140.&&mj14<=400."+met+tkveto, false);
  Bin r4("r4", "mt>140.&&mj14>400."+met+njets_nbm+tkveto, false);

  set<Block> blocks = {{"all", {{r1, r2}, {r3, r4}}}};

  WorkspaceGenerator wg(baseline, blocks, {ttbar, other}, signal, data,
                        syst_file, true, 0., 1.);
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
      {"syst_file", required_argument, 0, 0},
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
      }else if(optname == "syst_file"){
	syst_file = optarg;
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
