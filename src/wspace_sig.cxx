#include "make_workspace.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <initializer_list>
#include <vector>
#include <string>
#include <stdlib.h>
#include <ctime>
#include <sys/stat.h>

#include <unistd.h>
#include <getopt.h>

#include "TString.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "bin.hpp"
#include "process.hpp"
#include "utilities.hpp"
#include "systematic.hpp"
#include "cut.hpp"
#include "cross_sections.hpp"

#include "workspace_generator.hpp"

using namespace std;

namespace{
  double lumi = 2.138;
  double sig_strength = 0.;
  BlindLevel blind_level = BlindLevel::unblinded;
  bool no_kappa = false;
  bool do_syst = true;
  bool use_r4 = true;
  string minjets("6");
  string hijets("9");
  string himet("400");
  string mjthresh("400");
  unsigned n_toys = 0;
  string sigfile = "";
  string outfolder = "out/";
}

int main(int argc, char *argv[]){
  time_t begtime, endtime;
  time(&begtime);
  cout << fixed << setprecision(2);
  GetOptions(argc, argv);
  if(sigfile==""){
    cout<<endl<<"You need to specify the input file with -f. Exiting"<<endl<<endl;
    return 1;
  }
  string midjets(""); midjets += to_string(atoi(hijets.c_str())-1);
  string minjets2l(""); minjets2l += to_string(atoi(minjets.c_str())-1);
  string midjets2l(""); midjets2l += to_string(atoi(midjets.c_str())-1);

  string hostname = execute("echo $HOSTNAME");
  string basefolder("/net/cms2/cms2r0/babymaker/");
  if(Contains(hostname, "lxplus")) basefolder = "/afs/cern.ch/user/m/manuelf/work/";
  string skim("skim_abcd/");
  string foldermc(basefolder+"babies/2015_10_19/mc/merged_abcd/");
  string folderdata(basefolder+"babies/2015_11_20/data/singlelep/combined/"+skim);

  //Define processes. Try to minimize splitting
  Process ttbar{"ttbar", {
      {foldermc+"/*TTJets*.root/tree"}
    }};
  Process other{"other", {
      {foldermc+"/*other*.root/tree"}
    }};
  Process signal{"signal", {
      {sigfile+"/tree"}
    }, Cut(), false, true};

  string data_cuts("(trig[4]||trig[8])&&pass");

  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline1b{"1"};

  //Declare bins 
  Bin r1_lowmet_allnb{"r1_lowmet_allnb", "mt<=140&&mj<="+mjthresh+"&&met<="+himet,
      blind_level>=BlindLevel::blinded};
  Bin r1_highmet_allnb{"r1_highmet_allnb", "mt<=140&&mj<="+mjthresh+"&&met>"+himet,
      blind_level>=BlindLevel::blinded};

  Bin r2_lowmet_lownj_1b{"r2_lowmet_lownj_1b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==1",
      blind_level>=BlindLevel::blinded};
  Bin r2_lowmet_highnj_1b{"r2_lowmet_highnj_1b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==1",
      blind_level>=BlindLevel::blinded};
  Bin r2_lowmet_lownj_2b{"r2_lowmet_lownj_2b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==2",
      blind_level>=BlindLevel::blinded};
  Bin r2_lowmet_lownj_3b{"r2_lowmet_lownj_3b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm>2",
      blind_level>=BlindLevel::blinded};
  Bin r2_lowmet_highnj_2b{"r2_lowmet_highnj_2b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==2",
      blind_level>=BlindLevel::blinded};
  Bin r2_lowmet_highnj_3b{"r2_lowmet_highnj_3b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm>2",
      blind_level>=BlindLevel::blinded};

  Bin r2_highmet_lownj_1b{"r2_highmet_lownj_1b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm==1",
      blind_level>=BlindLevel::blinded};
  Bin r2_highmet_highnj_1b{"r2_highmet_highnj_1b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm==1",
      blind_level>=BlindLevel::blinded};
  Bin r2_highmet_lownj_2b{"r2_highmet_lownj_2b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm>=2",
      blind_level>=BlindLevel::blinded};
  Bin r2_highmet_highnj_2b{"r2_highmet_highnj_2b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm>=2",
      blind_level>=BlindLevel::blinded};

  Bin r3_lowmet_allnb{"r3_lowmet_allnb", "mt>140&&mj<="+mjthresh+"&&met<="+himet,
      blind_level>=BlindLevel::blinded};
  Bin r3_highmet_allnb{"r3_highmet_allnb", "mt>140&&mj<="+mjthresh+"&&met>"+himet,
      blind_level>=BlindLevel::blinded};

  Bin r4_lowmet_lownj_1b{"r4_lowmet_lownj_1b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==1",
      blind_level>BlindLevel::unblind_1b};
  Bin r4_lowmet_highnj_1b{"r4_lowmet_highnj_1b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==1",
      blind_level>BlindLevel::unblind_1b};
  Bin r4_lowmet_lownj_2b{"r4_lowmet_lownj_2b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==2",
      blind_level>BlindLevel::unblinded};
  Bin r4_lowmet_lownj_3b{"r4_lowmet_lownj_3b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm>2",
      blind_level>BlindLevel::unblinded};
  Bin r4_lowmet_highnj_2b{"r4_lowmet_highnj_2b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==2",
      blind_level>BlindLevel::unblinded};
  Bin r4_lowmet_highnj_3b{"r4_lowmet_highnj_3b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm>2",
      blind_level>BlindLevel::unblinded};

  Bin r4_highmet_lownj_1b{"r4_highmet_lownj_1b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm==1",
      blind_level>BlindLevel::unblind_1b};
  Bin r4_highmet_highnj_1b{"r4_highmet_highnj_1b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm==1",
      blind_level>BlindLevel::unblind_1b};
  Bin r4_highmet_lownj_2b{"r4_highmet_lownj_2b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm>=2",
      blind_level>BlindLevel::unblinded};
  Bin r4_highmet_highnj_2b{"r4_highmet_highnj_2b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm>=2",
      blind_level>BlindLevel::unblinded};

  //// METHOD 1BK: Adding 1b, fat R1/R3 integrated over njets, nb, but not MET
  set<Block> blocks_1bk{
    {"lowmet", {{r1_lowmet_allnb, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, r2_lowmet_lownj_2b, r2_lowmet_highnj_2b,
            r2_lowmet_lownj_3b, r2_lowmet_highnj_3b},
          {r3_lowmet_allnb, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, r4_lowmet_lownj_2b, r4_lowmet_highnj_2b,
              r4_lowmet_lownj_3b, r4_lowmet_highnj_3b}}},
      {"highmet", {{r1_highmet_allnb, r2_highmet_lownj_1b, r2_highmet_highnj_1b, r2_highmet_lownj_2b, r2_highmet_highnj_2b},
            {r3_highmet_allnb, r4_highmet_lownj_1b, r4_highmet_highnj_1b, r4_highmet_lownj_2b, r4_highmet_highnj_2b}}}
  };

  //// Parsing the gluino and LSP masses
  int mglu, mlsp;
  parseMasses(sigfile, mglu, mlsp);
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));
  double xsec, xsec_unc;
  xsec::signalCrossSection(mglu, xsec, xsec_unc);
  double rmax = 20.;
  if(mglu <= 1100. && mlsp <= 450.){
    rmax = 5.;
    if(mglu <= 625 && mlsp <= 375.){
      rmax = 1.5;
    }
  }

  //// Creating workspaces for the Nominal, uncert Up, and uncert Down signal cross sections
  Cut *pbaseline(&baseline1b);
  set<Block> *pblocks(&blocks_1bk);
  string sysfolder("/net/cms2/cms2r0/babymaker/sys/2015_11_28/scan/");
  if(Contains(hostname, "lxplus")) sysfolder = "txt/systematics/";
  string sysfile(sysfolder+"sys_SMS-T1tttt_"+glu_lsp+".txt");
  // If systematic file does not exist, use m1bk_nc for tests
  struct stat buffer;   
  if(stat (sysfile.c_str(), &buffer) != 0) {
    cout<<endl<<"WARNING: "<<sysfile<<" does not exist. Using ";
    sysfile = "txt/systematics/m1bk_nc.txt";
    cout<<sysfile<<" instead"<<endl<<endl;
  }

  gSystem->mkdir(outfolder.c_str(), kTRUE);
  string outname(outfolder+"/wspace_"+glu_lsp+"_xsecNom.root");

  WorkspaceGenerator wgNom(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1.);
  wgNom.SetRMax(rmax);
  wgNom.SetKappaCorrected(!no_kappa);
  wgNom.SetDoSystematics(do_syst);
  wgNom.SetLuminosity(lumi);
  wgNom.SetDoSystematics(do_syst);
  wgNom.AddToys(n_toys);
  wgNom.WriteToFile(outname);

  ReplaceAll(outname, "Nom", "Up");
  WorkspaceGenerator wgUp(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1+xsec_unc);
  wgUp.SetRMax(rmax);
  wgUp.SetKappaCorrected(!no_kappa);
  wgUp.SetDoSystematics(do_syst);
  wgUp.SetLuminosity(lumi);
  wgUp.SetDoSystematics(do_syst);
  wgUp.AddToys(n_toys);
  wgUp.WriteToFile(outname);

  ReplaceAll(outname, "Up", "Down");
  WorkspaceGenerator wgDown(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1-xsec_unc);
  wgDown.SetRMax(rmax);
  wgDown.SetKappaCorrected(!no_kappa);
  wgDown.SetDoSystematics(do_syst);
  wgDown.SetLuminosity(lumi);
  wgDown.SetDoSystematics(do_syst);
  wgDown.AddToys(n_toys);
  wgDown.WriteToFile(outname);

  time(&endtime); 
  cout<<"Finding workspaces took "<<fixed<<setprecision(0)<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  
}



void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sigfile", required_argument, 0, 'f'},
      {"lumi", required_argument, 0, 'l'},
      {"unblind", required_argument, 0, 'u'},
      {"no_syst", no_argument, 0, 0},
      {"lowj", required_argument, 0, 'j'},
      {"hij", required_argument, 0, 'h'},
      {"himet", required_argument, 0, 'm'},
      {"mj", required_argument, 0, 's'},
      {"nokappa", no_argument, 0, 'k'},
      {"use_r4", no_argument, 0, '4'},
      {"toys", required_argument, 0, 0},
      {"sig_strength", required_argument, 0, 'g'},
      {"outfolder", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:u:j:h:m:s:k4g:f:o:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 'g':
      sig_strength = atof(optarg);
      break;
    case 'u':
      if(string(optarg)=="all"){
        blind_level = BlindLevel::unblinded;
      }else if(string(optarg)=="sideband"){
        blind_level = BlindLevel::r4_blinded;
      }else if(string(optarg)=="1b"){
        blind_level = BlindLevel::unblind_1b;
      }else{
        blind_level = BlindLevel::blinded;
      }
      break;
    case 'o':
      outfolder = optarg;
      break;
    case 'f':
      sigfile = optarg;
      break;
    case 'j':
      minjets = optarg;
      break;
    case 'h':
      hijets = optarg;
      break;
    case 'm':
      himet = optarg;
      break;
    case 'k':
      no_kappa = true;
      break;
    case '4':
      use_r4 = true;
      break;
    case 's':
      mjthresh = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "no_syst"){
        do_syst = false;
      }else if(optname == "toys"){
        n_toys = atoi(optarg);
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
