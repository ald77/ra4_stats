#include "make_workspace.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <initializer_list>
#include <vector>
#include <string>
#include <stdlib.h>

#include <unistd.h>
#include <getopt.h>

#include "TString.h"

#include "bin.hpp"
#include "process.hpp"
#include "utilities.hpp"
#include "systematic.hpp"
#include "cut.hpp"

#include "workspace_generator.hpp"

using namespace std;

namespace{
  double lumi = 2.5;
  bool blinded = true;
  bool do_syst = true;
  string minjets("7");
  string hijets("9");
  string himet("400"); 
  bool do_compressed = false;
}

int main(int argc, char *argv[]){
  cout << fixed << setprecision(2);
  GetOptions(argc, argv);
  string midjets(""); midjets += to_string(atoi(hijets.c_str())-1); 

  //string foldermc("archive/2015_09_28_ana/skim/");
  string foldermc("/afs/cern.ch/user/m/manuelf/work/babies/2015_10_19/mc/skim_1lht500met200/");
  string folderdata("/afs/cern.ch/user/m/manuelf/work/babies/2015_10_19/data/singlelep/combined/skim_1lht500met200/");

  //Define processes. Try to minimize splitting
  Process ttbar{"ttbar", {
      {foldermc+"/*TTJets*.root/tree"}
    }};
  Process other{"other", {
      {foldermc+"/*DYJetsToLL*.root/tree"},
        {foldermc+"/*QCD_Pt*.root/tree"},
          {foldermc+"/*_TTWJets*.root/tree"},
	    {foldermc+"/*_TTZTo*.root/tree"},
	      {foldermc+"/*_ST_*.root/tree"},
		{foldermc+"/*_WJetsToLNu*.root/tree"},
		  {foldermc+"/*_WWTo*.root/tree"},
		    {foldermc+"/*ggZH_HToBB*.root/tree"},
		      {foldermc+"/*ttHJetTobb*.root/tree"}
    }};
  Process signal_nc{"signal", {
      {foldermc+"/*T1tttt*1500*100*.root/tree"}
    }};
  Process signal_c{"signal", {
      {foldermc+"/*T1tttt*1200*800*.root/tree"}
    }};
  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, "(trig[4]||trig[8])&&pass", true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline{"ht>500&&met>200&njets>="+minjets+"&&nbm>=2&&(nels+nmus)==1"};
  Cut baseline_david{"ht>450&&met>150&njets>=6&&nbm>=1&&(nels+nmus)==1"};

  //Declare bins
  //Method 2
  Bin r1_lowmet_lownb{"r1_lowmet_lownb", "mt<=140&&mj<=400&&met<="+himet+"&&nbm<=2"};
  Bin r1_lowmet_highnb{"r1_lowmet_highnb", "mt<=140&&mj<=400&&met<="+himet+"&&nbm>2"};
  Bin r1_highmet{"r1_highmet", "mt<=140&&mj<=400&&met>"+himet};

  Bin r2_lowmet_lownj_lownb{"r2_lowmet_lownj_lownb", "mt<=140&&mj>400&&met<="+himet+"&&njets<="+midjets+"&&nbm<=2"};
  Bin r2_lowmet_lownj_highnb{"r2_lowmet_lownj_highnb", "mt<=140&&mj>400&&met<="+himet+"&&njets<="+midjets+"&&nbm>2"};
  Bin r2_lowmet_highnj_lownb{"r2_lowmet_highnj_lownb", "mt<=140&&mj>400&&met<="+himet+"&&njets>"+midjets+"&&nbm<=2"};
  Bin r2_lowmet_highnj_highnb{"r2_lowmet_highnj_highnb", "mt<=140&&mj>400&&met<="+himet+"&&njets>"+midjets+"&&nbm>2"};
  Bin r2_highmet_lownj{"r2_highmet_lownj", "mt<=140&&mj>400&&met>"+himet+"&&njets<="+midjets};
  Bin r2_highmet_highnj{"r2_highmet_highnj", "mt<=140&&mj>400&&met>"+himet+"&&njets>"+midjets};

  Bin r3_lowmet_lownb{"r3_lowmet_lownb", "mt>140&&mj<=400&&met<="+himet+"&&nbm<=2"};
  Bin r3_lowmet_highnb{"r3_lowmet_highnb", "mt>140&&mj<=400&&met<="+himet+"&&nbm>2"};
  Bin r3_highmet{"r3_highmet", "mt>140&&mj<=400&&met>"+himet};

  Bin r4_lowmet_lownj_lownb{"r4_lowmet_lownj_lownb", "mt>140&&mj>400&&met<="+himet+"&&njets<="+midjets+"&&nbm<=2"};
  Bin r4_lowmet_lownj_highnb{"r4_lowmet_lownj_highnb", "mt>140&&mj>400&&met<="+himet+"&&njets<="+midjets+"&&nbm>2"};
  Bin r4_lowmet_highnj_lownb{"r4_lowmet_highnj_lownb", "mt>140&&mj>400&&met<="+himet+"&&njets>"+midjets+"&&nbm<=2"};
  Bin r4_lowmet_highnj_highnb{"r4_lowmet_highnj_highnb", "mt>140&&mj>400&&met<="+himet+"&&njets>"+midjets+"&&nbm>2"};
  Bin r4_highmet_lownj{"r4_highmet_lownj", "mt>140&&mj>400&&met>"+himet+"&&njets<="+midjets};
  Bin r4_highmet_highnj{"r4_highmet_highnj", "mt>140&&mj>400&&met>"+himet+"&&njets>"+midjets};

  //Method 1
  Bin m1_r1_lowmet_lownj{"m1_r1_lowmet_lownj", "mt<=140&&mj<=600&&met<="+himet+"&&njets<="+midjets};
  Bin m1_r1_lowmet_highnj{"m1_r1_lowmet_highnj", "mt<=140&&mj<=600&&met<="+himet+"&&njets>"+midjets};
  Bin m1_r1_highmet_lownj{"m1_r1_highmet_lownj", "mt<=140&&mj<=600&&met>"+himet+"&&njets<="+midjets};
  Bin m1_r1_highmet_highnj{"m1_r1_highmet_highnj", "mt<=140&&mj<=600&&met>"+himet+"&&njets>"+midjets};

  Bin m1_r2_lowmet_lownj{"m1_r2_lowmet_lownj", "mt<=140&&mj>600&&met<="+himet+"&&njets<="+midjets};
  Bin m1_r2_lowmet_highnj{"m1_r2_lowmet_highnj", "mt<=140&&mj>600&&met<="+himet+"&&njets>"+midjets};
  Bin m1_r2_highmet_lownj{"m1_r2_highmet_lownj", "mt<=140&&mj>600&&met>"+himet+"&&njets<="+midjets};
  Bin m1_r2_highmet_highnj{"m1_r2_highmet_highnj", "mt<=140&&mj>600&&met>"+himet+"&&njets>"+midjets};

  Bin m1_r3_lowmet_lownj{"m1_r3_lowmet_lownj", "mt>140&&mj<=600&&met<="+himet+"&&njets<="+midjets};
  Bin m1_r3_lowmet_highnj{"m1_r3_lowmet_highnj", "mt>140&&mj<=600&&met<="+himet+"&&njets>"+midjets};
  Bin m1_r3_highmet_lownj{"m1_r3_highmet_lownj", "mt>140&&mj<=600&&met>"+himet+"&&njets<="+midjets};
  Bin m1_r3_highmet_highnj{"m1_r3_highmet_highnj", "mt>140&&mj<=600&&met>"+himet+"&&njets>"+midjets};

  Bin m1_r4_lowmet_lownj{"m1_r4_lowmet_lownj", "mt>140&&mj>600&&met<="+himet+"&&njets<="+midjets};
  Bin m1_r4_lowmet_highnj{"m1_r4_lowmet_highnj", "mt>140&&mj>600&&met<="+himet+"&&njets>"+midjets};
  Bin m1_r4_highmet_lownj{"m1_r4_highmet_lownj", "mt>140&&mj>600&&met>"+himet+"&&njets<="+midjets};
  Bin m1_r4_highmet_highnj{"m1_r4_highmet_highnj", "mt>140&&mj>600&&met>"+himet+"&&njets>"+midjets};

  //Method David
  Bin r1{"r1", "mt<=140&&mj<=400"};
  Bin r2{"r2", "mt<=140&&mj>400"};
  Bin r3{"r3", "mt>140&&mj<=400"};
  Bin r4{"r4", "mt>140&&mj>400"};

  //Specify ABCD constraints
  set<Block> blocks_m2{
    {"lowmet_lownb", {{r1_lowmet_lownb, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb},
          {r3_lowmet_lownb, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb}}},
      {"lowmet_highnb", {{r1_lowmet_highnb, r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb},
            {r3_lowmet_highnb, r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb}}},
        {"highmet", {{r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
              {r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  set<Block> blocks_m1{
    {"lowmet_lownj", {{m1_r1_lowmet_lownj, m1_r2_lowmet_lownj},
          {m1_r3_lowmet_lownj, m1_r4_lowmet_lownj}}},
      {"lowmet_highnj", {{m1_r1_lowmet_highnj, m1_r2_lowmet_highnj},
            {m1_r3_lowmet_highnj, m1_r4_lowmet_highnj}}},
        {"highmet_lownj", {{m1_r1_highmet_lownj, m1_r2_highmet_lownj},
              {m1_r3_highmet_lownj, m1_r4_highmet_lownj}}},
          {"highmet_highnj", {{m1_r1_highmet_highnj, m1_r2_highmet_highnj},
                {m1_r3_highmet_highnj, m1_r4_highmet_highnj}}}
  };

  set<Block> blocks_david{
    {"all", {{r1, r2}, {r3, r4}}}
  };

  Process *psignal = &signal_nc;
  if(do_compressed) psignal = &signal_c;
  WorkspaceGenerator wg(baseline, blocks_m2, backgrounds, *psignal, data, "txt/systematics/method2.txt");
  if(!blinded){
    wg.SetBlindLevel(WorkspaceGenerator::BlindLevel::unblinded);
  }
  wg.SetLuminosity(lumi);
  wg.SetDoDilepton(true);
  wg.SetDoSystematics(true);

  TString lumi_s(""); lumi_s+=lumi; lumi_s.ReplaceAll(".","p");
  string signal = (do_compressed?"c":"nc");
  string outname("method2_"+signal+"_met"+himet+"_nj"+minjets+hijets+"_lumi"+lumi_s.Data()+".root");
  wg.WriteToFile(outname);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},
      {"unblind", no_argument, 0, 'u'},
      {"no_syst", no_argument, 0, 0},
      {"lowj", required_argument, 0, 'j'},
      {"hij", required_argument, 0, 'h'},
      {"himet", required_argument, 0, 'm'},
      {"compressed", required_argument, 0, 'c'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:uj:h:m:c", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 'u':
      blinded = false;
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
    case 'c':
      do_compressed = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "no_syst"){
        do_syst = false;
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
