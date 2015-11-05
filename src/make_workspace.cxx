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
  double lumi = 0.135;
  double sig_strength = 0.;
  bool blinded = true;
  bool no_kappa = false;
  bool do_syst = true;
  bool use_r4 = false;
  string method("m135");
  string minjets("6");
  string hijets("9");
  string himet("400");
  string mjthresh("400");
  unsigned n_toys = 0;
}

int main(int argc, char *argv[]){
  cout << fixed << setprecision(2);
  GetOptions(argc, argv);
  string midjets(""); midjets += to_string(atoi(hijets.c_str())-1);
  string minjets2l(""); minjets2l += to_string(atoi(minjets.c_str())-1);
  string midjets2l(""); midjets2l += to_string(atoi(midjets.c_str())-1);

  string skim("skim_1lht500met200/");
  if(Contains(method, "m135")) skim = "skim_1lht400/";
  string foldermc("/afs/cern.ch/user/m/manuelf/work/babies/2015_10_19/mc/"+skim);
  string folderdata("/afs/cern.ch/user/m/manuelf/work/babies/2015_10_25/data/singlelep/combined/"+skim);

  //Define processes. Try to minimize splitting
  Process ttbar{"ttbar", {
      {foldermc+"/*TTJets*Lept*.root/tree"},
        {foldermc+"/*TTJets_HT*.root/tree"}
    }};
  Process other{"other", {
      {foldermc+"/*_WJetsToLNu*.root/tree"},
        {foldermc+"/*_TTWJets*.root/tree"},
          {foldermc+"/*_TTZTo*.root/tree"},
            {foldermc+"/*_ST_*.root/tree"},
              {foldermc+"/*DYJetsToLL*.root/tree"},
                {foldermc+"/*QCD_HT*.root/tree"},
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

  string data_cuts("(trig[4]||trig[8])&&pass&&nonblind");
  if(method=="m2l") {
    data_cuts = "(trig[4]||trig[8])&&pass";
  }
  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline{"mj>250&&ht>500&&met>200&&njets>="+minjets+"&&nbm>=2&&nleps==1"};
  Cut baseline1b{"mj>250&&ht>500&&met>200&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Cut baseline2l{"mj>250&&ht>500&&met>200&&met<=400&&nbm<=2"};
  Cut baseline_135{"mj>250&&ht>450&&met>150"};

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

  //Declare bins
  //Method 2, m1b, and m1bk
  Bin r1_lowmet_1b{"r1_lowmet_1b", "mt<=140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==1"};
  Bin r1_highmet_1b{"r1_highmet_1b", "mt<=140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm==1"};
  Bin r1_lowmet_2b{"r1_lowmet_2b", "mt<=140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==2"};
  Bin r1_lowmet_3b{"r1_lowmet_3b", "mt<=140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm>2"};
  Bin r1_highmet{"r1_highmet", "mt<=140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm>=2"};

  Bin r1_lowmet_allnb{"r1_lowmet_allnb", "mt<=140&&mj<="+mjthresh+"&&met<="+himet};
  Bin r1_highmet_allnb{"r1_highmet_allnb", "mt<=140&&mj<="+mjthresh+"&&met>"+himet};
  Bin r1_allnb{"r1_allnb", "mt<=140&&mj<="+mjthresh};

  Bin r2_lowmet_lownj_1b{"r2_lowmet_lownj_1b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r2_highmet_lownj_1b{"r2_highmet_lownj_1b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r2_lowmet_highnj_1b{"r2_lowmet_highnj_1b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r2_highmet_highnj_1b{"r2_highmet_highnj_1b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r2_lowmet_lownj_2b{"r2_lowmet_lownj_2b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==2"};
  Bin r2_lowmet_lownj_3b{"r2_lowmet_lownj_3b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm>2"};
  Bin r2_lowmet_highnj_2b{"r2_lowmet_highnj_2b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==2"};
  Bin r2_lowmet_highnj_3b{"r2_lowmet_highnj_3b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm>2"};
  Bin r2_highmet_lownj{"r2_highmet_lownj", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm>=2"};
  Bin r2_highmet_highnj{"r2_highmet_highnj", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm>=2"};

  Bin r3_lowmet_1b{"r3_lowmet_1b", "mt>140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==1"};
  Bin r3_highmet_1b{"r3_highmet_1b", "mt>140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm==1"};
  Bin r3_lowmet_2b{"r3_lowmet_2b", "mt>140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==2"};
  Bin r3_lowmet_3b{"r3_lowmet_3b", "mt>140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm>2"};
  Bin r3_highmet{"r3_highmet", "mt>140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm>=2"};

  Bin r3_lowmet_allnb{"r3_lowmet_allnb", "mt>140&&mj<="+mjthresh+"&&met<="+himet};
  Bin r3_highmet_allnb{"r3_highmet_allnb", "mt>140&&mj<="+mjthresh+"&&met>"+himet};
  Bin r3_allnb{"r3_allnb", "mt>140&&mj<="+mjthresh};

  Bin r4_lowmet_lownj_1b{"r4_lowmet_lownj_1b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r4_highmet_lownj_1b{"r4_highmet_lownj_1b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r4_lowmet_highnj_1b{"r4_lowmet_highnj_1b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r4_highmet_highnj_1b{"r4_highmet_highnj_1b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r4_lowmet_lownj_2b{"r4_lowmet_lownj_2b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==2"};
  Bin r4_lowmet_lownj_3b{"r4_lowmet_lownj_3b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm>2"};
  Bin r4_lowmet_highnj_2b{"r4_lowmet_highnj_2b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==2"};
  Bin r4_lowmet_highnj_3b{"r4_lowmet_highnj_3b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm>2"};
  Bin r4_highmet_lownj{"r4_highmet_lownj", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm>=2"};
  Bin r4_highmet_highnj{"r4_highmet_highnj", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm>=2"};

  // Dilepton blocks
  Bin r1c_allnj{"r1c_allnj", "mt<=140&&mj<="+mjthresh+"&&nbm>=1&&njets>="+minjets+"&&nleps==1"};

  Bin r2c_lownj{"r2c_lownj", "mt<=140&&mj>"+mjthresh+"&&nbm>=1&&njets>="+minjets+"&&njets<="+midjets+"&&nleps==1"};
  Bin r2c_highnj{"r2c_highnj", "mt<=140&&mj>"+mjthresh+"&&nbm>=1&&njets>"+midjets+"&&nleps==1"};

  Bin d3_allnj{"d3_allnj", "mj<="+mjthresh+"&&njets>="+minjets2l+"&&nbm<=2&&nleps==2"};

  Bin d4_lownj{"d4_lownj", "mj>"+mjthresh+"&&njets>="+minjets2l+"&&njets<="+midjets2l+"&&nbm<=2&&nleps==2"};
  Bin d4_highnj{"d4_highnj", "mj>"+mjthresh+"&&njets>="+minjets2l+"&&njets>"+midjets2l+"&&nbm<=2&&nleps==2"};

  //Method 135
  Bin r1{"r1", "mt<=140&&mj<="+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin r2{"r2", "mt<=140&&mj>"+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin r3{"r3", "mt>140&&mj<="+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin r4{"r4", "mt>140&&mj>"+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin d3{"d3", "mj<="+mjthresh+"&&njets>="+minjets2l+"&&nbm>=0&&nleps==2"};
  Bin d4{"d4", "mj>"+mjthresh+"&&njets>="+minjets2l+"&&nbm>=0&&nleps==2"};

  //// METHOD 1BK: Adding 1b, fat R1/R3 integrated over njets, nb, but not MET
  set<Block> blocks_1bk{
    {"lowmet", {{r1_lowmet_allnb, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, r2_lowmet_lownj_2b, r2_lowmet_highnj_2b,
            r2_lowmet_lownj_3b, r2_lowmet_highnj_3b},
          {r3_lowmet_allnb, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, r4_lowmet_lownj_2b, r4_lowmet_highnj_2b,
              r4_lowmet_lownj_3b, r4_lowmet_highnj_3b}}},
      {"highmet", {{r1_highmet_allnb, r2_highmet_lownj_1b, r2_highmet_highnj_1b, r2_highmet_lownj, r2_highmet_highnj},
            {r3_highmet_allnb, r4_highmet_lownj_1b, r4_highmet_highnj_1b, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  //// METHOD 1
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

  //// METHOD 2
  set<Block> blocks_m2{
    {"lowmet_2b", {{r1_lowmet_2b, r2_lowmet_lownj_2b, r2_lowmet_highnj_2b},
          {r3_lowmet_2b, r4_lowmet_lownj_2b, r4_lowmet_highnj_2b}}},
      {"lowmet_3b", {{r1_lowmet_3b, r2_lowmet_lownj_3b, r2_lowmet_highnj_3b},
            {r3_lowmet_3b, r4_lowmet_lownj_3b, r4_lowmet_highnj_3b}}},
        {"highmet", {{r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
              {r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  //// METHOD 1BKALL: Adding 1b, fat R1/R3 integrated over njets, nb, and MET
  set<Block> blocks_1bkall{
    {"all", {{r1_allnb, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, r2_lowmet_lownj_2b, r2_lowmet_highnj_2b,
            r2_lowmet_lownj_3b, r2_lowmet_highnj_3b, r2_highmet_lownj_1b, r2_highmet_highnj_1b,
            r2_highmet_lownj, r2_highmet_highnj},
          {r3_allnb, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, r4_lowmet_lownj_2b, r4_lowmet_highnj_2b,
              r4_lowmet_lownj_3b, r4_lowmet_highnj_3b, r4_highmet_lownj_1b, r4_highmet_highnj_1b,
              r4_highmet_lownj, r4_highmet_highnj}}}
  };

  //// METHOD 1BK: Adding 1b, sharing RmT across all dimensions
  set<Block> blocks_1b{
    {"all", {{r1_lowmet_1b, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b,
            r1_lowmet_2b, r2_lowmet_lownj_2b, r2_lowmet_highnj_2b,
            r1_lowmet_3b, r2_lowmet_lownj_3b, r2_lowmet_highnj_3b,
            r1_highmet_1b, r2_highmet_lownj_1b, r2_highmet_highnj_1b,
            r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
          {r3_lowmet_1b, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b,
              r3_lowmet_2b, r4_lowmet_lownj_2b, r4_lowmet_highnj_2b,
              r3_lowmet_3b, r4_lowmet_lownj_3b, r4_lowmet_highnj_3b,
              r3_highmet_1b, r4_highmet_lownj_1b, r4_highmet_highnj_1b,
              r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  //// METHOD 135: simple ABCD for 135 ipb
  set<Block> blocks_135{
    {"all", {{r1, r2}, {r3, r4}}}
  };

  set<Block> blocks_135_2l{
    {"all", {{r1, r2}, {d3, d4}}}
  };

  //// METHOD 2l: Dilepton closure
  set<Block> blocks_2l{
    {"all", {{r1c_allnj, r2c_lownj, r2c_highnj},
          {d3_allnj, d4_lownj, d4_highnj}}}
  };

  Cut *pbaseline(&baseline);
  set<Block> *pblocks(&blocks_m2);
  string sysfile("txt/systematics/method2.txt");
  if(method == "m1b"){
    pbaseline = &baseline1b;
    pblocks = &blocks_1b;
    sysfile = "txt/systematics/m1b.txt";
  } else if(method == "m1bk"){
    pbaseline = &baseline1b;
    pblocks = &blocks_1bk;
    sysfile = "txt/systematics/m1bk.txt";
  } else if(method == "m1bk_nodilep"){
    pbaseline = &baseline1b;
    pblocks = &blocks_1bk;
    sysfile = "txt/systematics/m1bk_nodilep.txt";
  } else if(method == "m1bkall"){
    pbaseline = &baseline1b;
    pblocks = &blocks_1bkall;
    sysfile = "txt/systematics/m1bkall.txt";
  }else if(method == "m2l"){
    pbaseline = &baseline2l;
    pblocks = &blocks_2l;
    sysfile = "txt/systematics/m2l.txt";
  }else if(method == "m135"){
    pbaseline = &baseline_135;
    pblocks = &blocks_135;
    sysfile = "txt/systematics/m135.txt";
  }else if(method == "m135_2l"){
    pbaseline = &baseline_135;
    pblocks = &blocks_135_2l;
    sysfile = "txt/systematics/m135_2l.txt";
  }

  TString lumi_s("_lumi"); lumi_s += lumi; lumi_s.ReplaceAll(".","p"); lumi_s.ReplaceAll("00000000000001","");
  TString sig_s("_sig"); sig_s += sig_strength; sig_s.ReplaceAll(".","p"); sig_s.ReplaceAll("00000000000001","");
  string outname(method+(do_syst ? "" : "_nosys")+(use_r4 ? "" : "_nor4")
                 +(no_kappa ? "_nokappa" : "")+string("_c_met")
                 +himet+"_mj"+mjthresh+"_nj"+minjets+hijets
                 +sig_s+lumi_s.Data()+".root");

  for(unsigned itoy = 0; itoy <= n_toys; ++itoy){
    // Compressed SUSY
    WorkspaceGenerator wgc(*pbaseline, *pblocks, backgrounds, signal_c, data, sysfile, use_r4, sig_strength);
    if(!blinded) wgc.SetBlindLevel(WorkspaceGenerator::BlindLevel::unblinded);
    wgc.SetKappaCorrected(!no_kappa);
    wgc.SetDoSystematics(do_syst);
    wgc.SetToyNum(itoy);
    wgc.SetLuminosity(lumi);
    wgc.SetDoDilepton(false); // Applying dilep syst in text file
    wgc.SetDoSystematics(do_syst);
    ReplaceAll(outname, "_nc_", "_c_");
    wgc.WriteToFile(outname);

    // Non-compressed SUSY
    WorkspaceGenerator wgnc(*pbaseline, *pblocks, backgrounds, signal_nc, data, sysfile, use_r4, sig_strength);
    if(!blinded) wgnc.SetBlindLevel(WorkspaceGenerator::BlindLevel::unblinded);
    wgnc.SetKappaCorrected(!no_kappa);
    wgnc.SetDoSystematics(do_syst);
    wgnc.SetToyNum(itoy);
    wgnc.SetLuminosity(lumi);
    wgnc.SetDoDilepton(false); // Applying dilep syst in text file
    wgnc.SetDoSystematics(do_syst);
    ReplaceAll(outname, "_c_", "_nc_");
    wgnc.WriteToFile(outname);
  }
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
      {"mj", required_argument, 0, 's'},
      {"nokappa", no_argument, 0, 'k'},
      {"method", required_argument, 0, 't'},
      {"use_r4", no_argument, 0, '4'},
      {"toys", required_argument, 0, 0},
      {"sig_strength", required_argument, 0, 'g'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:uj:h:m:s:kt:4g:", long_options, &option_index);
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
    case 'k':
      no_kappa = true;
      break;
    case '4':
      use_r4 = true;
      break;
    case 's':
      mjthresh = optarg;
      break;
    case 't':
      method = optarg;
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
