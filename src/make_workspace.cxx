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
  bool blinded = true;
  bool no_kappa = false;
  bool do_syst = true;
  string method("m135");
  string minjets("6");
  string hijets("9");
  string himet("400"); 
  string mjthresh("400");
}

int main(int argc, char *argv[]){
  cout << fixed << setprecision(2);
  GetOptions(argc, argv);
  string midjets(""); midjets += to_string(atoi(hijets.c_str())-1); 
  string minjets2l(""); minjets2l += to_string(atoi(minjets.c_str())-1); 
  string midjets2l(""); midjets2l += to_string(atoi(midjets.c_str())-1); 

  string foldermc("/afs/cern.ch/user/m/manuelf/work/babies/2015_10_19/mc/skim_1lht400/");
  string folderdata("/afs/cern.ch/user/m/manuelf/work/babies/2015_10_25/data/singlelep/combined/skim_1lht400/");

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
    lumi = 1.26;
  }
  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline{"ht>500&&met>200&&njets>="+minjets+"&&nbm>=2&&(nels+nmus)==1"};
  Cut baseline1b{"ht>500&&met>200&&njets>="+minjets+"&&nbm>=1&&(nels+nmus)==1"};
  Cut baseline2l{"ht>500&&met>200&&met<=400&&njets>="+minjets2l+"&&nbm>=1&&nbm<=2"};
  Cut baseline2l0{"ht>450&&met>150&&met<=400&&njets>="+minjets2l+"&&nbm>=0&&nbm<=2"};
  Cut baseline_135{"ht>450&&met>150"};

  //Declare bins
  //Method 2, m1b, and m1bk
  Bin r1_lowmet_1b{"r1_lowmet_1b", "mt<=140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==1"};
  Bin r1_highmet_1b{"r1_highmet_1b", "mt<=140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm==1"};
  Bin r1_lowmet_lownb{"r1_lowmet_lownb", "mt<=140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==2"};
  Bin r1_lowmet_highnb{"r1_lowmet_highnb", "mt<=140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm>2"};
  Bin r1_highmet{"r1_highmet", "mt<=140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm>=2"};

  Bin r1_lowmet_allnb{"r1_lowmet_allnb", "mt<=140&&mj<="+mjthresh+"&&met<="+himet};
  Bin r1_highmet_allnb{"r1_highmet_allnb", "mt<=140&&mj<="+mjthresh+"&&met>"+himet};
  Bin r1_allnb{"r1_allnb", "mt<=140&&mj<="+mjthresh};

  Bin r2_lowmet_lownj_1b{"r2_lowmet_lownj_1b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r2_highmet_lownj_1b{"r2_highmet_lownj_1b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r2_lowmet_highnj_1b{"r2_lowmet_highnj_1b", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r2_highmet_highnj_1b{"r2_highmet_highnj_1b", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r2_lowmet_lownj_lownb{"r2_lowmet_lownj_lownb", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==2"};
  Bin r2_lowmet_lownj_highnb{"r2_lowmet_lownj_highnb", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm>2"};
  Bin r2_lowmet_highnj_lownb{"r2_lowmet_highnj_lownb", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==2"};
  Bin r2_lowmet_highnj_highnb{"r2_lowmet_highnj_highnb", "mt<=140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm>2"};
  Bin r2_highmet_lownj{"r2_highmet_lownj", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm>=2"};
  Bin r2_highmet_highnj{"r2_highmet_highnj", "mt<=140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm>=2"};

  Bin r3_lowmet_1b{"r3_lowmet_1b", "mt>140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==1"};
  Bin r3_highmet_1b{"r3_highmet_1b", "mt>140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm==1"};
  Bin r3_lowmet_lownb{"r3_lowmet_lownb", "mt>140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm==2"};
  Bin r3_lowmet_highnb{"r3_lowmet_highnb", "mt>140&&mj<="+mjthresh+"&&met<="+himet+"&&nbm>2"};
  Bin r3_highmet{"r3_highmet", "mt>140&&mj<="+mjthresh+"&&met>"+himet+"&&nbm>=2"};

  Bin r3_lowmet_allnb{"r3_lowmet_allnb", "mt>140&&mj<="+mjthresh+"&&met<="+himet};
  Bin r3_highmet_allnb{"r3_highmet_allnb", "mt>140&&mj<="+mjthresh+"&&met>"+himet};
  Bin r3_allnb{"r3_allnb", "mt>140&&mj<="+mjthresh};

  Bin r4_lowmet_lownj_1b{"r4_lowmet_lownj_1b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r4_highmet_lownj_1b{"r4_highmet_lownj_1b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm==1"};
  Bin r4_lowmet_highnj_1b{"r4_lowmet_highnj_1b", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r4_highmet_highnj_1b{"r4_highmet_highnj_1b", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm==1"};
  Bin r4_lowmet_lownj_lownb{"r4_lowmet_lownj_lownb", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm==2"};
  Bin r4_lowmet_lownj_highnb{"r4_lowmet_lownj_highnb", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&nbm>2"};
  Bin r4_lowmet_highnj_lownb{"r4_lowmet_highnj_lownb", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm==2"};
  Bin r4_lowmet_highnj_highnb{"r4_lowmet_highnj_highnb", "mt>140&&mj>"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&nbm>2"};
  Bin r4_highmet_lownj{"r4_highmet_lownj", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&nbm>=2"};
  Bin r4_highmet_highnj{"r4_highmet_highnj", "mt>140&&mj>"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&nbm>=2"};


  // Dilepton blocks
  Bin r1c_allnb{"r1c_allnb", "mt<=140&&mj<="+mjthresh+"&&njets>="+minjets+"&&(nels+nmus)==1"};

  Bin r2c_lownj_1b{"r2c_lownj_1b", "mt<=140&&mj>"+mjthresh+"&&njets<="+midjets+"&&nbm==1&&njets>="+minjets+"&&(nels+nmus)==1"};
  Bin r2c_highnj_1b{"r2c_highnj_1b", "mt<=140&&mj>"+mjthresh+"&&njets>"+midjets+"&&nbm==1&&njets>="+minjets+"&&(nels+nmus)==1"};
  Bin r2c_lownj_lownb{"r2c_lownj_lownb", "mt<=140&&mj>"+mjthresh+"&&njets<="+midjets+"&&nbm==2&&njets>="+minjets+"&&(nels+nmus)==1"};
  Bin r2c_highnj_lownb{"r2c_highnj_lownb", "mt<=140&&mj>"+mjthresh+"&&njets>"+midjets+"&&nbm==2&&njets>="+minjets+"&&(nels+nmus)==1"};

  Bin d3_allnb{"d3_allnb", "mj<="+mjthresh+"&&(nels+nmus)==2"};

  Bin d4_lownj_1b{"d4_lownj_1b", "mj>"+mjthresh+"&&njets<="+midjets2l+"&&nbm==1&&(nels+nmus)==2"};
  Bin d4_highnj_1b{"d4_highnj_1b", "mj>"+mjthresh+"&&njets>"+midjets2l+"&&nbm==1&&(nels+nmus)==2"};
  Bin d4_lownj_lownb{"d4_lownj_lownb", "mj>"+mjthresh+"&&njets<="+midjets2l+"&&nbm==2&&(nels+nmus)==2"};
  Bin d4_highnj_lownb{"d4_highnj_lownb", "mj>"+mjthresh+"&&njets>"+midjets2l+"&&nbm==2&&(nels+nmus)==2"};

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

  //Method 135
  Bin r1{"r1", "mt<=140&&mj<="+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin r2{"r2", "mt<=140&&mj>"+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin r3{"r3", "mt>140&&mj<="+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin r4{"r4", "mt>140&&mj>"+mjthresh+"&&njets>="+minjets+"&&nbm>=1&&nleps==1"};
  Bin d3{"d3", "mj<="+mjthresh+"&&njets>="+minjets2l+"&&nbm>=0&&nleps==2"};
  Bin d4{"d4", "mj>"+mjthresh+"&&njets>="+minjets2l+"&&nbm>=0&&nleps==2"};

  //Specify ABCD constraints
  set<Block> blocks_m2{
    {"lowmet_lownb", {{r1_lowmet_lownb, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb},
          {r3_lowmet_lownb, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb}}},
      {"lowmet_highnb", {{r1_lowmet_highnb, r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb},
            {r3_lowmet_highnb, r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb}}},
        {"highmet", {{r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
              {r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  set<Block> blocks_1bkall{
    {"all", {{r1_allnb, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb, 
	    r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb, r2_highmet_lownj_1b, r2_highmet_highnj_1b, 
	    r2_highmet_lownj, r2_highmet_highnj},
          {r3_allnb, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb, 
	      r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb, r4_highmet_lownj_1b, r4_highmet_highnj_1b, 
	      r4_highmet_lownj, r4_highmet_highnj}}}
  };

  set<Block> blocks_2l{
    {"all", {{r1c_allnb, r2c_lownj_1b, r2c_highnj_1b, r2c_lownj_lownb, r2c_highnj_lownb},
          {d3_allnb, d4_lownj_1b, d4_highnj_1b, d4_lownj_lownb, d4_highnj_lownb}}}
  };

  set<Block> blocks_1bk{
    {"lowmet", {{r1_lowmet_allnb, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb, 
	    r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb},
          {r3_lowmet_allnb, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb, 
	      r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb}}},
      {"highmet", {{r1_highmet_allnb, r2_highmet_lownj_1b, r2_highmet_highnj_1b, r2_highmet_lownj, r2_highmet_highnj},
	    {r3_highmet_allnb, r4_highmet_lownj_1b, r4_highmet_highnj_1b, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  set<Block> blocks_mfs{
    {"all", {{r1_lowmet_lownb, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb,
	    r1_lowmet_highnb, r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb, 
	    r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
          {r3_lowmet_lownb, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb, 
	      r3_lowmet_highnb, r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb, 
	      r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
  };

  set<Block> blocks_1b{
    {"all", {{r1_lowmet_1b, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, 
	    r1_lowmet_lownb, r2_lowmet_lownj_lownb, r2_lowmet_highnj_lownb,
	    r1_lowmet_highnb, r2_lowmet_lownj_highnb, r2_lowmet_highnj_highnb, 
	    r1_highmet_1b, r2_highmet_lownj_1b, r2_highmet_highnj_1b,
	    r1_highmet, r2_highmet_lownj, r2_highmet_highnj},
          {r3_lowmet_1b, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, 
	      r3_lowmet_lownb, r4_lowmet_lownj_lownb, r4_lowmet_highnj_lownb, 
	      r3_lowmet_highnb, r4_lowmet_lownj_highnb, r4_lowmet_highnj_highnb, 
	      r3_highmet_1b, r4_highmet_lownj_1b, r4_highmet_highnj_1b,
	      r3_highmet, r4_highmet_lownj, r4_highmet_highnj}}}
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

  set<Block> blocks_135{
    {"all", {{r1, r2}, {r3, r4}}}
  };

  set<Block> blocks_135_2l{
    {"all", {{r1, r2}, {d3, d4}}}
  };

  Cut *pbaseline(&baseline);
  set<Block> *pblocks(&blocks_m2);
  string sysfile("txt/systematics/method2.txt");
  if(method == "m1b"){
    pbaseline = &baseline1b;
    pblocks = &blocks_1b;
    sysfile = "txt/systematics/m1b.txt";
  } else if(method == "mfs"){
    pblocks = &blocks_mfs;
  } else if(method == "m1bk"){
    pbaseline = &baseline1b;
    pblocks = &blocks_1bk;
    sysfile = "txt/systematics/m1bk.txt";
  } else if(method == "m1bkall"){
    pbaseline = &baseline1b;
    pblocks = &blocks_1bkall;
    sysfile = "txt/systematics/m1bkall.txt";
  }else if(method == "m2l"){
    pbaseline = &baseline2l;
    pblocks = &blocks_2l;
    sysfile = "txt/systematics/m2l.txt";
    //do_syst = false; 
  }else if(method == "m135"){
    pbaseline = &baseline_135;
    pblocks = &blocks_135;
    sysfile = "txt/systematics/m135.txt";
  }else if(method == "m135_2l"){
    pbaseline = &baseline_135;
    pblocks = &blocks_135_2l;
    sysfile = "txt/systematics/m135_2l.txt";
  }
  WorkspaceGenerator wgnc(*pbaseline, *pblocks, backgrounds, signal_nc, data, sysfile);
  if(!blinded){
    wgnc.SetBlindLevel(WorkspaceGenerator::BlindLevel::unblinded);
  }
  if(no_kappa){
    wgnc.SetKappaCorrected(false);
  }
  wgnc.SetLuminosity(lumi);
  wgnc.SetDoDilepton(false); // Applying dilep syst in text file
  wgnc.SetDoSystematics(do_syst);

  TString lumi_s(""); lumi_s+=lumi; lumi_s.ReplaceAll(".","p");
  string outname(method+(no_kappa ? "_nokappa" : "")+string("_nc_met")+himet+"_mj"+mjthresh+"_nj"+minjets+hijets
		 +"_lumi"+lumi_s.Data()+".root");
  wgnc.WriteToFile(outname);

  WorkspaceGenerator wgc(*pbaseline, *pblocks, backgrounds, signal_c, data, sysfile);
  if(!blinded){
    wgc.SetBlindLevel(WorkspaceGenerator::BlindLevel::unblinded);
  }
  if(no_kappa){
    wgc.SetKappaCorrected(false);
  }
  wgc.SetLuminosity(lumi);
  wgc.SetDoDilepton(false); // Applying dilep syst in text file
  wgc.SetDoSystematics(do_syst);
  ReplaceAll(outname, "_nc_", "_c_");
  wgc.WriteToFile(outname);

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
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:uj:h:m:s:kt:", long_options, &option_index);
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
    case 'k':
      no_kappa = false;
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
