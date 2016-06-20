#include "make_workspace.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <initializer_list>
#include <vector>
#include <string>
#include <stdlib.h>
#include <ctime>

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
  double lumi = 0.815;
  double sig_strength = 0.;
  BlindLevel blind_level = BlindLevel::unblinded;
  bool no_kappa = false;
  bool do_syst = false;
  bool use_r4 = false;
  TString method("m2lveto");
  string minjets("6");
  string hijets("9");
  string himet("400");
  string mjthresh("400");
  unsigned n_toys = 0;
  string identifier = "";
}

int main(int argc, char *argv[]){
  time_t begtime, endtime;
  time(&begtime);
  cout << fixed << setprecision(2);
  GetOptions(argc, argv);
  string midjets(""); midjets += to_string(atoi(hijets.c_str())-1);
  string minjets2l(""); minjets2l += to_string(atoi(minjets.c_str())-1);
  string midjets2l(""); midjets2l += to_string(atoi(midjets.c_str())-1);

  string hostname = execute("echo $HOSTNAME");
  string bfolder("/net/cms2/cms2r0/babymaker/");
  if(Contains(hostname, "lxplus")) bfolder = "/afs/cern.ch/user/m/manuelf/work/";
  TString foldermc(bfolder+"babies/2016_06_14/mc/merged_standard/");
  TString folderdata(bfolder+"babies/2016_06_14/data/skim_standard/");
  if(method.Contains("met150")){
    foldermc = bfolder+"babies/2016_06_14/mc/merged_1lht500met150nj5/";
    folderdata = bfolder+"babies/2016_06_14/data/merged_1lht500met150nj5/";
  }

  cout<<"folderdata is "<<folderdata<<endl;
  cout<<"foldermc is "<<foldermc<<endl;
  //Define processes. Try to minimize splitting
  string stitch_cuts("stitch");
  Process ttbar{"ttbar", {
      {foldermc+"/*TTJets*Lept*.root/tree"},
        {foldermc+"/*TTJets_HT*.root/tree"}
    },stitch_cuts};
  Process other{"other", {
      {foldermc+"/*_WJetsToLNu*.root/tree"},
        {foldermc+"/*_TTWJets*.root/tree"},
          {foldermc+"/*_TTZTo*.root/tree"},
            {foldermc+"/*_ST_*.root/tree"},
              {foldermc+"/*DYJetsToLL*.root/tree"},
                {foldermc+"/*QCD_HT*.root/tree"},
                  {foldermc+"/*_WWTo*.root/tree"},
                    {foldermc+"/*_TTGJets*.root/tree"},
		      {foldermc+"/*_TTTT*.root/tree"},
			{foldermc+"/*_WZ*.root/tree"},
			  {foldermc+"/*_ZZ*.root/tree"},
			   {foldermc+"/*_ZJet*.root/tree"},
			     {foldermc+"/*_WH_HToBB*.root/tree"},
			       {foldermc+"/*_ZH_HToBB*.root/tree"},
			       //	 {foldermc+"/*ggZH_HToBB*.root/tree"},
				 {foldermc+"/*ttHJetTobb*.root/tree"}
    },stitch_cuts};
  Process signal_nc{"signal", {
      {foldermc+"/*T1tttt*1500*-100_*.root/tree"}
    }, Cut(), false, true};
  Process signal_c{"signal", {
      {foldermc+"/*T1tttt*1200*800*.root/tree"}
    }, Cut(), false, true};

  string data_cuts("(trig[4]||trig[8]||trig[13]||trig[33])&&nonblind");
  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline2l{"pass&&ht>500&&met>150&&met<=500"};

  // Cuts regions
  string c_r1="mt<=140 && mj14>250&&mj14<=400";
  string c_r2="mt<=140 && mj14>400";
  string c_r3="mj14>250&&mj14<=400"; // mt>140 cut only applied to nveto==1
  string c_r4="mj14>400"; // mt>140 cut only applied to nveto==1

  // Cuts MET
  string c_lowmet="&&met>200&&met<=350";
  string c_midmet="&&met>350&&met<=500";
  if (method.Contains("met150")){
      c_lowmet="&&met>150&&met<=200"; 
    }
  cout<<"c_lowmet is "<<c_lowmet<<endl;
  // Cuts njets
  string c_allnj="&&njets>="+minjets;
  string c_lownj="&&njets>="+minjets+"&&njets<="+midjets;
  string c_hignj="&&njets>"+midjets;
  string c_allnj2l="&&njets>="+minjets2l;
  string c_lownj2l="&&njets>="+minjets2l+"&&njets<="+midjets2l;
  string c_hignj2l="&&njets>"+midjets2l;

  // Cuts leptons, veto, nbm
  string c_1l="&&nleps==1&&nveto==0&&nbm>=1";
  string c_2l="nleps==2&&nbm<=2";
  string c_veto="nleps==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140";

  // Cuts combination 2l + veto
  string c_2lvetoallnj="&&(("+c_2l+c_allnj2l+")||("+c_veto+c_allnj+"))";
  string c_2lvetolownj="&&(("+c_2l+c_lownj2l+")||("+c_veto+c_lownj+"))";
  string c_2lvetohignj="&&(("+c_2l+c_hignj2l+")||("+c_veto+c_hignj+"))";
  c_2l = "&&"+c_2l; c_veto = "&&"+c_veto;

  // Methods m2l, mveto, m2lveto
  string c_2ltotallnj = c_2l+c_allnj2l, c_2ltotlownj = c_2l+c_lownj2l, c_2ltothignj = c_2l+c_hignj2l;
  if(method.Contains("mveto")) {
    c_2ltotallnj = c_veto+c_allnj;
    c_2ltotlownj = c_veto+c_lownj;
    c_2ltothignj = c_veto+c_hignj;
  }
  if(method.Contains("m2lveto")) {
    c_2ltotallnj = c_2lvetoallnj;
    c_2ltotlownj = c_2lvetolownj;
    c_2ltothignj = c_2lvetohignj;
  }

  cout<<c_r1 + c_1l + c_lowmet + c_allnj<<endl;
  cout<<c_r2 + c_1l + c_lowmet + c_lownj<<endl;
  cout<<c_r2 + c_1l + c_lowmet + c_hignj<<endl;
  cout<< c_r3 + c_2ltotallnj + c_lowmet<<endl;
  cout<< c_r4 + c_2ltotlownj + c_lowmet<<endl;
  cout<< c_r4 + c_2ltothignj + c_lowmet<<endl;


  //////////////////////////////// Dilepton+veto blocks ////////////////////////////////
  // Low MET: 200 < MET <= 350
  Bin r1_lowmet_allnj{"r1_lowmet_allnj", c_r1 + c_1l + c_lowmet + c_allnj, blind_level>=BlindLevel::blinded};

  Bin r2_lowmet_lownj{"r2_lowmet_lownj", c_r2 + c_1l + c_lowmet + c_lownj, blind_level>=BlindLevel::blinded};
  Bin r2_lowmet_hignj{"r2_lowmet_hignj", c_r2 + c_1l + c_lowmet + c_hignj, blind_level>=BlindLevel::blinded};

  Bin d3_lowmet_allnj{"d3_lowmet_allnj", c_r3 + c_2ltotallnj + c_lowmet, blind_level>=BlindLevel::blinded};

  Bin d4_lowmet_lownj{"d4_lowmet_lownj", c_r4 + c_2ltotlownj + c_lowmet, blind_level>=BlindLevel::blinded};
  Bin d4_lowmet_hignj{"d4_lowmet_hignj", c_r4 + c_2ltothignj + c_lowmet, blind_level>=BlindLevel::blinded};


  // Mid MET: 350 < MET <= 500
  Bin r1_midmet_allnj{"r1_midmet_allnj", c_r1 + c_1l + c_midmet + c_allnj, blind_level>=BlindLevel::blinded};

  Bin r2_midmet_lownj{"r2_midmet_lownj", c_r2 + c_1l + c_midmet + c_lownj, blind_level>=BlindLevel::blinded};
  Bin r2_midmet_hignj{"r2_midmet_hignj", c_r2 + c_1l + c_midmet + c_hignj, blind_level>=BlindLevel::blinded};

  Bin d3_midmet_allnj{"d3_midmet_allnj", c_r3 + c_2ltotallnj + c_midmet, blind_level>=BlindLevel::blinded};

  Bin d4_midmet_lownj{"d4_midmet_lownj", c_r4 + c_2ltotlownj + c_midmet, blind_level>=BlindLevel::blinded};
  Bin d4_midmet_hignj{"d4_midmet_hignj", c_r4 + c_2ltothignj + c_midmet, blind_level>=BlindLevel::blinded};




  //// METHOD 2l: Dilepton closure
  set<Block> blocks_2l{
    {"lowmet", {{r1_lowmet_allnj, r2_lowmet_lownj, r2_lowmet_hignj},
          {d3_lowmet_allnj, d4_lowmet_lownj, d4_lowmet_hignj}}},
    {"midmet", {{r1_midmet_allnj, r2_midmet_lownj, r2_midmet_hignj},
          {d3_midmet_allnj, d4_midmet_lownj, d4_midmet_hignj}}}
  };

  if(method.Contains("met150")){
    blocks_2l = {
    {"lowmet", {{r1_lowmet_allnj, r2_lowmet_lownj, r2_lowmet_hignj},
          {d3_lowmet_allnj, d4_lowmet_lownj, d4_lowmet_hignj}}}
    };
  }

  Cut *pbaseline(&baseline2l);
  set<Block> *pblocks(&blocks_2l);
  string sysfile("txt/systematics/m1bk.txt");
  if(method == "m2l"){
    pbaseline = &baseline2l;
    pblocks = &blocks_2l;
    sysfile = "txt/systematics/m2l.txt";
  }

  TString lumi_s("_lumi"); lumi_s += lumi; lumi_s.ReplaceAll(".","p"); 
  lumi_s.ReplaceAll("00000000000001",""); lumi_s.ReplaceAll("39999999999999","4"); lumi_s.ReplaceAll("499999999999995","5");
  TString sig_s("_sig"); sig_s += sig_strength; sig_s.ReplaceAll(".","p"); 
  sig_s.ReplaceAll("00000000000001","");
  string blind_name = "";
  if(blind_level == BlindLevel::unblinded){
    blind_name = "_unblinded";
  }else if(blind_level == BlindLevel::unblind_1b){
    blind_name = "_1bunblinded";
  }else if(blind_level == BlindLevel::r4_blinded){
    blind_name = "_r4blinded";
  }
  if(identifier != "") identifier = "_" + identifier;
  string outname(method+(do_syst ? "" : "_nosys")+(use_r4 ? "" : "_nor4")
                 +(no_kappa ? "_nokappa" : "")+string("_c_met")
                 +himet+"_mj"+mjthresh+"_nj"+minjets+hijets
                 +sig_s+lumi_s.Data()+blind_name+identifier+".root");

  // Compressed SUSY
  WorkspaceGenerator wgc(*pbaseline, *pblocks, backgrounds, signal_c, data, sysfile, use_r4, sig_strength);
  wgc.SetKappaCorrected(!no_kappa);
  wgc.SetDoSystematics(do_syst);
  wgc.SetLuminosity(lumi);
  wgc.SetDoDilepton(false); // Applying dilep syst in text file
  wgc.SetDoSystematics(do_syst);
  wgc.AddToys(n_toys);
  ReplaceAll(outname, "_nc_", "_c_");
  wgc.WriteToFile(outname);

  // Non-compressed SUSY
  ReplaceAll(sysfile, "_c", "_nc");
  WorkspaceGenerator wgnc(*pbaseline, *pblocks, backgrounds, signal_nc, data, sysfile, use_r4, sig_strength);
  wgnc.SetKappaCorrected(!no_kappa);
  wgnc.SetDoSystematics(do_syst);
  wgnc.SetLuminosity(lumi);
  wgnc.SetDoDilepton(false); // Applying dilep syst in text file
  wgnc.SetDoSystematics(do_syst);
  wgnc.AddToys(n_toys);
  ReplaceAll(outname, "_c_", "_nc_");
  wgnc.WriteToFile(outname);

  time(&endtime); 
  cout<<"Finding workspaces took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},
      {"unblind", required_argument, 0, 'u'},
      {"do_syst", no_argument, 0, 0},
      {"lowj", required_argument, 0, 'j'},
      {"hij", required_argument, 0, 'h'},
      {"himet", required_argument, 0, 'm'},
      {"mj", required_argument, 0, 's'},
      {"nokappa", no_argument, 0, 'k'},
      {"method", required_argument, 0, 't'},
      {"use_r4", no_argument, 0, '4'},
      {"toys", required_argument, 0, 0},
      {"sig_strength", required_argument, 0, 'g'},
      {"identifier", required_argument, 0, 'i'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:u:j:h:m:s:kt:4g:i:", long_options, &option_index);
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
    case 'i':
      identifier = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "do_syst"){
        do_syst = true;
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
