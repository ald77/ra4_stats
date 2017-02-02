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
  double lumi = 36.8;//12.9;
  double sig_strength = 0.;
  BlindLevel blind_level = BlindLevel::unblinded;
  bool no_kappa = false;
  bool do_syst = true;
  bool use_r4 = true;
  bool applyVeto = true;
  string binning = "alternate";
  string minjets("6");
  string hijets("9");
  string medmet("350");
  string himet("400");
  string vhimet("500");
  string mjdef("mj14");
  string lowmjthresh("250");
  string mjthresh("400");
  unsigned n_toys = 0;
  string sigfile = "";
  string injfile = "";
  bool inject_other_model = false;
  bool dummy_syst = false;
  string dummy_syst_file = "";
  string outfolder = "out/";
  bool nom_only = false;
  string CSVM = "0.8484";
  bool use_pois = false;
}
//nbm = Sum$(jets_csv>CSVM&&jets_pt>30&&!jets_islep)
int main(int argc, char *argv[]){
  time_t begtime, endtime;
  time(&begtime);
  cout << fixed << setprecision(2);
  GetOptions(argc, argv);
  if(sigfile==""){
    cout<<endl<<"You need to specify the input file with -f. Exiting"<<endl<<endl;
    return 1;
  }
  string midjets = to_string(atoi(hijets.c_str())-1);
  string minjets2l = to_string(atoi(minjets.c_str())-1);
  string midjets2l = to_string(atoi(midjets.c_str())-1);

  string hostname = execute("echo $HOSTNAME");
  string basefolder("/net/cms2/cms2r0/babymaker/");
  if(Contains(hostname, "lxplus")) basefolder = "/afs/cern.ch/user/m/manuelf/work/";
  string foldermc(basefolder+"babies/2016_08_10/mc/merged_mcbase_abcd/"); 
  string folderdata(basefolder+"babies/2016_11_08/data/merged_database_abcd/");

  cout<<"binning is "<<binning<<endl;
  if(applyVeto) cout<<"apply veto"<<endl;
  else cout<<"no veto"<<endl;

  cout<<"mj is "<<mjdef<<endl;
  cout<<"lumi is "<<to_string(lumi)<<endl;
  cout<<"outfolder is "<<outfolder<<endl;

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
  Process signal{"signal", {
      {sigfile+"/tree"}
    },"stitch", false, true};
  Process injection{"injection", {
      {injfile+"/tree"}
    },"stitch", false, true};

  string data_cuts("trig_ra4&&pass");

  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline1b{"met/met_calo<5.&&pass_ra2_badmu&&st>500&&met>200&&nleps==1&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>=1&&njets>="+minjets+"&&"+mjdef+">"+lowmjthresh};
  string veto("");
  if(applyVeto)
    veto = "&&nveto==0";
 
  set<Block> blocks_1bk;

  //Declare bins 
  if(binning=="debug"){
    Bin r1{"r1", "mt<=140&&"+mjdef+"<="+mjthresh+veto+"&&met>350&&met<=500",
	blind_level>=BlindLevel::blinded};
    Bin r2{"r2", "mt<=140&&"+mjdef+">"+mjthresh+veto+"&&met>350&&met<=500&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2&&njets>=6&&njets<=8",
	blind_level>=BlindLevel::blinded};
    Bin r3{"r3", "mt>140&&"+mjdef+"<="+mjthresh+veto+"&&met>350&&met<=500",
	blind_level>=BlindLevel::blinded};
    Bin r4{"r4", "mt>140&&"+mjdef+">"+mjthresh+veto+"&&met>350&&met<=500&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2&&njets>=6&&njets<=8",
	blind_level>BlindLevel::blinded};
    blocks_1bk = {
      {"inclusive", {{r1,r2},{r3,r4}}}
    };
  }else if(binning=="nominal"){
    Bin r1_lowmet_allnb{"r1_lowmet_allnb", "mt<=140&&"+mjdef+"<="+mjthresh+"&&met<="+himet,
	blind_level>=BlindLevel::blinded};
    Bin r1_highmet_allnb{"r1_highmet_allnb", "mt<=140&&"+mjdef+"<="+mjthresh+"&&met>"+himet,
	blind_level>=BlindLevel::blinded};

    Bin r2_lowmet_lownj_1b{"r2_lowmet_lownj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1",
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_highnj_1b{"r2_lowmet_highnj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1",
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_lownj_2b{"r2_lowmet_lownj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2",
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_lownj_3b{"r2_lowmet_lownj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2",
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_highnj_2b{"r2_lowmet_highnj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2",
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_highnj_3b{"r2_lowmet_highnj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2",
	blind_level>=BlindLevel::blinded};

    Bin r2_highmet_lownj_1b{"r2_highmet_lownj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1",
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_highnj_1b{"r2_highmet_highnj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1",
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_lownj_2b{"r2_highmet_lownj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>=2",
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_highnj_2b{"r2_highmet_highnj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>=2",
	blind_level>=BlindLevel::blinded};

    Bin r3_lowmet_allnb{"r3_lowmet_allnb", "mt>140&&"+mjdef+"<="+mjthresh+"&&met<="+himet+veto,
	blind_level>=BlindLevel::blinded};
    Bin r3_highmet_allnb{"r3_highmet_allnb", "mt>140&&"+mjdef+"<="+mjthresh+"&&met>"+himet+veto,
	blind_level>=BlindLevel::blinded};

    Bin r4_lowmet_lownj_1b{"r4_lowmet_lownj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_lowmet_highnj_1b{"r4_lowmet_highnj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_lowmet_lownj_2b{"r4_lowmet_lownj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_lowmet_lownj_3b{"r4_lowmet_lownj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_lowmet_highnj_2b{"r4_lowmet_highnj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_lowmet_highnj_3b{"r4_lowmet_highnj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};

    Bin r4_highmet_lownj_1b{"r4_highmet_lownj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_highmet_highnj_1b{"r4_highmet_highnj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_highmet_lownj_2b{"r4_highmet_lownj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>=2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_highmet_highnj_2b{"r4_highmet_highnj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+himet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>=2"+veto,
	blind_level>BlindLevel::unblinded};

    //// METHOD 1BK: Adding 1b, fat R1/R3 integrated over njets, nb, but not MET
    blocks_1bk = {
      {"lowmet", {{r1_lowmet_allnb, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, r2_lowmet_lownj_2b, r2_lowmet_highnj_2b,
		   r2_lowmet_lownj_3b, r2_lowmet_highnj_3b},
		  {r3_lowmet_allnb, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, r4_lowmet_lownj_2b, r4_lowmet_highnj_2b,
		   r4_lowmet_lownj_3b, r4_lowmet_highnj_3b}}},
      {"highmet", {{r1_highmet_allnb, r2_highmet_lownj_1b, r2_highmet_highnj_1b, r2_highmet_lownj_2b, r2_highmet_highnj_2b},
		   {r3_highmet_allnb, r4_highmet_lownj_1b, r4_highmet_highnj_1b, r4_highmet_lownj_2b, r4_highmet_highnj_2b}}}
    };
  }
  else if(binning=="alternate"){
    
    Bin r1_lowmet_allnb{"r1_lowmet_allnb", "mt<=140&&"+mjdef+"<="+mjthresh+"&&met<="+medmet+veto,
	blind_level>=BlindLevel::blinded};
    Bin r1_medmet_allnb{"r1_medmet_allnb", "mt<=140&&"+mjdef+"<="+mjthresh+"&&met>"+medmet+"&&met<="+vhimet+veto,
	blind_level>=BlindLevel::blinded};
    Bin r1_highmet_allnb{"r1_highmet_allnb", "mt<=140&&"+mjdef+"<="+mjthresh+"&&met>"+vhimet+veto,
	blind_level>=BlindLevel::blinded};

    Bin r2_lowmet_lownj_1b{"r2_lowmet_lownj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_highnj_1b{"r2_lowmet_highnj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_lownj_2b{"r2_lowmet_lownj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_lownj_3b{"r2_lowmet_lownj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_highnj_2b{"r2_lowmet_highnj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_lowmet_highnj_3b{"r2_lowmet_highnj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>=BlindLevel::blinded};


    Bin r2_medmet_lownj_1b{"r2_medmet_lownj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_medmet_highnj_1b{"r2_medmet_highnj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_medmet_lownj_2b{"r2_medmet_lownj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_medmet_lownj_3b{"r2_medmet_lownj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_medmet_highnj_2b{"r2_medmet_highnj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_medmet_highnj_3b{"r2_medmet_highnj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>=BlindLevel::blinded};

    Bin r2_highmet_lownj_1b{"r2_highmet_lownj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_highnj_1b{"r2_highmet_highnj_1b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_lownj_2b{"r2_highmet_lownj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_lownj_3b{"r2_highmet_lownj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_highnj_2b{"r2_highmet_highnj_2b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>=BlindLevel::blinded};
    Bin r2_highmet_highnj_3b{"r2_highmet_highnj_3b", "mt<=140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>=BlindLevel::blinded};

    Bin r3_lowmet_allnb{"r3_lowmet_allnb", "mt>140&&"+mjdef+"<="+mjthresh+"&&met<="+medmet+veto,
	blind_level>=BlindLevel::blinded};
    Bin r3_medmet_allnb{"r3_medmet_allnb", "mt>140&&"+mjdef+"<="+mjthresh+"&&met>"+medmet+"&&met<="+vhimet+veto,
	blind_level>=BlindLevel::blinded};
    Bin r3_highmet_allnb{"r3_highmet_allnb", "mt>140&&"+mjdef+"<="+mjthresh+"&&met>"+vhimet+veto,
	blind_level>=BlindLevel::blinded};

    Bin r4_lowmet_lownj_1b{"r4_lowmet_lownj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_lowmet_highnj_1b{"r4_lowmet_highnj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_lowmet_lownj_2b{"r4_lowmet_lownj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_lowmet_lownj_3b{"r4_lowmet_lownj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_lowmet_highnj_2b{"r4_lowmet_highnj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_lowmet_highnj_3b{"r4_lowmet_highnj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};


    Bin r4_medmet_lownj_1b{"r4_medmet_lownj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_medmet_highnj_1b{"r4_medmet_highnj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_medmet_lownj_2b{"r4_medmet_lownj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_medmet_lownj_3b{"r4_medmet_lownj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_medmet_highnj_2b{"r4_medmet_highnj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_medmet_highnj_3b{"r4_medmet_highnj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met<="+vhimet+"&&met>"+medmet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};

    Bin r4_highmet_lownj_1b{"r4_highmet_lownj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_highmet_highnj_1b{"r4_highmet_highnj_1b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==1"+veto,
	blind_level>BlindLevel::unblind_1b};
    Bin r4_highmet_lownj_2b{"r4_highmet_lownj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_highmet_lownj_3b{"r4_highmet_lownj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets<="+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_highmet_highnj_2b{"r4_highmet_highnj_2b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)==2"+veto,
	blind_level>BlindLevel::unblinded};
    Bin r4_highmet_highnj_3b{"r4_highmet_highnj_3b", "mt>140&&"+mjdef+">"+mjthresh+"&&met>"+vhimet+"&&njets>"+midjets+"&&Sum$(jets_csv>"+CSVM+"&&jets_pt>30&&!jets_islep)>2"+veto,
	blind_level>BlindLevel::unblinded};
    
    blocks_1bk = {
      {"lowmet", {{r1_lowmet_allnb, r2_lowmet_lownj_1b, r2_lowmet_highnj_1b, r2_lowmet_lownj_2b, r2_lowmet_highnj_2b,
		   r2_lowmet_lownj_3b, r2_lowmet_highnj_3b},
		  {r3_lowmet_allnb, r4_lowmet_lownj_1b, r4_lowmet_highnj_1b, r4_lowmet_lownj_2b, r4_lowmet_highnj_2b,
		   r4_lowmet_lownj_3b, r4_lowmet_highnj_3b}}},
      {"medmet", {{r1_medmet_allnb, r2_medmet_lownj_1b, r2_medmet_highnj_1b, r2_medmet_lownj_2b, r2_medmet_highnj_2b,
		   r2_medmet_lownj_3b, r2_medmet_highnj_3b},
		  {r3_medmet_allnb, r4_medmet_lownj_1b, r4_medmet_highnj_1b, r4_medmet_lownj_2b, r4_medmet_highnj_2b,
		   r4_medmet_lownj_3b, r4_medmet_highnj_3b}}},
      {"highmet", {{r1_highmet_allnb, r2_highmet_lownj_1b, r2_highmet_highnj_1b, r2_highmet_lownj_2b, r2_highmet_highnj_2b,
		    r2_highmet_lownj_3b, r2_highmet_highnj_3b},
		   {r3_highmet_allnb, r4_highmet_lownj_1b, r4_highmet_highnj_1b, r4_highmet_lownj_2b, r4_highmet_highnj_2b,
		    r4_highmet_lownj_3b, r4_highmet_highnj_3b}}}
    };
  }


  //// Parsing the gluino and LSP masses
  int mglu, mlsp;
  parseMasses(sigfile, mglu, mlsp);
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));
  
  //// Creating workspaces for the Nominal, uncert Up, and uncert Down signal cross sections
  Cut *pbaseline(&baseline1b);
  set<Block> *pblocks(&blocks_1bk);
  string model = "T1tttt";
  
  string sysfolder = "/net/cms2/cms2r0/babymaker/sys/2017_02_01/T1tttt/";
  //Protect default
  if(binning=="nominal" && lumi < 3) sysfolder = "/net/cms2/cms2r0/babymaker/sys/2016_01_11/scan/";
  
  if(Contains(hostname, "lxplus")) sysfolder = "txt/systematics/";
  if(Contains(sigfile, "T5tttt")) {
    sysfolder = "/net/cms2/cms2r0/babymaker/sys/2016_02_09/T5tttt/";
    model = "T5tttt";
  }
  if(Contains(sigfile, "T2tt")) {
    sysfolder = "/net/cms2/cms2r0/babymaker/sys/2016_02_09/T2tt/";
    model = "T2tt";
  }
  if(Contains(sigfile, "T6ttWW")) {
    sysfolder = "/net/cms2/cms2r0/babymaker/sys/2016_02_09/T6ttWW/";
    model = "T6ttWW";
  }
  cout<<"sysfolder is "<<sysfolder<<endl;
  
  //  string sysfile(sysfolder+"sys_SMS-"+model+"_"+glu_lsp+"_"+to_string(lumi)+"ifb");
  string sysfile(sysfolder+"sys_SMS-"+model+"_"+glu_lsp+"_12.9ifb");
  if(binning=="alternate") sysfile+="_altbins.txt";
  else sysfile+="_nominal.txt";
  
  if(binning=="nominal" && lumi < 3) sysfile = sysfolder+"sys_SMS-"+model+"_"+glu_lsp+".txt";
  if(dummy_syst) sysfile = dummy_syst_file;
  cout<<"sysfile is "<<sysfile<<endl;
  // If systematic file does not exist, use m1bk_nc for tests
  struct stat buffer;   
  if(stat (sysfile.c_str(), &buffer) != 0) {
    cout<<endl<<"WARNING: "<<sysfile<<" does not exist. Using ";
    sysfile = "txt/systematics/m1bk_nc.txt";
    cout<<sysfile<<" instead"<<endl<<endl;
  }
  
  // Cross sections
  float xsec, xsec_unc;
  if(model=="T1tttt" || model=="T5tttt") xsec::signalCrossSection(mglu, xsec, xsec_unc);
  else xsec::stopCrossSection(mglu, xsec, xsec_unc);
  double rmax = 20.;
  if(mglu <= 1500 && mlsp <= 800){
    rmax = 5.;
    if(mglu <= 1200 && mlsp <= 550){
      rmax = 1.25;
      if(mglu <= 900 && mlsp <= 350){
	rmax = 0.5;
      }
    }
  }
  
  gSystem->mkdir(outfolder.c_str(), kTRUE);
  string outname(outfolder+"/wspace_"+model+"_"+glu_lsp+"_xsecNom.root");
  if(!use_r4) ReplaceAll(outname, "wspace_","wspace_nor4_");
  if(no_kappa) ReplaceAll(outname, "wspace_","wspace_nokappa_");
  if(!do_syst) ReplaceAll(outname, "wspace_","wspace_nosyst_");
  ReplaceAll(CSVM,"0.","");
  ReplaceAll(outname,"wspace_","wspace_"+CSVM+"_");

  WorkspaceGenerator wgNom(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1.);
  wgNom.UseGausApprox(!use_pois);
  wgNom.SetRMax(rmax);
  wgNom.SetKappaCorrected(!no_kappa);
  wgNom.SetLuminosity(lumi);
  wgNom.SetDoSystematics(do_syst);
  if(inject_other_model){
    wgNom.SetInjectionModel(injection);
  }
  wgNom.AddToys(n_toys);
  wgNom.WriteToFile(outname);

  if(!nom_only){
    ReplaceAll(outname, "Nom", "Up");
    WorkspaceGenerator wgUp(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1+xsec_unc);
    wgUp.UseGausApprox(!use_pois);
    wgUp.SetRMax(rmax);
    wgUp.SetKappaCorrected(!no_kappa);
    wgUp.SetLuminosity(lumi);
    wgUp.SetDoSystematics(do_syst);
    if(inject_other_model){
      wgUp.SetInjectionModel(injection);
    }
    wgUp.AddToys(n_toys);
    wgUp.WriteToFile(outname);

    ReplaceAll(outname, "Up", "Down");
    WorkspaceGenerator wgDown(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1-xsec_unc);
    wgDown.UseGausApprox(!use_pois);
    wgDown.SetRMax(rmax);
    wgDown.SetKappaCorrected(!no_kappa);
    wgDown.SetLuminosity(lumi);
    wgDown.SetDoSystematics(do_syst);
    if(inject_other_model){
      wgDown.SetInjectionModel(injection);
    }
    wgDown.AddToys(n_toys);
    wgDown.WriteToFile(outname);
  }

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
      {"lowmjthresh", required_argument, 0, 'a'},
      {"mj_var", required_argument, 0, 'd'},
      {"nokappa", no_argument, 0, 'k'},
      {"no_r4", no_argument, 0, '4'},
      {"useVeto", required_argument, 0, 'v'},
      {"alt_binning", required_argument, 0, 'b'},
      {"toys", required_argument, 0, 0},
      {"sig_strength", required_argument, 0, 'g'},
      {"dummy_syst", required_argument, 0, 0},
      {"outfolder", required_argument, 0, 'o'},
      {"inject", required_argument, 0, 'i'},
      {"nominal", no_argument, 0, 'n'},
      {"csvthres", required_argument, 0, 'c'},
      {"poisson", no_argument, 0, 'p'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:u:j:h:m:s:a:d:k4b:v:g:f:o:i:nc:p", long_options, &option_index);

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
      use_r4 = false;
      break;
    case 'v':
      if(string(optarg) =="on") applyVeto = true;
      break;
    case 'b':
      binning = optarg;
      break;
    case 's':
      mjthresh = optarg;
      break;
    case 'a':
      lowmjthresh = optarg;
      break;
    case 'd':
      mjdef = optarg;
      break;
    case 'i':
      injfile = optarg;
      inject_other_model = true;
      break;
    case 'n':
      nom_only = true;
      break;
    case 'c':
      CSVM = optarg;
      break;
    case 'p':
      use_pois = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "no_syst"){
        do_syst = false;
      }else if(optname == "toys"){
        n_toys = atoi(optarg);
      }else if(optname == "dummy_syst"){
	dummy_syst = true;
	dummy_syst_file = optarg;
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
