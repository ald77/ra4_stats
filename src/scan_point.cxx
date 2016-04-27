#include "scan_point.hpp"

#include <cstdlib>

#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>

#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "utilities.hpp"
#include "cross_sections.hpp"

using namespace std;

namespace{
  string file_name = "";
  bool do_signif = false;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  if(file_name == "") ERROR("Must supply an input file name");

  TFile file(file_name.c_str(), "read");
  if(!file.IsOpen()) ERROR("Could not open "+file_name);

  string model = "T1tttt";
  if(Contains(file_name, "T5tttt")) model = "T5tttt";
  if(Contains(file_name, "T2tt")) model = "T2tt";
  if(Contains(file_name, "T6ttWW")) model = "T6ttWW";

  //// Parsing the gluino and LSP masses
  int mglu, mlsp;
  parseMasses(file_name, mglu, mlsp);
  float xsec, xsec_unc;
  if(model=="T1tttt" || model=="T5tttt") xsec::signalCrossSection(mglu, xsec, xsec_unc);
  else xsec::stopCrossSection(mglu, xsec, xsec_unc);
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));

  //string workdir = MakeDir("scan_point_"+glu_lsp);
  string workdir = "scan_point_"+model+"_"+glu_lsp+"/";
  gSystem->mkdir(workdir.c_str(), kTRUE);
 
  ostringstream command;
  string done = " < /dev/null &> /dev/null; ";
  done = "; ";
  //Need to get modify these file names
  string up_file_name = file_name;   ReplaceAll(up_file_name, "xsecNom", "xsecUp");
  string down_file_name = file_name; ReplaceAll(down_file_name, "xsecNom", "xsecDown");
  command
    << "export origdir=$(pwd); "
    << "ln -s $(readlink -f " << file_name << ") " << workdir << done
    << "ln -s $(readlink -f " << up_file_name << ") " << workdir << done
    << "ln -s $(readlink -f " << down_file_name << ") " << workdir << done
    << "cd " << workdir << done
    << "combine -M Asymptotic " << GetBaseName(file_name) << done
    << "combine -M Asymptotic --run observed --name Up " << GetBaseName(up_file_name) << done
    << "combine -M Asymptotic --run observed --name Down " << GetBaseName(down_file_name) << done;
  if(do_signif){
    command
      << "combine -M ProfileLikelihood --significance --expectSignal=1 --verbose=999999 " << GetBaseName(file_name)
      << " < /dev/null &> signif_obs.log; "
      << "combine -M ProfileLikelihood --significance --expectSignal=1 -t -1 --verbose=999999 " << GetBaseName(file_name)
      << " < /dev/null &> signif_exp.log; ";
  }
  command << flush;
  execute(command.str());
  
  string limits_file_name = workdir+"/higgsCombineTest.Asymptotic.mH120.root";
  TFile limits_file(limits_file_name.c_str(), "read");
  if(!limits_file.IsOpen()) ERROR("Could not open limits file "+limits_file_name);
  TTree *tree = static_cast<TTree*>(limits_file.Get("limit"));
  if(tree == nullptr) ERROR("Could not get limits tree");
  double limit;
  tree->SetBranchAddress("limit", &limit);
  int num_entries = tree->GetEntries();
  if(num_entries != 6) ERROR("Expected 6 tree entries. Saw "+to_string(num_entries));
  tree->GetEntry(1);
  double exp_down = limit;
  tree->GetEntry(2);
  double exp = limit;
  tree->GetEntry(3);
  double exp_up = limit;
  tree->GetEntry(5);
  double obs = limit;
  limits_file.Close();

  string up_limits_file_name = workdir+"/higgsCombineUp.Asymptotic.mH120.root";
  TFile up_limits_file(up_limits_file_name.c_str(), "read");
  if(!up_limits_file.IsOpen()) ERROR("No \"up\" file "+up_limits_file_name);
  tree = static_cast<TTree*>(up_limits_file.Get("limit"));
  if(tree == nullptr) ERROR("Could not get \"up\" limits tree");
  tree->SetBranchAddress("limit", &limit);
  num_entries = tree->GetEntries();
  if(num_entries != 1) ERROR("Expected 1 \"up\" tree entry. Saw "+to_string(num_entries));
  tree->GetEntry(0);
  double obs_up = limit;
  up_limits_file.Close();

  string down_limits_file_name = workdir+"/higgsCombineDown.Asymptotic.mH120.root";
  TFile down_limits_file(down_limits_file_name.c_str(), "read");
  if(!down_limits_file.IsOpen()) ERROR("No \"down\" file "+down_limits_file_name);
  tree = static_cast<TTree*>(down_limits_file.Get("limit"));
  if(tree == nullptr) ERROR("Could not get \"down\" limits tree");
  tree->SetBranchAddress("limit", &limit);
  num_entries = tree->GetEntries();
  if(num_entries != 1) ERROR("Expected 1 \"down\" tree entry. Saw "+to_string(num_entries));
  tree->GetEntry(0);
  double obs_down = limit;
  down_limits_file.Close();

  double sig_obs, sig_exp;
  if(do_signif){
    sig_obs = GetSignif(workdir+"/signif_obs.log");
    sig_exp = GetSignif(workdir+"/signif_exp.log");
  }

  //execute("rm -rf "+workdir);

  cout
    << setprecision(numeric_limits<double>::max_digits10)
    << ' ' << mglu
    << ' ' << mlsp
    << ' ' << xsec
    << ' ' << obs
    << ' ' << obs_up
    << ' ' << obs_down
    << ' ' << exp
    << ' ' << exp_up
    << ' ' << exp_down;
  if(do_signif){
    cout
      << ' ' << sig_obs
      << ' ' << sig_exp;
  }
  cout << endl;

  string txtname(workdir+"/limits_"+model+"_"+glu_lsp+".txt");
  ofstream txtfile(txtname);
  txtfile
    << setprecision(numeric_limits<double>::max_digits10)
    << ' ' << mglu
    << ' ' << mlsp
    << ' ' << xsec
    << ' ' << obs
    << ' ' << obs_up
    << ' ' << obs_down
    << ' ' << exp
    << ' ' << exp_up
    << ' ' << exp_down;
  if(do_signif){
    txtfile
      << ' ' << sig_obs
      << ' ' << sig_exp;
  }
  txtfile << endl;
}

double GetSignif(const string &filename){
  double signif = 0.;
  ifstream file(filename);
  string line;
  while(getline(file, line)){
    auto pos = line.find("Significance: ");
    if(pos != 0) continue;
    string val = line.substr(14);
    signif = stod(val);
  }
  return signif;
}

string GetBaseName(const string &path){
  auto pos = path.rfind("/");
  if(pos == string::npos){
    return path;
  }else{
    return path.substr(pos+1);
  }
}

double ExtractNumber(const string &results, const string &key){
  auto pos = results.find(key);
  if(pos != string::npos){
    pos += key.size();
    istringstream iss(results.substr(pos));
    double result;
    iss >> result;
    return result;
  }else{
    return -1.;
  }
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"filename", required_argument, 0, 'f'},
      {"signif", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:s", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      file_name = optarg;
      break;
    case 's':
      do_signif = true;
      break;
    default:
      cerr << "Bad option! getopt_long returned character code " << static_cast<int>(opt) << endl;
      break;
    }
  }
}
