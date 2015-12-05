#include "scan_point.hpp"

#include <cstdlib>

#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <limits>

#include <getopt.h>

#include "TFile.h"
#include "TTree.h"

#include "utilities.hpp"

using namespace std;

namespace{
  string file_name = "";
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  if(file_name == "") throw runtime_error("Must supply an input file name");

  string workdir = MakeDir("scan_point_");

  TFile file(file_name.c_str(), "read");
  if(!file.IsOpen()) throw runtime_error("Could not open "+file_name);
  //Need to read these from workspace file
  double x = 1500.;
  double y = 100.;
  double xsec = 0.0141903;
  file.Close();
  
  ostringstream command;
  string done = " < /dev/null &> /dev/null; ";
  //Need to get modify these file names
  string up_file_name = file_name;
  string down_file_name = file_name;
  command
    << "export origdir=$(pwd); "
    << "cp " << file_name << ' ' << workdir << done
    << "cp " << up_file_name << ' ' << workdir << done
    << "cp " << down_file_name << ' ' << workdir << done
    << "cd " << workdir << done
    << "combine -M Asymptotic " << file_name << done
    << "combine -M Asymptotic --run observed --name Up " << up_file_name << done
    << "combine -M Asymptotic --run observed --name Down " << down_file_name << done
    << flush;
  execute(command.str());
  
  string limits_file_name = workdir+"/higgsCombineTest.Asymptotic.mH120.root";
  TFile limits_file(limits_file_name.c_str(), "read");
  if(!limits_file.IsOpen()) throw runtime_error("Could not open limits file "+limits_file_name);
  TTree *tree = static_cast<TTree*>(limits_file.Get("limit"));
  if(tree == nullptr) throw runtime_error("Could not get limits tree");
  double limit;
  tree->SetBranchAddress("limit", &limit);
  int num_entries = tree->GetEntries();
  if(num_entries != 6) throw runtime_error("Expected 6 tree entries. Saw "+to_string(num_entries));
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
  if(!up_limits_file.IsOpen()) throw runtime_error("No \"up\" file "+up_limits_file_name);
  tree = static_cast<TTree*>(up_limits_file.Get("limit"));
  if(tree == nullptr) throw runtime_error("Could not get \"up\" limits tree");
  tree->SetBranchAddress("limit", &limit);
  num_entries = tree->GetEntries();
  if(num_entries != 1) throw runtime_error("Expected 1 \"up\" tree entry. Saw "+to_string(num_entries));
  tree->GetEntry(0);
  double obs_up = limit;
  up_limits_file.Close();

  string down_limits_file_name = workdir+"/higgsCombineDown.Asymptotic.mH120.root";
  TFile down_limits_file(down_limits_file_name.c_str(), "read");
  if(!down_limits_file.IsOpen()) throw runtime_error("No \"down\" file "+down_limits_file_name);
  tree = static_cast<TTree*>(down_limits_file.Get("limit"));
  if(tree == nullptr) throw runtime_error("Could not get \"down\" limits tree");
  tree->SetBranchAddress("limit", &limit);
  num_entries = tree->GetEntries();
  if(num_entries != 1) throw runtime_error("Expected 1 \"down\" tree entry. Saw "+to_string(num_entries));
  tree->GetEntry(0);
  double obs_down = limit;
  down_limits_file.Close();

  execute("rm -rf "+workdir);

  cout
    << setprecision(numeric_limits<double>::max_digits10)
    << ' ' << x
    << ' ' << y
    << ' ' << xsec
    << ' ' << obs
    << ' ' << obs_up
    << ' ' << obs_down
    << ' ' << exp
    << ' ' << exp_up
    << ' ' << exp_down
    << endl;
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
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      file_name = optarg;
      break;
    default:
      cerr << "Bad option! getopt_long returned character code " << static_cast<int>(opt) << endl;
      break;
    }
  }
}
