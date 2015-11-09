#include "sensitivity.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <initializer_list>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "TIterator.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"

#include "utilities.hpp"
#include "styles.hpp"

using namespace std;

namespace{
  string temp_name = "my_temp_name.root";
  double max_lim = 2.;
  double sig_scale = 1.1;
  double min_lumi = 0.;
  double lumi_increment = 0.5;
  double max_lumi = 6.;
  double lumi_in_file = 1.264;
  bool do_sys = false;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  styles style("LargeLabels");
  style.setDefaultStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.12);

  vector<string> files, names;
  if(!do_sys){
    files.push_back("m1bk_nc_met400_mj400_nj69_sig0_lumi1p264_unblinded.root"); names.push_back("T1tttt(1500,100)");
    files.push_back("m1bk_c_met400_mj400_nj69_sig0_lumi1p264_unblinded.root");  names.push_back("T1tttt(1200,800)");
  } else {
    files.push_back("m1bk_nc_met400_mj400_nj69_sig0_lumi1p264.root"); names.push_back("All systs");
    files.push_back("m1bk_nodilep_nc_met400_mj400_nj69_sig0_lumi1p264.root"); names.push_back("No dilep");
    files.push_back("m1bk_nosys_nc_met400_mj400_nj69_sig0_lumi1p264.root"); names.push_back("No systs");
  }

  vector<vector<double> > signif(files.size()), limit(files.size());
  vector<double> lumis;

  double max_sig = 0.;
  for(size_t ifile = 0; ifile < files.size(); ++ifile){
    for(double lumi = min_lumi; lumi <= max_lumi; lumi += lumi_increment){
      if(lumi <= 0) continue;
      if(ifile == 0) lumis.push_back(lumi);
      signif.at(ifile).push_back(GetSignificance(files.at(ifile), lumi));
      limit.at(ifile).push_back(GetLimit(files.at(ifile), lumi));
      cout << files.at(ifile) << ", " << lumi << ": "
	   << signif.at(ifile).back() << ", " << limit.at(ifile).back() << endl;
      if(signif.at(ifile).back() > max_sig){
	max_sig = signif.at(ifile).back();
      }
    }
  }
  max_sig = 5.0/sig_scale;
  for(size_t ifile = 0; ifile < files.size(); ++ifile){
    for(size_t ilumi = 0; ilumi < signif.at(ifile).size(); ++ilumi){
      signif.at(ifile).at(ilumi) *= max_lim/(sig_scale*max_sig);
      if(signif.at(ifile).at(ilumi) < 0.) signif.at(ifile).at(ilumi) = 0.;
      if(limit.at(ifile).at(ilumi) < 0.) limit.at(ifile).at(ilumi) = max_lim;
      if(limit.at(ifile).at(ilumi) > max_lim) limit.at(ifile).at(ilumi) = max_lim;
    }
  }

  vector<TGraph> signif_graphs, limit_graphs, dumb_graphs;
  for(size_t ifile = 0; ifile < files.size(); ++ifile){
    signif_graphs.push_back(TGraph(lumis.size(), &lumis.at(0), &signif.at(ifile).at(0)));
    limit_graphs.push_back(TGraph(lumis.size(), &lumis.at(0), &limit.at(ifile).at(0)));
    dumb_graphs.push_back(TGraph(0));

    signif_graphs.back().SetLineStyle(ifile+1);
    signif_graphs.back().SetLineWidth(5);
    signif_graphs.back().SetLineColor(kRed);

    limit_graphs.back().SetLineStyle(ifile+1);
    limit_graphs.back().SetLineWidth(5);
    limit_graphs.back().SetLineColor(kBlue);

    dumb_graphs.back().SetLineStyle(ifile+1);
    dumb_graphs.back().SetLineWidth(5);
    dumb_graphs.back().SetLineColor(kBlack);
  }

  TH1D h("h", ";Luminosity [fb^{-1}];Expected Limit/X-Section", 1, min_lumi, max_lumi);
  h.SetMaximum(max_lim);
  TAxis &yaxis = *h.GetYaxis();
  yaxis.SetLabelColor(kBlue);
  yaxis.SetTitleColor(kBlue);
  yaxis.SetTitleOffset(.8);
  TCanvas c;
  c.SetTicks(1,0);
  TGaxis *raxis = new TGaxis(max_lumi, 0.,
                             max_lumi, max_lim,
                             0., max_sig*sig_scale, 510, "+L");
  raxis->SetLabelColor(kRed);
  raxis->SetTitleColor(kRed);
  raxis->SetTitle("Expected Significance");
  raxis->SetTitleOffset(yaxis.GetTitleOffset());
  raxis->SetTitleSize(yaxis.GetTitleSize());
  //raxis->SetTextSize(yaxis.GetTextSize());
  raxis->SetLabelSize(yaxis.GetLabelSize());
  raxis->SetLabelFont(yaxis.GetLabelFont());
  raxis->SetTitleFont(yaxis.GetTitleFont());

  h.Draw("hist");
  TLegend l(1.0-gStyle->GetPadRightMargin()-0.4, 1.0-gStyle->GetPadTopMargin()-0.25,
	    1.0-gStyle->GetPadRightMargin(), 1.0-gStyle->GetPadTopMargin());
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.SetFillStyle(4000);
  for(size_t ifile = files.size() - 1; ifile < files.size(); --ifile){
    signif_graphs.at(ifile).Draw("samel");
    limit_graphs.at(ifile).Draw("samel");
    dumb_graphs.at(ifile).Draw("samel");
    l.AddEntry(&dumb_graphs.at(ifile), names.at(ifile).c_str(), "l");
  }
  l.Draw("same");
  raxis->Draw("same");

  string pname("sensitivity");  if(do_sys) pname += "_sys";
  pname += ".pdf";
  c.Print(pname.c_str());
}

double GetSignificance(const string &file_name, double lumi){
  ModifyLumi(file_name, lumi);
  string results = execute("combine -M ProfileLikelihood --significance --expectSignal=1 -t -1 "+temp_name);
  return ExtractNumber(results, "Significance: ");
}

double GetLimit(const string &file_name, double lumi){
  ModifyLumi(file_name, lumi);
  string results = execute("combine -M Asymptotic -t -1 "+temp_name);
  return ExtractNumber(results, "Median for expected limits: ");
}

void ModifyLumi(const string &file_name, double lumi){
  execute("cp "+file_name+" "+temp_name);
  TFile file(temp_name.c_str(), "read");
  RooWorkspace *w = static_cast<RooWorkspace*>(file.Get("w"));
  if(w == nullptr) return;
  const RooArgSet &vars = w->allVars();
  TIterator *iter_ptr = vars.createIterator();
  for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
    RooRealVar *var = static_cast<RooRealVar*>(*(*iter_ptr));
    if(var == nullptr) continue;
    if(lumi < 0.) continue;
    double ratio = lumi/lumi_in_file;
    string name = var->GetName();
    if(Contains(name, "norm_")
       || Contains(name, "nobs_")
       || Contains(name, "wmc_")){
      var->setVal(ratio*var->getVal());
      if(var->hasMin() && var->hasMax() && var->getMin()>=0. && var->getMax()>=var->getMax()){
	var->setMax(ratio*var->getMax());
      }
    }else if(Contains(name, "strength") && Contains(name, "dilep")){
      if(ratio>0.) var->setVal(var->getVal()/sqrt(ratio));
    }
  }
  w->writeToFile(temp_name.c_str());
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
      {"sys", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      do_sys = true;
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
