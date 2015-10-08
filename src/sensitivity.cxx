#include "sensitivity.hpp"

#include <vector>
#include <string>
#include <sstream>

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

string temp_name = "my_temp_name.root";
double min_lumi = 0.;
double lumi_increment = 0.5;
double max_lumi = 6.;
double min_sig = 0.;
double max_sig = 4.;
double min_lim = 0.;
double max_lim = 2.;
string file_nc = "methoddavidnc_30.root";
string file_c = "methoddavidc_30.root";

int main(){
  styles style("LargeLabels");
  style.setDefaultStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.12);

  vector<double> sig_nc, sig_c, lim_nc, lim_c, lumis;

  for(double lumi = min_lumi; lumi <= max_lumi; lumi += lumi_increment){
    lumis.push_back(lumi);
    sig_nc.push_back(GetSignificance(file_nc, lumi));
    lim_nc.push_back(GetLimit(file_nc, lumi));
    sig_c.push_back(GetSignificance(file_c, lumi));
    lim_c.push_back(GetLimit(file_c, lumi));
    cout << lumi << " " << sig_nc.back() << " " << lim_nc.back() << " " << sig_c.back() << " " << lim_c.back() << endl;
    FixSignificance(sig_nc.back());
    FixSignificance(sig_c.back());
    FixLimit(lim_nc.back());
    FixLimit(lim_c.back());
  }

  TGraph snc(lumis.size(), &lumis.at(0), &sig_nc.at(0));
  TGraph sc(lumis.size(), &lumis.at(0), &sig_c.at(0));
  TGraph lnc(lumis.size(), &lumis.at(0), &lim_nc.at(0));
  TGraph lc(lumis.size(), &lumis.at(0), &lim_c.at(0));
  snc.SetLineStyle(1);
  snc.SetLineWidth(5);
  TGraph dnc = snc;
  snc.SetLineColor(kRed);
  sc.SetLineStyle(2);
  sc.SetLineWidth(5);
  TGraph dc = sc;
  sc.SetLineColor(kRed);
  lnc.SetLineColor(kBlue);
  lnc.SetLineStyle(1);
  lnc.SetLineWidth(5);
  lc.SetLineColor(kBlue);
  lc.SetLineStyle(2);
  lc.SetLineWidth(5);
  TH1D h("h", ";Luminosity [fb^{-1}];Expected Limit/X-Section", 1, min_lumi, max_lumi);
  h.SetMaximum(max_lim);
  TAxis &yaxis = *h.GetYaxis();
  yaxis.SetLabelColor(kBlue);
  yaxis.SetTitleColor(kBlue);
  yaxis.SetTitleOffset(.8);
  TCanvas c;
  c.SetTicks(1,0);
  TGaxis *raxis = new TGaxis(max_lumi, min_lim,
                             max_lumi, max_lim,
                             min_sig, max_sig, 510, "+L");
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
  snc.Draw("samel");
  sc.Draw("samel");
  lnc.Draw("samel");
  lc.Draw("samel");
  raxis->Draw("same");
  TLegend l(1.0-gStyle->GetPadRightMargin()-0.4, 1.0-gStyle->GetPadTopMargin()-0.25,
	    1.0-gStyle->GetPadRightMargin(), 1.0-gStyle->GetPadTopMargin());
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.SetFillStyle(4000);
  l.AddEntry(&dnc, "T1tttt(1500,100)", "l");
  l.AddEntry(&dc, "T1tttt(1200,800)", "l");
  l.Draw("same");

  c.Print("sensitivity.pdf");
}

double GetSignificance(const string &file_name, double lumi){
  ModifyLumi(file_name, lumi);
  string results = execute("combine -M ProfileLikelihood --significance --expectSignal=1 -t -1 "+temp_name);
  return ExtractNumber(results, "Significance: ");
}

double GetLimit(const string &file_name, double lumi){
  ModifyLumi(file_name, lumi);
  string results = execute("combine -M Asymptotic "+temp_name);
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
    double ratio = lumi/3.0;
    string name = var->GetName();
    if(Contains(name, "norm_")
       || Contains(name, "nobs_")
       || (Contains(name, "rate_BLK_") && Contains(name, "_PRC_signal"))){
      var->setVal(ratio*var->getVal());
      if(var->hasMin() && var->hasMax() && var->getMin()>=0. && var->getMax()>=var->getMax()){
	var->setMax(ratio*var->getMax());
      }
    }else if(Contains(name, "strength_dilep_")){
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

void FixSignificance(double &x){
  if(x<min_sig) x = min_sig;
  if(x>max_sig) x = max_sig;
  double sig_scale = max_lim/max_sig;
  x *= sig_scale;
}

void FixLimit(double &x){
  if(x<=0. || x>max_lim) x = max_lim;
  if(x>0. && x<min_lim) x = min_lim;
}
