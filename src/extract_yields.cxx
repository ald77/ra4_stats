#include "extract_yields.hpp"

#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <string>

#include "TFile.h"
#include "TIterator.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLegend.h"

#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"

#include "utilities.hpp"
#include "styles.hpp"

using namespace std;

int main(int argc, char *argv[]){
  if(argc < 2) return 1;
  styles style("RA4");
  string init_file = argv[1];

  string command = "combine -M MultiDimFit -H ProfileLikelihood --algo cross --saveWorkspace --floatAllNuisances 1 --robustFit 1 --profilingMode all --floatOtherPOIs 1 "+init_file;
  cout << command << endl;
  string output = execute(command);
  cout << output << endl;

  TFile w_file("higgsCombineTest.MultiDimFit.mH120.root", "read");
  if(!w_file.IsOpen()) return 1;
  RooWorkspace *w = static_cast<RooWorkspace*>(w_file.Get("w"));
  if(w == nullptr) return 1;
  if(!w->loadSnapshot("MultiDimFit")) return 1;
  RooRealVar *vr = w->var("r");
  double r = (vr != nullptr) ? vr->getVal() : -999.;
  ostringstream r_string;
  r_string << "r = " << r << flush;
  RooArgSet args = w->allVars();
  TIterator *iter_ptr = args.createIterator();
  if(iter_ptr == nullptr) return 1;

  map<string, Yields> ymap;
  map<string, double> kmap;

  for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
    RooAbsArg *var = static_cast<RooAbsArg*>(*(*iter_ptr));
    if(var == nullptr) continue;
    string name = var->GetName();
    if(StartsWith(name, "nobs_BLK_")
       || StartsWith(name, "nbkg_BLK_")
       || StartsWith(name, "nexp_BLK_")
       || StartsWith(name, "ymc_BLK_")
       || StartsWith(name, "kappamc_BLK_")){
      string base_name = name.substr(5);
      if(StartsWith(name, "ymc_BLK_")) base_name = name.substr(4);
      if(StartsWith(name, "kappamc_BLK_")) base_name = name.substr(8);

      RooAbsReal *vnobs = static_cast<RooAbsReal*>(w->arg(("nobs_"+base_name).c_str()));
      RooAbsReal *vnbkg = static_cast<RooAbsReal*>(w->arg(("nbkg_"+base_name).c_str()));
      RooAbsReal *vnexp = static_cast<RooAbsReal*>(w->arg(("nexp_"+base_name).c_str()));
      RooAbsReal *vnmc = static_cast<RooAbsReal*>(w->arg(("ymc_"+base_name).c_str()));
      RooAbsReal *vkappa = static_cast<RooAbsReal*>(w->arg(("kappamc_"+base_name).c_str()));
      if(vnobs == nullptr
	 || vnbkg == nullptr
	 || vnexp == nullptr
	 || vnmc == nullptr) continue;
      Yields yields;
      yields.nobs = vnobs->getVal();
      yields.nbkg = vnbkg->getVal();
      yields.nexp = vnexp->getVal();
      yields.nmc = vnmc->getVal();
      ymap[base_name] = yields;
      kmap[base_name] = (vkappa != nullptr) ? vkappa->getVal() : 1.;
    }
  }

  TH1D hobs("hobs", ";;Yield", ymap.size(), 0.5, ymap.size()+0.5);
  TH1D hbkg("hbkg", ";;Yield", ymap.size(), 0.5, ymap.size()+0.5);
  TH1D hexp("hexp", ";;Yield", ymap.size(), 0.5, ymap.size()+0.5);
  TH1D hmc("hmc", ";;Yield", ymap.size(), 0.5, ymap.size()+0.5);
  TH1D hkappa("hkappa", ";;#kappa", ymap.size(), 0.5, ymap.size()+0.5);

  size_t bin = 1;
  for(auto yield = ymap.cbegin(); yield != ymap.cend(); ++yield){
    hobs.SetBinContent(bin, yield->second.nobs);
    hbkg.SetBinContent(bin, yield->second.nbkg);
    hexp.SetBinContent(bin, yield->second.nexp);
    hmc.SetBinContent(bin, yield->second.nmc);

    ++bin;
  }
  bin = 1;
  for(auto kappa = kmap.cbegin(); kappa != kmap.cend(); ++kappa){
    hkappa.SetBinContent(bin, kappa->second);
    hkappa.GetXaxis()->SetBinLabel(bin, (kappa->first.substr(4)).c_str());
    ++bin;
  }

  double font_size = 0.1;
  double offset = 0.5;

  hkappa.SetLineColor(kBlue);
  hkappa.SetLineWidth(6);
  hkappa.SetLineStyle(1);
  hkappa.SetFillColor(0);
  hkappa.SetFillStyle(4000);
  hkappa.SetStats(0);
  hkappa.SetLabelSize(font_size, "XY");
  hkappa.SetTitleSize(font_size, "XY");
  hkappa.SetTitleOffset(offset, "Y");
  hkappa.GetXaxis()->LabelsOption("U");

  hobs.SetLineColor(kBlack);
  hobs.SetLineWidth(6);
  hobs.SetLineStyle(1);
  hobs.SetFillColor(0);
  hobs.SetFillStyle(4000);
  hobs.SetStats(false);
  hobs.SetLabelSize(font_size, "Y");
  hobs.SetTitleSize(font_size, "Y");
  hobs.SetTitleOffset(offset, "Y");

  hmc.SetLineColor(kBlue);
  hmc.SetLineWidth(3);
  hmc.SetLineStyle(2);
  hmc.SetFillColor(0);
  hmc.SetFillStyle(4000);

  hbkg.SetLineColor(0);
  hbkg.SetLineWidth(0);
  hbkg.SetLineStyle(0);
  hbkg.SetFillColor(kGreen+2);

  hexp.SetLineColor(0);
  hexp.SetLineWidth(0);
  hexp.SetLineStyle(0);
  hexp.SetFillColor(kRed+2);

  TCanvas c;
  c.cd();
  TPad bot_pad("bot_pad", "bot_pad", 0., 0., 1., 0.4);
  bot_pad.SetFillColor(0);
  bot_pad.SetFillStyle(4000);
  bot_pad.SetMargin(0.1, 0., 0.4, 0.);
  bot_pad.Draw();
  c.cd();
  TPad mid_pad("mid_pad", "mid_pad", 0., 0.4, 1., 0.9);
  mid_pad.SetFillColor(0);
  mid_pad.SetFillStyle(4000);
  mid_pad.SetMargin(0.1, 0., 0., 0.);
  mid_pad.Draw();
  c.cd();
  TPad top_pad("top_pad", "top_pad", 0., 0.9, 1., 1.);
  top_pad.SetFillColor(0);
  top_pad.SetFillStyle(4000);
  top_pad.SetMargin(0.1, 0., 0., 0.);
  top_pad.Draw();

  mid_pad.cd();
  mid_pad.SetLogy();
  hobs.Draw();
  hexp.Draw("same");
  hbkg.Draw("same");
  hobs.Draw("same");
  hmc.Draw("same");

  top_pad.cd();
  TLegend l(0.1, 0., 1., 1.);
  l.SetNColumns(3);
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.SetFillStyle(4000);
  l.AddEntry(&hobs, "Observed", "l");
  l.AddEntry(&hbkg, "Background Fit", "f");
  l.AddEntry(&hobs, r_string.str().c_str(), "");
  l.AddEntry(&hmc, "MC Yield", "l");
  l.AddEntry(&hexp, "Signal Fit", "f");
  l.Draw("same");

  bot_pad.cd();
  hkappa.Draw();
  TLine line(0.5, 1., kmap.size()+0.5, 1.);
  line.SetLineColor(kBlack);
  line.SetLineStyle(2);
  line.SetLineWidth(1);
  line.Draw("same");

  c.Print((init_file+"_yields.pdf").c_str());

  iter_ptr->Reset();
  delete iter_ptr;
  iter_ptr = nullptr;
}
