#include "limit_scan.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"

#include "utilities.hpp"
#include "styles.hpp"

using namespace std;

namespace{
  string filename = "";
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  //styles style("RA4");
  //style.setDefaultStyle();
  SetupColors();
  
  if(filename == "") throw runtime_error("No input file provided");
  vector<double> vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown;

  ifstream infile(filename);
  string line;
  bool start_read = false;
  while(getline(infile, line)){
    if(Contains(line, "---")){
      start_read = true;
    }else if(start_read){
      istringstream iss(line);
      double pmx, pmy, pxsec, pobs, pobsup, pobsdown, pexp, pup, pdown;
      iss >> pmx >> pmy >> pxsec >> pobs >> pobsup >> pobsdown >> pexp >> pup >> pdown;
      vmx.push_back(pmx);
      vmy.push_back(pmy);
      vxsec.push_back(pxsec);
      vobs.push_back(pobs);
      vobsup.push_back(pobsup);
      vobsdown.push_back(pobsdown);
      vexp.push_back(pexp);
      vup.push_back(pup);
      vdown.push_back(pdown);
    }
  }
  infile.close();

  if(vmx.size() <= 2) throw runtime_error("Need at least 3 models to draw scan");
  if(vmx.size() != vmy.size()
     || vmx.size() != vxsec.size()
     || vmx.size() != vobs.size()
     || vmx.size() != vobsup.size()
     || vmx.size() != vobsdown.size()
     || vmx.size() != vexp.size()
     || vmx.size() != vup.size()
     || vmx.size() != vdown.size()) throw runtime_error("Error parsing text file. Model point not fully specified");
  
  vector<double> vlim(vxsec.size());
  for(size_t i = 0; i < vxsec.size(); ++i){
    vlim.at(i) = vxsec.at(i) * vobs.at(i);
  }

  TGraph2D glim("glim", "Cross-Section Limit", vlim.size(), &vmx.at(0), &vmy.at(0), &vlim.at(0));
  TGraph2D gobs("gobs", "Observed Limit", vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));
  TGraph2D gobsup("gobsup", "Observed +1#sigma Limit", vobsup.size(), &vmx.at(0), &vmy.at(0), &vobsup.at(0));
  TGraph2D gobsdown("gobsdown", "Observed -1#sigma Limit", vobsdown.size(), &vmx.at(0), &vmy.at(0), &vobsdown.at(0));
  TGraph2D gexp("gexp", "Expected Limit", vexp.size(), &vmx.at(0), &vmy.at(0), &vexp.at(0));
  TGraph2D gup("gup", "Expected +1#sigma Limit", vup.size(), &vmx.at(0), &vmy.at(0), &vup.at(0));
  TGraph2D gdown("gdown", "Expected -1#sigma Limit", vdown.size(), &vmx.at(0), &vmy.at(0), &vdown.at(0));
  TGraph dots(vmx.size(), &vmx.at(0), &vmy.at(0));

  double xmin = *min_element(vmx.cbegin(), vmx.cend());
  double xmax = *max_element(vmx.cbegin(), vmx.cend());
  double ymin = *min_element(vmy.cbegin(), vmy.cend());
  double ymax = *max_element(vmy.cbegin(), vmy.cend());
  double bin_size = 12.5;
  int nxbins = max(1, min(500, static_cast<int>(ceil((xmax-xmin)/bin_size))));
  int nybins = max(1, min(500, static_cast<int>(ceil((ymax-ymin)/bin_size))));
  glim.SetNpx(nxbins);
  glim.SetNpy(nybins);

  TH2D *hlim = glim.GetHistogram();
  if(hlim == nullptr) throw runtime_error("Could not retrieve histogram");
  hlim->SetTitle(";m_{gluino} [GeV];m_{LSP} [GeV]");
  
  TCanvas c("","",800,800);
  c.SetLogz();
  hlim->SetMinimum(*min_element(vlim.cbegin(), vlim.cend()));
  hlim->SetMaximum(*max_element(vlim.cbegin(), vlim.cend()));
  glim.Draw("colz");
  TLegend l(gStyle->GetPadLeftMargin(), 1.-gStyle->GetPadTopMargin(),
            1.-gStyle->GetPadRightMargin(), 1.);
  l.SetNColumns(2);
  l.SetBorderSize(0);
  TGraph cup = DrawContours(gup, 2, 2);
  TGraph cdown = DrawContours(gdown, 2, 2);
  TGraph cexp = DrawContours(gexp, 2, 1, &l, "Expected");
  TGraph cobsup = DrawContours(gobsup, 1, 2);
  TGraph cobsdown = DrawContours(gobsdown, 1, 2);
  TGraph cobs = DrawContours(gobs, 1, 1, &l, "Observed");
  l.Draw("same");
  dots.Draw("p same");
  c.Print("limit_scan.pdf");

  TFile file("limit_scan.root","recreate");
  hlim->Write("hXsec_exp_corr");
  cobs.Write("graph_smoothed_Obs");
  cobsup.Write("graph_smoothed_ObsP");
  cobsdown.Write("graph_smoothed_ObsM");
  cexp.Write("graph_smoothed_Exp");
  cup.Write("graph_smoothed_ExpP");
  cdown.Write("graph_smoothed_ExpM");
  file.Close();
}

TGraph DrawContours(TGraph2D &g2, int color, int style,
                    TLegend *leg, const string &name){
  TGraph out;
  TList *l = g2.GetContourList(1.);
  if(l == nullptr) return out;
  bool added = false;
  int max_points = -1;
  for(int i = 0; i < l->GetSize(); ++i){
    TGraph *g = static_cast<TGraph*>(l->At(i));
    if(g == nullptr) continue;
    int n_points = g->GetN();
    if(n_points > max_points){
      out = *g;
      max_points = n_points;
    }
    g->SetLineColor(color);
    g->SetLineStyle(style);
    g->SetLineWidth(5);
    g->Draw("L same");
    if(!added && leg != nullptr && name != ""){
      leg->AddEntry(g, name.c_str(), "l");
      added = true;
    }
  }
  return out;
}
  
void SetupColors(){
  const unsigned num = 5;
  const int bands = 255;
  int colors[bands];
  double stops[num] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[num] = {0.50, 0.50, 1.00, 1.00, 1.00};
  double green[num] = {0.50, 1.00, 1.00, 0.60, 0.50};
  double blue[num] = {1.00, 1.00, 0.50, 0.40, 0.50};
  /*const unsigned num = 6;
  double red[num] =   {1.,0.,0.,0.,1.,1.};
  double green[num] = {0.,0.,1.,1.,1.,0.};
  double blue[num] =  {1.,1.,1.,0.,0.,0.};
  double stops[num] = {0.,0.2,0.4,0.6,0.8,1.};*/
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    colors[i] = fi+i;
  }
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      filename = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == ""){
        printf("Bad option! Found option name %s\n", optname.c_str());
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
