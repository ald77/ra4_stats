#include "limit_scan.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLegend.h"

#include "utilities.hpp"
#include "styles.hpp"

using namespace std;

namespace{
  string filename = "";
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  styles style("RA4");
  style.setDefaultStyle();
  SetupColors();
  
  if(filename == "") throw runtime_error("No input file provided");
  vector<double> vmx, vmy, vxsec, vobs, vexp, vup, vdown;

  ifstream infile(filename);
  string line;
  bool start_read = false;
  while(getline(infile, line)){
    if(Contains(line, "---")){
      start_read = true;
    }else if(start_read){
      istringstream iss(line);
      double pmx, pmy, pxsec, pobs, pexp, pup, pdown;
      iss >> pmx >> pmy >> pxsec >> pobs >> pexp >> pup >> pdown;
      vmx.push_back(pmx);
      vmy.push_back(pmy);
      vxsec.push_back(pxsec);
      vobs.push_back(pobs);
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
     || vmx.size() != vexp.size()
     || vmx.size() != vup.size()
     || vmx.size() != vdown.size()) throw runtime_error("Error parsing text file. Model point not fully specified");
  
  vector<double> vlim(vxsec.size());
  for(size_t i = 0; i < vxsec.size(); ++i){
    vlim.at(i) = vxsec.at(i) * vobs.at(i);
  }

  TGraph2D glim(vlim.size(), &vmx.at(0), &vmy.at(0), &vlim.at(0));
  TGraph2D gobs(vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));
  TGraph2D gexp(vexp.size(), &vmx.at(0), &vmy.at(0), &vexp.at(0));
  TGraph2D gup(vup.size(), &vmx.at(0), &vmy.at(0), &vup.at(0));
  TGraph2D gdown(vdown.size(), &vmx.at(0), &vmy.at(0), &vdown.at(0));

  glim.SetNpx(500);
  glim.SetNpy(500);

  TH2D *hlim = glim.GetHistogram();
  if(hlim == nullptr) throw runtime_error("Could not retrieve histogram");
  hlim->SetTitle(";m_{gluino} [GeV];m_{LSP} [GeV]");
  
  TCanvas c;
  glim.Draw("colz");
  TLegend l(gStyle->GetPadLeftMargin(), 1.-gStyle->GetPadTopMargin(),
            1.-gStyle->GetPadRightMargin(), 1.);
  l.SetNColumns(2);
  l.SetBorderSize(0);
  DrawContours(gobs, 1, 1, &l, "Observed");
  DrawContours(gexp, 2, 1, &l, "Expected");
  DrawContours(gup, 2, 2);
  DrawContours(gdown, 2, 2);
  l.Draw("same");
  c.Print("limit_scan.pdf");
}

void DrawContours(TGraph2D &g2, int color, int style,
                  TLegend *leg, const string &name){
  TList *l = g2.GetContourList(1.);
  if(l == nullptr) return;
  bool added = false;
  for(int i = 0; i < l->GetSize(); ++i){
    TGraph *g = static_cast<TGraph*>(l->At(i));
    if(g == nullptr) continue;
    g->SetLineColor(color);
    g->SetLineStyle(style);
    g->SetLineWidth(5);
    g->Draw("L same");
    if(!added && leg != nullptr && name != ""){
      leg->AddEntry(g, name.c_str(), "l");
      added = true;
    }
  }
}
  
void SetupColors(){
  const int bands = 999;
  int colors[bands];
  const unsigned num = 6;
  double red[num] =   {1.,0.,0.,0.,1.,1.};
  double green[num] = {0.,0.,1.,1.,1.,0.};
  double blue[num] =  {1.,1.,1.,0.,0.,0.};
  double stops[num] = {0.,0.2,0.4,0.6,0.8,1.};
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
