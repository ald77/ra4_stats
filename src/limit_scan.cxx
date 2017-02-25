#include "limit_scan.hpp"

#include <fstream>
#include <sstream>
#include <string>
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
#include "cross_sections.hpp"

using namespace std;

namespace{
  int Nsmooth = 4; // Number of times to smooth TH2D
  string filename = "txt/t1tttt_limit_scan.txt";
  string model = "T1tttt";
}
void reverseGraph(TGraph &graph);
void fixGraph(TGraph &graph);

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  styles style("2Dscan"); style.setDefaultStyle();
  SetupColors();
  
  if(filename == "") ERROR("No input file provided");
  vector<double> vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown, vsigobs, vsigexp;

  ifstream infile(filename);
  string line;

  while(getline(infile, line)){
    istringstream iss(line);
    double pmx, pmy, pxsec, pobs, pobsup, pobsdown, pexp, pup, pdown, sigobs, sigexp;
    iss >> pmx >> pmy >> pxsec >> pobs >> pobsup >> pobsdown >> pexp >> pup >> pdown >> sigobs >> sigexp;
    int mglu(pmx);
    float xsec, exsec;
    xsec::signalCrossSection(mglu, xsec, exsec);
    // int factor(50), mlsp(pmy);
    // if((mglu%factor!=0 || mlsp%factor!=0) && mglu-mlsp!=225 && mlsp!=1450) continue;
    // if(mglu-mlsp==225 && mglu%factor!=0) continue;
    vmx.push_back(pmx);
    vmy.push_back(pmy);
    if(Contains(model, "T1") || Contains(model, "T5")) vxsec.push_back(xsec);
    else vxsec.push_back(pxsec);
    vobs.push_back(pobs);
    vobsup.push_back(pobsup);
    vobsdown.push_back(pobsdown);
    vexp.push_back(pexp);
    vup.push_back(pup);
    vdown.push_back(pdown);
    vsigobs.push_back(sigobs);
    vsigexp.push_back(sigexp);
  }
  infile.close();

  if(vmx.size() <= 2) ERROR("Need at least 3 models to draw scan");
  if(vmx.size() != vmy.size()
     || vmx.size() != vxsec.size()
     || vmx.size() != vobs.size()
     || vmx.size() != vobsup.size()
     || vmx.size() != vobsdown.size()
     || vmx.size() != vexp.size()
     || vmx.size() != vup.size()
     || vmx.size() != vdown.size()
     || vmx.size() != vsigobs.size()
     || vmx.size() != vsigexp.size()) ERROR("Error parsing text file. Model point not fully specified");
  
  vector<double> vlim(vxsec.size());
  for(size_t i = 0; i < vxsec.size(); ++i){
    //vobs.at(i) /= 5;
    vlim.at(i) = vxsec.at(i) * vobs.at(i);
  }

  TGraph2D glim("glim", "Cross-Section Limit", vlim.size(), &vmx.at(0), &vmy.at(0), &vlim.at(0));
  TGraph2D gobs("gobs", "Observed Limit", vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));
  TGraph2D gobsup("gobsup", "Observed +1#sigma Limit", vobsup.size(), &vmx.at(0), &vmy.at(0), &vobsup.at(0));
  TGraph2D gobsdown("gobsdown", "Observed -1#sigma Limit", vobsdown.size(), &vmx.at(0), &vmy.at(0), &vobsdown.at(0));
  TGraph2D gexp("gexp", "Expected Limit", vexp.size(), &vmx.at(0), &vmy.at(0), &vexp.at(0));
  TGraph2D gup("gup", "Expected +1#sigma Limit", vup.size(), &vmx.at(0), &vmy.at(0), &vup.at(0));
  TGraph2D gdown("gdown", "Expected -1#sigma Limit", vdown.size(), &vmx.at(0), &vmy.at(0), &vdown.at(0));
  TGraph2D gsigobs("gsigobs", "Observed Significance", vsigobs.size(), &vmx.at(0), &vmy.at(0), &vsigobs.at(0));
  TGraph2D gsigexp("gsigexp", "Expected Significance", vsigexp.size(), &vmx.at(0), &vmy.at(0), &vsigexp.at(0));
  vector<double> vmx_excl, vmy_excl, vmx_incl, vmy_incl;
  for(size_t i = 0; i < vobs.size(); ++i){
    if(vobs.at(i) < 1.){
      vmx_excl.push_back(vmx.at(i));
      vmy_excl.push_back(vmy.at(i));
    }else{
      vmx_incl.push_back(vmx.at(i));
      vmy_incl.push_back(vmy.at(i));
    }
  }
  TGraph dots_excl(vmx_excl.size(), &vmx_excl.at(0), &vmy_excl.at(0));
  TGraph dots_incl(vmx_incl.size(), &vmx_incl.at(0), &vmy_incl.at(0));

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
  if(hlim == nullptr) ERROR("Could not retrieve histogram");
  TString xparticle("gluino"), yparticle("LSP");
  if(model=="T2tt") xparticle = "stop";
  if(model=="T6ttWW") {
    xparticle = "sbottom";
    yparticle = "chargino";
  }
  hlim->SetTitle(";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV]; 95% C.L. upper limit on cross section [pb]");
  hlim->SetMinimum(*min_element(vlim.cbegin(), vlim.cend()));
  hlim->SetMaximum(*max_element(vlim.cbegin(), vlim.cend()));
  
  TH2D *hsigobs = gsigobs.GetHistogram();
  if(hsigobs == nullptr) ERROR("Could not retrieve histogram");
  hsigobs->SetTitle(";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV]; Observed Significance");
  //hsigobs->SetMinimum(-2.);
  //hsigobs->SetMaximum(2.);
  
  TH2D *hsigexp = gsigexp.GetHistogram();
  if(hsigexp == nullptr) ERROR("Could not retrieve histogram");
  hsigexp->SetTitle(";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV]; Expected Significance");
  hsigexp->SetMinimum(0.);
  hsigexp->SetMaximum(6.);
  
  TCanvas c("","");
  c.SetLogz();
  //hlim->Draw("colz");
  glim.Draw("colz");
  //  glim.Draw("triw same");
  TLegend l(gStyle->GetPadLeftMargin(), 1.-gStyle->GetPadTopMargin(),
            1.-gStyle->GetPadRightMargin(), 1.);
  l.SetNColumns(2); l.SetTextSize(0.05);
  l.SetBorderSize(0);
  TGraph cup = DrawContours(gup, 2, 2);
  TGraph cdown = DrawContours(gdown, 2, 2);
  TGraph cexp = DrawContours(gexp, 2, 1, &l, "Expected");
  TGraph cobsup = DrawContours(gobsup, 1, 2);
  TGraph cobsdown = DrawContours(gobsdown, 1, 2);
  TGraph cobs = DrawContours(gobs, 1, 1, &l, "Observed");
  l.Draw("same");
  dots_excl.SetMarkerStyle(23); dots_excl.SetMarkerColor(28); dots_excl.SetMarkerSize(0.6);
  dots_incl.SetMarkerStyle(22); dots_incl.SetMarkerColor(28); dots_incl.SetMarkerSize(0.6);
  dots_excl.Draw("p same");
  dots_incl.Draw("p same");

  TString filebase = model+"_limit_scan";
  if(Nsmooth>0) {filebase += "_smooth"; filebase += Nsmooth;}
  c.Print(filebase+".pdf");

  c.SetLogz(false);

  gsigobs.Draw("colz");
  cup = DrawContours(gup, 2, 2);
  cdown = DrawContours(gdown, 2, 2);
  cexp = DrawContours(gexp, 2, 1);
  cobsup = DrawContours(gobsup, 1, 2);
  cobsdown = DrawContours(gobsdown, 1, 2);
  cobs = DrawContours(gobs, 1, 1);
  l.Draw("same");
  dots_excl.Draw("p same");
  dots_incl.Draw("p same");
  c.Print((model+"_sigobs.pdf").c_str());

  gsigexp.Draw("colz");
  cup = DrawContours(gup, 2, 2);
  cdown = DrawContours(gdown, 2, 2);
  cexp = DrawContours(gexp, 2, 1);
  cobsup = DrawContours(gobsup, 1, 2);
  cobsdown = DrawContours(gobsdown, 1, 2);
  cobs = DrawContours(gobs, 1, 1);
  l.Draw("same");
  dots_excl.Draw("p same");
  dots_incl.Draw("p same");
  c.Print((model+"_sigexp.pdf").c_str());

  TFile file(filebase+".root","recreate");

  cout<<endl<<"Saved limit curves in "<<filebase<<".root"<<endl<<endl;
  hlim->SetZTitle("");
  hlim->Write("hXsec_exp_corr");
  cobs.Write("graph_smoothed_Obs");
  cobsup.Write("graph_smoothed_ObsP");
  cobsdown.Write("graph_smoothed_ObsM");
  cexp.Write("graph_smoothed_Exp");
  cup.Write("graph_smoothed_ExpP");
  cdown.Write("graph_smoothed_ExpM");
  hsigobs->Write("hsig_obs_corr");
  hsigexp->Write("hsig_exp_corr");
  file.Close();
}

void Style(TGraph *g, int color, int style){
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(5);
}

TGraph DrawContours(TGraph2D &g2, int color, int style,
                    TLegend *leg, const string &name){
  TGraph graph;
  
  TList *l;
  //// Finding the TH2D, smoothing it, and creating a TGraph2D to get a new Delauny interpolation
  if(Nsmooth>0){
    g2.SetNpx(100);
    g2.SetNpy(100);
    TH2D *histo2d = g2.GetHistogram();
    TH2D *hclone = static_cast<TH2D*>(histo2d->Clone("clone"));
    for(int ind=0; ind<Nsmooth; ind++) histo2d->Smooth(1,"k5b");
    vector<double> vx, vy, vz;
    double glu_lsp = 225;
    for(int binx=1; binx<=histo2d->GetNbinsX(); binx++){
      for(int biny=1; biny<=histo2d->GetNbinsY(); biny++){
	double x = histo2d->GetXaxis()->GetBinCenter(binx);
	double y = histo2d->GetYaxis()->GetBinCenter(biny);
	double z = histo2d->GetBinContent(histo2d->GetBin(binx,biny));
	vx.push_back(x);
	vy.push_back(y);
	if(x-y>glu_lsp+85) vz.push_back(z);
	else vz.push_back(hclone->GetBinContent(hclone->GetBin(binx,biny)));
      }
    }
    TGraph2D gsmooth("gsmooth", "Cross-Section Limit", vx.size(), &vx.at(0), &vy.at(0), &vz.at(0));
    gsmooth.GetHistogram();
    l = gsmooth.GetContourList(1.);
  } else {
    g2.GetHistogram();
    l = g2.GetContourList(1.);
  }
  if(l == nullptr) return graph;
  bool added = false;
  int max_points = -1;
  for(int i = 0; i < l->GetSize(); ++i){
    TGraph *g = static_cast<TGraph*>(l->At(i));
    if(g == nullptr) continue;
    int n_points = g->GetN();
    if(n_points > max_points){
      if(Nsmooth>0) fixGraph(*g);
      graph = *g;
      max_points = n_points;
    }
    Style(g, color, style);
    g->Draw("L same");
    if(!added && leg != nullptr && name != ""){
      leg->AddEntry(g, name.c_str(), "l");
      added = true;
    }
  }


  return graph;
}

void fixGraph(TGraph &graph){
  double glu_lsp = 225;
  if(model=="T5tttt") glu_lsp = 265.;
  int np(graph.GetN());
  double mglu, iniglu, endglu, mlsp, inilsp, endlsp;
  //cout<<endl<<endl;
  for(int point(0); point < np; point++){
    graph.GetPoint(point, mglu, mlsp);
    //cout<<mglu<<", "<<mlsp<<endl;
  }
  graph.GetPoint(0, iniglu, inilsp);
  graph.GetPoint(np-1, endglu, endlsp);
  // Reversing graph if printed towards decreasing mgluino
  if(inilsp < endlsp) {
    reverseGraph(graph);
    endglu = iniglu;
    endlsp = inilsp;
  }
  // Adding a point so that it goes down to mLSP = 0, but not for WZ,SOS
  //cout<<"endlsp "<<endlsp<<", endglu "<<endglu<<endl;
  if(endlsp<30){
    graph.SetPoint(graph.GetN(), endglu, 0);
    np++;
  }

  reverseGraph(graph);
  // Adding a point at mLSP = 0, and removing points beyond the diagonal
  for(int point(0); point < np; point++){
    graph.GetPoint(point, mglu, mlsp);
    if(mlsp > mglu-glu_lsp-5){
      //cout<<"point "<<point<<", np "<<np<<endl;
      while(point <= graph.GetN() && point!=0) {
	graph.RemovePoint(graph.GetN()-1);
	np--;
      }
      break;
    }
  }
}
void reverseGraph(TGraph &graph){
  int np(graph.GetN());
  double mglu, mlsp;
  vector<double> mglus, mlsps;
  for(int point(np-1); point >= 0; point--){
    graph.GetPoint(point, mglu, mlsp);
    mglus.push_back(mglu);
    mlsps.push_back(mlsp);
  }
  for(int point(0); point < np; point++)
    graph.SetPoint(point, mglus[point], mlsps[point]);
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
      {"Nsmooth", required_argument, 0, 's'},
      {"model", required_argument, 0, 'm'},
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:m:s:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      Nsmooth = atoi(optarg);
      break;
    case 'm':
      model = optarg;
      break;
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
