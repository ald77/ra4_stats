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
#include "TLatex.h"
#include "TLine.h"

#include "utilities.hpp"
#include "styles.hpp"
#include "cross_sections.hpp"

using namespace std;

namespace{
  int num_smooth_ = 4; // Number of times to smooth TH2D
  string filename_ = "txt/t1tttt_limit_scan.txt";
  string model_ = "T1tttt";
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  styles style("2Dscan"); style.setDefaultStyle();

  if(filename_ == "") ERROR("No input file provided");

  vector<double> vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown, vsigobs, vsigexp;
  ReadPoints(vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown, vsigobs, vsigexp);

  vector<double> vlim(vxsec.size());
  for(size_t i = 0; i < vxsec.size(); ++i){
    vlim.at(i) = vxsec.at(i) * vobs.at(i);
  }

  TH2D hsigobs = MakeObservedSignificancePlot(vmx, vmy, vsigobs);
  TH2D hsigexp = MakeExpectedSignificancePlot(vmx, vmy, vsigexp);

  MakeLimitPlot(vmx, vmy, vlim,
                vobs, vobsup, vobsdown,
                vexp, vup, vdown,
                hsigobs, hsigexp);
}

void ReadPoints(vector<double> &vmx,
		vector<double> &vmy,
		vector<double> &vxsec,
		vector<double> &vobs,
		vector<double> &vobsup,
		vector<double> &vobsdown,
		vector<double> &vexp,
		vector<double> &vup,
		vector<double> &vdown,
		vector<double> &vsigobs,
		vector<double> &vsigexp){
  ifstream infile(filename_);
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
    if(Contains(model_, "T1") || Contains(model_, "T5")) vxsec.push_back(xsec);
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

  if(vmx.size() <= 2) ERROR("Need at least 3 model_s to draw scan");
  if(vmx.size() != vmy.size()
     || vmx.size() != vxsec.size()
     || vmx.size() != vobs.size()
     || vmx.size() != vobsup.size()
     || vmx.size() != vobsdown.size()
     || vmx.size() != vexp.size()
     || vmx.size() != vup.size()
     || vmx.size() != vdown.size()
     || vmx.size() != vsigobs.size()
     || vmx.size() != vsigexp.size()) ERROR("Error parsing text file. Model_ point not fully specified");
}

TH2D MakeObservedSignificancePlot(vector<double> vmx,
                                  vector<double> vmy,
                                  vector<double> vobs){
  SetupSignedColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"};Observed Significance";

  TGraph2D g("", title.c_str(), vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));

  double the_max = 0.;
  for(int i = 0; i < g.GetN(); ++i){
    double z = fabs(g.GetZ()[i]);
    if(z>the_max) the_max = z;
  }
  g.SetMinimum(-the_max);
  g.SetMaximum(the_max);

  g.SetNpx(GetNumBins(vmx, 12.5));
  g.SetNpy(GetNumBins(vmy, 12.5));

  g.GetHistogram()->SetTitle(title.c_str());
  g.GetHistogram()->SetTickLength(0., "Z");

  TCanvas c;
  c.cd();

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
	      "#font[62]{CMS}#scale[0.76]{#font[52]{ Preliminary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
	       "#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);

  g.Draw("colz");

  vector<TLine> lines;
  for(double z = 0.; z < the_max; z+=0.5){
    int style = min(5, static_cast<int>(2.*fabs(z))+1);
    double width = max(1., 3.-fabs(z));
    DrawContours(g, 1, style, width, 0, z);
    double x1 = 1.-c.GetRightMargin()+0.0047;
    double x2 = 1.-c.GetRightMargin()+0.05;
    double ybot = c.GetBottomMargin();
    double ytop = 1.-c.GetTopMargin();
    double zpos = ybot+(ytop-ybot)*(the_max+z)/(2.*the_max);
    lines.emplace_back(x1, zpos, x2, zpos);
    lines.back().SetLineColor(1);
    lines.back().SetLineStyle(style);
    lines.back().SetLineWidth(width);
    lines.back().SetNDC(true);
    if(z != 0.){
      DrawContours(g, 1, style, width, 0, -z);
      double zneg = ybot+(ytop-ybot)*(the_max-z)/(2.*the_max);
      lines.emplace_back(x1, zneg, x2, zneg);
      lines.back().SetLineColor(1);
      lines.back().SetLineStyle(style);
      lines.back().SetLineWidth(width);
      lines.back().SetNDC(true);
    }
  }
  for(auto &l: lines){
    l.Draw("same");
  }
  
  ltitle.Draw("same");
  rtitle.Draw("same");

  c.Print((model_+"_sigobs.pdf").c_str());

  TH2D h = *g.GetHistogram();
  h.SetTitle("Observed Significance");
  return h;
}

TH2D MakeExpectedSignificancePlot(vector<double> vmx,
                                  vector<double> vmy,
                                  vector<double> vobs){
  SetupColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"}; Expected Significance";

  TGraph2D g("", title.c_str(), vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));

  double the_max = 0.;
  for(int i = 0; i < g.GetN(); ++i){
    double z = g.GetZ()[i];
    if(z>the_max) the_max = z;
  }
  if(the_max > 6.) the_max = 6.;
  g.SetMinimum(0.);
  g.SetMaximum(the_max);

  g.SetNpx(GetNumBins(vmx, 12.5));
  g.SetNpy(GetNumBins(vmy, 12.5));

  g.GetHistogram()->SetTitle(title.c_str());

  TCanvas c;
  c.cd();

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
	      "#font[62]{CMS}#scale[0.76]{#font[52]{ Preliminary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
	       "#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);

  g.Draw("colz");

  for(double z = 0.; z < the_max; z+=1.){
    int style = (z == 5. ? 1 : 2);
    double width = (z == 5. ? 2. : 1.);
    DrawContours(g, 1, style, width, 0, z);
  }
  
  ltitle.Draw("same");
  rtitle.Draw("same");
  
  c.Print((model_+"_sigexp.pdf").c_str());

  TH2D h = *g.GetHistogram();
  h.SetTitle("Expected Significance");
  return h;
}

void MakeLimitPlot(vector<double> vmx,
                   vector<double> vmy,
                   vector<double> vlim,
                   vector<double> vobs,
                   vector<double> vobsup,
                   vector<double> vobsdown,
                   vector<double> vexp,
                   vector<double> vup,
                   vector<double> vdown,
                   const TH2D &hsigobs,
                   const TH2D &hsigexp){
  SetupColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"};95% CL upper limit on cross section [pb]";
  
  TGraph2D glim("", title.c_str(), vlim.size(), &vmx.at(0), &vmy.at(0), &vlim.at(0));
  TGraph2D gobs("", "Observed Limit", vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));
  TGraph2D gobsup("", "Observed +1#sigma Limit", vobsup.size(), &vmx.at(0), &vmy.at(0), &vobsup.at(0));
  TGraph2D gobsdown("", "Observed -1#sigma Limit", vobsdown.size(), &vmx.at(0), &vmy.at(0), &vobsdown.at(0));
  TGraph2D gexp("", "Expected Limit", vexp.size(), &vmx.at(0), &vmy.at(0), &vexp.at(0));
  TGraph2D gup("", "Expected +1#sigma Limit", vup.size(), &vmx.at(0), &vmy.at(0), &vup.at(0));
  TGraph2D gdown("", "Expected -1#sigma Limit", vdown.size(), &vmx.at(0), &vmy.at(0), &vdown.at(0));

  glim.SetMinimum(0.001);
  glim.SetMaximum(2.);

  glim.SetNpx(GetNumBins(vmx, 12.5));
  glim.SetNpy(GetNumBins(vmy, 12.5));

  glim.SetTitle(title.c_str());

  TLegend l(gStyle->GetPadLeftMargin(), 1.-2.*gStyle->GetPadTopMargin(),
            1.-gStyle->GetPadRightMargin(), 1.-gStyle->GetPadTopMargin());
  l.SetNColumns(2);
  l.SetTextSize(0.05);
  l.SetBorderSize(0);
  
  TCanvas c;
  c.cd();

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
                "#font[62]{CMS}#scale[0.76]{#font[52]{ Preliminary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
                "#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);

  c.SetTopMargin(2.*c.GetTopMargin());
  c.SetLogz();
  glim.Draw("colz");

  TGraph cup = DrawContours(gup, 2, 2, 5, num_smooth_);
  TGraph cdown = DrawContours(gdown, 2, 2, 5, num_smooth_);
  TGraph cexp = DrawContours(gexp, 2, 1, 5, num_smooth_, 1.);
  TGraph cobsup = DrawContours(gobsup, 1, 2, 5, num_smooth_);
  TGraph cobsdown = DrawContours(gobsdown, 1, 2, 5, num_smooth_);
  TGraph cobs = DrawContours(gobs, 1, 1, 5, num_smooth_, 1.);

  l.AddEntry(&cexp, "Expected", "l");
  l.AddEntry(&cobs, "Observed", "l");

  l.Draw("same");

  ltitle.Draw("same");
  rtitle.Draw("same");
  
  string filebase = model_+"_limit_scan";
  if(num_smooth_>0){
    filebase += "_smooth";
    filebase += to_string(num_smooth_);
  }

  c.Print((filebase+".pdf").c_str());
  
  TFile file((filebase+".root").c_str(), "recreate");
  glim.GetHistogram()->Write("hXsec_exp_corr");
  cobs.Write("graph_smoothed_Obs");
  cobsup.Write("graph_smoothed_ObsP");
  cobsdown.Write("graph_smoothed_ObsM");
  cexp.Write("graph_smoothed_Exp");
  cup.Write("graph_smoothed_ExpP");
  cdown.Write("graph_smoothed_ExpM");
  hsigobs.Write("hsig_obs_corr");
  hsigexp.Write("hsig_exp_corr");
  file.Close();
  cout << "\nSaved limit curves in " << filebase << ".root\n" << endl;
}

int GetNumBins(const vector<double> &pts, double width){
  double pmin = *min_element(pts.cbegin(), pts.cend());
  double pmax = *max_element(pts.cbegin(), pts.cend());
  return max(1, min(500, static_cast<int>(ceil((pmax-pmin)/width))));
}

void GetParticleNames(string &xparticle, string &yparticle){
  if(model_=="T1tttt"){
    xparticle = "gluino";
    yparticle = "LSP";
  }else if(model_=="T5tttt"){
    xparticle = "gluino";
    yparticle = "LSP";
  }else if(model_=="T2tt"){
    xparticle = "gluino";
    yparticle = "LSP";
  }else if(model_=="T6ttWW"){
    xparticle = "sbottom";
    yparticle = "chargino";
  }else{
    DBG(("Unknown model: "+model_));
    xparticle = "gluino";
    yparticle = "LSP";
  }
}

void Style(TGraph *g, int color, int style, float width){
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
}

TGraph DrawContours(TGraph2D &g2, int color, int style, double width,
		    int n_smooth, double val){
  TGraph graph;

  TList *l;
  //// Finding the TH2D, smoothing it, and creating a TGraph2D to get a new Delauny interpolation
  if(n_smooth>0){
    TH2D *histo2d = g2.GetHistogram();
    TH2D htemp("", "",
               100, histo2d->GetXaxis()->GetXmin(), histo2d->GetXaxis()->GetXmax(),
               100, histo2d->GetYaxis()->GetXmin(), histo2d->GetYaxis()->GetXmax());
    for(int binx=1; binx<=htemp.GetNbinsX(); ++binx){
      double x = htemp.GetXaxis()->GetBinCenter(binx);
      for(int biny=1; biny<=htemp.GetNbinsY(); ++biny){
        double y = htemp.GetYaxis()->GetBinCenter(biny);
        double z = g2.Interpolate(x,y);
        if(z!=0.){
          htemp.SetBinContent(htemp.GetBin(binx, biny), z);
        }
      }
    }
    
    for(int ind=0; ind<n_smooth; ++ind){
      htemp.Smooth(1,"k5b");
    }
    
    vector<double> vx, vy, vz;
    double glu_lsp = 225;
    for(int binx=1; binx<=htemp.GetNbinsX(); ++binx){
      double x = htemp.GetXaxis()->GetBinCenter(binx);
      for(int biny=1; biny<=htemp.GetNbinsY(); ++biny){
	double y = htemp.GetYaxis()->GetBinCenter(biny);
	double z = htemp.GetBinContent(htemp.GetBin(binx,biny));
        
	vx.push_back(x);
	vy.push_back(y);
	if(x-y>glu_lsp+85){
          vz.push_back(z);
	}else{
          vz.push_back(g2.Interpolate(x,y));
        }
      }
    }
    
    TGraph2D gsmooth("gsmooth", "Cross-Section Limit", vx.size(), &vx.at(0), &vy.at(0), &vz.at(0));
    gsmooth.GetHistogram();
    l = gsmooth.GetContourList(val);
  } else {
    g2.GetHistogram();
    l = g2.GetContourList(val);
  }
  if(l == nullptr) return graph;
  int max_points = -1;
  for(int i = 0; i < l->GetSize(); ++i){
    TGraph *g = static_cast<TGraph*>(l->At(i));
    Style(g, color, style, width);
    if(g == nullptr) continue;
    int n_points = g->GetN();
    if(n_points > max_points){
      if(n_smooth>0) FixGraph(*g);
      graph = *g;
      max_points = n_points;
    }
    g->Draw("L same");
  }

  return graph;
}

void FixGraph(TGraph &graph){
  double glu_lsp = 225;
  if(model_=="T5tttt") glu_lsp = 265.;
  int np(graph.GetN());
  double iniglu, endglu, inilsp, endlsp;

  graph.GetPoint(0, iniglu, inilsp);
  graph.GetPoint(np-1, endglu, endlsp);

  // Reversing graph if printed towards decreasing mgluino
  if(inilsp < endlsp) {
    ReverseGraph(graph);
    endglu = iniglu;
    endlsp = inilsp;
  }

  // Adding a point so that it goes down to mLSP = 0, but not for WZ,SOS
  if(endlsp<30){
    graph.SetPoint(graph.GetN(), endglu, 0);
    np++;
  }

  ReverseGraph(graph);
  // Adding a point at mLSP = 0, and removing points beyond the diagonal
  for(int point(0); point < np; point++){
    double mglu, mlsp;
    graph.GetPoint(point, mglu, mlsp);
    if(mlsp > mglu-glu_lsp-5){
      while(point <= graph.GetN() && point!=0) {
	graph.RemovePoint(graph.GetN()-1);
	np--;
      }
      break;
    }
  }
}

void ReverseGraph(TGraph &graph){
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

void SetupSignedColors(){
  const unsigned num = 3;
  const int bands = 255;
  int colors[bands];
  double stops[num] = {0.0, 0.5, 1.0};
  double red[num]   = {1.0, 1.0, 0.0};
  double green[num] = {0.0, 1.0, 0.0};
  double blue[num]  = {0.0, 1.0, 1.0};
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
      {"num_smooth", required_argument, 0, 's'},
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
      num_smooth_ = atoi(optarg);
      break;
    case 'm':
      model_ = optarg;
      break;
    case 'f':
      filename_ = optarg;
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
