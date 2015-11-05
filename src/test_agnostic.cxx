#include <string>

#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"

using namespace std;

namespace{
  double a = 120;
  double b = 80;
  double c = 40;
  double d = 30;
  double e = 20;
  double f = 10;
}

int main(){
  double N = a+b+c+d+e+f;
  double rx1 = (b+e)/(a+d);
  double rx2 = (c+f)/(a+d);
  double ry1 = (d+e+f)/(a+b+c);
  
  RooWorkspace w("w");
  w.factory(("N["+to_string(N)+",0.001,"+to_string(max(10.*N,20.))+"]").c_str());
  w.factory(("rx1["+to_string(rx1)+",0.001,10.]").c_str());
  w.factory(("rx2["+to_string(rx2)+",0.001,10.]").c_str());
  w.factory(("ry1["+to_string(ry1)+",0.001,10.]").c_str());
  w.factory("s1_1[0.,-10.,10.]");
  w.factory("s2_1[0.,-10.,10.]");
  w.factory("sum::rxnorm(1.,rx1,rx2)");
  w.factory("sum::rynorm(1.,ry1)");
  w.factory("prod::rnorm(rxnorm,rynorm)");
  w.factory("expr::rscale('N/rnorm',N,rnorm)");
  w.factory("prod::bkg_a(rscale)");
  w.factory("prod::bkg_b(rscale,rx1)");
  w.factory("prod::bkg_c(rscale,rx2)");
  w.factory("prod::bkg_d(rscale,ry1)");
  w.factory("prod::bkg_e(rscale,rx1,ry1)");
  w.factory("prod::bkg_f(rscale,rx2,ry1)");
  w.factory("prod::sig_a(bkg_a)");
  w.factory("prod::sig_b(bkg_b)");
  w.factory("prod::sig_c(bkg_c)");
  w.factory("prod::sig_d(bkg_d)");
  w.factory("expr::ms1_1('exp(s1_1)',s1_1)");
  w.factory("expr::ms2_1('exp(s2_1)',s2_1)");
  w.factory("prod::sig_e(bkg_e,ms1_1)");
  w.factory("prod::sig_f(bkg_f,ms2_1)");
  w.factory(("obsa["+to_string(a)+"]").c_str());
  w.factory(("obsb["+to_string(b)+"]").c_str());
  w.factory(("obsc["+to_string(c)+"]").c_str());
  w.factory(("obsd["+to_string(d)+"]").c_str());
  w.factory(("obse["+to_string(e)+"]").c_str());
  w.factory(("obsf["+to_string(f)+"]").c_str());
  w.factory("RooPoisson::nulla(obsa,bkg_a)");
  (static_cast<RooPoisson*>(w.pdf("nulla")))->setNoRounding();
  w.factory("RooPoisson::nullb(obsb,bkg_b)");
  (static_cast<RooPoisson*>(w.pdf("nullb")))->setNoRounding();
  w.factory("RooPoisson::nullc(obsc,bkg_c)");
  (static_cast<RooPoisson*>(w.pdf("nullc")))->setNoRounding();
  w.factory("RooPoisson::nulld(obsd,bkg_d)");
  (static_cast<RooPoisson*>(w.pdf("nulld")))->setNoRounding();
  w.factory("RooPoisson::nulle(obse,bkg_e)");
  (static_cast<RooPoisson*>(w.pdf("nulle")))->setNoRounding();
  w.factory("RooPoisson::nullf(obsf,bkg_f)");
  (static_cast<RooPoisson*>(w.pdf("nullf")))->setNoRounding();
  w.factory("RooPoisson::alta(obsa,sig_a)");
  (static_cast<RooPoisson*>(w.pdf("alta")))->setNoRounding();
  w.factory("RooPoisson::altb(obsb,sig_b)");
  (static_cast<RooPoisson*>(w.pdf("altb")))->setNoRounding();
  w.factory("RooPoisson::altc(obsc,sig_c)");
  (static_cast<RooPoisson*>(w.pdf("altc")))->setNoRounding();
  w.factory("RooPoisson::altd(obsd,sig_d)");
  (static_cast<RooPoisson*>(w.pdf("altd")))->setNoRounding();
  w.factory("RooPoisson::alte(obse,sig_e)");
  (static_cast<RooPoisson*>(w.pdf("alte")))->setNoRounding();
  w.factory("RooPoisson::altf(obsf,sig_f)");
  (static_cast<RooPoisson*>(w.pdf("altf")))->setNoRounding();
  w.factory("RooGaussian::pdfx(x[0.,-10.,10.],0.,1.)");
  w.factory("PROD::model_b(nulla,nullb,nullc,nulld,nulle,nullf,pdfx)");
  w.factory("PROD::model_s(alta,altb,altc,altd,alte,altf,pdfx)");
  w.defineSet("POI","s1_1,s2_1");
  w.defineSet("nuisances","N,rx1,rx2,ry1,x");
  w.defineSet("observables","obsa,obsb,obsc,obsd,obse,obsf");
  w.defineSet("globalObservables", "");
  RooDataSet data("data_obs", "data_obs", *w.set("observables"));
  data.add(*w.set("observables"));
  w.import(data);

  RooStats::ModelConfig model_config("ModelConfig", &w);
  model_config.SetPdf(*w.pdf("model_s"));
  model_config.SetParametersOfInterest(*w.set("POI"));
  model_config.SetObservables(*w.set("observables"));
  model_config.SetNuisanceParameters(*w.set("nuisances"));
  model_config.SetGlobalObservables(*w.set("globalObservables"));

  RooStats::ModelConfig model_config_bonly("ModelConfig_bonly", &w);
  model_config_bonly.SetPdf(*w.pdf("model_b"));
  model_config_bonly.SetParametersOfInterest(*w.set("POI"));
  model_config_bonly.SetObservables(*w.set("observables"));
  model_config_bonly.SetNuisanceParameters(*w.set("nuisances"));
  model_config_bonly.SetGlobalObservables(*w.set("globalObservables"));

  w.import(model_config);
  w.import(model_config_bonly);

  w.writeToFile("agnostic.root");
}
