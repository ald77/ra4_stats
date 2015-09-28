#include <cmath>

#include <iostream>
#include <limits>

#include "TMath.h"

using namespace std;

int main(int argc, char *argv[]){
  double inf = numeric_limits<double>::infinity();
  double sum= 0.;
  for(int iarg = 1; iarg < argc; ++iarg){
    sum += -2.*log(erfc(atof(argv[iarg])/sqrt(2.)));
  }
  double p = TMath::Prob(sum, 2*(argc-1));
  double cp = 1.-0.5*p;
  double z = cp <= 0. ? -inf
    : cp >= 1. ? inf
    : TMath::NormQuantile(cp);
  cout << z << endl;
}
