#include <iostream>

#include "block_yields.hpp"

#include "utilities.hpp"

using namespace std;

BlockYields::BlockYields(const Block &block,
			 const set<Process> &processes,
			 const Cut &cut,
			 const YieldManager &yields):
  gps_(block.Bins().size(), 
       block.Bins().size()
       ? vector<GammaParams>(block.Bins().at(0).size())
       : vector<GammaParams>(0)){
  size_t irow = 0;
  for(const auto &vbin: block.Bins()){
    size_t icol = 0;
    for(const auto &bin: vbin){
      GammaParams &gps = gps_.at(irow).at(icol);
      gps = GammaParams(0., 0.);
      for(const auto &process: processes){
	YieldKey key(bin, process, cut);
	gps += yields.GetYield(key);
      }
      ++icol;
    }
    ++irow;
  }
}

vector<GammaParams> BlockYields::RowSums() const{
  vector<GammaParams> sums(gps_.size(), GammaParams(0., 0.));
  for(size_t irow = 0; irow < gps_.size(); ++irow){
    for(size_t icol = 0; icol < gps_.at(irow).size(); ++icol){
      sums.at(irow) += gps_.at(irow).at(icol);
    }
  }
  return sums;
}

vector<GammaParams> BlockYields::ColSums() const{
  vector<GammaParams> sums(gps_.size() ? gps_.at(0).size() : 0, GammaParams(0., 0.));
  for(size_t irow = 0; irow < gps_.size(); ++irow){
    for(size_t icol = 0; icol < gps_.at(irow).size(); ++icol){
      sums.at(icol) += gps_.at(irow).at(icol);
    }
  }
  return sums;
}

size_t BlockYields::MaxRow() const{
  vector<GammaParams> gps = RowSums();
  if(gps.size() == 0) return -1;
  size_t imax = 0;
  double yield_max = gps.at(0).Yield();
  for(size_t i = 1; i < gps.size(); ++i){
    if(gps.at(i).Yield() > yield_max){
      yield_max = gps.at(i).Yield();
      imax = i;
    }
  }
  return imax;
}

size_t BlockYields::MaxCol() const{
  vector<GammaParams> gps = ColSums();
  if(gps.size() == 0) return -1;
  size_t imax = 0;
  double yield_max = gps.at(0).Yield();
  for(size_t i = 1; i < gps.size(); ++i){
    if(gps.at(i).Yield() > yield_max){
      yield_max = gps.at(i).Yield();
      imax = i;
    }
  }
  return imax;
}

GammaParams BlockYields::Total() const{
  GammaParams gps(0., 0.);
  for(const auto &vbin: gps_){
    for(const auto &bin: vbin){
      gps += bin;
    }
  }
  return gps;
}
