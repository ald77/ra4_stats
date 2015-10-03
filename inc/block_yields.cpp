#include "block_yields.hpp"

#include "utility.hpp"

BlockYields(const Block &block,
	    const set<Process> &processes,
	    const Cut &cut,
	    const map<YieldKey, GammaParams> &yields):
  gps_(block.Bins().size(), 
       block.Bins().size() ? block.Bins().at(0).size() : 0){
  size_t irow = 0, icol = 0;
  for(const auto &vbin: block.Bins()){
    for(const auto &bin: vbin){
      GammaParams &gps = gps_.at(irow).at(icol);
      gps = GammaParams(0., 0.);
      for(const auto &process: processes){
	YieldKey key(bin, process, cut);
	gps += yields.at(key);
      }
      ++icol;
    }
    ++irow;
  }
}

vector<GammaParams> BlockYields::RowSums() const{
  vector<GammaParams> sums(gps_.size(), GammaParams(0., 0.));
  for(size_t irow = 0; irow < gps_.size(); ++irow){
    for(size_t icol = 0; icol < gps_.size(); ++icol){
      sums.at(irow) += gps_.at(irow).at(icol);
    }
  }
  return sums;
}

vector<GammaParams> BlockYields::ColSums() const{
  vector<GammaParams> sums(gps_.size() ? gps_.at(0).size() : 0, GammaParams(0., 0.));
  for(size_t irow = 0; irow < gps_.size(); ++irow){
    for(size_t icol = 0; icol < gps_.size(); ++icol){
      sums.at(icol) += gps_.at(irow).at(icol);
    }
  }
  return sums;
}

size_t BlockYields::MaxRow() const{
  return MaxIndex(RowSums());
}

size_t BlockYields::MaxCol() const{
  return MaxIndex(ColSums());
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
