#ifndef H_BLOCK_YIELDS
#define H_BLOCK_YIELDS

#include <set>
#include <vector>
#include <map>

#include "block.hpp"
#include "process.hpp"
#include "cut.hpp"
#include "gamma_params.hpp"
#include "yield_key.hpp"
#include "yield_manager.hpp"

class BlockYields{
public:
  BlockYields(const Block &block,
	      const std::set<Process> &processes,
	      const Cut &cut,
	      const YieldManager &yields);

  std::vector<GammaParams> RowSums() const;
  std::vector<GammaParams> ColSums() const;

  size_t MaxRow() const;
  size_t MaxCol() const;

  GammaParams Total() const;

private:
  std::vector<std::vector<GammaParams> > gps_;
};

#endif
