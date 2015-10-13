#ifndef H_YIELD_MANAGER
#define H_YIELD_MANAGER

#include <map>

#include "yield_key.hpp"
#include "gamma_params.hpp"
#include "bin.hpp"
#include "process.hpp"
#include "cut.hpp"

class YieldManager{
public:
  explicit YieldManager(double lumi = 4.);

  GammaParams GetYield(const YieldKey &key) const;
  GammaParams GetYield(const Bin &bin,
		       const Process &process,
		       const Cut &cut) const;

  const double & Luminosity() const;
  double & Luminosity();

private:
  static std::map<YieldKey, GammaParams> yields_;
  static const double store_lumi_;
  double local_lumi_;
  bool verbose_;

  bool HaveYield(const YieldKey &key) const;
  void ComputeYield(const YieldKey &key) const;
};

#endif
