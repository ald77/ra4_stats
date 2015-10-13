#ifndef H_YIELD_MANAGER
#define H_YIELD_MANAGER

#include <map>

#include "yield_key.hpp"
#include "gamma_params.hpp"

class YieldManager{
public:
  explicit YieldManager(double lumi = 4.);

  GammaParams GetYield(const YieldKey &key) const;

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
