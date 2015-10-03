#ifndef H_SYSTEMATIC
#define H_SYSTEMATIC

#include <string>
#include <tuple>
#include <type_traits>
#include <functional>
#include <utility>

class Systematic{
public:
  Systematic(const std::string &name,
	     double strength);

  const std::string & Name() const;
  Systematic & Name(const std::string &name);

  double Strength() const;
  Systematic & Strength(double strength);

  bool operator<(const Systematic &systematic) const;
  bool operator==(const Systematic &systematic) const;

private:
  std::string name_;
  double strength_;

  auto ComparisonTuple() const{
    return make_tuple(name_);
  }
};

#endif
