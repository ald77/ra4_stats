#ifndef H_SYSTEMATIC
#define H_SYSTEMATIC

#include <string>
#include <tuple>
#include <ostream>

class Systematic{
public:
  Systematic(const std::string &name,
	     double strength);

  const std::string & Name() const;
  Systematic & Name(const std::string &name);

  const double & Strength() const;
  double & Strength();

  bool operator<(const Systematic &systematic) const;
  bool operator==(const Systematic &systematic) const;

private:
  std::string name_;
  double strength_;
};

std::ostream & operator<<(std::ostream &stream, const Systematic &syst);

#endif
