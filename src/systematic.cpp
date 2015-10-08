#include "systematic.hpp"

#include <string>

using namespace std;

Systematic::Systematic(const string &name,
		       double strength):
  name_(name),
  strength_(strength){
}

const std::string & Systematic::Name() const{
  return name_;
}

Systematic & Systematic::Name(const std::string &name){
  name_ = name;
  return *this;
}

const double & Systematic::Strength() const{
  return strength_;
}

double & Systematic::Strength(){
  return strength_;
}

bool Systematic::operator<(const Systematic &systematic) const{
  return name_ < systematic.name_;
}

bool Systematic::operator==(const Systematic &systematic) const{
  return name_ == systematic.name_;
}

ostream & operator<<(ostream &stream, const Systematic &syst){
  stream << "Systematic::" << syst.Name()
	 << "(strength=" << syst.Strength() << ")";
  return stream;
}
