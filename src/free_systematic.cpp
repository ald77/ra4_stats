#include "free_systematic.hpp"

using namespace std;

FreeSystematic::FreeSystematic(const string &name):
  name_(name),
  strengths_(){
}

const string & FreeSystematic::Name() const{
  return name_;
}

string & FreeSystematic::Name(){
  return name_;
}

bool FreeSystematic::HasEntry(const Bin &bin, const Process &process) const{
  pair<Bin, Process> bp(bin, process);
  return strengths_.find(bp) != strengths_.end()
    && strengths_.at(bp) != 0.;
}

double FreeSystematic::Strength(const Bin &bin, const Process &process) const{
  if(HasEntry(bin, process)){
    pair<Bin, Process> bp(bin, process);
    return strengths_.at(bp);
  }else{
    return 0.;
  }
}

double & FreeSystematic::Strength(const Bin &bin, const Process &process){
  pair<Bin, Process> bp(bin, process);
  return strengths_[bp];
}

bool FreeSystematic::operator<(const FreeSystematic &syst) const{
  return tie(name_, strengths_) < tie(syst.name_, syst.strengths_);
}

bool FreeSystematic::operator==(const FreeSystematic &syst) const{
  return tie(name_, strengths_) == tie(syst.name_, syst.strengths_);
}

ostream & operator<<(ostream &stream, const FreeSystematic &syst){
  stream << "FreeSystematic::" << syst.Name();
  return stream;
}
