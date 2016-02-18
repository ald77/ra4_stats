#include "bin.hpp"

#include <algorithm>

#include "systematic.hpp"
#include "utilities.hpp"

using namespace std;

Bin::Bin(const string &name, const class Cut &cut,
	 bool is_blind, const SystCollection &systematics):
  cut_(cut),
  name_(name),
  systematics_(systematics),
  is_blind_(is_blind){
  ReplaceAll(name_, " ", "");
  }

const string & Bin::Name() const{
  return name_;
}

Bin & Bin::Name(const string &name){
  name_ = name;
  return *this;
}

const class Cut & Bin::Cut() const{
  return cut_;
}

class Cut & Bin::Cut(){
  return cut_;
}

bool Bin::Blind() const{
  return is_blind_;
}

bool & Bin::Blind(){
  return is_blind_;
}

const Bin::SystCollection & Bin::Systematics() const{
  return systematics_;
}

Bin & Bin::Systematics(const SystCollection &systematics){
  systematics_ = systematics;
  return *this;
}

Bin & Bin::AddSystematic(const Systematic &systematic){
  if(!HasSystematic(systematic)){
    systematics_.insert(systematics_.end(), systematic);
  }
  return *this;
}

Bin & Bin::AddSystematics(const SystCollection &systematics){
  for(const auto& systematic: systematics){
    AddSystematic(systematic);
  }
  return *this;
}

bool Bin::HasSystematic(const Systematic &systematic) const{
  return find(systematics_.cbegin(), systematics_.cend(), systematic) != systematics_.cend();
}

Bin & Bin::RemoveSystematic(const Systematic &systematic){
  try{
    systematics_.erase(find(systematics_.begin(), systematics_.end(), systematic));
  }catch(const out_of_range &e){
    ERROR(string(e.what())+": bin "+name_+" does not contain systematic "+systematic.Name()+".");
  }
  return *this;
}

Bin & Bin::RemoveSystematics(){
  systematics_.clear();
  return *this;
}

Bin & Bin::SetSystematicStrength(const std::string &name, double strength){
  bool found_it = false;
  for(auto systematic = systematics_.cbegin(); systematic != systematics_.cend(); ++systematic){
    if(systematic->Name() == name){
      Systematic new_syst = *systematic;
      new_syst.Strength() = strength;
      found_it = true;
      systematics_.erase(systematic);
      systematics_.insert(systematics_.end(), new_syst);
    }
  }
  if(!found_it){
    ERROR("Bin "+name_+" does not contain systematic "+name);
  }
  return *this;
}

bool Bin::operator<(const Bin &b) const{
  return tie(cut_, systematics_, is_blind_) < tie(b.cut_, b.systematics_, b.is_blind_);
}

bool Bin::operator==(const Bin &b) const{
  return tie(cut_, systematics_, is_blind_) == tie(b.cut_, b.systematics_, b.is_blind_);
}

ostream & operator<<(ostream &stream, const Bin &bin){
  stream << "Bin::" << bin.Name()
	 << "(cut=" << bin.Cut()
	 << ",blind=" << bin.Blind()
	 << ")";
  return stream;
}
