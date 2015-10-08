#include "cut.hpp"

#include "utilities.hpp"

using namespace std;

Cut::Cut(const string &cut):
  cut_(cut){
  Clean();
}

Cut::Cut(const char *cut):
  cut_(cut){
  Clean();
}

Cut & Cut::Replace(const Cut &orig, const Cut &rep){
  ReplaceAll(cut_, orig.cut_, rep.cut_);
  return *this;
}

Cut & Cut::RmCutOn(const Cut &to_rm, const Cut &rep){
  ::RmCutOn(cut_, to_rm.cut_, rep.cut_);
  return *this;
}

Cut & Cut::operator &= (const Cut &cut){
  cut_ = "("+cut_+")&&("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator |= (const Cut &cut){
  cut_ = "("+cut_+")||("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator += (const Cut &cut){
  cut_ = "("+cut_+")+("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator -= (const Cut &cut){
  cut_ = "("+cut_+")-("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator *= (const Cut &cut){
  cut_ = "("+cut_+")*("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator /= (const Cut &cut){
  cut_ = "("+cut_+")/("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator %= (const Cut &cut){
  cut_ = "("+cut_+")%("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator ^= (const Cut &cut){
  cut_ = "("+cut_+")^("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator <<= (const Cut &cut){
  cut_ = "("+cut_+")<<("+cut.cut_+")";
  Clean();
  return *this;
}

Cut & Cut::operator >>= (const Cut &cut){
  cut_ = "("+cut_+")>>("+cut.cut_+")";
  Clean();
  return *this;
}

Cut::operator string() const{
  return cut_;
}

Cut::operator const char*() const{
  return cut_.c_str();
}

bool Cut::operator<(const Cut &cut) const{
  return cut_ < cut.cut_;
}

bool Cut::operator==(const Cut &cut) const{
  return cut_ == cut.cut_;
}

void Cut::Clean(){
  ReplaceAll(cut_, " ", "");
}

Cut operator& (Cut a, Cut b){
  return (a&=b);
}

Cut operator&& (Cut a, Cut b){
  return (a&=b);
}

Cut operator| (Cut a, Cut b){
  return (a|=b);
}

Cut operator|| (Cut a, Cut b){
  return (a|=b);
}

Cut operator+ (Cut a, Cut b){
  return (a+=b);
}

Cut operator- (Cut a, Cut b){
  return (a-=b);
}

Cut operator* (Cut a, Cut b){
  return (a*=b);
}

Cut operator/ (Cut a, Cut b){
  return (a/=b);
}

Cut operator% (Cut a, Cut b){
  return (a%=b);
}

Cut operator^ (Cut a, Cut b){
  return (a^=b);
}

Cut operator<< (Cut a, Cut b){
  return (a<<=b);
}

Cut operator>> (Cut a, Cut b){
  return (a>>=b);
}

ostream & operator<<(ostream &stream, const Cut &cut){
  stream << static_cast<string>(cut);
  return stream;
}
