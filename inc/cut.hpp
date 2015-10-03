#ifndef H_CUT
#define H_CUT

#include <string>

class Cut{
public:
  Cut(const std::string &cut = "1");
  Cut(const char *cut);

  Cut & Replace(const Cut &orig, const Cut &rep);
  Cut & RmCutOn(const Cut &to_rm, const Cut &rep = Cut());

  Cut & operator &= (const Cut &cut);
  Cut & operator |= (const Cut &cut);
  Cut & operator += (const Cut &cut);
  Cut & operator -= (const Cut &cut);
  Cut & operator *= (const Cut &cut);
  Cut & operator /= (const Cut &cut);
  Cut & operator %= (const Cut &cut);
  Cut & operator ^= (const Cut &cut);
  Cut & operator <<= (const Cut &cut);
  Cut & operator >>= (const Cut &cut);

  operator std::string() const;
  operator const char *() const;

private:
  std::string cut_;

  void Clean();
};

Cut operator&(Cut a, Cut b);
Cut operator&&(Cut a, Cut b);
Cut operator|(Cut a, Cut b);
Cut operator||(Cut a, Cut b);
Cut operator+(Cut a, Cut b);
Cut operator-(Cut a, Cut b);
Cut operator*(Cut a, Cut b);
Cut operator/(Cut a, Cut b);
Cut operator%(Cut a, Cut b);
Cut operator^(Cut a, Cut b);
Cut operator<<(Cut a, Cut b);
Cut operator>>(Cut a, Cut b);

#endif
