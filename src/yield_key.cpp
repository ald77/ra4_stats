#include "yield_key.hpp"

using namespace std;

const Bin & GetBin(const YieldKey &yk){
  return get<0>(yk);
}

const Process & GetProcess(const YieldKey &yk){
  return get<1>(yk);
}

const Cut & GetCut(const YieldKey &yk){
  return get<2>(yk);
}

ostream & operator<<(ostream &stream, const YieldKey &key){
  stream << "YieldKey(" << GetBin(key)
	 << "," << GetProcess(key)
	 << "," << GetCut(key)
	 << ")";
  return stream;
}
