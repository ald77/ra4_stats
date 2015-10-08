#ifndef H_YIELD_KEY
#define H_YIELD_KEY

#include <tuple>
#include <ostream>

#include "bin.hpp"
#include "process.hpp"
#include "cut.hpp"

using YieldKey = std::tuple<Bin, Process, Cut>;

const Bin & GetBin(const YieldKey &yk);
const Process & GetProcess(const YieldKey &yk);
const Cut & GetCut(const YieldKey &yk);

std::ostream & operator<<(std::ostream &stream, const YieldKey &key);

#endif
