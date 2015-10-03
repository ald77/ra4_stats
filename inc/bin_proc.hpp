#ifndef H_BIN_PROC
#define H_BIN_PROC

#include "bin.hpp"
#include "process.hpp"

struct BinProc{
  BinProc(const Bin &bin, const Process &process);

  bool operator<(const BinProc &bp) const;

  const Process &process_;
  const Bin &bin_;
};

#endif
