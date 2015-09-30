#ifndef H_BIN_PROC
#define H_BIN_PROC

#include "bin.hpp"
#include "process.hpp"

struct BinProc{
  BinProc(const Bin &bin, Process &process);

  bool operator<(const BinProc &bp) const;

  Process &process_;
  Bin bin_;
};

#endif
