#include "bin_proc.hpp"

#include "bin.hpp"
#include "process.hpp"

BinProc::BinProc(const Bin &bin, const Process &process):
  process_(process),
  bin_(bin){
  }

bool BinProc::operator<(const BinProc &bp) const{
  return bin_ < bp.bin_
    || (!(bp.bin_ < bin_)
        && process_ < bp.process_);
}

