#ifndef H_BLOCK
#define H_BLOCK

#include <string>
#include <vector>
#include <initializer_list>

#include "bin.hpp"

struct Block{
  Block(const std::string &name, const std::vector<std::vector<Bin> > &bins);
  Block(const std::string &name, std::initializer_list<std::vector<Bin> > bins);

  std::vector<std::vector<Bin> > bins_;
  std::string name_;
};

#endif
