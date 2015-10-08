#ifndef H_BLOCK
#define H_BLOCK

#include <string>
#include <vector>
#include <initializer_list>
#include <tuple>
#include <ostream>

#include "bin.hpp"

class Block{
public:
  Block(const std::string &name, const std::vector<std::vector<Bin> > &bins);
  Block(const std::string &name, std::initializer_list<std::vector<Bin> > bins);

  const std::string & Name() const;
  Block & Name(const std::string &name);

  const std::vector<std::vector<Bin> > & Bins() const;
  std::vector<std::vector<Bin> > & Bins();

  bool operator<(const Block &b) const;
  bool operator==(const Block &b) const;

private:
  std::vector<std::vector<Bin> > bins_;
  std::string name_;
};

std::ostream & operator<<(std::ostream &stream, const Block &block);

#endif
