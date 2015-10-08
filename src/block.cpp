#include "block.hpp"

#include <string>
#include <vector>
#include <initializer_list>

using namespace std;

Block::Block(const string &name, const vector<vector<Bin> > &bins):
  bins_(bins),
  name_(name){
  }

Block::Block(const string &name, initializer_list<vector<Bin> > bins):
  bins_(bins),
  name_(name){
  }

const string & Block::Name() const{
  return name_;
}

Block & Block::Name(const std::string &name){
  name_ = name;
  return *this;
}

const vector<vector<Bin> > & Block::Bins() const{
  return bins_;
}

vector<vector<Bin> > & Block::Bins(){
  return bins_;
}

bool Block::operator<(const Block &b) const{
  return bins_ < b.bins_;
}

bool Block::operator==(const Block &b) const{
  return bins_ == b.bins_;
}

ostream & operator<<(ostream &stream, const Block &block){
  stream << "Block::" << block.Name();
  return stream;
}
