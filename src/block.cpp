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
