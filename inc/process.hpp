#ifndef H_PROCESS
#define H_PROCESS

#include <string>
#include <vector>
#include <initializer_list>

#include "TChain.h"

struct Process{
  Process(const std::string &name,
          const std::vector<std::string> &file_names,
          const std::string &cut = "1",
          bool count_zeros = true);
  Process(const std::string &name,
          std::initializer_list<std::string> file_names,
          const std::string &cut = "1",
          bool count_zeros = true);

  Process(Process&& p) = delete;
  Process(const Process &p) = delete;
  Process& operator=(const Process &p) = delete;

  bool operator<(const Process &p) const;

  TChain chain_;
  std::string name_, cut_;
  bool count_zeros_;
};

#endif
