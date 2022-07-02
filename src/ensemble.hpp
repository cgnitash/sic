
#pragma once

#include <chrono>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sic
{
struct EnsembleError
{
};

struct Sequence
{
  std::string sequence;
  std::string label;
  double      weight;
};

struct Summary
{
  int                                 D;
  int                                 L;
  double                              total_weight = 0;
  std::map<int, int>                  length_counts;
  std::map<char, std::pair<int, int>> symbol_counts;
  std::vector<char>                   symbols;

  void print() const;
};

class Ensemble
{
  friend class PWM_1;
  friend class PWM_2;
  friend class PWM_3;
  friend class PWM_4;

private:
  std::vector<Sequence> sequences;
  Summary               summary;

public:
  Ensemble(std::vector<Sequence> const &seqs);
  void print_summary() const;
  bool
      lengthsAligned() const
  {
    return summary.length_counts.size() == 1;
  }
  void
      verify() const
  {
    if (not lengthsAligned())
    {
      std::cout << "Lengths must be aligned in order to generate PWM. Print "
                   "Summary for details.\n";
      throw EnsembleError{};
    }
  }
};

void removeLowerCaseResidues(std::vector<Sequence> &sequences,
                             std::string const     &true_target);

std::string extractSingleA2Msequence(std::istream &is);

std::vector<Sequence> extractA2MSequencesFromFile(std::string file);

std::vector<Sequence> extractSequencesFromFile(std::string file,
                                               std::string delimiter,
                                               std::string sequence,
                                               std::string label,
                                               std::string weight);

void filter(std::vector<Sequence> &sequences, std::string value);

std::vector<Sequence>
    sample(std::vector<Sequence> const &sequences, int fraction, int replicate);

std::vector<std::string> split(std::string const &, char delim);

}   // namespace sic
