
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>

#include "ensemble.hpp"

namespace sic
{
std::vector<std::string>
    split(std::string const &line, char delim)
{
  std::istringstream       iss{ line };
  std::vector<std::string> v;
  std::string              s;
  while (std::getline(iss, s, delim))
    v.push_back(s);
  std::getline(iss, s);
  v.push_back(s);
  return v;
}

std::vector<Sequence>
    extractSequencesFromFile(std::string file,
                             std::string delimiter,
                             std::string sequence,
                             std::string label,
                             std::string weight)
{

  std::ifstream ifs{ file };
  if (not ifs.is_open())
  {
    std::cout << "Error: file " << file << " not found";
    throw EnsembleError{};
  }
  std::cout << "Loading file " << file << " ...\n";

  std::string line;
  char const  delim = delimiter[0];   // multi-char delimiters not supported
  std::getline(ifs, line);
  std::vector<std::string> columns = split(line, delim);

  auto column_index = [&](std::string field) -> int
  {
    if (field == "__")
      return -1;

    auto f = std::find(std::begin(columns), std::end(columns), field);
    if (f == std::end(columns))
    {
      std::cout << "Error: " << field << " is not a column in " << file << "\n";
      throw EnsembleError{};
    }
    return std::distance(std::begin(columns), f);
  };

  auto const sequence_index = column_index(sequence);
  auto const label_index    = column_index(label);
  auto const weight_index   = column_index(weight);

  if (sequence_index == -1)
  {
    std::cout << "Error: " << sequence << " is not a column in " << file
              << " A valid column must be specified for sequences\n";
    throw EnsembleError{};
  }

  std::vector<Sequence> all_sequences;
  while (std::getline(ifs, line))
  {
    auto const row = split(line, delim);

    all_sequences.push_back(
        { row[sequence_index],
          label == "__" ? label : row[label_index],
          weight_index == -1 ? 1 : std::stod(row[weight_index]) });
  }

  std::cout << "File " << file << " succesfully loaded.\n";

  return all_sequences;
}

void
    filter(std::vector<Sequence> &sequences, std::string value)
{
  sequences.erase(std::remove_if(std::begin(sequences),
                                 std::end(sequences),
                                 [&](auto const &sequence)
                                 { return sequence.label != value; }),
                  std::end(sequences));
}

std::vector<Sequence>
    sample(std::vector<Sequence> const &all_sequences,
           int                          fraction,
           int                          replicate)
{

  std::vector<Sequence> sequences;
  std::mt19937          gen;
  gen.seed(replicate);
  std::sample(std::begin(all_sequences),
              std::end(all_sequences),
              std::back_inserter(sequences),
              static_cast<int>(fraction / 100.0 * all_sequences.size()) +
                  1,   // guarantee at least one sequence
              gen);
  return sequences;
}

void
    removeLowerCaseResidues(std::vector<Sequence> &sequences,
                            std::string const     &true_target)
{
  for (auto &[sequence, label, weight] : sequences)
  {
    std::string fixed_sequence;
    for (int i = 0; i < static_cast<int>(sequence.length()); ++i)
      if (not std::islower(true_target[i]))
        fixed_sequence.push_back(sequence[i]);
    sequence = fixed_sequence;
  }
}

Ensemble::Ensemble(std::vector<Sequence> const &seqs)
{
  sequences = seqs;

  if (sequences.empty())
  {
    std::cout << "Error: no sequences provided\n";
    throw EnsembleError{};
  }

  for (auto const &[sequence, label, weight] : sequences)
  {

    summary.total_weight += weight;

    std::set<char> uniq_symbols;
    for (auto const &c : sequence)
    {
      uniq_symbols.insert(c);
      summary.symbol_counts[c].first++;
    }
    for (auto const &c : uniq_symbols)
      summary.symbol_counts[c].second++;

    summary.length_counts[sequence.size()]++;
  }
  for (auto const &symbol : summary.symbol_counts)
    summary.symbols.push_back(symbol.first);

  summary.D = static_cast<int>(summary.symbol_counts.size());
  summary.L = static_cast<int>(sequences[0].sequence.size());
}

void
    Ensemble::print_summary() const
{
  std::cout << "------\nSummary report for ensemble"
            << "\n------\n";
  std::cout << "Ensemble contains " << sequences.size() << " sequences\n";
  summary.print();
}

void
    Summary::print() const
{
  if (length_counts.size() == 1)
    std::cout << "All sequences are of length "
              << std::begin(length_counts)->first << "\n";
  else
  {
    for (auto const &[len, cnt] : length_counts)
      std::cout << std::setw(5) << len << " sequences found of length " << cnt
                << "\n";
  }
  std::cout << "Total Weight of sequences in ensemble: " << total_weight
            << "\n";
  std::cout << "Ensemble contains " << symbol_counts.size()
            << " different characters\n";
  std::cout << "Symbol:              ";
  auto const padding = 8;
  for (auto const &mon : symbol_counts)
    std::cout << std::setw(padding) << mon.first;
  std::cout << "\nFrequency (Global):  ";
  for (auto const &mon : symbol_counts)
    std::cout << std::setw(padding) << mon.second.first;
  std::cout << "\nFrequency (Per-site):";
  for (auto const &mon : symbol_counts)
    std::cout << std::setw(padding) << mon.second.second;
  std::cout << std::endl;
}

std::string
    extractSingleA2Msequence(std::istream &is)
{
  std::string res, line;
  while (std::getline(is, line) and line[0] != '>')
  {
    res += line;
  }
  return res;
}

std::pair<std::vector<Sequence>, int>
    extractA2MSequencesFromFile(std::string file)
{
  std::ifstream ifs{ file };
  if (not ifs.is_open())
  {
    std::cout << "Error: file " << file << " not found";
    throw EnsembleError{};
  }

  std::string line;
  std::getline(ifs, line);
  std::regex  r{ R"(.*/(\d+)-\d+)" };
  std::smatch m;
  if (not std::regex_match(line, m, r))
  {
    std::cout << "Error: Cannot compute true offset. a2m file header is: "
              << line << "\n";
    throw EnsembleError{};
  }
  auto true_offset = std::stoi(m[1].str());

  std::vector<Sequence> result;
  auto                  target = extractSingleA2Msequence(ifs);
  result.push_back({ target, "__", 1.0 });

  std::string seq;
  while (not(seq = extractSingleA2Msequence(ifs)).empty())
    result.push_back({ seq, "__", 1.0 });

  return { result, true_offset };
}
}   // namespace sic
