
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>

#include "ensemble.hpp"

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

void
    Ensemble::load_ensemble(std::string sf,
                            std::string tf,
                            std::string tv,
                            std::string wf,
                            std::string fn)
{
  sequence_field = sf;
  train_field    = tf;
  train_value    = tv;
  weight_field   = wf;
  file_name      = fn;
  std::ifstream ifs{ file_name };
  if (not ifs.is_open())
  {
    std::cout << "Error: ensemble file " << file_name << " not found";
    throw EnsembleError{};
  }
  std::cout << "Loading ensemble file " << file_name << " ...\n";

  std::string line;
  std::getline(ifs, line);
  std::vector<std::string> columns      = split(line);
  auto                     column_index = [&](std::string field) -> int
  {
    if (field == " ")
      return -1;

    auto f = std::find(std::begin(columns), std::end(columns), field);
    if (f == std::end(columns))
    {
      std::cout << "Error: " << field << " is not a column in " << file_name
                << "\n";
      throw EnsembleError{};
    }
    return std::distance(std::begin(columns), f);
  };

  train_field_index    = column_index(train_field);
  weight_field_index   = column_index(weight_field);
  sequence_field_index = column_index(sequence_field);
  if (sequence_field_index == -1)
  {
    std::cout << "Error: " << sequence_field << " is not a column in "
              << file_name << ". A valid column must be specified\n";
    throw EnsembleError{};
  }

  if (train_value == " ")
  {
    if (train_field_index != -1)
    {
      std::cout
          << "Error: Training column specified without value. If all sequences "
             "should be used to train, the column must not be specified.\n";
      throw EnsembleError{};
    }
    std::cout << "No training column specified. All sequences will be used.\n";
  }
  else if (train_field_index == -1)
  {
    std::cout
        << "Error: Training value specified without column. If all sequences "
           "should be used to train, the value must not be specified.\n";
    throw EnsembleError{};
  }

  while (std::getline(ifs, line))
  {
    auto const row = split(line);

    if (train_field_index != -1 and row[train_field_index] != train_value)
      continue;

    auto sequence = row[sequence_field_index];
    sequences.push_back(
        { sequence,
          weight_field_index == -1 ? 1 : std::stod(row[weight_field_index]) });

    total_weight += sequences.back().weight;

    std::set<char> uniq_symbols;
    for (auto const &c : sequence)
    {
      uniq_symbols.insert(c);
      symbol_counts[c].first++;
    }
    for (auto const &c : uniq_symbols)
      symbol_counts[c].second++;

    length_counts[sequence.size()]++;
  }
  for (auto const &symbol : symbol_counts)
    symbols.push_back(symbol.first);

  if (sequences.empty())
  {
    std::cout << "Error: ensemble file " << file_name << " is empty\n";
    throw EnsembleError{};
  }
  D = static_cast<int>(symbol_counts.size());
  L = static_cast<int>(sequences[0].sequence.size());
  std::cout << "Ensemble file " << file_name << " successfully loaded\n";
}

void
    Ensemble::summary() const
{
  std::cout << "------\nSummary report for ensemble file " << file_name
            << "\n------\n";
  std::cout << "Ensemble contains " << sequences.size() << " sequences\n";
  if (length_counts.size() == 1)
    std::cout << "All sequences are of length " << sequences[0].sequence.size()
              << "\n";
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

void
    Ensemble::generate_pwm_1()
{
  std::cout << "Generating PWM of order 1 ...\n" << std::flush;

  /*
    std::ofstream ofs{ file_name + "." + train_field + "." + train_value +
                     ".pwm_1" };

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  */

  pwm_1.clear();   // map had better be empty
  for (auto const &sequence : sequences)
    for (int i = 0; i < L; ++i)
      pwm_1[{ i, sequence.sequence[i] }] += sequence.weight;

  for (auto &[key, val] : pwm_1)
  {
    val /= total_weight;
    // auto const [i, a] = key;
    // ofs << i << " " << a << " " << val << "\n";
  }
  std::cout << " PWM of order 1 generated.\n" << std::flush;
}

void
    Ensemble::generate_pwm_2()
{
  std::cout << "Generating PWM of order 2 ...\n" << std::flush;

  /*
   std::ofstream ofs{ file_name + "." + train_field + "." + train_value +
                     ".pwm_2" };

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  */

  pwm_2.clear();   // map had better be empty
  for (auto const &sequence : sequences)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        pwm_2[{ i, j, sequence.sequence[i], sequence.sequence[j] }] +=
            sequence.weight;
  }

  for (auto &[key, val] : pwm_2)
  {
    val /= total_weight;
    // auto const [i, j, a, b] = key;
    // ofs << i << " " << j << " " << a << " " << b << " " << val << "\n";
  }
  std::cout << " PWM of order 2 generated.\n" << std::flush;
}

void
    Ensemble::generate_pwm_3()
{
  std::cout << "Generating PWM of order 3 ...\n" << std::flush;

  /*
  std::ofstream ofs{ file_name + "." + train_field + "." + train_value +
                     ".pwm_3" };

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  */

  pwm_3.clear();   // map has to be empty
  for (auto const &sequence : sequences)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          pwm_3[{ i,
                  j,
                  k,
                  sequence.sequence[i],
                  sequence.sequence[j],
                  sequence.sequence[k] }] += sequence.weight;
  }

  for (auto &[key, val] : pwm_3)
  {
    val /= total_weight;
    // auto const [i, j, k, a, b, c] = key;
    // ofs << i << " " << j << " " << k << " " << a << " " << b << " " << c << "
    // "
    //     << val << "\n";
  }
  std::cout << " PWM of order 3 generated.\n" << std::flush;
}

double
    Ensemble::calculate_individual_score_1(std::string const &sequence) const
{

  auto score = 0.;
  for (int i = 0; i < L; ++i)
  {
    auto val = pwm_1.find({ i, sequence[i] });

    score += std::log(D * ((val != pwm_1.end() ? val->second : 0.) + c)) /
             std::log(D);
  }
  return score;
}

double
    Ensemble::calculate_individual_score_2(std::string const &sequence) const
{

  auto score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
    {
      auto val = pwm_2.find({ i, j, sequence[i], sequence[j] });

      score += std::log(D * D * ((val != pwm_2.end() ? val->second : 0.) + c)) /
               std::log(D);
    }
  return score;
}

double
    Ensemble::calculate_individual_score_3(std::string const &sequence) const
{
  auto score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
      for (int k = j + 1; k < L; ++k)
      {
        auto val =
            pwm_3.find({ i, j, k, sequence[i], sequence[j], sequence[k] });

        score += std::log(D * D * D *
                          ((val != pwm_3.end() ? val->second : 0.) + c)) /
                 std::log(D);
      }
  return score;
}

void
    Ensemble::run_tests() const
{
  std::ifstream ifs{ file_name };
  if (not ifs.is_open())
  {
    std::cout << "Error: test file " << file_name << " not found";
    throw EnsembleError{};
  }
  std::cout << "Loading test file " << file_name << " ...\n" << std::flush;

  std::ofstream ofs{ file_name + "." + train_field + "." + train_value +
                     ".scores" };

  std::string line;
  std::getline(ifs, line);
  split(line);   // read header and ignore since this was used in training
  ofs << line;
  switch (pwm_order)
  {
    case 3:
      ofs << ",score_3";
      [[fallthrough]];
    case 2:
      ofs << ",score_2";
      [[fallthrough]];
    case 1:
      ofs << ",score_1";
  }

  std::vector<std::string> row;

  ofs << "\n";
  while (std::getline(ifs, line))
  {
    auto const row      = split(line);
    auto       sequence = row[sequence_field_index];
    if (not check_validity(sequence))
      continue;

    ofs << line;
    switch (pwm_order)
    {
      case 3:
        ofs << "," << calculate_individual_score_3(sequence);
        [[fallthrough]];
      case 2:
        ofs << "," << calculate_individual_score_2(sequence);
        [[fallthrough]];
      case 1:
        ofs << "," << calculate_individual_score_1(sequence);
    }
    ofs << "\n";
  }
  std::cout << "All test sequences are scored.\n" << std::flush;
}

void
    Ensemble::generate_pwms(int n)
{
  if (length_counts.size() != 1)
  {
    std::cout << "Error: PWM can only be generated for aligned ensembles\n";
    throw EnsembleError{};
  }
  pwm_order = n;
  switch (pwm_order)
  {
    case 3:
      timer([this] { generate_pwm_3(); });
      [[fallthrough]];
    case 2:
      timer([this] { generate_pwm_2(); });
      [[fallthrough]];
    case 1:
      timer([this] { generate_pwm_1(); });
  }
}

bool
    Ensemble::check_validity(std::string const &sequence) const
{
  if (sequence.length() != sequences[0].sequence.length())
  {
    std::cout << "Warning: sequence " << sequence << " has length "
              << sequence.length() << " but ensemble sequences have length "
              << sequences[0].sequence.length()
              << ". Skipping this sequence ...\n";
    return false;
  }

  if (auto f = std::find_if(sequence.begin(),
                            sequence.end(),
                            [&](auto c) {
                              return std::find(symbols.begin(),
                                               symbols.end(),
                                               c) == symbols.end();
                            });
      f != sequence.end())
  {
    std::cout << "Warning: sequence " << sequence << " contains a symbol '"
              << *f
              << "' that was never seen in the ensemble. Skipping this "
                 "sequence ...\n";
    return false;
  }
  return true;
}

