
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <tuple>

#include "ensemble.hpp"

void
    Ensemble::load(std::string const &fn)
{
  ensemble_file_name = fn;
  std::ifstream ifs{ ensemble_file_name };
  if (not ifs.is_open())
  {
    std::cout << "Error: ensemble file " << ensemble_file_name << " not found";
    throw EnsembleError{};
  }
  std::cout << "Loading ensemble file " << ensemble_file_name << " ...\n";
  std::string sequence;
  std::string weight;
  std::string id;
  while (std::getline(ifs, sequence, ',') and std::getline(ifs, weight, ',') and
         std::getline(ifs, id))
  {
    id.pop_back();   // clear control character ^M (not sure why it's there in
                     // the first place
    sequences.push_back({ sequence, std::stod(weight), id });

    total_weight += std::stod(weight);
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
    std::cout << "Error: ensemble file " << ensemble_file_name << " is empty\n";
    throw EnsembleError{};
  }
  std::cout << "Ensemble file " << ensemble_file_name
            << " successfully loaded\n";
}

void
    Ensemble::summary()
{
  std::cout << "------\nSummary report for ensemble file " << ensemble_file_name
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

  if (length_counts.size() != 1)
  {
    std::cout << "Error: PWM can only be generated for aligned ensembles\n";
    throw EnsembleError{};
  }
  std::ofstream ofs{ ensemble_file_name + ".pwm_1" };
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  for (auto const &sequence : sequences)
    for (int i = 0; i < L; ++i)
      pwm_1[{ i, sequence.sequence[i] }] += use_weights ? sequence.weight : 1.0;

  for (auto &[key, val] : pwm_1)
  {
    auto const [i, a] = key;
    val /= use_weights ? total_weight : sequences.size();
    ofs << i << " " << a << " " << val << "\n";
  }
}

void
    Ensemble::generate_pwm_2()
{

  if (length_counts.size() != 1)
  {
    std::cout << "Error: PWM can only be generated for aligned ensembles\n";
    throw EnsembleError{};
  }
  std::ofstream ofs{ ensemble_file_name + ".pwm_2" };
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  pwm_2.clear();   // map had better be empty
  for (auto const &sequence : sequences)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        pwm_2[{ i, j, sequence.sequence[i], sequence.sequence[j] }] +=
            use_weights ? sequence.weight : 1.0;
  }

  for (auto &[key, val] : pwm_2)
  {
    val /= use_weights ? total_weight : sequences.size();
    auto const [i, j, a, b] = key;
    ofs << i << " " << j << " " << a << " " << b << " " << val << "\n";
  }
}

void
    Ensemble::generate_pwm_3()
{

  if (length_counts.size() != 1)
  {
    std::cout << "Error: PWM can only be generated for aligned ensembles\n";
    throw EnsembleError{};
  }
  std::ofstream ofs{ ensemble_file_name + ".pwm_3" };
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  pwm_3.clear();
  for (auto const &sequence : sequences)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          pwm_3[{ i,
                  j,
                  k,
                  sequence.sequence[i],
                  sequence.sequence[i],
                  sequence.sequence[i] }] +=
              use_weights ? sequence.weight : 1.0;
  }

  for (auto &[key, val] : pwm_3)
  {
    val /= use_weights ? total_weight : sequences.size();
    auto const [i, j, k, a, b, c] = key;
    ofs << i << " " << j << " " << k << " " << a << " " << b << " " << c << " "
        << val << "\n";
  }
}

double
    Ensemble::calculate_individual_score_1(Sequence const &sequence)
{

  auto const c     = 0.000001;
  auto const D     = static_cast<int>(symbol_counts.size());
  auto const L     = static_cast<int>(sequences[0].sequence.size());
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
  {
    auto val = pwm_1.find({ i, sequence.sequence[i] });

    score += std::log(D * ((val != pwm_1.end() ? val->second : 0.) + c)) /
             std::log(D);
  }
  return score;
}

void
    Ensemble::calculate_true_scores_1()
{
  std::ofstream ofs{ ensemble_file_name + ".ts_1" };
  for (auto const &sequence : sequences)
  {
    auto const score = calculate_individual_score_1(sequence);
    // true_scores_1.push_back(score);
    ofs << score << "\n";
  }
}

double
    Ensemble::calculate_individual_score_2(Sequence const &sequence)
{

  auto const c     = 0.000001;
  auto const D     = static_cast<int>(symbol_counts.size());
  auto const L     = static_cast<int>(sequences[0].sequence.size());
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
    {
      auto val =
          pwm_2.find({ i, j, sequence.sequence[i], sequence.sequence[j] });

      score += std::log(D * D * ((val != pwm_2.end() ? val->second : 0.) + c)) /
               std::log(D);
    }
  return score;
}

void
    Ensemble::calculate_true_scores_2()
{
  std::ofstream ofs{ ensemble_file_name + ".ts_2" };
  for (auto const &sequence : sequences)
  {
    auto const score = calculate_individual_score_2(sequence);
    // true_scores_2.push_back(score);
    ofs << score << "\n";
  }
}

double
    Ensemble::calculate_individual_score_3(Sequence const &sequence)
{
  auto const c     = 0.000001;
  auto const D     = static_cast<int>(symbol_counts.size());
  auto const L     = static_cast<int>(sequences[0].sequence.size());
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
      for (int k = j + 1; k < L; ++k)
      {
        auto val = pwm_3.find({ i,
                                j,
                                k,
                                sequence.sequence[i],
                                sequence.sequence[j],
                                sequence.sequence[k] });

        score += std::log(D * D * D *
                          ((val != pwm_3.end() ? val->second : 0.) + c)) /
                 std::log(D);
      }
  return score;
}

void
    Ensemble::calculate_true_scores_3()
{
  std::ofstream ofs{ ensemble_file_name + ".ts_3" };
  for (auto const &sequence : sequences)
  {
    auto const score = calculate_individual_score_3(sequence);
    // true_scores_3.push_back(score);
    ofs << score << "\n";
  }
}

void
    Ensemble::load_tests(std::string const &test_file_name)
{
  std::ifstream ifs{ test_file_name };
  if (not ifs.is_open())
  {
    std::cout << "Error: test file " << test_file_name << " not found";
    throw EnsembleError{};
  }
  std::cout << "Loading test file " << test_file_name << " ...\n";
  std::string sequence;
  std::string weight;
  std::string id;

  std::ofstream ofs1{ ensemble_file_name + "." + test_file_name + ".scores_1" };
  std::ofstream ofs2{ ensemble_file_name + "." + test_file_name + ".scores_2" };
  std::ofstream ofs_sparse_2{ ensemble_file_name + "." + test_file_name +
                              ".scores_sparse_2" };
  std::ofstream ofs3{ ensemble_file_name + "." + test_file_name + ".scores_3" };
  while (std::getline(ifs, sequence, ',') and std::getline(ifs, weight, ',') and
         std::getline(ifs, id))
  {
    id.pop_back();
    // std::cout << sequence << ":" << weight << ":" << id << "\n";
    auto sequence_with_data = Sequence({ sequence, std::stod(weight), id });
    if (not check_validity(sequence_with_data))
      continue;

    switch (pwm_order)
    {
      case 3:
        ofs3 << calculate_individual_score_3(sequence_with_data) << "," << id
             << "\n";
        [[fallthrough]];
      case 2:
        ofs2 << calculate_individual_score_2(sequence_with_data) << "," << id
             << "\n";
        [[fallthrough]];
      case 1:
        ofs1 << calculate_individual_score_1(sequence_with_data) << "," << id
             << "\n";
    }
  }
}

void
    Ensemble::generate_pwms_and_true_scores(int n)
{
  pwm_order = n;
  switch (pwm_order)
  {
    case 3:
      generate_pwm_3();
      calculate_true_scores_3();
      [[fallthrough]];
    case 2:
      generate_pwm_2();
      calculate_true_scores_2();
      [[fallthrough]];
    case 1:
      generate_pwm_1();
      calculate_true_scores_1();
  }
}

bool
    Ensemble::check_validity(Sequence const &sequence_data)
{
  if (sequence_data.sequence.length() != sequences[0].sequence.length())
  {
    std::cout << "Warning: sequence with id: " << sequence_data.id
              << " has length " << sequence_data.sequence.length()
              << " but ensemble sequences have length "
              << sequences[0].sequence.length()
              << ". Skipping this sequence ...\n";
    return false;
  }

  if (auto f = std::find_if(sequence_data.sequence.begin(),
                            sequence_data.sequence.end(),
                            [&](auto c) {
                              return std::find(symbols.begin(),
                                               symbols.end(),
                                               c) == symbols.end();
                            });
      f != sequence_data.sequence.end())
  {
    std::cout << "Warning: sequence with id: " << sequence_data.id
              << " contains a symbol '" << *f
              << "' that was never seen in the ensemble. Skipping this "
                 "sequence ...\n";
    return false;
  }
  return true;
}
