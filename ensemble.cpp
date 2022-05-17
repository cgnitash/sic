
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "ensemble.hpp"

void
    Ensemble::load(std::string const &fn)
{
  file_name = fn;
  std::ifstream ifs{ file_name };
  if (not ifs.is_open())
  {
    std::cout << "Error: ensemble file " << file_name << " not found";
    std::exit(1);
  }
  std::cout << "Loading ensemble file " << file_name << " ...\n";
  std::string sequence;
  std::string weight;
  std::string id;
  while (std::getline(ifs, sequence, ',') and std::getline(ifs, weight, ',') and
         std::getline(ifs, id))
  {
    // std::cout << sequence << ":" << weight << ":" << id << "\n";
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
    std::cout << "Error: ensemble file " << file_name << " is empty\n";
    std::exit(1);
  }
  std::cout << "Ensemble file " << file_name << " successfully loaded\n";
}

void
    Ensemble::summary()
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

int
    Ensemble::symbol_index(char c) const
{
  return std::distance(std::begin(symbols),
                       std::find(std::begin(symbols), std::end(symbols), c));
};

void
    Ensemble::generate_pwm_1()
{

  if (length_counts.size() != 1)
  {
    std::cout << "Error: PWM can only be generated for aligned ensembles\n";
    std::exit(1);
  }
  std::ofstream ofs{ file_name + ".pwm_1" };
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  pwm_1 = std::vector<std::vector<double>>(L, std::vector<double>(D, 0));
  for (auto const &sequence : sequences)
    for (int i = 0; i < L; ++i)
      pwm_1[i][symbol_index(sequence.sequence[i])]++;

  for (auto &row : pwm_1)
  {
    for (auto &val : row)
    {
      val /= sequences.size();
      ofs << val << " ";
    }
    ofs << "\n";
  }
}

int
    Ensemble::C(int n, int r)
{
  int c = 1;
  for (int i = n - r + 1; i <= n; ++i)
    c *= i;
  for (int i = 2; i <= r; ++i)
    c /= i;
  return c;
}

void
    Ensemble::generate_pwm_2()
{

  if (length_counts.size() != 1)
  {
    std::cout << "Error: PWM can only be generated for aligned ensembles\n";
    std::exit(1);
  }
  std::ofstream ofs{ file_name + ".pwm_2" };
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  pwm_2 =
      std::vector<std::vector<double>>(C(L, 2), std::vector<double>(D * D, 0));
  for (auto const &sequence : sequences)
  {
    int pos = 0;
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        pwm_2[pos++][symbol_index(sequence.sequence[i]) * D +
                     symbol_index(sequence.sequence[j])]++;
  }

  for (auto &row : pwm_2)
  {
    for (auto &val : row)
    {
      val /= sequences.size();
      ofs << val << " ";
    }
    ofs << "\n";
  }
}

void
    Ensemble::generate_pwm_3()
{

  if (length_counts.size() != 1)
  {
    std::cout << "Error: PWM can only be generated for aligned ensembles\n";
    std::exit(1);
  }
  std::ofstream ofs{ file_name + ".pwm_3" };
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());

  ofs << D << "\n" << L << "\n";
  for (auto const &symbol : symbols)
  {
    ofs << symbol;
  }
  ofs << "\n";
  pwm_3 = std::vector<std::vector<double>>(C(L, 3),
                                           std::vector<double>(D * D * D, 0));
  for (auto const &sequence : sequences)
  {
    int pos = 0;
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          pwm_3[pos++][symbol_index(sequence.sequence[i]) * D * D +
                       symbol_index(sequence.sequence[j]) * D +
                       symbol_index(sequence.sequence[k])]++;
  }

  for (auto &row : pwm_3)
  {
    for (auto &val : row)
    {
      val /= sequences.size();
      ofs << val << " ";
    }
    ofs << "\n";
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
    score += std::log(D * (pwm_1[i][symbol_index(sequence.sequence[i])] + c)) /
             std::log(D);
  return score;
}

void
    Ensemble::calculate_true_scores_1()
{
  std::ofstream ofs{ file_name + ".ts_1" };
  for (auto const &sequence : sequences)
  {
    auto const score = calculate_individual_score_1(sequence);
    true_scores_1.push_back(score);
    ofs << score << "\n";
  }
}

void
    Ensemble::calculate_true_scores_2()
{
  auto const    c = 0.000001;
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());
  std::ofstream ofs{ file_name + ".ts_2" };
  for (auto const &sequence : sequences)
  {
    auto score = 0.;
    int  pos   = 0;
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        score += std::log(D * D *
                          (pwm_2[pos++][symbol_index(sequence.sequence[i]) * D +
                                        symbol_index(sequence.sequence[j])] +
                           c)) /
                 std::log(D);
    true_scores_2.push_back(score);
    ofs << score << "\n";
  }
}

void
    Ensemble::calculate_true_scores_3()
{
  auto const    c = 0.000001;
  auto const    D = static_cast<int>(symbol_counts.size());
  auto const    L = static_cast<int>(sequences[0].sequence.size());
  std::ofstream ofs{ file_name + ".ts_3" };
  for (auto const &sequence : sequences)
  {
    auto score = 0.;
    int  pos   = 0;
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          score +=
              std::log(
                  D * D * D *
                  (pwm_3[pos++][symbol_index(sequence.sequence[i]) * D * D +
                                symbol_index(sequence.sequence[j]) * D +
                                symbol_index(sequence.sequence[k])] +
                   c)) /
              std::log(D);
    true_scores_3.push_back(score);
    ofs << score << "\n";
  }
}

void
    Ensemble::load_tests(std::string const &fn)
{
  file_name = fn;
  std::ifstream ifs{ file_name };
  if (not ifs.is_open())
  {
    std::cout << "Error: test file " << file_name << " not found";
    std::exit(1);
  }
  std::cout << "Loading test file " << file_name << " ...\n";
  std::string sequence;
  std::string weight;
  std::string id;

  std::ofstream ofs{ file_name + ".scores" };
  while (std::getline(ifs, sequence, ',') and std::getline(ifs, weight, ',') and
         std::getline(ifs, id))
  {
    id.pop_back();
    // std::cout << sequence << ":" << weight << ":" << id << "\n";
    auto sequence_with_data = Sequence({ sequence, std::stod(weight), id });

    if (sequence.length() != sequences[0].sequence.length())
    {
      std::cout << "Warning: sequence with id: " << id << " has length "
                << sequence.length() << " but ensemble sequences have length "
                << sequences[0].sequence.length()
                << ". Skipping this sequence ...\n";
      continue;
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
      std::cout << "Warning: sequence with id: " << id << " contains a symbol '"
                << *f
                << "' that was never seen in the ensemble. Skipping this "
                   "sequence ...\n";
      continue;
    }
    ofs << id << "," << calculate_individual_score_1(sequence_with_data)
        << "\n";
  }
}
