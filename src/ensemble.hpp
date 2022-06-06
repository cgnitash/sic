
#include <chrono>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

struct EnsembleError
{
};

struct SequenceData
{
  std::string sequence;
  double      weight;
};

class Ensemble
{
public:
  void load_ensemble(std::string sequence_field,
                     std::string train_field,
                     std::string train_value,
                     int         fraction,
                     std::string weight_field,
                     int         replicate,
                     std::string file_name);
  void run_tests() const;
  void summary() const;
  void generate_pwms(int);

private:
  std::vector<SequenceData> sequences;

  // initialized in load_ensemble
  std::string                         file_name;
  std::string                         sequence_field;
  int                                 sequence_field_index;
  std::string                         train_field;
  int                                 train_field_index;
  std::string                         weight_field;
  int                                 weight_field_index;
  std::string                         train_value;
  int                                 fraction;
  int                                 replicate;
  int                                 D;
  int                                 L;
  double                              total_weight = 0;
  std::map<int, int>                  length_counts;
  std::map<char, std::pair<int, int>> symbol_counts;
  std::vector<char>                   symbols;

  const double c =
      0.000001;   // should be user-provided eventually, currently hard-coded

  // constructed in generate_pwm*
  int                                                           pwm_order = 0;
  std::map<std::tuple<int, char>, double>                       pwm_1;
  std::map<std::tuple<int, int, char, char>, double>            pwm_2;
  std::map<std::tuple<int, int, int, char, char, char>, double> pwm_3;
  std::map<std::tuple<int, int, int, int, char, char, char, char>, double>
      pwm_4;

  void   generate_pwm_1();
  double calculate_individual_score_1(std::string const &) const;

  void   generate_pwm_2();
  double calculate_individual_score_2(std::string const &) const;

  void   generate_pwm_3();
  double calculate_individual_score_3(std::string const &) const;

  void   generate_pwm_4();
  double calculate_individual_score_4(std::string const &) const;

  bool check_validity(std::string const &) const;
};

std::vector<std::string> split(std::string const &, char delim = ',');

template <typename Func>
void
    timer(Func func)
// should also support functions taking arguments (maybe make it a lambda)
{
  // from https://en.cppreference.com/w/cpp/chrono
  auto start = std::chrono::steady_clock::now();
  func();
  auto                          end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}
