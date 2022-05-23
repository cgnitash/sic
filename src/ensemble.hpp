
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

class Sequence
{
public:
  std::string sequence;
  double      weight = 0.;
  std::string id;
};

class Ensemble
{
public:
  void load(std::string const &);
  void load_tests(std::string const &fn);

  void
      enable_weights()
  {
    use_weights = true;
  }

  void summary();

  void generate_pwms_and_true_scores(int);

private:
  std::string                         ensemble_file_name;
  std::vector<Sequence>               sequences;
  std::map<int, int>                  length_counts;
  std::map<char, std::pair<int, int>> symbol_counts;
  std::vector<char>                   symbols;
  int                                 pwm_order    = 0;
  double                              total_weight = 0;
  bool                                use_weights  = false;

  std::vector<std::vector<double>>                              pwm_1;
  std::map<std::tuple<int, int, char, char>, double>            pwm_2;
  std::map<std::tuple<int, int, int, char, char, char>, double> pwm_3;

  void   generate_pwm_1();
  void   calculate_true_scores_1();
  double calculate_individual_score_1(Sequence const &);

  void   generate_pwm_2();
  void   calculate_true_scores_2();
  double calculate_individual_score_2(Sequence const &);

  void   generate_pwm_3();
  void   calculate_true_scores_3();
  double calculate_individual_score_3(Sequence const &);

  int symbol_index(char) const;
  int C(int, int);

  bool check_validity(Sequence const &);
};

