
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

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
  void   load(std::string const &);
  void   summary();
  void   generate_pwm_1();
  void   generate_pwm_2();
  void   generate_pwm_3();
  void   calculate_true_scores_1();
  double calculate_individual_score_1(Sequence const &);
  void   calculate_true_scores_2();
  void   calculate_true_scores_3();
  int    symbol_index(char) const;
  int    C(int, int);

 void load_tests(std::string const &fn);
private:
  std::string                         file_name;
  std::vector<Sequence>               sequences;
  std::map<int, int>                  length_counts;
  std::map<char, std::pair<int, int>> symbol_counts;
  std::vector<std::vector<double>>    pwm_1;
  std::vector<std::vector<double>>    pwm_2;
  std::vector<std::vector<double>>    pwm_3;
  std::vector<char>                   symbols;
  std::vector<double>                 true_scores_1;
  std::vector<double>                 true_scores_2;
  std::vector<double>                 true_scores_3;
  double total_weight = 0;
};


