
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

#include "ensemble.hpp"

namespace sic
{
class PWM_1
{
private:
  std::map<std::tuple<int, char>, double> pwm;
  std::vector<std::map<char, double>>     pwm_t;
  Summary                                 summary;

public:
  double
      evaluate(std::string const &sequence, bool use_threads, double c, bool use_bias) const;
  PWM_1(Ensemble const &ensemble, bool use_threads);
  PWM_1() = default;
};

class WT_PWM_1
{
private:
  std::map<std::tuple<int, char>, double> wt_pwm;
  double                                  wt_score;
  Summary                                 summary;

public:
  double
      evaluate(Ensemble const &ensemble, Mutant const &mutant, double c, bool use_bias) const;

  WT_PWM_1(Ensemble const &ensemble, double c, bool use_bias);
  WT_PWM_1() = default;
};

class PWM_2
{
private:
  std::map<std::tuple<int, int, char, char>, double>         pwm;
  std::vector<std::map<std::tuple<int, char, char>, double>> pwm_t;
  Summary                                                    summary;

public:
  double
      evaluate(std::string const &sequence, bool use_threads, double c, bool use_bias) const;
  PWM_2(Ensemble const &ensemble, bool use_threads);
  PWM_2() = default;
};

class WT_PWM_2
{
private:
  std::map<std::tuple<int, int, char, char>, double> wt_pwm;
  double                                             wt_score;
  Summary                                            summary;

public:
  double
      evaluate(Ensemble const &ensemble, Mutant const &mutant, double c, bool use_bias) const;

  WT_PWM_2(Ensemble const &ensemble, double c, bool use_bias);
  WT_PWM_2() = default;
};

class PWM_3
{
private:
  std::map<std::tuple<int, int, int, char, char, char>, double>         pwm;
  std::vector<std::map<std::tuple<int, int, char, char, char>, double>> pwm_t;
  Summary                                                               summary;

public:
  double
      evaluate(std::string const &sequence, bool use_threads, double c, bool use_bias) const;
  PWM_3(Ensemble const &ensemble, bool use_threads);
  PWM_3() = default;
};

class WT_PWM_3
{
private:
  std::map<std::tuple<int, int, int, char, char, char>, double> wt_pwm;
  double                                                        wt_score;
  Summary                                                       summary;

public:
  double
      evaluate(Ensemble const &ensemble, Mutant const &mutant, double c, bool use_bias) const;

  WT_PWM_3(Ensemble const &ensemble, double c, bool use_bias);
  WT_PWM_3() = default;
};

class PWM_4
{
private:
  std::map<std::tuple<int, int, int, int, char, char, char, char>, double> pwm;
  Summary summary;

public:
  double
      evaluate(std::string const &sequence, bool use_threads, double c, bool use_bias) const;
  PWM_4(Ensemble const &ensemble, bool use_threads);
  PWM_4() = default;
};

struct Mutant
{
  std::string                        descriptor;
  std::vector<std::tuple<int, char>> mutations;
  bool                               valid_mutation;
};

std::tuple<PWM_1, PWM_2, PWM_3, PWM_4>
    generatePWMs(Ensemble const &ensemble, int order, bool use_threads);

void test(std::string const                            &out_file_name,
          std::string const                            &train_column,
          std::vector<Sequence> const                  &sequences,
          std::tuple<PWM_1, PWM_2, PWM_3, PWM_4> const &pwms,
          int                                           order,
          double                                        c);

void testA2M(std::string const                            &out_file_name,
             std::string const                            &train_file,
             std::string const                            &true_wild_type,
             std::tuple<PWM_1, PWM_2, PWM_3, PWM_4> const &pwms,
             int                                           order,
             int                                           true_offset,
             bool                                          use_threads,
             double                                        c, bool use_bias);

void testA2MWithoutPWMs(std::string const &out_file_name,
                        std::string const &train_file,
                        std::string const &true_wild_type,
                        Ensemble const    &ensemble,
                        int                order,
                        int                true_offset,
                        bool               use_threads,
                        double             c, bool use_bias);

bool mutateSequence(std::string            &sequence,
                    std::string const      &col,
                    std::string const      &true_wild_type,
                    std::vector<int> const &valid_positions,
                    int                     true_offset,
                    std::ofstream          &fails);

std::vector<Mutant> generateMutants(std::string const      &train_file,
                                    std::string const      &true_wild_type,
                                    std::vector<int> const &valid_positions,
                                    int                     true_offset,
                                    std::ofstream          &fails);

}   // namespace sic
