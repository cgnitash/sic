
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
  const double c =
      0.000001;   // should be user-provided eventually, currently hard-coded
  std::map<std::tuple<int, char>, double> pwm;
  Summary                                 summary;

public:
  double evaluate(std::string const &sequence) const;
  PWM_1(Ensemble const &ensemble);
  PWM_1() = default;
};

class PWM_2
{
  const double c =
      0.000001;   // should be user-provided eventually, currently hard-coded
  std::map<std::tuple<int, int, char, char>, double> pwm;
  Summary                                            summary;

public:
  double evaluate(std::string const &sequence) const;
  PWM_2(Ensemble const &ensemble);
  PWM_2() = default;
};

class PWM_3
{
  const double c =
      0.000001;   // should be user-provided eventually, currently hard-coded
  std::map<std::tuple<int, int, int, char, char, char>, double> pwm;
  Summary                                                       summary;

public:
  double evaluate(std::string const &sequence) const;
  PWM_3(Ensemble const &ensemble);
  PWM_3() = default;
};

class PWM_4
{
  const double c =
      0.000001;   // should be user-provided eventually, currently hard-coded
  std::map<std::tuple<int, int, int, int, char, char, char, char>, double> pwm;
  Summary summary;

public:
  double evaluate(std::string const &sequence) const;
  PWM_4(Ensemble const &ensemble);
  PWM_4() = default;
};

std::tuple<PWM_1, PWM_2, PWM_3, PWM_4> generatePWMs(Ensemble const &ensemble,
                                                    int             order);

void test(std::string const                            &out_file_name,
          std::string const                            &train_column,
          std::vector<Sequence> const                  &sequences,
          std::tuple<PWM_1, PWM_2, PWM_3, PWM_4> const &pwms,
          int                                           order);
}   // namespace sic
