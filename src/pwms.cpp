
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pwms.hpp"

namespace sic
{

PWM_1::PWM_1(Ensemble const &ensemble):summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;
  pwm.clear();   // map had better be empty
  for (auto const &sequence : ensemble.sequences)
    for (int i = 0; i < L; ++i)
      pwm[{ i, sequence.sequence[i] }] += sequence.weight;

  for (auto &[key, val] : pwm)
  {
    val /= ensemble.summary.total_weight;
  }
}

PWM_2::PWM_2(Ensemble const &ensemble):summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;
  pwm.clear();   // map had better be empty
  for (auto const &sequence : ensemble.sequences)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        pwm[{ i, j, sequence.sequence[i], sequence.sequence[j] }] +=
            sequence.weight;
  }

  for (auto &[key, val] : pwm)
  {
    val /= ensemble.summary.total_weight;
  }
}

PWM_3::PWM_3(Ensemble const &ensemble):summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;
  pwm.clear();   // map has to be empty
  for (auto const &sequence : ensemble.sequences)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          pwm[{ i,
                j,
                k,
                sequence.sequence[i],
                sequence.sequence[j],
                sequence.sequence[k] }] += sequence.weight;
  }

  for (auto &[key, val] : pwm)
  {
    val /= ensemble.summary.total_weight;
  }
}

PWM_4::PWM_4(Ensemble const &ensemble):summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;
  pwm.clear();   // map has to be empty
  for (auto const &sequence : ensemble.sequences)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          for (int l = k + 1; l < L; ++l)
            pwm[{ i,
                  j,
                  k,
                  l,
                  sequence.sequence[i],
                  sequence.sequence[j],
                  sequence.sequence[k],
                  sequence.sequence[l] }] += sequence.weight;
  }

  for (auto &[key, val] : pwm)
  {
    val /= ensemble.summary.total_weight;
  }
}

double
    PWM_1::evaluate(std::string const &sequence) const
{
  auto const L     = summary.L;
  auto const D     = summary.D;
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
  {
    auto val = pwm.find({ i, sequence[i] });

    score +=
        std::log(D * ((val != pwm.end() ? val->second : 0.) + c)) / std::log(D);
  }
  return score;
}

double
    PWM_2::evaluate(std::string const &sequence) const
{
  auto const L     = summary.L;
  auto const D     = summary.D;
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
    {
      auto val = pwm.find({ i, j, sequence[i], sequence[j] });

      score += std::log(D * D * ((val != pwm.end() ? val->second : 0.) + c)) /
               std::log(D);
    }
  return score;
}

double
    PWM_3::evaluate(std::string const &sequence) const
{
  auto const L     = summary.L;
  auto const D     = summary.D;
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
      for (int k = j + 1; k < L; ++k)
      {
        auto val = pwm.find({ i, j, k, sequence[i], sequence[j], sequence[k] });

        score +=
            std::log(D * D * D * ((val != pwm.end() ? val->second : 0.) + c)) /
            std::log(D);
      }
  return score;
}

double
    PWM_4::evaluate(std::string const &sequence) const
{
  auto const L     = summary.L;
  auto const D     = summary.D;
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
      for (int k = j + 1; k < L; ++k)
        for (int l = k + 1; l < L; ++l)
        {
          auto val = pwm.find({ i,
                                j,
                                k,
                                l,
                                sequence[i],
                                sequence[j],
                                sequence[l],
                                sequence[k] });

          score += std::log(D * D * D *
                            ((val != pwm.end() ? val->second : 0.) + c)) /
                   std::log(D);
        }
  return score;
}

std::tuple<PWM_1, PWM_2, PWM_3, PWM_4>
    generatePWMs(Ensemble const &ensemble, int order)
{
  switch (order)
  {
    case 1:
      return { PWM_1{ ensemble }, {}, {}, {} };
    case 2:
      return { PWM_1{ ensemble }, PWM_2{ ensemble }, {}, {} };
    case 3:
      return { PWM_1{ ensemble }, PWM_2{ ensemble }, PWM_3{ ensemble }, {} };
    case 4:
      return { PWM_1{ ensemble },
               PWM_2{ ensemble },
               PWM_3{ ensemble },
               PWM_4{ ensemble } };
    default:
      std::cout << "Error: PWM order must be between 1 and 4\n";
      throw EnsembleError{};
  }
}

void
    test(std::string const                            &out_file_name,
         std::string const                            &train_column,
         std::vector<Sequence> const                  &sequences,
         std::tuple<PWM_1, PWM_2, PWM_3, PWM_4> const &pwms,
         int                                           order)
{
  assert(order < 5 and order > 0);
  std::ofstream ofs{ out_file_name + ".scores" };

  ofs << train_column;
  switch (order)
  {
    case 4:
      ofs << ",score_4";
      [[fallthrough]];
    case 3:
      ofs << ",score_3";
      [[fallthrough]];
    case 2:
      ofs << ",score_2";
      [[fallthrough]];
    case 1:
      ofs << ",score_1";
  }

  ofs << "\n";
  for (auto const &sequence : sequences)
  {
    ofs << sequence.label;
    switch (order)
    {
      case 4:
        ofs << "," << std::get<3>(pwms).evaluate(sequence.sequence);
        [[fallthrough]];
      case 3:
        ofs << "," << std::get<2>(pwms).evaluate(sequence.sequence);
        [[fallthrough]];
      case 2:
        ofs << "," << std::get<1>(pwms).evaluate(sequence.sequence);
        [[fallthrough]];
      case 1:
        ofs << "," << std::get<0>(pwms).evaluate(sequence.sequence);
    }
    ofs << "\n";
  }
  std::cout << "All test sequences are scored.\n" << std::flush;
}
}   // namespace sic

