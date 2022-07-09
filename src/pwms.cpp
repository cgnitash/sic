
#include <algorithm>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pwms.hpp"

namespace sic
{

PWM_1::PWM_1(Ensemble const &ensemble, bool use_threads)
    : summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;
  if (not use_threads)
  {
    pwm.clear();   // map had better be empty
    for (auto const &sequence : ensemble.sequences)
      for (int i = 0; i < L; ++i)
        pwm[{ i, sequence.sequence[i] }] += sequence.weight;

    for (auto &[key, val] : pwm)
    {
      val /= ensemble.summary.total_weight;
    }
  }
  else
  {
    pwm_t = std::vector(L, std::map<char, double>{});

    std::vector<std::thread> v;
    for (int i = 0; i < L; ++i)
      v.emplace_back(
          [&, i_t = i]
          {
            for (auto const &sequence : ensemble.sequences)
              pwm_t[i_t][sequence.sequence[i_t]] += sequence.weight;

            for (auto &[key, val] : pwm_t[i_t])
              val /= ensemble.summary.total_weight;
          });

    for (int i = 0; i < L; ++i)
      v[i].join();
  }
}

PWM_2::PWM_2(Ensemble const &ensemble, bool use_threads)
    : summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;
  if (not use_threads)
  {
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
  else
  {
    pwm_t = std::vector(L, std::map<std::tuple<int, char, char>, double>{});

    std::vector<std::thread> v;
    for (int i = 0; i < L; ++i)
      v.emplace_back(
          [&, i_t = i]
          {
            for (auto const &sequence : ensemble.sequences)
              for (int j = i_t + 1; j < L; ++j)
                pwm_t[i_t]
                     [{ j, sequence.sequence[i_t], sequence.sequence[j] }] +=
                    sequence.weight;

            for (auto &[key, val] : pwm_t[i_t])
              val /= ensemble.summary.total_weight;
          });

    for (int i = 0; i < L; ++i)
      v[i].join();
  }
}

PWM_3::PWM_3(Ensemble const &ensemble, bool use_threads)
    : summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;

  if (not use_threads)
  {
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
  else
  {

    pwm_t = std::vector(
        L, std::map<std::tuple<int, int, char, char, char>, double>{});

    std::vector<std::thread> v;
    for (int i = 0; i < L; ++i)
      v.emplace_back(
          [&, i_t = i]
          {
            for (auto const &sequence : ensemble.sequences)
              for (int j = i_t + 1; j < L; ++j)
                for (int k = j + 1; k < L; ++k)
                  pwm_t[i_t][{ j,
                               k,
                               sequence.sequence[i_t],
                               sequence.sequence[j],
                               sequence.sequence[k] }] += sequence.weight;

            for (auto &[key, val] : pwm_t[i_t])
              val /= ensemble.summary.total_weight;
          });

    for (int i = 0; i < L; ++i)
      v[i].join();
  }
}

PWM_4::PWM_4(Ensemble const &ensemble, bool use_threads)
    : summary(ensemble.summary)
{
  ensemble.verify();
  auto const L = summary.L;
  if (not use_threads)
  {
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
  else
  {
    std::cout << "Generating PWMs of order 4 with threads not implemented. "
                 "PWMS of this size are not expected to fit in memory.\n";
    throw EnsembleError{};
  }
}

double
    PWM_1::evaluate(std::string const &sequence, bool use_threads) const
{
  auto const L     = summary.L;
  auto const D     = summary.D;
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
  {
    if (use_threads)
    {
      auto val = pwm_t[i].find(sequence[i]);

      score += std::log(D * ((val != pwm_t[i].end() ? val->second : 0.) + c)) /
               std::log(D);
    }
    else
    {
      auto val = pwm.find({ i, sequence[i] });

      score += std::log(D * ((val != pwm.end() ? val->second : 0.) + c)) /
               std::log(D);
    }
  }
  return score;
}

double
    PWM_2::evaluate(std::string const &sequence, bool use_threads) const
{
  auto const L     = summary.L;
  auto const D     = summary.D;
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
    {
      if (use_threads)
      {
        auto val = pwm_t[i].find({ j, sequence[i], sequence[j] });

        score +=
            std::log(D * D * ((val != pwm_t[i].end() ? val->second : 0.) + c)) /
            std::log(D);
      }
      else
      {
        auto val = pwm.find({ i, j, sequence[i], sequence[j] });

        score += std::log(D * D * ((val != pwm.end() ? val->second : 0.) + c)) /
                 std::log(D);
      }
    }
  return score;
}

double
    PWM_3::evaluate(std::string const &sequence, bool use_threads) const
{
  auto const L     = summary.L;
  auto const D     = summary.D;
  auto       score = 0.;
  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
      for (int k = j + 1; k < L; ++k)
      {
        if (use_threads)
        {
          auto val =
              pwm_t[i].find({ j, k, sequence[i], sequence[j], sequence[k] });

          score += std::log(D * D * D *
                            ((val != pwm_t[i].end() ? val->second : 0.) + c)) /
                   std::log(D);
        }
        else
        {
          auto val =
              pwm.find({ i, j, k, sequence[i], sequence[j], sequence[k] });

          score += std::log(D * D * D *
                            ((val != pwm.end() ? val->second : 0.) + c)) /
                   std::log(D);
        }
      }
  return score;
}

double
    PWM_4::evaluate(std::string const &sequence, bool use_threads) const
{
  if (use_threads)
  {
    std::cout << "Error: Can't evaluate with order 4 PWM.\n";
    throw EnsembleError{};
  }
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
    generatePWMs(Ensemble const &ensemble, int order, bool use_threads)
{
  switch (order)
  {
    case 1:
      return { PWM_1{ ensemble, use_threads }, {}, {}, {} };
    case 2:
      return {
        PWM_1{ ensemble, use_threads }, PWM_2{ ensemble, use_threads }, {}, {}
      };
    case 3:
      return { PWM_1{ ensemble, use_threads },
               PWM_2{ ensemble, use_threads },
               PWM_3{ ensemble, use_threads },
               {} };
    case 4:
      return { PWM_1{ ensemble, use_threads },
               PWM_2{ ensemble, use_threads },
               PWM_3{ ensemble, use_threads },
               PWM_4{ ensemble, use_threads } };
    default:
      std::cout << "Error: PWM order must be between 1 and 4\n";
      throw EnsembleError{};
  }
}

bool
    mutateSequence(std::string            &sequence,
                   std::string const      &mutation,
                   std::string const      &true_wild_type,
                   bool                    ignore_lower,
                   std::vector<int> const &valid_positions,
                   int                     true_offset,
                   std::ofstream          &fails)
{
  bool        valid_mutation = true;
  std::regex  r{ R"((\w)(\d+)(\w))" };
  std::smatch m;
  if (not std::regex_match(mutation, m, r))
  {
    std::cout << "Error: Not able to match " << mutation << "\n";
    throw EnsembleError{};
  }

  auto index = std::stoi(m[2].str()) - true_offset;

  if (index >= static_cast<int>(true_wild_type.length()) or index < 0)
  {
    fails << "Skipping: Can't test mutation at position " << index + 1
          << " when target only contains " << true_wild_type.length()
          << " positions" << std::endl;
    valid_mutation = false;
  }

  if (true_wild_type[index] != m[1].str()[0])
  {
    fails << "Target sequence does not work for " << mutation
          << ", target sequences contains '" << true_wild_type[index]
          << "' at that position" << std::endl;
    valid_mutation = false;
  }
  if (ignore_lower and std::islower(true_wild_type[index]))
    valid_mutation = false;

  if (valid_mutation)
  {
    if (ignore_lower)
      index = valid_positions[index];

    sequence[index] = m[3].str()[0];
  }
  return valid_mutation;
}

void
    testA2M(std::string const                            &out_file_name,
            std::string const                            &train_file,
            std::string const                            &true_wild_type,
            std::tuple<PWM_1, PWM_2, PWM_3, PWM_4> const &pwms,
            int                                           order,
            int                                           true_offset,
            bool                                          ignore_lower,
            bool                                          use_threads)
{
  assert(order < 5 and order > 0);
  std::ofstream ofs{ out_file_name + ".scores" };
  std::ofstream fails{ out_file_name + ".fails" };

  ofs << "label";
  switch (order)
  {
    case 4:
      ofs << ";score_4";
      [[fallthrough]];
    case 3:
      ofs << ";score_3";
      [[fallthrough]];
    case 2:
      ofs << ";score_2";
      [[fallthrough]];
    case 1:
      ofs << ";score_1";
  }

  ofs << "\n";
  std::ifstream ifs{ train_file };
  if (not ifs.is_open())
  {
    std::cout << "Error: file " << train_file << " not found";
    throw EnsembleError{};
  }

  std::vector<int> valid_positions;
  int              counter = 0;
  for (unsigned char c : true_wild_type)
    valid_positions.push_back(std::islower(c) ? counter : counter++);

  auto wild_type = true_wild_type;
  if (ignore_lower)
    wild_type.erase(std::remove_if(std::begin(wild_type),
                                   std::end(wild_type),
                                   [](unsigned char c)
                                   { return std::islower(c); }),
                    std::end(wild_type));

  std::string line;
  while (std::getline(ifs, line))
  {
    if (line[0] == '#')   // skip comments
      continue;
    auto col = line.substr(0, line.find(';'));
    if (col == "mutant")   // skip header
      continue;
    auto sequence       = wild_type;
    bool valid_mutation = true;
    if (col != "WT" and col != "wt")
    {
      auto const mutations = split(col, ',');
      for (auto const &mutation : mutations)
        valid_mutation &= mutateSequence(sequence,
                                         mutation,
                                         true_wild_type,
                                         ignore_lower,
                                         valid_positions,
                                         true_offset,
                                         fails);
    }
    ofs << col;
    switch (order)
    {
      case 4:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";" << std::get<3>(pwms).evaluate(sequence, use_threads);
        [[fallthrough]];
      case 3:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";" << std::get<2>(pwms).evaluate(sequence, use_threads);
        [[fallthrough]];
      case 2:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";" << std::get<1>(pwms).evaluate(sequence, use_threads);
        [[fallthrough]];
      case 1:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";" << std::get<0>(pwms).evaluate(sequence, use_threads);
    }
    ofs << "\n";
  }
  std::cout << "All test sequences are scored.\n" << std::flush;
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
        ofs << "," << std::get<3>(pwms).evaluate(sequence.sequence, false);
        [[fallthrough]];
      case 3:
        ofs << "," << std::get<2>(pwms).evaluate(sequence.sequence, false);
        [[fallthrough]];
      case 2:
        ofs << "," << std::get<1>(pwms).evaluate(sequence.sequence, false);
        [[fallthrough]];
      case 1:
        ofs << "," << std::get<0>(pwms).evaluate(sequence.sequence, false);
    }
    ofs << "\n";
  }
  std::cout << "All test sequences are scored.\n" << std::flush;
}
}   // namespace sic

