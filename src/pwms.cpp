
#include <algorithm>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
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

#include "ensemble.hpp"
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
    PWM_1::evaluate(std::string const &sequence,
                    bool               use_threads,
                    double             c,
                    bool               use_bias) const
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

      auto const biased_D =
          use_bias ? static_cast<double>(summary.N * summary.L) /
                         summary.symbol_counts.at(sequence[i]).first
                   : D;

      score +=
          std::log(biased_D * ((val != pwm.end() ? val->second : 0.) + c)) /
          std::log(D);
    }
  }
  return score;
}

double
    PWM_2::evaluate(std::string const &sequence,
                    bool               use_threads,
                    double             c,
                    bool               use_bias) const
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

        auto const biased_D_i =
            use_bias ? static_cast<double>(summary.N * summary.L) /
                           summary.symbol_counts.at(sequence[i]).first
                     : D;

        auto const biased_D_j =
            use_bias ? static_cast<double>(summary.N * summary.L) /
                           summary.symbol_counts.at(sequence[j]).first
                     : D;

        score += std::log(biased_D_i * biased_D_j *
                          ((val != pwm.end() ? val->second : 0.) + c)) /
                 std::log(D);
      }
    }
  return score;
}

double
    PWM_3::evaluate(std::string const &sequence,
                    bool               use_threads,
                    double             c,
                    bool) const
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
    PWM_4::evaluate(std::string const &sequence,
                    bool               use_threads,
                    double             c,
                    bool) const
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

std::tuple<int, char, bool>
    checkValidMutation(std::string const      &mutation,
                       std::string const      &true_wild_type,
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
  if (std::islower(true_wild_type[index]))
    valid_mutation = false;

  if (valid_mutation)
  {
    index = valid_positions[index];
  }
  return { index, m[3].str()[0], valid_mutation };
}

std::vector<Mutant>
    generateMutants(std::string const      &train_file,
                    std::string const      &true_wild_type,
                    std::vector<int> const &valid_positions,
                    int                     true_offset,
                    std::ofstream          &fails)
{
  std::vector<Mutant> mutants;

  std::ifstream ifs{ train_file };
  if (not ifs.is_open())
  {
    std::cout << "Error: file " << train_file << " not found";
    throw EnsembleError{};
  }

  std::string line;
  while (std::getline(ifs, line))
  {
    if (line[0] == '#')   // skip comments
      continue;
    auto col = line.substr(0, line.find(';'));
    if (col == "mutant")   // skip header
      continue;
    Mutant mutant;
    mutant.valid_mutation = true;
    mutant.descriptor     = col;
    if (col != "WT" and col != "wt")
    {
      auto const mutations = split(col, ',');
      for (auto const &mutation : mutations)
      {
        auto const [position, replacement, valid_mutation] = checkValidMutation(
            mutation, true_wild_type, valid_positions, true_offset, fails);

        mutant.valid_mutation &= valid_mutation;
        mutant.mutations.push_back({ position, replacement });
      }
    }
    mutants.push_back(mutant);
  }
  return mutants;
}

WT_PWM_1::WT_PWM_1(Ensemble const &ensemble, double c, bool use_bias)
    : summary(ensemble.summary)
{
  wt_score = 0.0;

  ensemble.verify();
  auto const L = ensemble.summary.L;
  auto const D = ensemble.summary.D;

  auto const wild_type = ensemble.sequences[0].sequence;

  wt_pwm.clear();
  for (auto const &sequence : ensemble.sequences)
    for (int i = 0; i < L; ++i)
      if (sequence.sequence[i] == wild_type[i])
        wt_pwm[{ i, wild_type[i] }] +=
            sequence.weight / ensemble.summary.total_weight;

  for (int i = 0; i < L; ++i)
  {
    if (summary.symbol_counts.find(wild_type[i]) == summary.symbol_counts.end())
    {
      std::cout << "err: " << wild_type[i] << "not found\n";
      ensemble.print_summary();
    }
    auto const biased_D = use_bias
                              ? static_cast<double>(summary.N * summary.L) /
                                    summary.symbol_counts.at(wild_type[i]).first
                              : D;

    wt_score +=
        std::log(biased_D * (wt_pwm[{ i, wild_type[i] }] + c)) / std::log(D);
  }
}

double
    WT_PWM_1::evaluate(Ensemble const &ensemble,
                       Mutant const   &mutant,
                       double          c,
                       bool            use_bias) const
{
  auto const D = ensemble.summary.D;

  auto const wild_type = ensemble.sequences[0].sequence;

  auto score = wt_score;

  for (auto const &[pos, rep] : mutant.mutations)
  {
    auto const wild_type_biased_D =
        use_bias ? static_cast<double>(summary.N * summary.L) /
                       summary.symbol_counts.at(wild_type[pos]).first
                 : D;

    score -= std::log(wild_type_biased_D *
                      (wt_pwm.at({ pos, wild_type[pos] }) + c)) /
             std::log(D);

    auto combo_score = 0.0;
    for (auto const &sequence : ensemble.sequences)
      if (sequence.sequence[pos] == rep)
        combo_score += sequence.weight / ensemble.summary.total_weight;

    auto const biased_D = use_bias
                              ? static_cast<double>(summary.N * summary.L) /
                                    summary.symbol_counts.at(rep).first
                              : D;

    score += std::log(biased_D * (combo_score + c)) / std::log(D);
  }

  return score;
}

WT_PWM_2::WT_PWM_2(Ensemble const &ensemble, double c, bool use_bias)
    : summary(ensemble.summary)
{
  wt_score = 0.0;

  ensemble.verify();
  auto const L = ensemble.summary.L;
  auto const D = ensemble.summary.D;

  auto const wild_type = ensemble.sequences[0].sequence;

  wt_pwm.clear();
  for (auto const &sequence : ensemble.sequences)
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        if (sequence.sequence[i] == wild_type[i] and
            sequence.sequence[j] == wild_type[j])
          wt_pwm[{ i, j, wild_type[i], wild_type[j] }] +=
              sequence.weight / ensemble.summary.total_weight;

  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
    {
      auto const biased_D_i =
          use_bias ? static_cast<double>(summary.N * summary.L) /
                         summary.symbol_counts.at(wild_type[i]).first
                   : D;

      auto const biased_D_j =
          use_bias ? static_cast<double>(summary.N * summary.L) /
                         summary.symbol_counts.at(wild_type[j]).first
                   : D;

      wt_score += std::log(biased_D_i * biased_D_j *
                           (wt_pwm[{ i, j, wild_type[i], wild_type[j] }] + c)) /
                  std::log(D);
    }
}

double
    WT_PWM_2::evaluate(Ensemble const &ensemble,
                       Mutant const   &mutant,
                       double          c,
                       bool            use_bias) const
{
  auto const D = ensemble.summary.D;
  auto const L = ensemble.summary.L;

  auto const wild_type = ensemble.sequences[0].sequence;

  auto score = wt_score;

  for (auto const &[pos, rep] : mutant.mutations)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        if (i == pos or j == pos)
        {
          auto const wild_type_biased_D_i =
              use_bias ? static_cast<double>(summary.N * summary.L) /
                             summary.symbol_counts.at(wild_type[i]).first
                       : D;

          auto const wild_type_biased_D_j =
              use_bias ? static_cast<double>(summary.N * summary.L) /
                             summary.symbol_counts.at(wild_type[j]).first
                       : D;
          score -=
              std::log(wild_type_biased_D_i * wild_type_biased_D_j *
                       (wt_pwm.at({ i, j, wild_type[i], wild_type[j] }) + c)) /
              std::log(D);
        }

    std::vector<double> adjusted_scores(L, 0.0);

    for (auto const &sequence : ensemble.sequences)
      for (int i = 0; i < L; ++i)
        if (i != pos and sequence.sequence[i] == wild_type[i] and
            sequence.sequence[pos] == rep)
          adjusted_scores[i] += sequence.weight / ensemble.summary.total_weight;

    for (int i = 0; i < L; ++i)
      if (i != pos)
      {
        auto const biased_D_i =
            use_bias ? static_cast<double>(summary.N * summary.L) /
                           summary.symbol_counts.at(wild_type[i]).first
                     : D;

        auto const biased_D_rep =
            use_bias ? static_cast<double>(summary.N * summary.L) /
                           summary.symbol_counts.at(rep).first
                     : D;
        score +=
            std::log(biased_D_i * biased_D_rep * (adjusted_scores[i] + c)) /
            std::log(D);
      }
  }

  return score;
}

WT_PWM_3::WT_PWM_3(Ensemble const &ensemble, double c, bool)
    : summary(ensemble.summary)
{
  wt_score = 0.0;

  ensemble.verify();
  auto const L = ensemble.summary.L;
  auto const D = ensemble.summary.D;

  auto const wild_type = ensemble.sequences[0].sequence;

  wt_pwm.clear();
  for (auto const &sequence : ensemble.sequences)
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          if (sequence.sequence[i] == wild_type[i] and
              sequence.sequence[j] == wild_type[j] and
              sequence.sequence[k] == wild_type[k])
            wt_pwm[{ i, j, k, wild_type[i], wild_type[j], wild_type[k] }] +=
                sequence.weight / ensemble.summary.total_weight;

  for (int i = 0; i < L; ++i)
    for (int j = i + 1; j < L; ++j)
      for (int k = j + 1; k < L; ++k)
        wt_score +=
            std::log(
                D * D * D *
                (wt_pwm[{ i, j, k, wild_type[i], wild_type[j], wild_type[k] }] +
                 c)) /
            std::log(D);
}

double
    WT_PWM_3::evaluate(Ensemble const &ensemble,
                       Mutant const   &mutant,
                       double          c,
                       bool) const
{
  auto const D = ensemble.summary.D;
  auto const L = ensemble.summary.L;

  auto const wild_type = ensemble.sequences[0].sequence;

  auto score = wt_score;

  for (auto const &[pos, rep] : mutant.mutations)
  {
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        for (int k = j + 1; k < L; ++k)
          if (i == pos or j == pos or k == pos)
            score -= std::log(D * D * D *
                              (wt_pwm.at({ i,
                                           j,
                                           k,
                                           wild_type[i],
                                           wild_type[j],
                                           wild_type[k] }) +
                               c)) /
                     std::log(D);

    std::vector<std::vector<double>> adjusted_scores(
        L, std::vector<double>(L, 0.0));

    for (auto const &sequence : ensemble.sequences)
      for (int i = 0; i < L; ++i)
        for (int j = i + 1; j < L; ++j)
          if (i != pos and j != pos and sequence.sequence[i] == wild_type[i] and
              sequence.sequence[j] == wild_type[j] and
              sequence.sequence[pos] == rep)
            adjusted_scores[i][j] +=
                sequence.weight / ensemble.summary.total_weight;

    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j)
        if (i != pos and j != pos)
          score +=
              std::log(D * D * D * (adjusted_scores[i][j] + c)) / std::log(D);
  }

  return score;
}

void
    testA2MWithoutPWMs(std::string const &out_file_name,
                       std::string const &train_file,
                       std::string const &true_wild_type,
                       Ensemble const    &ensemble,
                       int                order,
                       int                true_offset,
                       bool,
                       double c,
                       bool   use_bias)
{
  assert(order < 5 and order > 0);
  std::ofstream ofs{ out_file_name + ".scores" };

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

  std::vector<int> valid_positions;
  int              counter = 0;
  for (unsigned char c : true_wild_type)
    valid_positions.push_back(std::islower(c) ? counter : counter++);

  auto wild_type = true_wild_type;
  wild_type.erase(std::remove_if(std::begin(wild_type),
                                 std::end(wild_type),
                                 [](unsigned char c)
                                 { return std::islower(c); }),
                  std::end(wild_type));

  std::ofstream fails{ out_file_name + ".fails" };
  auto const   &mutants = generateMutants(
      train_file, true_wild_type, valid_positions, true_offset, fails);

  auto const wt_pwm_1 = WT_PWM_1{ ensemble, c, use_bias };
  WT_PWM_2   wt_pwm_2;
  if (order > 1)
    wt_pwm_2 = WT_PWM_2{ ensemble, c, use_bias };
  WT_PWM_3 wt_pwm_3;
  if (order > 2)
    wt_pwm_3 = WT_PWM_3{ ensemble, c, use_bias };

  for (auto const &mutant : mutants)
  {
    ofs << mutant.descriptor;
    switch (order)
    {
        /*
      case 4:
        if (not mutant.valid_mutation)
          ofs << ";";
        else
           ofs << ";" << std::get<3>(pwms).evaluate(sequence, use_threads);
        [[fallthrough]];
          */
      case 3:
        if (not mutant.valid_mutation)
          ofs << ";";
        else
          // ofs << ";" << std::get<2>(pwms).evaluate(sequence, use_threads);
          ofs << ";" << wt_pwm_3.evaluate(ensemble, mutant, c, use_bias);
        [[fallthrough]];
      case 2:
        if (not mutant.valid_mutation)
          ofs << ";";
        else
          // ofs << ";" << std::get<1>(pwms).evaluate(sequence, use_threads);
          ofs << ";" << wt_pwm_2.evaluate(ensemble, mutant, c, use_bias);
        [[fallthrough]];
      case 1:
        if (not mutant.valid_mutation)
          ofs << ";";
        else
          //  ofs << ";" << std::get<0>(pwms).evaluate(sequence, use_threads);
          ofs << ";" << wt_pwm_1.evaluate(ensemble, mutant, c, use_bias);
    }
    ofs << "\n";
  }
  std::cout << "All test sequences are scored.\n" << std::flush;
}

void
    testA2M(std::string const                            &out_file_name,
            std::string const                            &train_file,
            std::string const                            &true_wild_type,
            std::tuple<PWM_1, PWM_2, PWM_3, PWM_4> const &pwms,
            int                                           order,
            int                                           true_offset,
            bool                                          use_threads,
            double                                        c,
            bool                                          use_bias)
{
  assert(order < 5 and order > 0);
  std::ofstream ofs{ out_file_name + ".scores" };

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

  std::vector<int> valid_positions;
  int              counter = 0;
  for (unsigned char c : true_wild_type)
    valid_positions.push_back(std::islower(c) ? counter : counter++);

  auto wild_type = true_wild_type;
  wild_type.erase(std::remove_if(std::begin(wild_type),
                                 std::end(wild_type),
                                 [](unsigned char c)
                                 { return std::islower(c); }),
                  std::end(wild_type));

  std::ofstream fails{ out_file_name + ".fails" };
  auto const   &mutants = generateMutants(
      train_file, true_wild_type, valid_positions, true_offset, fails);

  for (auto const &[descriptor, mutations, valid_mutation] : mutants)
  {
    ofs << descriptor;
    auto sequence = wild_type;
    for (auto [pos, rep] : mutations)
      sequence[pos] = rep;
    switch (order)
    {
      case 4:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";"
              << std::get<3>(pwms).evaluate(sequence, use_threads, c, use_bias);
        [[fallthrough]];
      case 3:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";"
              << std::get<2>(pwms).evaluate(sequence, use_threads, c, use_bias);
        [[fallthrough]];
      case 2:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";"
              << std::get<1>(pwms).evaluate(sequence, use_threads, c, use_bias);
        [[fallthrough]];
      case 1:
        if (not valid_mutation)
          ofs << ";";
        else
          ofs << ";"
              << std::get<0>(pwms).evaluate(sequence, use_threads, c, use_bias);
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
         int                                           order,
         double                                        c,
         bool                                          use_bias)
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
        ofs << ","
            << std::get<3>(pwms).evaluate(
                   sequence.sequence, false, c, use_bias);
        [[fallthrough]];
      case 3:
        ofs << ","
            << std::get<2>(pwms).evaluate(
                   sequence.sequence, false, c, use_bias);
        [[fallthrough]];
      case 2:
        ofs << ","
            << std::get<1>(pwms).evaluate(
                   sequence.sequence, false, c, use_bias);
        [[fallthrough]];
      case 1:
        ofs << ","
            << std::get<0>(pwms).evaluate(
                   sequence.sequence, false, c, use_bias);
    }
    ofs << "\n";
  }
  std::cout << "All test sequences are scored.\n" << std::flush;
}
}   // namespace sic

