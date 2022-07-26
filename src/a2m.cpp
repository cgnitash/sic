
#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "clap.hpp"
#include "ensemble.hpp"
#include "pwms.hpp"

template <typename Time>
void
    printTime(Time t)
{
  auto hrs  = std::chrono::duration_cast<std::chrono::hours>(t);
  auto mins = std::chrono::duration_cast<std::chrono::minutes>(t - hrs);
  auto secs = std::chrono::duration_cast<std::chrono::seconds>(t - hrs - mins);
  auto mils = std::chrono::duration_cast<std::chrono::milliseconds>(
      t - hrs - mins - secs);
  std::cout << hrs.count() << ":" << mins.count() << ":" << secs.count() << ":"
            << mils.count() << std::endl;
}

int
    main(int argc, char **argv)
try
{
  CommandLineArgParser c;
  c.add_argument("Testing File",
                 "Testing File containing mutations",
                 { "-of", "--testing-file" },
                 {},
                 "");

  c.add_argument("Training File",
                 "Training File containing sequences, labels, and weights",
                 { "-if", "--training-file" },
                 {},
                 "");

  c.add_argument("Train Fraction",
                 "Fraction (per 100) of valid sequences to be used in training",
                 { "-f", "--fraction" },
                 {},
                 "100");

  c.add_argument("Replicate",
                 "Replicate number used as seed for training fraction",
                 { "-r", "--replicate" },
                 {},
                 "1");

  c.add_argument("Summarize",
                 "Display summary of training data (Y/N)",
                 { "-s", "--summarize" },
                 { "Y", "yes", "N", "no" },
                 "N");

  c.add_argument("PWMSize",
                 "Correlation order for PWM",
                 { "-o", "--order" },
                 { "1", "2", "3", "4" },
                 "1");

  c.add_argument("Multi Threaded",
                 "Use multiple threads (Y/N)",
                 { "-t", "--threads" },
                 { "Y", "yes", "N", "no" },
                 "N");

  c.add_argument("Use PWMs",
                 "Use PWMs generated from ensemble (Y/N)",
                 { "-u", "--use-pwms" },
                 { "Y", "yes", "N", "no" },
                 "Y");

  c.add_argument(
      "Adjust Weights",
      "Adjust weights according to similarity (Y/N/U) Y requires percentage",
      { "-a", "--adjust-weights" },
      { "U", "uniform", "Y", "yes", "N", "no" },
      "N");

  c.add_argument("Similarity Percentage",
                 "Percentage similarity threshold for weights",
                 { "-sim", "--similarity" },
                 {},
                 "100");

  c.add_argument("Pseudo Count",
                 "Pseudo-count value N -> 1/10^N",
                 { "-p", "--pseudo-count" },
                 {},
                 "6");

  auto const args = c.parse_arguments(argc, argv);

  try
  {
    std::stoi(args.at("Replicate"));
  }
  catch (std::invalid_argument const &)
  {
    std::cout << "Error: Replicate must be an integer\n";
    throw sic::EnsembleError{};
  }
  try
  {
    int fraction = std::stoi(args.at("Train Fraction"));
    if (fraction < 1 or fraction > 100)
      throw std::invalid_argument("");
  }
  catch (std::invalid_argument const &)
  {
    std::cout << "Error: Fraction must be an integer between 1 and 100\n";
    throw sic::EnsembleError{};
  }

  auto start = std::chrono::system_clock::now();
  auto end   = std::chrono::system_clock::now();

  auto [all_seqs, true_offset] =
      sic::extractA2MSequencesFromFile(args.at("Training File"));

  std::cout << "\n ---> Training file : " << args.at("Training File") << "\n";

  auto const true_target = all_seqs[0].sequence;
  sic::removeLowerCaseResidues(all_seqs, true_target);

  end = std::chrono::system_clock::now();
  std::cout << "time to extract and clean ";
  printTime(end - start);
  start = std::chrono::system_clock::now();

  auto const adjust_arg = args.at("Adjust Weights");

  if (adjust_arg == "U" or adjust_arg == "uniform")
  {
    start = std::chrono::system_clock::now();
    std::cout << "Similarity percentage will be ignored if provided...\n";
    sic::adjustWeightsUniformly(all_seqs);

    end = std::chrono::system_clock::now();
    std::cout << "time to adjust weights (uniform) ";
    printTime(end - start);
  }

  if (adjust_arg == "Y" or adjust_arg == "yes")
  {
    auto const sim_perc = std::stoi(args.at("Similarity Percentage"));
    start               = std::chrono::system_clock::now();
    sic::adjustWeights(all_seqs, sim_perc);

    end = std::chrono::system_clock::now();
    std::cout << "time to adjust weights (by similarity) ";
    printTime(end - start);
  }

  auto ensemble = sic::Ensemble(all_seqs);

  if (auto const summary = args.at("Summarize");
      summary == "Y" or summary == "yes")
    ensemble.print_summary();

  auto const thread_arg  = args.at("Multi Threaded");
  auto const use_threads = thread_arg == "Y" or thread_arg == "yes";

  auto const use_pwms_arg = args.at("Use PWMs");
  auto const use_pwms     = use_pwms_arg == "Y" or use_pwms_arg == "yes";

  std::tuple<sic::PWM_1, sic::PWM_2, sic::PWM_3, sic::PWM_4> all_pwms;
  if (use_pwms)
  {
    start = std::chrono::system_clock::now();
    all_pwms =
        sic::generatePWMs(ensemble, std::stoi(args.at("PWMSize")), use_threads);
    end = std::chrono::system_clock::now();
    std::cout << "time to generate ";
    printTime(end - start);
  }

  std::cout << "\n ---> Testing file : " << args.at("Testing File") << "\n";

  auto pseudo_count_arg = args.at("Pseudo Count");

  auto out_file_name = args.at("Testing File");
  if (auto slash = out_file_name.find_last_of('/'); slash != std::string::npos)
    out_file_name = out_file_name.substr(slash + 1);

  out_file_name = out_file_name.substr(0, out_file_name.find('.'));
  if (adjust_arg == "Y" or adjust_arg == "yes")
    out_file_name +=
        "_" + args.at("Similarity Percentage") + "_" + pseudo_count_arg;

  auto const pseudo_count = 1.0 / std::pow(10.0, std::stod(pseudo_count_arg));

  start = std::chrono::system_clock::now();
  if (use_pwms)
    sic::testA2M(out_file_name,
                 args.at("Testing File"),
                 true_target,
                 all_pwms,
                 std::stoi(args.at("PWMSize")),
                 true_offset,
                 use_threads,
                 pseudo_count);
  else
    sic::testA2MWithoutPWMs(out_file_name,
                            args.at("Testing File"),
                            true_target,
                            ensemble,
                            std::stoi(args.at("PWMSize")),
                            true_offset,
                            use_threads,
                            pseudo_count);

  end = std::chrono::system_clock::now();
  std::cout << "time to test ";
  printTime(end - start);

  return 0;
}
catch (RuntimeError const &)
{
}
catch (sic::EnsembleError const &)
{
}
catch (...)
{
  std::cout << "Internal bug: Unknown exception\n";
}
