
#include <chrono>
#include <fstream>
#include <ranges>
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

  auto all_seqs = sic::extractA2MSequencesFromFile(args.at("Training File"));

  auto const ensemble = sic::Ensemble(all_seqs, true);

  if (auto const summary = args.at("Summarize");
      summary == "Y" or summary == "yes")
    ensemble.print_summary();

  auto start = std::chrono::system_clock::now();

  auto const all_pwms =
      sic::generatePWMs(ensemble, std::stoi(args.at("PWMSize")));

  auto end = std::chrono::system_clock::now();
  std::cout << "time to generate ";
  printTime(end - start);

  auto out_file_name = args.at("Testing File");
  if (auto slash = out_file_name.find_last_of('/'); slash != std::string::npos)
    out_file_name = out_file_name.substr(slash + 1);

  out_file_name = out_file_name.substr(0, out_file_name.find('.'));

  start = std::chrono::system_clock::now();
  sic::testA2M(out_file_name,
               args.at("Testing File"),
               all_seqs[0].sequence,
               all_pwms,
               std::stoi(args.at("PWMSize")),
               true);
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
