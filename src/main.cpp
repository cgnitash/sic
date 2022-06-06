
#include <fstream>
#include <stdexcept>

#include "clap.hpp"
#include "ensemble.hpp"

int
    main(int argc, char **argv)
try
{
  CommandLineArgParser c;
  c.add_argument("Ensemble",
                 "Ensemble file containing sequences, and possibly metadata",
                 { "-e", "--ensemble" },
                 {},
                 "");
  c.add_argument("Sequence Field",
                 "Column in file containing sequences",
                 { "-s", "--sequence" },
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
  c.add_argument(
      "Train Field",
      "Select field (column name) used to generate PWMs (no argument will use "
      "entire ensemble for training)",
      { "-t", "--train-column" },
      {},
      " ");   // space is sentinel to indicate no argument
  c.add_argument("Train Value",
                 "Select value used to generate PWMs (no argument will use "
                 "entire ensemble for training)",
                 { "-v", "--train-value" },
                 {},
                 " ");   // space is sentinel to indicate no argument
  c.add_argument(
      "Weight",
      "Field used for weighting sequences (no argument uses unit weight)",
      { "-w", "--weight" },
      {},
      " ");   // space is sentinel to indicate no argument
  c.add_argument("Summarize",
                 "Display summary of training data (Y/N)",
                 { "-k", "--summarize" },   // k is summarize for now
                 { "Y", "yes", "N", "no" },
                 "N");
  c.add_argument("PWMSize",
                 "Correlation order for PWM",
                 { "-o", "--order" },
                 { "1", "2", "3", "4" },
                 "1");
  auto const args = c.parse_arguments(argc, argv);

  Ensemble e;
  try
  {
    std::stoi(args.at("Replicate"));
  }
  catch (std::invalid_argument const &)
  {
    std::cout << "Error: Replicate must be an integer\n";
    throw EnsembleError{};
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
    throw EnsembleError{};
  }

  e.load_ensemble(args.at("Sequence Field"),
                  args.at("Train Field"),
                  args.at("Train Value"),
                  std::stoi(args.at("Train Fraction")),
                  args.at("Weight"),
                  std::stoi(args.at("Replicate")),
                  args.at("Ensemble"));

  if (auto const summary = args.at("Summarize");
      summary == "Y" or summary == "yes")
    e.summary();

  e.generate_pwms(std::stoi(args.at("PWMSize")));

  timer([&] { e.run_tests(); });
}
catch (RuntimeError const &)
{
}
catch (EnsembleError const &)
{
}
catch (...)
{
  std::cout << "Internal bug: Unknown exception\n";
}
