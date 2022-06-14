
#include <fstream>
#include <ranges>
#include <stdexcept>

#include "clap.hpp"
#include "ensemble.hpp"
#include "pwms.hpp"

int
    main(int argc, char **argv)
try
{
  CommandLineArgParser c;
  c.add_argument("Testing File",
                 "Testing File containing sequences",
                 { "-of", "--testing-file" },
                 {},
                 "");
  c.add_argument("Test Sequence Column",
                 "Column in testing file containing sequences",
                 { "-osc", "--test-sequence-column" },
                 {},
                 "");
  c.add_argument("Test Label Column",
                 "Column in test file containing true labels (optional)",
                 { "-olc", "--test-label-column" },
                 {},
                 "__");

  c.add_argument("Training File",
                 "Training File containing sequences, labels, and weights",
                 { "-if", "--training-file" },
                 {},
                 "");
  c.add_argument("Train Sequence Column",
                 "Column in train file containing sequences",
                 { "-isc", "--train-sequence-column" },
                 {},
                 "");
  c.add_argument("Train Label Column",
                 "Column in file containing labels (If none is provided, all "
                 "sequences will be used)",
                 { "-ilc", "--train-label-column" },
                 {},
                 "__");
  c.add_argument("Train Label Value",
                 "Select label used to generate PWMs (no argument will use "
                 "entire ensemble for training)",
                 { "-lv", "--train-label-value" },
                 {},
                 "__");   // space is sentinel to indicate no argument
  c.add_argument(
      "Weight Column",
      "Field used for weighting sequences (no argument uses unit weight)",
      { "-wc", "--weight-column" },
      {},
      "__");

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

  auto all_seqs =
      sic::extractSequencesFromFile(args.at("Training File"),
                                    args.at("Train Sequence Column"),
                                    args.at("Train Label Column"),
                                    args.at("Weight Column"));

  if (args.at("Train Label Value") == "__")
  {
    if (args.at("Train Label Column") != "__")
    {
      std::cout
          << "Error: Training column specified without value. If all sequences "
             "should be used to train, the column must not be specified.\n";
      throw sic::EnsembleError{};
    }
    std::cout << "No training column specified. All sequences will be used.\n";
  }

  sic::filter(all_seqs, args.at("Train Label Value"));

  auto const seqs = sic::sample(all_seqs,
                                std::stoi(args.at("Train Fraction")),
                                std::stoi(args.at("Replicate")));

  auto const ensemble = sic::Ensemble(seqs);

  if (auto const summary = args.at("Summarize");
      summary == "Y" or summary == "yes")
    ensemble.print_summary();

  auto const all_pwms =
      sic::generatePWMs(ensemble, std::stoi(args.at("PWMSize")));

  auto test_seqs =
      sic::extractSequencesFromFile(args.at("Testing File"),
                                    args.at("Test Sequence Column"),
                                    args.at("Test Label Column"),
                                    "__");

  auto const out_file_name =
      args.at("Training File") + "-" + args.at("Train Label Column") + "-" +
      args.at("Train Label Value") + "-" + args.at("Testing File") + "-" +
      args.at("Test Label Column") + "-" + args.at("Train Fraction") + "-" +
      args.at("Replicate");

  sic::test(out_file_name,
            args.at("Test Label Column"),
            test_seqs,
            all_pwms,
            std::stoi(args.at("PWMSize")));

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
