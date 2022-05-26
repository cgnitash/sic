
#include <fstream>

#include "clap.hpp"
#include "ensemble.hpp"

int
    main(int argc, char **argv)
try
{
  CommandLineArgParser c;
  c.add_argument("Ensemble",
                 "Ensemble file used to generate PWM",
                 { "-e", "--ensemble" },
                 {},
                 "");
  c.add_argument("Test",
                 "Test file containing sequences to score based on PWM",
                 { "-t", "--test" },
                 {},
                 "");
  c.add_argument("UseWeights",
                 "Should sequences be weighted? (Y/N)",
                 { "-w", "--use-weights" },
                 { "Y", "yes", "N", "no" },
                 "N");
  c.add_argument("Summarize",
                 "Display summary of ensemble data (Y/N)",
                 { "-s", "--summarize" },
                 { "Y", "yes", "N", "no" },
                 "N");
  c.add_argument("PWMSize",
                 "Correlation order for PWM",
                 { "-o", "--order" },
                 { "1", "2", "3" },
                 "1");
  auto const args = c.parse_arguments(argc, argv);

  Ensemble e;
  e.load_ensemble(args.at("Ensemble"));

  if (auto const summary = args.at("UseWeights");
      summary == "Y" or summary == "yes")
    e.enable_weights();

  if (auto const summary = args.at("Summarize");
      summary == "Y" or summary == "yes")
    e.summary();

  e.generate_pwms_and_true_scores(std::stoi(args.at("PWMSize")));

  e.load_tests(args.at("Test"));
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
