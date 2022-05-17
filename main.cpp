
#include "ensemble.hpp"
#include<fstream>

int
    main(int argc, char **argv)
{
  if (argc == 1)
  {
    std::cout << "Error: no ensemble file provided\n";
    std::exit(1);
  }
  Ensemble e;
  e.load(argv[1]);
  e.summary();
  e.generate_pwm_1();
  e.calculate_true_scores_1();
  //e.generate_pwm_2();
  //e.calculate_true_scores_2();
  //e.generate_pwm_3();
  //e.calculate_true_scores_3();
  if (argc == 2)
  {
    std::cout << "Warning: no test file provided, so no tests will be run\n";
    std::exit(1);
  }
  e.load_tests(argv[2]);
}
