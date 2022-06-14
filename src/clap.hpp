
#pragma once

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

class Argument
{
public:
  std::string              name;
  std::string              description;
  std::vector<std::string> aliases;
  std::vector<std::string> valid_values;
  std::string              default_value;
};

struct RuntimeError
{
};

class CommandLineArgParser
{
private:
  std::set<std::string> all_aliases;
  std::set<std::string> all_names;
  std::vector<Argument> arguments;

  void add_single_argument(std::map<std::string, std::string> &,
                           std::string const &,
                           std::string const &);

  Argument argument_for_alias(std::string const &);
  void     issue_runtime_diagnostic(std::string const &);

  void check_provided_values_are_valid(
      std::map<std::string, std::string> const &);
  void check_and_update_defaults(std::map<std::string, std::string> &);

public:
  // should also have a fluent version, or make this fluent, i.e. return *this
  void add_argument(std::string const &             name,
                    std::string const &             description,
                    std::vector<std::string> const &aliases,
                    std::vector<std::string> const &valid_values,
                    std::string const &             default_value);

  std::map<std::string, std::string> parse_arguments(int, char **);

};
