
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "clap.hpp"

void
    CommandLineArgParser::issue_runtime_diagnostic(std::string const &message)
{
  std::cout << "Error: " << message << "\nUsage:\n\t" << std::setw(15)
            << std::left << "Argument Name "
            << ":\t" << std::setw(20) << std::left << "Aliases"
            << ":\t" << std::setw(60) << std::left << "Description"
            << ":\tRequired / Defaults\n"
            << std::string(150, '-') << "\n";
  for (auto const &argument : arguments)
  {
    std::cout << "\t";
    std::ostringstream os;
    for (auto const &alias : argument.aliases)
      os << alias << "  ";
    std::cout << std::setw(15) << std::left << argument.name << ":\t"
              << std::setw(20) << std::left << os.str() << ":\t"
              << std::setw(60) << std::left << argument.description << ":\t("
              << (argument.default_value.empty()
                      ? std::string{ "Required" }
                      : std::string{ "Default - " } + argument.default_value)
              << ")\n";
  }
  std::cout << std::flush;
  throw RuntimeError{};
}
void
    CommandLineArgParser::add_argument(
        std::string const              &name,
        std::string const              &description,
        std::vector<std::string> const &aliases,
        std::vector<std::string> const &valid_values,
        std::string const              &error_message)
{
  if (all_names.find(name) != all_names.end())
  {
    std::cout << "Error: flag name : " << name << " is already used\n";
    throw RuntimeError{};
  }
  all_names.insert(name);
  for (auto const &alias : aliases)
    if (all_aliases.find(alias) != all_aliases.end())
    {
      std::cout << "Error: flag alias : " << alias << " is already used\n";
      throw RuntimeError{};
    }
    else
      all_aliases.insert(alias);

  arguments.push_back(
      { name, description, aliases, valid_values, error_message });
}

Argument
    CommandLineArgParser::argument_for_alias(std::string const &user_alias)
{

  for (auto const &argument : arguments)
    for (auto const &alias : argument.aliases)
      if (alias == user_alias)
      {
        return argument;
      }
  issue_runtime_diagnostic(std::string{ "flag " } + user_alias + " not valid");
  return {};   // will never execute, silences warning for control flowing off
               // the end. There must be a better solution
}

std::map<std::string, std::string>
    CommandLineArgParser::parse_arguments(int argc, char **argv)
{
  std::map<std::string, std::string> args;

  for (int i = 1; i < argc; ++i)
  {
    std::string arg(argv[i]);
    if (auto eq_pos = arg.find('='); eq_pos != std::string::npos)
    {
      auto alias          = arg.substr(0, eq_pos);
      auto value          = arg.substr(eq_pos + 1);
      auto argument       = argument_for_alias(alias);
      args[argument.name] = value;
    }
    else if (i == argc - 1)
    {
      issue_runtime_diagnostic(std::string{ "Error: trailing flag " } + arg +
                               " without argument\n");
    }
    else
    {
      auto alias          = arg;
      auto value          = std::string{ argv[i + 1] };
      auto argument       = argument_for_alias(alias);
      args[argument.name] = value;
      ++i;
    }
  }
  check_provided_values_are_valid(args);
  check_and_update_defaults(args);
  return args;
}

void
    CommandLineArgParser::check_provided_values_are_valid(
        std::map<std::string, std::string> const &args)
{
  for (auto const &[name, value] : args)
  {
    auto f = std::find_if(std::begin(arguments),
                          std::end(arguments),
                          [name = name](auto const &arg)
                          { return arg.name == name; });
    if (not f->valid_values.empty() and
        std::find(std::begin(f->valid_values),
                  std::end(f->valid_values),
                  value) == std::end(f->valid_values))
      issue_runtime_diagnostic(std::string{ "value '" } + value +
                               "' not valid for flag '" + name + "'\n");
  }
}
void
    CommandLineArgParser::check_and_update_defaults(
        std::map<std::string, std::string> &args)
{
  for (auto const &argument : arguments)
  {
    if (args.find(argument.name) != args.end())
      continue;

    if (argument.default_value.empty())
      issue_runtime_diagnostic(std::string{ "flag '" } + argument.name +
                               "' has no default value and must be specified");

    args[argument.name] = argument.default_value;
  }
}

