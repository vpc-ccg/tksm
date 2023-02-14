#ifndef _MODULE_H_
#define _MODULE_H_

#include <cxxopts/cxxopts.hpp>

#include "util.h"

class tksm_module {
    // Pure abstract class for the modules
    // The modules are the classes that will be used to
    // implement the different functionalities of the tksm
    // Each module will have private functions to parse the parameters
    // using the cxxopts library and a function to validate the given parameters
    // modules will be constructed using argc and argv from the main function
    // and will be called using the run function which will be implemented
    // in each module

protected:
    string program_name;
    string program_description;
    cxxopts::Options options;

    virtual int validate_arguments() = 0;

public:
    tksm_module(string program_name, string program_description)
        : program_name(program_name),
          program_description(program_description),
          options{program_name, program_description} {
        // clang-format off
        options.add_options("utility")
            (
                "h,help",
                "Print help"
            )(
                "version",
                "Print version"
            )(
                "s,seed",
                "random seed",
                cxxopts::value<int>()->default_value("42")
            )(
                "v,verbosity",
                fmt::format("Verbosity level to choosen from [{}]",LogLevels::log_choices()),
                cxxopts::value<string>()->default_value("INFO")
            );
        // clang-format on
    };
    int help_or_version_is_used(cxxopts::ParseResult &args) {
        // Used fmt to print
        if (args.count("help") > 0) {
            fmt::print("{}\n", options.help());
            return 1;
        }
        if (args.count("version") > 0) {
            fmt::print("version: {}\n", VERSION);
            return 1;
        }
        return 0;
    }

    virtual int run()               = 0;
    virtual void describe_program() = 0;
};

#endif
