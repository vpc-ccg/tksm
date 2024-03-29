#ifndef _MODULE_H_
#define _MODULE_H_

#include <cxxopts.hpp>

#include "util.h"

template <class RAND_GENERATOR>
class tksm_submodule_base {
    string submodule_name;
    string submodule_description;

public:
    RAND_GENERATOR &rand_gen;
    enum class submodule_status { RUN, DONT_RUN, ERROR } status;

    submodule_status update_status(submodule_status new_status) {
        status = new_status;
        return status;
    }

    tksm_submodule_base(string submodule_name, string submodule_description, RAND_GENERATOR &rand_gen)
        : submodule_name(submodule_name),
          submodule_description(submodule_description),
          rand_gen(rand_gen),
          status{submodule_status::DONT_RUN} {}

    virtual void add_options(cxxopts::Options &options) = 0;

    virtual submodule_status receive_arguments(const cxxopts::ParseResult &args) = 0;

    virtual ~tksm_submodule_base()                                  = default;
    virtual void describe_program(const cxxopts::ParseResult &args) = 0;
};

using tksm_submodule = tksm_submodule_base<std::mt19937>;

template <class RAND_GENERATOR>
class tksm_module_base {
protected:
    vector<string> added;
    struct OptionAdder : public cxxopts::OptionAdder {
        tksm_module_base<RAND_GENERATOR> *owner;
        OptionAdder(cxxopts::OptionAdder &&opt, tksm_module_base<RAND_GENERATOR> *owner) : cxxopts::OptionAdder(opt), owner{owner} {}
        OptionAdder &operator () (const string &opt, const string &help, const std::shared_ptr<const cxxopts::Value>& value 
      = ::cxxopts::value<bool>(),
    std::string arg_help = ""){
            owner->added.push_back(opt);
            cxxopts::OptionAdder::operator() (opt, help, value, arg_help);
            return *this;
        }
    };
    struct Options : public cxxopts::Options {
        tksm_module_base<RAND_GENERATOR> *owner;
        Options (const string &name, const string &desc, tksm_module_base<RAND_GENERATOR> *owner) : cxxopts::Options{name ,desc }, owner{owner} {}
        OptionAdder add_options(const std::string &group=""){
            return OptionAdder{cxxopts::Options::add_options(group), owner};
        }

    };
    string program_name;
    string program_description;
    Options options;
    RAND_GENERATOR rand_gen;

    virtual int validate_arguments() = 0;

public:
    void print_option_names(const cxxopts::ParseResult &result) const{
        for(const auto&lines : added){
            fmt::print("{}\n",rsplit(lines,",").back());
        }
    }

    tksm_module_base(const string &program_name, const string &program_description)
        : program_name(program_name),
          program_description(program_description),
          options{program_name, program_description, this} {
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
                "list",
                "list available options"
            )(
                "verbosity",
                fmt::format("Verbosity level to choosen from [{}]",LogLevels::log_choices()),
                cxxopts::value<string>()->default_value("INFO")
            )(
                "log-file",
                "Log file to write the logs (or stderr/stdout)",
                cxxopts::value<string>()->default_value("stderr")
            );
        // clang-format on
    };
    virtual ~tksm_module_base() = default;
    int process_utility_arguments(const cxxopts::ParseResult &args) {
        int seed = args["seed"].as<int>();
        rand_gen.seed(seed);

        string verbosity = args["verbosity"].as<string>();
        fmtlog::setLogLevel(LogLevels::parse_loglevel(verbosity));

        string log_file = args["log-file"].as<string>();
        if (log_file == "stderr") {
            fmtlog::setLogFile(stderr, false);
        }
        else if (log_file == "stdout") {
            fmtlog::setLogFile(stdout, false);
        }
        else {
            fmtlog::setLogFile(log_file.c_str(), true);
        }

        return help_or_version_is_used(args);
    }
    int list(const cxxopts::ParseResult &args) const{
        if(args.count("list") > 0){
            print_option_names(args);
            return 1;
        }
        return 0;
    }
    int help_or_version_is_used(const cxxopts::ParseResult &args) const {
        // Used fmt to print
        if (args.count("help") > 0) {
            fmt::print("{}\n", options.help());
            return 1;
        }
        if (args.count("version") > 0) {
            fmt::print("version: {}\n", VERSION);
            return 1;
        }
        if (args.count("list") > 0){
            print_option_names(args);
            return 1;
        }
        return 0;
    }

    virtual int run()               = 0;
    virtual void describe_program() = 0;
};

using tksm_module = tksm_module_base<std::mt19937>;

#endif
