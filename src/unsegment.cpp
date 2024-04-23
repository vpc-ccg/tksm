#include "unsegment.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using namespace std::string_literals;

#include "pimpl.h"

class Unsegment_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "i,input",
                "input mdf file",
                cxxopts::value<string>()
            )(
                "o,output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "p,probability",
                "probability of concatenation",
                cxxopts::value<double>()
            ) 
            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv)
        : tksm_module{"Unsegment", "Concatenate adjacent molecules with a probability"}, args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "probability"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                loge("{} is required!", param);
                ++missing_parameters;
            }
        }
        // Other parameter checks here

        if (missing_parameters > 0) {
            fmt::print(stderr, "{}\n", options.help());
            return 1;
        }
        return 0;
    }
    int run() {
        if (process_utility_arguments(args)) {
            return 0;
        }
        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        string input_file = args["input"].as<string>();

        string output_file = args["output"].as<string>();

        double probability = args["probability"].as<double>();

        ifstream input(input_file);

        ofstream output(output_file);

        molecule_descriptor current;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (auto &md : stream_mdf(input)) {
            if (current.get_id() == "") {
                current = md;
                continue;
            }
            if (dist(rand_gen) < probability) {
                current.concat(md);

                current.add_comment("Cat", md.get_id());
            }
            else {
                output << current;
                current = md;
            }
        }
        return 0;
    }

    void describe_program() {
        logi("Running {}", program_name);
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Unsegment_module);
