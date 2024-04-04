#include "strand_man.h"

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

#include "pimpl.h"

class StrandMan_module::impl : public tksm_module {
    ParseResult parse(int argc, char **argv) {
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
                "p,flip-probability",
                "probability of flipping the strand of a molecule",
                cxxopts::value<double>()
            )
            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    ParseResult args;

    std::uniform_real_distribution<double> dist{0.0, 1.0};

    auto strand_flip_transformer(double flip_probability = 0) {
        return std::ranges::views::transform([&, flip_probability](auto &md) {
            if (dist(rand_gen) < flip_probability) {
                return flip_molecule(md);
            }
            return md;
        });
    }

public:
    impl(int argc, char **argv) : tksm_module{"<MODULE>", "<MODULE> description"}, args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "flip-probability"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                report_missing_parameter(param);
                ++missing_parameters;
            }
        }

        if (missing_parameters > 0) {
            fmt::print(stderr, "{}\n", options.help());
            return 1;
        }

        // Other parameter checks here

        if (args["flip-probability"].as<double>() < 0.0 || args["flip-probability"].as<double>() > 1.0) {
            loge("Flip probability must be between 0 and 1");
            ++missing_parameters;
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

        ifstream input(input_file);

        ofstream output(output_file);

        double flip_probability = args["flip-probability"].as<double>();
        for (const auto &md :
             stream_mdf(args["input"].as<string>(), true) | strand_flip_transformer(flip_probability)) {
            output << md;
        }
        return 0;
    }

    void describe_program() {
        logi("Running [Flip module]");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());

        logi("Molecule flip probability: {}%", args["flip-probability"].as<double>() * 100.0);

        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(StrandMan_module);
