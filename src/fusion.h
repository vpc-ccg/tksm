#ifndef _FUSION_H_
#define _FUSION_H_
#include <random>
#include <string>
#include <vector>

#include "cxxopts/cxxopts.hpp"
#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

class Fusion_module : public tksm_module {
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
            );
        // clang-format on

        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;


public:
    Fusion_module(int argc, char **argv) : tksm_module{"fusion", "Fusion module"}, args(parse(argc, argv)) {}

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                std::cerr << param << " is required!\n";
                ++missing_parameters;
            }
        }

        if (missing_parameters > 0) {
            std::cerr << options.help() << std::endl;
            return 1;
        }
        return 0;
    }
    int run() {
        if(process_utility_arguments(args)) {
            return 0;
        }

        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        string mdf_file_path{args["input"].as<string>()};
        std::ifstream mdf_file{mdf_file_path};
        auto streamer = stream_mdf(mdf_file);

        return 0;
    }

    void describe_program() {
        logi("Running Fusion module");

        fmtlog::poll(true);
    }
};

#endif
