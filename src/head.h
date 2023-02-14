#ifndef HEAD_H
#define HEAD_H
#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

class Head_module : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "n,count",
                "Number of items to be printed",
                cxxopts::value<int>()->default_value("10")
            );
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;
    std::mt19937 rand_gen;

public:
    Head_module(int argc, char **argv) : tksm_module{"head", "Head module"}, args(parse(argc, argv)) {}

    int validate_arguments() { return 0; }
    int run() {
        if (process_utility_arguments(args)) {
            return 0;
        }

        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        auto streamer = stream_mdf(std::cin, true);

        int count = args["count"].as<int>();
        for (int i = 0; (i < count) && streamer; ++i) {
            molecule_descriptor md = streamer();
            std::cout << md;
        }
        return 0;
    }

    void describe_program() {
        logi("Running Head module");
        logi("Number of items to be printed: {}", args["count"].as<int>());
        fmtlog::poll(true);
    }
};

#endif
