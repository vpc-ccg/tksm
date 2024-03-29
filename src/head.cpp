#include "head.h"

#include <cxxopts.hpp>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <variant>

#include "interval.h"
#include "mdf.h"
#include "util.h"

using std::string;

std::mt19937 rand_gen{std::random_device{}()};

int
main(int argc, char **argv) {
    cxxopts::Options options("tksm head module", "head");

    options.add_options()("n,count", "Number of molecule description", cxxopts::value<int>()->default_value("10"))(
        "h,help", "Help screen");
    auto args = options.parse(argc, argv);

    if (args.count("help") > 0) {
        fmt::print("{}\n", options.help());
        return 0;
    }
    std::vector<string> mandatory = {};

    int missing_parameters = 0;
    for (string &param : mandatory) {
        if (args.count(param) == 0) {
            report_missing_parameter(param);
            ++missing_parameters;
        }
    }

    if (missing_parameters > 0) {
        std::cerr << options.help() << std::endl;
        return 1;
    }

    int cnt = args["count"].as<int>();
    for (const auto &md : stream_mdf(std::cin)) {
        std::cout << md;
        if (--cnt == 0) {
            break;
        }
    }
    return 0;
}
