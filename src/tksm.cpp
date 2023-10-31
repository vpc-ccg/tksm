#include <fmt/core.h>

#include <set>
#include <string>
#include <vector>

// Import kisims
#include "abundance.h"
#include "filter.h"
#include "head.h"
#include "model_truncation.h"
#include "pcr.h"
#include "polyA.h"
#include "random_wgs.h"
#include "scb.h"
#include "sequence.h"
#include "shuffle.h"
#include "strand_man.h"
#include "tag.h"
#include "transcribe.h"
#include "truncate.h"
#include "unsegment.h"
#include "util.h"
#include "mutate.h"

using std::set;
using std::string;
using std::vector;

#ifndef VERSION
#define VERSION "0.0.0"
#endif

// clang-format off
vector<std::pair<string, string>> kisims = {
    {"transcribe", "Builds RNA transcript molecules given abundances"},
    {"tag", "Adds a tag to each molecule given tag pattern"},
    {"polyA", "Adds polyA tail to each molecule"},
    {"scb", "Adds single cell barcode (from CB tags added by transcribe)"},
    {"pcr", "Simulates PCR amplification"},
    {"flip", "Simulates strand flipping"},
    {"truncate", "Simulates read truncation"},
    {"shuffle", "Shuffles an mdf file"},
    {"sequence", "Simulates reads given molecules"},
    {"random-wgs", "Simulates random WGS reads"},
    {"unsegment", "Concatenate adjacent molecules with random probability"},
    {"mutate", "Mutate molecules given mutations"},
};

vector<std::pair<string,string>> utility = {
    {"abundance", "Computes the abundance of a long read RNA-seq experiment"},
    {"model-truncation", "Kernel density estimation"},
    {"model-errors", "Models sequencing errors"},
    {"model-qscores", "Models sequencing quality scores"},
    {"head", "Prints the first n molecules of a file"},
    {"filter", "Filters a file based on a condition"},
};

vector<string> info = {
    "version",
    "help"
};
// clang-format on

void
help(char **argv, auto file) {
    fmt::print(file, "{}\n", ASCII_ART);
    fmt::print(file, "Usage: {} <module> [options]", argv[0]);
    fmt::print("\nAvailable modules: \n");
    for (const auto &[kisim, description] : kisims) {
        fmt::print("\t{}\t\t{}\n", kisim, description);
    }

    fmt::print("\nAvailable utilities: \n");
    for (auto [util, description] : utility) {
        fmt::print("\t{}\t\t{}\n", util, description);
    }

    fmt::print("----------------\n");
    for (auto info : info) {
        fmt::print("\t{}\t\n", info);
    }
}

int
main(int argc, char **argv) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);
    if (argc == 1) {
        help(argv, stderr);
        return 1;
    }

    set<string> all_kisims;
    for (const auto &[kisim, description] : kisims) {
        all_kisims.insert(kisim);
    }
    for (const auto &[util, description] : utility) {
        all_kisims.insert(util);
    }

    all_kisims.insert(info.begin(), info.end());

    if (all_kisims.find(argv[1]) == all_kisims.end()) {
        // check if --help is called
        for (int i = 1; i < argc; i++) {
            if (string(argv[i]) == "--help") {
                help(argv, stdout);
                return 0;
            }
        }
    }

    string kisim = argv[1];

    if (kisim == "version") {
        fmt::print("Version: {}\n", VERSION);
        return 0;
    }
    else if (kisim == "help") {
        help(argv, stdout);
        return 0;
    }
    else if (kisim == "abundance") {
        return Abundance_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "transcribe") {
        return Splicer_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "tag") {
        return TAG_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "polyA") {
        return PolyA_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "scb") {
        return SingleCellBarcoder_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "pcr") {
        return PCR_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "truncate") {
        return Truncate_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "flip") {
        return StrandMan_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "sequence") {
        return Sequencer_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "model-truncation") {
        return KDE_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "head") {
        return Head_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "filter") {
        return Filter_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "shuffle") {
        return Shuffle_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "unsegment") {
        return Unsegment_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "mutate"){
        return Mutate_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "model-errors") {
        fmt::print("Model errors\n");
    }
    else if (kisim == "model-qscores") {
        fmt::print("Model qscores\n");
    }
    else if (kisim == "random-wgs") {
        return RWGS_module{argc - 1, argv + 1}.run();
    }
    else {
        fmt::print(stderr, "Unknown kisim: {}\n", kisim);
        return 1;
    }

    return 0;
}
