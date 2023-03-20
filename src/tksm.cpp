#include <fmt/core.h>

#include <set>
#include <string>
#include <vector>

#include "util.h"

// Import kisims
#include "abundance.h"
#include "head.h"
#include "kde.h"
#include "pcr.h"
#include "polyA.h"
#include "sequencer.h"
#include "single-cell-barcoder.h"
#include "splicer.h"
#include "truncate.h"
#include "tag.h"
#include "strand_man.h"

using std::set;
using std::string;
using std::vector;

#ifndef VERSION
#define VERSION "0.0.0"
#endif

// clang-format off
vector<string> kisims = {
    "abundance",
    "splicer",
    "tag",
    "polyA",
    "single-cell-barcoder",
    "pcr",
    "flip",
    "truncate",
    "sequencer",
};

vector<string> utility = {
    "kde",
    "head",
    "model-errors",
    "model-qscores"
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
    for (auto kisim : kisims) {
        fmt::print("\t{}\t\t{}\n", kisim, "description");
    }

    fmt::print("\nAvailable utilities: \n");
    for (auto util : utility) {
        fmt::print("\t{}\t\t{}\n", util, "description");
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

    set<string> all_kisims = set<string>(kisims.begin(), kisims.end());
    all_kisims.insert(utility.begin(), utility.end());
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
    else if (kisim == "splicer") {
        return Splicer_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "tag") {
        return TAG_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "polyA") {
        return PolyA_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "single-cell-barcoder") {
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
    else if (kisim == "sequencer") {
        return Sequencer_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "kde") {
        return KDE_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "head") {
        return Head_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "model-errors") {
        fmt::print("Model errors\n");
    }
    else if (kisim == "model-qscores") {
        fmt::print("Model qscores\n");
    }
    else {
        fmt::print(stderr, "Unknown kisim: {}\n", kisim);
        return 1;
    }

    return 0;
}
