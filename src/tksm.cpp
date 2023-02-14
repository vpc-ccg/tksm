
#include <fmt/core.h>

#include <set>
#include <string>
#include <vector>

#include "util.h"

// Import kisims
#include "abundance.h"
// #include "fusion.h"
#include "head.h"
#include "kde.h"
#include "polyA.h"
#include "sequencer.h"
#include "single-cell-barcoder.h"
#include "splicer.h"
#include "truncate.h"
#include "umi.h"

using std::set;
using std::string;
using std::vector;

#ifndef VERSION
#define VERSION "0.0.0"
#endif

vector<string> kisims = {"abundance", "splicer", "umi",      "polyA", "single-cell-barcoder",
                         "truncate",  "fusion",  "sequencer"};

vector<string> utility = {"kde", "head", "model-errors", "model-qscores"};

vector<string> info = {"version", "help"};

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
        run_abundance(argc - 1, argv + 1);
    }
    else if (kisim == "splicer") {
        return splicer_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "umi") {
        return UMI_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "polyA") {
        return PolyA_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "single-cell-barcoder") {
        return SingleCellBarcoder_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "truncate") {
        return Truncate_module{argc - 1, argv + 1}.run();
    }
    else if (kisim == "fusion") {
        fmt::print("Fusion\n");
        return 1;
    }
    else if (kisim == "sequencer") {
        // Inject default badread model path into argv
        // If one provided by user, it will override this
        // This is a hack to avoid environment variables
        const char *path_cmd = "--badread-model-path=" INSTALL_PATH "/badread_models";
        argc--;
        argv = argv + 1;
        fmt::print("{}\n", path_cmd);
        char **argv_cpy = (char **)malloc((argc + 2) * sizeof(char *));
        argv_cpy[0]     = argv[0];
        for (int i = 1; i < argc; i++) {
            argv_cpy[i + 1] = argv[i];
        }
        argv_cpy[1] = (char *)path_cmd;
        return run_sequencer(argc + 1, argv_cpy);
    }
    else if (kisim == "kde") {
        return run_kde(argc - 1, argv + 1);
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
