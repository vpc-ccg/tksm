
#include "util.h"
#include <fmt/core.h>

#include <vector>
#include <string>
#include <set>

//Import modules
#include "kde.h"
#include "abundance.h"
#include "sequencer.h"

#include "umi.h"

using std::vector;
using std::string;
using std::set;

#ifndef VERSION
#define VERSION "0.0.0"
#endif

vector<string> modules = {
    "abundance",
    "splicer",
    "umi",
    "polya",
    "single-cell-barcoder",
    "truncate",
    "fusion",
    "sequencer"
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

void help(char **argv, auto file){
    fmt::print(file, "Usage: {} <module> [options]", argv[0]);
    fmt::print("\nAvailable modules: \n");
    for (auto module : modules){
        fmt::print("\t{}\t\t{}\n", module, "description");
    }

    fmt::print("\nAvailable utilities: \n");
    for (auto util : utility){
        fmt::print("\t{}\t\t{}\n", util,"description");
    }
}

int main(int argc, char**argv){


    if (argc == 1){
        help(argv, stderr);
        return 1;
    }
    
    set<string> all_modules = set<string>(modules.begin(), modules.end());
    all_modules.insert(utility.begin(), utility.end());
    all_modules.insert(info.begin(), info.end());

    if (all_modules.find(argv[1]) == all_modules.end()){
        //check if --help is called
        for(int i = 1; i < argc; i++){
            if (string(argv[i]) == "--help"){
                help(argv, stdout);
                return 0;
            }
        }
    }
    
    string module = argv[1];
    
    if(module == "version"){
        fmt::print("Version: {}\n", VERSION);
        return 0;
    }
    else if(module == "help"){
        help(argv, stdout);
        return 0;
    }
    else if(module == "abundance"){
        run_abundance(argc - 1, argv + 1);
    }
    else if(module == "splicer"){
        fmt::print("Splicer\n");
    }
    else if(module == "umi"){
        return UMI_module{argc - 1, argv + 1}.run();
    }
    else if(module == "polya"){
        return Polya_module{argc - 1, argv + 1}.run();
    }
    else if(module == "single-cell-barcoder"){
        fmt::print("Single cell barcoder\n");
    }
    else if(module == "truncate"){
        fmt::print("Truncate\n");
    }
    else if(module == "fusion"){
        fmt::print("Fusion\n");
    }
    else if(module == "sequencer"){

        //Inject default badread model path into argv
        //If one provided by user, it will override this
        //This is a hack to avoid environment variables
        const char *path_cmd = "--badread-model-path=" INSTALL_PATH "/badread_models";
        argc--;
        argv = argv + 1;
        fmt::print("{}\n", path_cmd);
        char **argv_cpy = (char **)malloc((argc+2)*sizeof(char *));
        argv_cpy[0] = argv[0];
        for(int i=1; i<argc; i++){
            argv_cpy[i + 1] = argv[i];
        } 
        argv_cpy[1] = (char *)path_cmd;
        run_sequencer(argc + 1, argv_cpy);
    }
    else if(module == "kde"){
        run_kde(argc - 1, argv + 1);
    }
    else if(module == "head"){
        fmt::print("Head\n");
    }
    else if(module == "model-errors"){
        fmt::print("Model errors\n");
    }
    else if(module == "model-qscores"){
        fmt::print("Model qscores\n");
    }
    else{
        fmt::print(stderr, "Unknown module: {}\n", module);
        return 1;
    }


    return 0;
}
