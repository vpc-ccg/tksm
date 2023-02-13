#ifndef HEAD_H
#define HEAD_H
#include <string>
#include <vector>
#include "interval.h"
#include "mdf.h"
#include "util.h"
#include "module.h"


#include <random>
#include "cxxopts/cxxopts.hpp"

class Head_module : public tksm_module{

    cxxopts::ParseResult parse(int argc, char **argv){
        options.add_options("main")
            ("n,count", "Number of items to be printed", cxxopts::value<int>()->default_value("10"))
            ;
        return  options.parse(argc, argv);
    }
    
    cxxopts::ParseResult args;
    std::mt19937 rand_gen;
    public:
    Head_module( int argc, char **argv) : tksm_module{"head", "Head module"}, args(parse(argc, argv)){
    }

    int  validate_arguments(){

        return 0;
    }
    int run(){
        fmtlog::setLogLevel(LogLevels::parse_loglevel(args["verbosity"].as<string>()));
        fmtlog::flushOn(fmtlog::DBG);

        if(help_or_version_is_used(args)){
            return 0;
        }

        if(validate_arguments()){
            return 1;
        }
        describe_program();

        auto streamer = stream_mdf(std::cin, true);
    
        int count = args["count"].as<int>();
        for(int i = 0; (i < count) && streamer; ++i){
            molecule_descriptor md = streamer();
            std::cout << md;
        } 
        return 0;   
    }

    void describe_program(){  
        logi("Running Head module");
        logi("Number of items to be printed: {}", args["count"].as<int>());
        fmtlog::poll(true);
    }
};

#endif
