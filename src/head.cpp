#include "umi.h"

#include <iostream>
#include <iterator>
#include <string>
#include <fstream>
#include <random>
#include <variant>

#include <cxxopts/cxxopts.hpp>

#include "interval.h"
#include "mdf.h"
#include "util.h"

using std::string;

std::mt19937 rand_gen{std::random_device{}()};

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser head module", "head");

    options.add_options()
        ("n,count", "Number of molecule description", cxxopts::value<int>()->default_value("10"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {};

    int missing_parameters = 0;
    for( string &param : mandatory){
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }

    if(missing_parameters  > 0){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    auto streamer = stream_mdf(std::cin, true);


    for( int i = 0; i < args["count"].as<int>() && streamer; ++i){
        std::cout << streamer();
    }
    return 0;
}