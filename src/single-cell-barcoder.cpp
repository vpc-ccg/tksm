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

    cxxopts::Options options("RNAInfuser Single-Cell barcoder", "Single-Cell barcoder module");

    options.add_options()
        ("i,input",  "input mdf file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("f,fasta", "FASTA file containing barcode sequences (should be given to sequencer module)", cxxopts::value<string>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("keep-meta-barcodes", "Keep the barcodes in the mdf metadata", cxxopts::value<bool>())
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"input","output", "fasta"};

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

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);

    string mdf_file_path {args["input"].as<string>()};
    std::ifstream mdf_file {mdf_file_path};
    
    //vector<molecule_descriptor> molecules = parse_mdf(mdf_file, true);

    auto streamer = stream_mdf(mdf_file, true);

    std::ofstream outfile{args["output"].as<string>()};
    std::ofstream fastafile{args["fasta"].as<string>()};
    std::unordered_set<string> used_barcodes;
    while(streamer){
        molecule_descriptor md = streamer();
        const string &barcode_str = md.get_comment("CB")[0];
        const string barcode_ctg_id = "CB_" + barcode_str;
        if( barcode_str != "."){
            md.prepend_segment(ginterval{barcode_ctg_id, 0, (int)barcode_str.size(), "+"});
            if( used_barcodes.count(barcode_str) == 0){
                fastafile << ">" << barcode_ctg_id << "\n" << barcode_str << "\n";
                used_barcodes.insert(barcode_str);
            }
        }
        if( !args["keep-meta-barcodes"].as<bool>()){
            md.drop_comment("CB");
        }
        outfile << md;
    }

    return 0;
}
