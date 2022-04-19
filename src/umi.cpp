
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <variant>

#include "extern/cxxopts.hpp"

#include "interval.h"
#include "mdf.h"

#define FMT_HEADER_ONLY 1
#include <fmt/format.h>

using std::string;

std::mt19937 rand_gen{std::random_device{}()};

class num2seq{
    public:
        std::array<char,4> table  {{'A','T','C','G'}};
        uint64_t max;
        num2seq(uint64_t max):max(max){

        }
        string operator [] (uint64_t num) const{
            std::stringstream st;
            int back = 0;
            for(uint64_t mask = 3; mask < max; mask <<=2){
                st << table[ (num & mask) >> back];
                back+=2;
            }
            return st.str();
        }
};


void add_UMIs(vector<pcr_copy> &copies, ostream &umifile, int len){

    uint64_t lim_from_bitc = std::pow(4, len);
    num2seq make_seq(lim_from_bitc);
    std::uniform_int_distribution<uint64_t> dist(0,lim_from_bitc);


    int index = 0;
    for(pcr_copy &pcp : copies){
        uint64_t umi_int = dist(rand_gen);
        string umi_seq = make_seq[umi_int];
        string umi_ctg_name = fmt::format("RI_umi_ctg_{}",index);

        umifile << fmt::format(">{}\n", umi_ctg_name);
        umifile << umi_seq << "\n";
        pcp.prepend(ginterval{umi_ctg_name, 0, len, "+"});
        ++index;
    }
}

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser UMI", "UMI module of RNAInfuser");

    options.add_options()
        ("i,input",  "input mdf file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("a,umi-reference", "Output umi reference file", cxxopts::value<string>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("length", "Length of the UMI sequences", cxxopts::value<int>()->default_value("0"))
        ("back-ratio", "Ratio of umi's inserted to the 3' of the molecule", cxxopts::value<double>()->default_value("0.0"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"input","output", "umi-reference", "length"};

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

    int length = args["length"].as<int>();
    string mdf_file_path {args["input"].as<string>()};
    std::ifstream mdf_file {mdf_file_path};
    
    vector<pcr_copy> molecules = parse_mdf(mdf_file);



    vector<pcr_copy> nm;
    for(const pcr_copy &pcp : molecules){
        for(int i = 0; i < pcp.depth; ++i){
            nm.push_back(pcp);
            nm.back().depth = 1;
            nm.back().id += ("_" + std::to_string(i));
        }
    }
    molecules = nm;

    std::string umi_ref_file = args["umi-reference"].as<string>();
    std::ofstream umi_file{umi_ref_file};

    add_UMIs(molecules, umi_file, length);

    string outfile_name = args["output"].as<string>();
    std::ofstream outfile{outfile_name};
    print_all_mdf(outfile, molecules);
    
    return 0;
}
