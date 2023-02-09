#include "umi.h"

#include <cassert>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>
#include <random>
#include <variant>

#include <cxxopts/cxxopts.hpp>

#include "interval.h"
#include "mdf.h"


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

class fmt2seq{
    public:
        std::array<char,4> table  {{'A','T','C','G'}};
        std::array<string, 128> lookup;
        string fmt;
        fmt2seq(const string &fmt): fmt(fmt){
#define set_l(A,B) lookup[A[0]]=B
#define set_sl(A,B) set_l(#A,#B)
            set_sl(A,A);
            set_sl(G,G);
            set_sl(T,T);
            set_sl(C,C);
            set_sl(U,U);
            set_sl(R,GA);
            set_sl(Y,TC);
            set_sl(K,GT);
            set_sl(M,AC);
            set_sl(S,GC);
            set_sl(W,AT);
            set_sl(B,GTC);
            set_sl(D,GAT);
            set_sl(H,ACT);
            set_sl(V,GCA);
            set_sl(N,AGCT);
#undef set_l
#undef set_sl

        }
        template< class RANDGEN>
        string operator [] (RANDGEN &gen) const{
            std::stringstream st;

            for( char c : fmt){
                string buffer;
                assert(lookup[c].size() > 0);
                std::sample(lookup[c].begin(),lookup[c].end(),std::back_inserter(buffer),1,gen);
                st << buffer;
            }
            return st.str();
        }
};

void add_UMIs(vector<molecule_descriptor> &copies, ostream &umifile, const string &format5, const string &format3){

    fmt2seq make_seq5(format5);
    fmt2seq make_seq3(format3);

    int index = 0;
    for(molecule_descriptor &pcp : copies){
        string umi_seq5 = make_seq5[rand_gen];
        string umi_seq3 = make_seq3[rand_gen];
        string umi_ctg_name = fmt::format("RI_umi_ctg_{}",index);

        umifile << fmt::format(">{}\n", umi_ctg_name);
        umifile << umi_seq5 << umi_seq3 << "\n";

        int len5 = static_cast<int>(umi_seq5.size());
        int len3 = static_cast<int>(umi_seq3.size());
        if(len5 > 0){
            pcp.prepend_segment(ginterval{umi_ctg_name, 0, len5, "+"});
        }
        if(len3 > 0){
            pcp.append_segment(ginterval{umi_ctg_name, len5, len5 + len3,  "+"});
        }
        ++index;
    }
}

void add_UMIs(vector<molecule_descriptor> &copies, ostream &umifile, int len){

    uint64_t lim_from_bitc = std::pow(4, len);
    num2seq make_seq(lim_from_bitc);
    std::uniform_int_distribution<uint64_t> dist(0,lim_from_bitc);


    int index = 0;
    for(molecule_descriptor &pcp : copies){
        uint64_t umi_int = dist(rand_gen);
        string umi_seq = make_seq[umi_int];
        string umi_ctg_name = fmt::format("RI_umi_ctg_{}",index);

        umifile << fmt::format(">{}\n", umi_ctg_name);
        umifile << umi_seq << "\n";
        pcp.prepend_segment(ginterval{umi_ctg_name, 0, len, "+"});
        ++index;
    }
}

void describe_program( const cxxopts::ParseResult &args){
    logi("Running UMI tagging module");
    logi("Input MDF: {}", args["input"].as<string>());
    logi("Output MDF: {}", args["output"].as<string>());
    logi("Output UMI FASTA: {}", args["umi-fasta"].as<string>());
    logi("5' UMI tag format: {}", args["format5"].as<string>());
    logi("3' UMI tag format: {}", args["format3"].as<string>());
    logi("Random seed: {}", args["seed"].as<int>());
    fmtlog::poll(true);
}

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser UMI", "UMI module of RNAInfuser");

    options.add_options()
        ("i,input",  "input mdf file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("f,umi-fasta", "Output umi reference file", cxxopts::value<string>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("format5", "UMI sequence format or length", cxxopts::value<string>()->default_value("16"))
        ("format3", "UMI sequence format or length to be attached to 3'", cxxopts::value<string>()->default_value(""))
        ("verbosity", "choose verbosity among [DEBUG, INFO, WARN, ERROR, OFF]", cxxopts::value<string>()->default_value("WARN"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);
    fmtlog::setLogLevel(parse_loglevel(args["verbosity"].as<string>()));
    fmtlog::flushOn(fmtlog::DBG);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"input","output", "umi-fasta"};

    int missing_parameters = 0;
    for( string &param : mandatory){
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }

    if(args.count("format5") + args.count("format3") < 1){
        missing_parameters++;
        std::cerr << "format5 and/or format3 is required!\n";
    }
    if(missing_parameters  > 0){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    describe_program(args);

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);

    string format5 = args["format5"].as<string>();
    string format3 = args["format3"].as<string>();

    string mdf_file_path {args["input"].as<string>()};
    std::ifstream mdf_file {mdf_file_path};
   
    logi("Parsing MDF file {}", mdf_file_path);

    auto streamer = stream_mdf(mdf_file);


    logi("Adding UMI tags");

    if(isdigit(format5[0])){
        int len = std::stoi(format5);
        format5 = std::string(len, 'N');
    }
    if(isdigit(format3[0])){
        int len = std::stoi(format3);
        format3 = std::string(len, 'N');
    }

    fmt2seq make_seq5(format5);
    fmt2seq make_seq3(format3);

    string outfile_name = args["output"].as<string>();
    logi("Adding UMIs and printing to: {}", outfile_name);
    fmtlog::poll(true);
    std::ofstream outfile{outfile_name};
    std::string umi_ref_file = args["umi-fasta"].as<string>();
    std::ofstream umifile{umi_ref_file};

    int index = 0;
    while(streamer){
        molecule_descriptor md = streamer();
        string umi_seq5 = make_seq5[rand_gen];
        string umi_seq3 = make_seq3[rand_gen];
        string umi_ctg_name = fmt::format("RI_umi_ctg_{}",index);

        umifile << fmt::format(">{}\n", umi_ctg_name);
        umifile << umi_seq5 << umi_seq3 << "\n";

        int len5 = static_cast<int>(umi_seq5.size());
        int len3 = static_cast<int>(umi_seq3.size());
        if(len5 > 0){
            md.prepend_segment(ginterval{umi_ctg_name, 0, len5, "+"});
        }
        if(len3 > 0){
            md.append_segment(ginterval{umi_ctg_name, len5, len5 + len3,  "+"});
        }
        ++index;
        outfile << md;
    }

    return 0;
}
