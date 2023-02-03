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

class fmt2seq{
    public:
        std::array<char,4> table  {{'A','T','C','G'}};
        std::array<string, 128> lookup;
        string fmt;
        fmt2seq(const string &fmt):fmt(fmt){
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
                std::sample(lookup[c].begin(),lookup[c].end(),std::back_inserter(buffer),1,gen);
                st << buffer;
            }
            return st.str();
        }
};

void add_UMIs(vector<molecule_descriptor> &copies, ostream &umifile, const string &format, const string &format_back){

    fmt2seq make_seq(format);
    fmt2seq make_back_seq(format_back);

    int index = 0;
    for(molecule_descriptor &pcp : copies){
        string umi_seq = make_seq[rand_gen];
        string umi_ctg_name = fmt::format("RI_umi_ctg_{}",index);

        string umi_back_seq = make_back_seq[rand_gen];

        umifile << fmt::format(">{}\n", umi_ctg_name);
        umifile << umi_seq << umi_back_seq << "\n";

        int len = static_cast<int>(umi_seq.size());
        int len_back = static_cast<int>(umi_back_seq.size());
        if(len > 0){
            pcp.prepend_segment(ginterval{umi_ctg_name, 0, len, "+"});
        }
        if(len_back > 0){
            pcp.append_segment(ginterval{umi_ctg_name, len, len + len_back,  "+"});
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

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser UMI", "UMI module of RNAInfuser");

    options.add_options()
        ("i,input",  "input mdf file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("f,umi-fasta", "Output umi reference file", cxxopts::value<string>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("format5", "UMI sequence format or length", cxxopts::value<string>()->default_value("16"))
        ("format3", "UMI sequence format or length to be attached to 3'", cxxopts::value<string>()->default_value(""))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

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

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);

    string umi_format = args["format5"].as<string>();

    string mdf_file_path {args["input"].as<string>()};
    std::ifstream mdf_file {mdf_file_path};
    
    vector<molecule_descriptor> molecules = parse_mdf(mdf_file);


//Unroll the depth
    vector<molecule_descriptor> nm;
    for(const molecule_descriptor &pcp : molecules){
        for(int i = 0; i < pcp.get_depth(); ++i){
            nm.push_back(pcp);
            nm.back().depth(1);
            nm.back().id(nm.back().get_id() +"_" + std::to_string(i));
        }
    }
    molecules = nm;

    std::string umi_ref_file = args["umi-fasta"].as<string>();
    std::ofstream umi_file{umi_ref_file};

    if(std::isdigit(umi_format[0])){
        int length = stoi(umi_format);
        add_UMIs(molecules, umi_file, length);
    }
    else{
        add_UMIs(molecules,umi_file, umi_format, args["format3"].as<string>());
    }

    string outfile_name = args["output"].as<string>();
    std::ofstream outfile{outfile_name};
    for(const molecule_descriptor &md : molecules){
        outfile << md;
    }
    
    return 0;
}
