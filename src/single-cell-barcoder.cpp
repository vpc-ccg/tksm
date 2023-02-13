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
/*
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
*/
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
