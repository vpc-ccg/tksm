#include <exception>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <fstream>
#include <set>

#include <climits>

#include <cgranges/IITree.h>
#include <cxxopts/cxxopts.hpp>

#include "interval.h"
#include "tree.h"
#include "gtf.h"
#include "mdf.h"

#include "umi.h"
#include <cstring>
#include <utility>

using std::pair;
using std::ofstream;
using std::vector;
using std::string;
using std::map;

std::mt19937 rand_gen{std::random_device{}()};


vector<string> read_vec(ifstream &file, size_t count){

    vector<string> vec;
    string buffer;

    for(size_t i = 0; i < count; ++i){
        std::getline(file, buffer);
        vec.push_back(buffer);
    }
    return vec;
}

class SparseMatrix{
    map<pair<int,int>, int> mtx;
    vector<string> features;
    vector<string> barcodes;
    public:
    SparseMatrix(const string &matrix_path)
//features{read_vec(feature_path)}, barcodes{read_vec(barcode_path)}
    {
        ifstream file(matrix_path);
        string buffer;
        std::getline(file, buffer);
        int feat_count, barcode_count, entry_count;
        string comment;
        std::stringstream(buffer) >> feat_count >> barcode_count >> entry_count >> comment;

        features = read_vec(file, feat_count);
        barcodes = read_vec(file, barcode_count);
        while(std::getline(file, buffer)){
            int row_id, col_id, count;
            sscanf(buffer.c_str(), "%d\t%d\t%d", &row_id, &col_id, &count);
            mtx[std::make_pair(row_id, col_id)] = count;
        }
    }

    friend ostream& operator<< (ostream &ost, const SparseMatrix &sparse){

        for(const string &bc: sparse.barcodes){
            ost << "\t" << bc ;
        }
        ost << "\n";
        for(size_t feat_index = 0; feat_index < sparse.features.size(); ++feat_index){
            ost << sparse.features[feat_index] << "\t";
            for(size_t bar_index = 0; bar_index < sparse.barcodes.size(); ++bar_index){
                auto ptr = sparse.mtx.find(std::make_pair(feat_index, bar_index));
                if(ptr != sparse.mtx.end()){
                    ost << ptr->second << "\t";
                }
                else{
                    ost << 0 << "\t";
                }
            }
            ost << "\n";
        }
        return ost;
    }
};


int single_cell_splice(int argc, char **argv){

    cxxopts::Options options("RNAInfuser Splicer", "Splicer module of RNAInfuser");

    options.add_options()
        ("g,gtf",  "Path to gtf annotation file", cxxopts::value<string>())
        ("m,count-matrix",  "Barcode to transcript matrix of counts", cxxopts::value<string>())
        ("use-whole-id", "Use whole transcript id instead of first 15 characters (ENSEMBL ids)", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("molecule-count", "Number of molecules to simulate, (requires abundance table", cxxopts::value<int>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"gtf", "output", "count-matrix"};

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

    //parameters will be moved to argument parser when done
    string path_to_gtf   {args["gtf"].as<string>()};
    string out_path {args["output"].as<string>()};
    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);
   


//    add_UMIs(std::vector<molecule_descriptor> &copies, ostream &umifile, const std::string &format);

    return 0;
}

int run_splice(const map<string, int> &counts, const vector<molecule_descriptor> &isoforms){

    for(const molecule_descriptor &cpy : isoforms){
        auto fptr = counts.find(cpy.get_id());
        if(fptr == counts.end()){
            continue;
        }

    }
    return 0;
}

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser Splicer", "Splicer module of RNAInfuser");

    options.add_options()
        ("g,gtf",  "Path to gtf annotation file", cxxopts::value<string>())
        ("a,abundance-table",  "Path to tab separated abundance table (Formatted as transcript_id\\tcount\\tpm)", cxxopts::value<string>())
        ("use-whole-id", "Use whole transcript id instead of first 15 characters (ENSEMBL ids)", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("molecule-count", "Number of molecules to simulate, (requires abundance table", cxxopts::value<int>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("non-coding", "Process non-coding genes/transcripts as well", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("default-depth", "Default depth for transcripts that are not in expression table", cxxopts::value<int>()->default_value("0"))
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"gtf", "output", "abundance-table", "molecule-count"};
    
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

    //parameters will be moved to argument parser when done
    string path_to_gtf   {args["gtf"].as<string>()};
    string out_path {args["output"].as<string>()};

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);
   

    auto isoforms = read_gtf_transcripts(path_to_gtf, args["default-depth"].as<int>());


    ifstream table_file(args["abundance-table"].as<string>());
    ofstream outfile( out_path);
    if(!table_file){
        std::cerr << "Cannot open " << args["abundance-table"].as<string>() << ". Terminating.\n";
        return 1;
    }
    double molecule_count = args["molecule-count"].as<int>();

    string tid = "BEG";
    double tpm;
    string comment;



    std::uniform_real_distribution<> dist(0, 1);

    std::string buffer;
    std::getline(table_file, buffer); //Read header

    while(std::getline(table_file, buffer)){
        std::istringstream(buffer) >> tid >> tpm >> comment;
        if(!args["use-whole-id"].as<bool>()){
            tid = tid.substr(0,15);
        }
        molecule_descriptor molecule = isoforms[tid];
        comment.append(";");
        double count = tpm*molecule_count/1'000'000;
        double carry = count - int(count);
        if( dist(rand_gen) > carry){
            count++;
        }
        molecule.comment(comment)->depth(count); //rounded down to int

        outfile << molecule << "\n";
    }
}


