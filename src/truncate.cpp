#include <cstddef>
#include <exception>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <fstream>
#include <set>
#include <sstream>
#include <climits>

#include "extern/IITree.h"
#include "extern/cxxopts.hpp"
#include "interval.h"
#include "tree.h"
#include "gtf.h"
#include "mdf.h"

using std::ofstream;
using std::vector;
using std::string;

std::mt19937 rand_gen{std::random_device{}()};

template <class Distribution>
void truncate(pcr_copy &pcp, Distribution &dist){
    int rand_val = dist(rand_gen);

    size_t i = 0;
    int len_so_far = 0;
    for( const ginterval &g : pcp.segments){
        len_so_far += (g.end - g.start);
        if(len_so_far >= rand_val){
            break;
        }
        ++i;
    }
    if( i != pcp.segments.size()){
        std::stringstream ss;
        ss << pcp.segments[i].chr << ':' << pcp.segments[i].end - (len_so_far - rand_val) << '-' << pcp.segments[i].end << ',';
        pcp.segments[i].end -= (len_so_far - rand_val); //Update last segment;
        for(size_t j = i + 1;j<pcp.segments.size();++j){
            ss << pcp.segments[j] << ',';
        }
        string trunc_str = ss.str();
        trunc_str.back() = ';';
        pcp.comment += "truncated:" + trunc_str;
        pcp.segments.resize(i+1);
    }

    i = 0;
    for( const auto &p : pcp.errors_so_far){
        if(p.first > rand_val){
            break;
        }
        ++i;
    }
    pcp.errors_so_far.resize(i);
}

int main(int argc, char **argv){
    cxxopts::Options options("RNAInfuser Splicer", "Splicer module of RNAInfuser");

    options.add_options()
        ("i,input", "Molecule description file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("normal", "Use Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
        ("lognormal", "Use Log-Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
        ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    vector<string> mandatory = {"input", "output"};

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
    vector<string> distributions = {"lognormal", "normal"};


    string chosen_dist;
    int arg_dist_count = 0;
    for( const string &dist_name : distributions){
        if( args.count(dist_name) > 0 ){
            ++arg_dist_count;
            chosen_dist = dist_name;
        }
    }

    if(arg_dist_count != 1){
        std::cerr << "Please provide one distribution type!\n";
        return 1;
    }
    vector<double> dist_params = args[chosen_dist].as<vector<double>>();

    if( dist_params.size() != 2){
        std::cerr << chosen_dist << " requires 2 parameters! Make sure there is no space between them.\n";
        return 1;
    }

    string mdf_file_path {args["input"].as<string>()};

    std::ifstream mdf_file {mdf_file_path};

    std::normal_distribution<> normal_dist(dist_params[0], dist_params[1]);
    std::lognormal_distribution<> lognormal_dist(dist_params[0], dist_params[1]);

    vector<pcr_copy> molecules = parse_mdf(mdf_file);
    for( pcr_copy &pcp : molecules){
        if(chosen_dist == "normal"){
            truncate(pcp, normal_dist);
        }
        else if(chosen_dist == "lognormal"){
            truncate(pcp, lognormal_dist);
        }
    }

    string outfile_name = args["output"].as<string>();
    std::ofstream outfile{outfile_name};
    print_all_mdf(outfile, molecules);
    
    return 0;
}
