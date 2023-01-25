
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <variant>

#include <cxxopts/cxxopts.hpp>

#include "interval.h"
#include "mdf.h"

#ifndef POLYA_CONTIG_NAME
#define POLYA_CONTIG_NAME "rf_polyA_contig"
#endif

using std::string;

std::mt19937 rand_gen{std::random_device{}()};

template<class Distribution>
int add_polyA(vector<molecule_descriptor> &copies, Distribution &dist, int min_polya_len = 0){

    int max_size = 0;
    for(molecule_descriptor &pcp : copies){
        int poly_a_len = static_cast<int>(dist(rand_gen));
        if(poly_a_len > max_size){
            max_size = poly_a_len;
        }
        pcp.prepend_segment(ginterval{POLYA_CONTIG_NAME, 0, poly_a_len, "+"});
    }

        
    return max_size;
}

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser polyA", "PolyA module of RNAInfuser");

    options.add_options()
        ("i,input",  "input mdf file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("a,polya-reference", "Output poly A reference file", cxxopts::value<string>())
        ("expand-isoforms", "Expand isoforms to molecules", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("gamma", "Use Gamma distribution [α,β]", cxxopts::value<vector<double>>())
        ("poisson", "Use Poisson distribution [λ]", cxxopts::value<vector<double>>())
        ("weibull", "Use Weibull distribution [α,β]", cxxopts::value<vector<double>>())
        ("normal", "Use Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
        ("min-length", "Minimum length of polyA", cxxopts::value<int>()->default_value("0"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"input","output", "polya-reference"};

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

    std::vector<string> distributions = {"gamma", "poisson", "weibull", "normal"};


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
    if ( chosen_dist != "poisson"){
        if( dist_params.size() != 2){
            std::cerr << chosen_dist << " requires 2 parameters! Make sure there is no space between them.\n";
            return 1;
        }
    }
    else{
        if( dist_params.size() != 1){
            std::cerr << chosen_dist << " requires 1 parameter! Make sure there is no space between them.\n";
            return 1;
        }
    }

    int min_polya_len = args["min-length"].as<int>();
    string mdf_file_path {args["input"].as<string>()};
    std::ifstream mdf_file {mdf_file_path};
    
    vector<molecule_descriptor> molecules = parse_mdf(mdf_file);
    int largest_poly_a = 0; 
    std::variant<
        std::weibull_distribution<>,
        std::poisson_distribution<>,
        std::normal_distribution<>,
        std::gamma_distribution<>
    > dist;
    if( chosen_dist == "poisson"){
        dist = std::poisson_distribution<>{dist_params[0]};
    }
    else if(chosen_dist == "gamma"){
        dist = std::gamma_distribution<> (dist_params[0], dist_params[1]);
    }
    else if(chosen_dist == "normal"){
        dist = std::normal_distribution<> (dist_params[0], dist_params[1]);
    }
    else if(chosen_dist == "weibull"){
        dist = std::weibull_distribution<> (dist_params[0], dist_params[1]);
    }
    else{
        std::cerr << "Distribution " << chosen_dist << " is not implemented!\n";
        return 1;
    }

    if(args["expand-isoforms"].as<bool>()){
        vector<molecule_descriptor> nm;
        for(const molecule_descriptor &pcp : molecules){
            for(int i = 0; i < pcp.get_depth(); ++i){
                nm.push_back(pcp);
                nm.back().depth(1)
                    ->id(nm.back().get_id() + "_" + std::to_string(i));
            }
        }
        molecules = nm;
    }
    std::visit([&molecules, min_polya_len, &largest_poly_a](auto dist){
        largest_poly_a = add_polyA(molecules, dist, min_polya_len);
    }, dist);
    string outfile_name = args["output"].as<string>();
    std::ofstream outfile{outfile_name};
    for(const molecule_descriptor& md : molecules){
        outfile << md;
    }
    
    std::string polya_ref_file = args["polya-reference"].as<string>();
    std::ofstream polya_file{polya_ref_file};

    polya_file << ">" << POLYA_CONTIG_NAME << "\n";
    string polya_tail_string(largest_poly_a, 'A');

    polya_file << polya_tail_string << "\n";
    return 0;
}
