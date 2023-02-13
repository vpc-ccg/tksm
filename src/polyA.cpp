
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
int add_polyA(molecule_descriptor &md, Distribution &dist, int min_polya_len, int max_polya_len){
    int poly_a_len = static_cast<int>(dist(rand_gen));
    if(poly_a_len < min_polya_len){
        poly_a_len = min_polya_len;
    }
    if(poly_a_len > max_polya_len){
        poly_a_len = max_polya_len;
    }
    md.prepend_segment(ginterval{POLYA_CONTIG_NAME, 0, poly_a_len, "+"});
    return poly_a_len;
}

template<class Distribution>
int add_polyA_vec(vector<molecule_descriptor> &copies, Distribution &dist, int min_polya_len = 0){

    int max_size = 0;
    for(molecule_descriptor &md : copies){
        int poly_a_len = add_polyA(md , dist, min_polya_len);
        if(poly_a_len > max_size){
            max_size = poly_a_len;
        }
    }   
    return max_size;
}

int main(int argc, char **argv){

    cxxopts::Options options("tksm polyA", "PolyA module of tksm");

    options.add_options()
        ("i,input",  "input mdf file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("f,polya-reference", "Output poly A reference file", cxxopts::value<string>())
        ("expand-isoforms", "Expand isoforms to molecules", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("gamma", "Use Gamma distribution [α,β]", cxxopts::value<vector<double>>())
        ("poisson", "Use Poisson distribution [λ]", cxxopts::value<vector<double>>())
        ("weibull", "Use Weibull distribution [α,β]", cxxopts::value<vector<double>>())
        ("normal", "Use Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
        ("min-length", "Minimum length of polyA", cxxopts::value<int>()->default_value("0"))
        ("max-length", "Maximum length of polyA", cxxopts::value<int>()->default_value("5000"))
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
    int max_polya_len = args["max-length"].as<int>();
    string mdf_file_path {args["input"].as<string>()};
    std::ifstream mdf_file {mdf_file_path};
   

//    vector<molecule_descriptor> molecules = parse_mdf(mdf_file);
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

    std::string polya_ref_file = args["polya-reference"].as<string>();
    do{
        std::ofstream polya_file{polya_ref_file};
        polya_file << ">" << POLYA_CONTIG_NAME << "\n";
        string polya_tail_string(args["max-length"].as<int>(), 'A');
        polya_file << polya_tail_string << "\n";
    }while(0);
    string outfile_name = args["output"].as<string>();
    std::ofstream outfile{outfile_name};
    auto mdf_stream = stream_mdf(mdf_file, args["expand-isoforms"].as<bool>());
    while(mdf_stream){
        molecule_descriptor md = mdf_stream();
        std::visit([&md, min_polya_len, max_polya_len](auto dist){
            add_polyA(md, dist, min_polya_len, max_polya_len);
        }, dist);
        outfile << md;
    }
    return 0;
}
