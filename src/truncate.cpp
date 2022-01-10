#include <cstddef>
#include <exception>
#include <numeric>
#include <stdexcept>
#include <string>
#include <random>
#include <tuple>
#include <vector>
#include <map>
#include <fstream>
#include <set>
#include <sstream>
#include <climits>
#include <algorithm>
#include <numbers>
#include <functional>

#include "extern/IITree.h"
#include "extern/cxxopts.hpp"
#include "interval.h"
#include "tree.h"
#include "gtf.h"
#include "mdf.h"
#include "paf.h"

using std::ofstream;
using std::vector;
using std::string;

std::mt19937 rand_gen{std::random_device{}()};

template <class Distribution, class MeanFoo>
void multi_truncate(pcr_copy &pcp, MeanFoo mu_func, double sigma, Distribution &dist){
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


struct identity{
    constexpr auto operator () (auto val) const {
        return val;
    }
};

template<class Accessor>
auto mean_std(const vector<int> &lengths, Accessor accessor){
    double mean = 0;
    int i = 0;

    for(int val : lengths){
        mean = mean + (accessor(val) - mean) / (i+1);
        ++i;
    }
    double varsum = 0;
    for(int val : lengths){
        varsum = varsum + std::pow( accessor(val)- mean,2) ;
    }
    
    double stdev = std::sqrt(varsum / (i-1));
    
    return std::make_tuple(mean, stdev);
}

auto mean_std(const vector<int> &lengths){
    return mean_std(lengths, identity{});
}

auto fit(const vector<int> &lengths){

    using namespace std::numbers;

    std::map<string, std::function<std::tuple<double, double, double> (const vector<int> &)>> ll_computers = {
        { 
            "normal" , [] (const vector<int> &ll){
                auto [_mean, std] = mean_std(ll);
                double mean = _mean; //This is needed because clangd complains 
                size_t n = ll.size();
                double var_cal = std::accumulate(ll.begin(), ll.end(), 0.0L, [mean] (double sum, int val) -> double {
                    return sum + (val - mean) * (val - mean); 
                        });
                double LL = -(n/2.0)*std::log(2 * pi) -(n/2.0) *std::log(std*std) - (1.0/(2*std*std)) * var_cal;
                return std::make_tuple(LL, mean, std);
            }
        },
        { 
            "lognormal" , [] (const vector<int> &ll){
                auto [_mean, std] = mean_std(ll, [](int val){ return std::log(val);});
                double mean = _mean; //This is needed because clangd complains 
                size_t n = ll.size();
                double var_cal = std::accumulate(ll.begin(), ll.end(), 0.0L, [mean] (double sum, int val) -> double {
                    return sum + (std::log(val) - mean) * (std::log(val) - mean); 
                        });
                double LL = -(n/2.0)*std::log(2 * pi) -(n/2.0) *std::log(std*std) - (1.0/(2*std*std)) * var_cal;
                return std::make_tuple(LL, mean, std);
            }
        }
    };

    string max_dist = "None";
    double max_ll = std::numeric_limits<double>::lowest();
    double chosen_mean = -1;
    double chosen_std  = -1;
    for( auto &p : ll_computers){
        auto [LL, mean, std] = p.second(lengths);
        std::cerr << p.first << ": ";
        std::cerr << "ll: " << LL << " "; 
        std::cerr << "mean: " << mean << " "; 
        std::cerr << "std: " << std << "\n";
        if(LL > max_ll){
            max_dist = p.first;
            max_ll   = LL;
            chosen_mean = mean;
            chosen_std = std;
        }
    }
    return make_tuple(max_dist, max_ll, chosen_mean, chosen_std);
}


auto multi_variate(const vector<int> &rlens, const vector<int> &tlens){
    auto [dist, ll, mean, std] = fit(rlens);
    double mu1 = mean;
    double sig1 = std;

    auto [tmean, tstd] = mean_std(tlens);
    double mu2 = tmean;
    double sig2 = tstd;

    auto tlen_iter = tlens.begin();


    double xy_sum = 0;
    double x_sum = 0;
    double y_sum = 0;
    for( auto rlen_iter = rlens.begin(); rlen_iter != rlens.end(); ++rlen_iter, ++tlen_iter){
        double xdif = (*rlen_iter - mean);
        double ydif = (*tlen_iter - tmean);
        xy_sum += xdif * ydif;
        x_sum  += xdif * xdif;
        y_sum  += ydif * ydif;
    }

    double ro = xy_sum / (std::sqrt(x_sum) * std::sqrt(y_sum));

    auto mean_adjuster = [ro, mu1, mu2, sig1, sig2] (double alpha){
        return mu1 + (sig1/sig2) * ro * (alpha - mu2);
    };
    double adjusted_sig = std::sqrt(( 1 - (ro*ro)) * sig1*sig1);
    return std::make_tuple(dist, mean_adjuster, adjusted_sig);
}

int main(int argc, char **argv){
    cxxopts::Options options("RNAInfuser Splicer", "Splicer module of RNAInfuser");

    options.add_options()
        ("i,input", "Molecule description file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("fit", "Use estimated distribution from the given fast{a,q}", cxxopts::value<string>())
        ("multivariate", "Use estimated distribution from the given paf", cxxopts::value<string>())
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
    vector<string> distributions = {"lognormal", "normal", "fit", "multivariate"};


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
    vector<double> dist_params;
    if( chosen_dist == "multivariate"){
        string dist_fit_file_path {args["multivariate"].as<string>()};
        std::ifstream dist_fit_file(dist_fit_file_path);

        string buffer;
        vector<int> lens;
        vector<int> tlens;

        
        while(std::getline(dist_fit_file, buffer)){ 
            paf pf{buffer};
            if(!pf.primary()){
                continue;
            }
            lens.push_back( pf.qend - pf.qstart);
            tlens.push_back( pf.tlen);
        }

        auto [distname, mu_func, sigma] = multi_variate(lens, tlens);

        string mdf_file_path {args["input"].as<string>()};

        std::ifstream mdf_file {mdf_file_path};

        vector<pcr_copy> molecules = parse_mdf(mdf_file);
        for( pcr_copy &pcp : molecules){

            std::normal_distribution<> normal_dist(mu_func(pcp.size()), sigma);
            std::lognormal_distribution<> lognormal_dist(mu_func(pcp.size()), sigma);
            if( distname == "normal"){
                truncate(pcp, normal_dist);
            }
            else if( distname == "lognormal"){
                truncate(pcp, lognormal_dist);
            }
            else{
                std::cerr << distname << " not implemented!\n";
                return -1;
            }
        }

        string outfile_name = args["output"].as<string>();
        std::ofstream outfile{outfile_name};
        print_all_mdf(outfile, molecules);
        return 0;   
    }
    else if( chosen_dist == "fit"){
        string dist_fit_file_path {args["fit"].as<string>()};
        std::ifstream dist_fit_file(dist_fit_file_path);

        string buffer;
        vector<int> lens;
        while(std::getline(dist_fit_file, buffer)){ //Subsamples fasta files 1/2 ratio. TODO: use a proper fast{a,q} reader
            std::getline(dist_fit_file, buffer);
            lens.push_back(buffer.size() - 1);
            std::getline(dist_fit_file, buffer);
            std::getline(dist_fit_file, buffer);
        }
        auto [name, ll, mean, std] = fit(lens);
        std::cerr << name << " is chosen with " << ll << " log-likelihood\n";
        chosen_dist = name;
        dist_params = {mean, std};

    }
    else{
        dist_params = args[chosen_dist].as<vector<double>>();
        if( dist_params.size() != 2){
            std::cerr << chosen_dist << " requires 2 parameters! Make sure there is no space between them.\n";
            return 1;
        }

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
