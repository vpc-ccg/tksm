#ifndef _STATS_H
#define _STATS_H
#include <random>
#include <vector>
#include <iostream>




template<class RESULT_TYPE=int>
struct dirichelet_distribution{
    using result_type = RESULT_TYPE;
    using param_type = std::vector<result_type>;
    
    std::gamma_distribution<double> dist;
    std::vector<std::gamma_distribution<double>::param_type> gamma_dists;
    double scale;
 
    dirichelet_distribution (const param_type &R, double scale = 1)  :scale{scale}{
        for (const auto &r: R){
            gamma_dists.emplace_back(r, 1);
        }
    }

    auto operator() (auto &rand_gen){
        std::vector<double> sample;
        for(auto &d : gamma_dists){
            sample.emplace_back(dist(rand_gen,d));
        }
        double sum = 0;
        for(const auto &s : sample){
            sum +=s;
        }
        for( auto &s : sample){
            s = scale * s / sum;
        }
        return sample;
    }

    //TODO implement std random dist requirements
};
#endif
