#ifndef _POLYA_H
#define _POLYA_H
#include "module.h"
#include "util.h"

class PolyA_module : public ri_module{
     //copy of the UMI_module from src/umi.h
    
    private:
    cxxopts::ParseResult parse(int argc, char **argv){
        options.add_options()
            ("i,input", "Input file", cxxopts::value<std::string>())
            ("o,output", "Output file", cxxopts::value<std::string>())
            ("f,polya-reference", "Output polyA reference file", cxxopts::value<std::string>())   
            ("gamma", "Use Gamma distribution [α,β]", cxxopts::value<vector<double>>())
            ("poisson", "Use Poisson distribution [λ]", cxxopts::value<vector<double>>())
            ("weibull", "Use Weibull distribution [α,β]", cxxopts::value<vector<double>>())
            ("normal", "Use Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
            ("min-length", "Minimum length of polyA", cxxopts::value<int>()->default_value("0"))
            ("max-length", "Maximum length of polyA", cxxopts::value<int>()->default_value("5000"))
            ("h,help", "Print help")
        ;
        auto result = options.parse(argc, argv);
        return result;
    }
    cxxopts::ParseResult args;
   
    int validate_arguments(){
        if (args.count("help")){
            std::cout << options.help() << std::endl;
            exit(0);
        }
        if (!args.count("input")){
            std::cerr << "Error: input file not specified" << std::endl;
            exit(1);
        }
        if (!args.count("output")){
            std::cerr << "Error: output file not specified" << std::endl;
            exit(1);
        }
        return 0;
    }
    public:
        PolyA_module( int argc, char **argv ):
            ri_module( "polyA module", "Adds polyA tails to molecules with given size distribution" ),
            args( parse( argc, argv ) )
    {}
      
    int run(){
        std::cout << "Hello world!" << std::endl;
        return 0;
    }
   
    void describe_program(){//Use logi 
        logi("Running polyA module");
        logi("Input file: {}", args["input"].as<std::string>());
        logi("Output file: {}", args["output"].as<std::string>());
        logi("Output polyA reference file: {}", args["polya-reference"].as<std::string>());
        logi("Minimum length of polyA: {}", std::to_string(args["min-length"].as<int>()));
        logi("Maximum length of polyA: {}", std::to_string(args["max-length"].as<int>()));

        logi("Distribution: ");
        if (args.count("gamma")){
            auto vec = args["gamma"].as<vector<double>>();
            logi("Gamma distribution [α,β]: {},{}", vec[0], vec[1]);
        }
        else if (args.count("poisson")){
            auto vec = args["poisson"].as<vector<double>>();
            logi("Poisson distribution [λ]: {}", vec[0]);
        }
        else if (args.count("weibull")){
            auto vec = args["weibull"].as<vector<double>>();
            logi("Weibull distribution [α,β]: {},{}", vec[0], vec[1]);
        }
        else if (args.count("normal")){
            auto vec = args["normal"].as<vector<double>>();
            logi("Normal distribution [μ,σ]: {},{}", vec[0], vec[1]);
        }
        else{
            logi("No distribution specified");
        }

    }
};
#endif
