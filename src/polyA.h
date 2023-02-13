#ifndef _POLYA_H
#define _POLYA_H
#include <variant>


#include "module.h"
#include "util.h"

#ifndef POLYA_REFERENCE_NAME
#define POLYA_REFERENCE_NAME "tksm_polyA_reference"
#endif

class PolyA_module : public tksm_module{
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
        ;
        auto result = options.parse(argc, argv);
        return result;
    }
    cxxopts::ParseResult args;
    std::mt19937 rand_gen;
   
    int validate_arguments(){
        if (args.count("help")){
            std::cout << options.help() << std::endl;
            exit(0);
        }
        if (!args.count("input")){
            loge("Error: input file not specified");
            exit(1);
        }
        if (!args.count("output")){
            loge("Error: output file not specified");
            exit(1);
        }
        if (!args.count("polya-reference")){
            loge("Error: output polyA reference file not specified");
            exit(1);
        }
        if (!args.count("gamma") && !args.count("poisson") && !args.count("weibull") && !args.count("normal")){
            loge("Error: no distribution specified");
            exit(1);
        }
        if (args.count("gamma") + args.count("poisson") + args.count("weibull") + args.count("normal") > 1){
            loge("Error: multiple distributions specified");
            exit(1);
        }
        if (args.count("gamma")){
            auto vec = args["gamma"].as<vector<double>>();
            if (vec.size() != 2){
                loge("Error: gamma distribution requires two parameters");
                exit(1);
            }
        }
        if (args.count("poisson")){
            auto vec = args["poisson"].as<vector<double>>();
            if (vec.size() != 1){
                loge("Error: poisson distribution requires one parameter");
                exit(1);
            }
        }
        if (args.count("weibull")){
            auto vec = args["weibull"].as<vector<double>>();
            if (vec.size() != 2){
                loge("Error: weibull distribution requires two parameters");
                exit(1);
            }
        }
        if (args.count("normal")){
            auto vec = args["normal"].as<vector<double>>();
            if (vec.size() != 2){
                loge("Error: normal distribution requires two parameters");
                exit(1);
            }
        }
        if (args["min-length"].as<int>() < 0){
            loge("Error: minimum length of polyA cannot be negative");
            exit(1);
        }
        if (args["max-length"].as<int>() < 0){
            loge("Error: maximum length of polyA cannot be negative");
            exit(1);
        }
        if (args["min-length"].as<int>() > args["max-length"].as<int>()){
            loge("Error: minimum length of polyA cannot be greater than maximum length of polyA");
            exit(1);
        }
            
        return 0;
    }
    public:
        PolyA_module( int argc, char **argv ):
            tksm_module( "polyA module", "Adds polyA tails to molecules with given size distribution" ),
            args( parse( argc, argv ) )
    {}

    template<class Distribution>
    int add_polyA(molecule_descriptor &md, Distribution &dist, int min_polya_len, int max_polya_len){
        int poly_a_len = static_cast<int>(dist(rand_gen));
        if(poly_a_len < min_polya_len){
            poly_a_len = min_polya_len;
        }
        if(poly_a_len > max_polya_len){
            poly_a_len = max_polya_len;
        }
        md.prepend_segment(ginterval{POLYA_REFERENCE_NAME, 0, poly_a_len, "+"});
        return poly_a_len;
    }

          
        int run(){
        fmtlog::setLogLevel(LogLevels::parse_loglevel(args["verbosity"].as<string>()));
        fmtlog::setFlushDelay(1e8);

        if(help_or_version_is_used(args)){
            return 0;
        }

        if(validate_arguments()){
            return 1;
        }
        describe_program();

        int seed = args["seed"].as<int>();;
        rand_gen.seed(seed);

        int min_length = args["min-length"].as<int>();
        int max_length = args["max-length"].as<int>();

        std::variant<std::gamma_distribution<double>, std::poisson_distribution<int>, std::weibull_distribution<double>, std::normal_distribution<double>> dist;

        if (args.count("gamma")){
            auto vec = args["gamma"].as<vector<double>>();
            dist = std::gamma_distribution<double>(vec[0], vec[1]);
        }
        if (args.count("poisson")){
            auto vec = args["poisson"].as<vector<double>>();
            dist = std::poisson_distribution<int>(vec[0]);
        }
        if (args.count("weibull")){
            auto vec = args["weibull"].as<vector<double>>();
            dist = std::weibull_distribution<double>(vec[0], vec[1]);
        }
        if (args.count("normal")){
            auto vec = args["normal"].as<vector<double>>();
            dist = std::normal_distribution<double>(vec[0], vec[1]);
        }

        std::string input_file = args["input"].as<std::string>();
        std::string output_file = args["output"].as<std::string>();

        std::string output_polya_reference_file = args["polya-reference"].as<std::string>();

        do{
            std::ofstream output_polya_reference(output_polya_reference_file);
            if (!output_polya_reference.is_open()){
                loge("Error: cannot open output polyA reference file");
                break;
            }

            output_polya_reference << ">" << POLYA_REFERENCE_NAME << "\n";
            string polya_tail_string (max_length, 'A');
            output_polya_reference << polya_tail_string << "\n";

        }while(0);

        std::ifstream input(input_file);
        if (!input.is_open()){
            loge("Error: cannot open input file: {}" , input_file);
            return 1;
        }

        std::ofstream output(output_file);
        if (!output.is_open()){
            loge("Error: cannot open output file: {}" , output_file);
            return 1;
        }

        auto mdf_stream = stream_mdf(input, true);

        while(mdf_stream){
            molecule_descriptor md = mdf_stream();

            std::visit( [&](auto dist){
                add_polyA(md, dist, min_length, max_length);
                    }, dist);

            output << md;
        }
        
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
        fmtlog::poll(true);
    }
};
#endif
