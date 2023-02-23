#include "pcr.h"
#include "pimpl.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::ifstream;
using std::map;
using std::ofstream;
using std::pair;
using std::string;
using std::vector;




class PCR{
    int     cycles;
    double  efficiency;
    double  error_rate;
    std::uniform_real_distribution<> z1dist  {0,1};
    std::uniform_int_distribution<> basedist {0,3};
    char bases[4] = {'A','C','T','G'};

    std::mt19937 &rand_gen;
    public:
    PCR( int cycles, double efficiency, double error_rate, std::mt19937 &rand_gen) :
        cycles{cycles},
        efficiency{efficiency},
        error_rate{(4*error_rate)/3}, //Our error definition is a base to be assigned to random base not changing to other bases so we are adjusting the probability
        rand_gen{rand_gen}
    {}

    void do_pcr (ostream &ost, const molecule_descriptor& md, int step, const vector<int> &positions, vector<std::pair<int, char>> mutations, double drop_ratio) {
            if( z1dist(rand_gen) > efficiency) { //Molecule is not captured by PCR
                return;
            }
            double expected_mutation_count = error_rate * positions.size();

            vector<int> mutation_pos;
            std::sample(positions.begin(), positions.end(), std::back_inserter(mutation_pos), (int)expected_mutation_count, rand_gen);
            

            molecule_descriptor mdc{md};
            for(int pos : mutation_pos){ //TODO: We can make this a normal distribution    
                mutations.push_back(std::make_pair(pos, bases[basedist(rand_gen)]));
                mdc.add_error(mutations.back());
            }

            mdc.id(mdc._id + "." + std::to_string(step));
            if ( z1dist(rand_gen) > drop_ratio){ // Molecule is captured by the sequencing
                ost << mdc;
            }
            for(int cycle = step+1; cycle < cycles; ++cycle){
                do_pcr(ost, mdc, cycle, positions, mutations, drop_ratio);
            }
    } 

    int perform(const vector<molecule_descriptor> &molecules, int number_of_target_reads, ostream &ost){
//Calculate number of molecules to be produced with pcr
        int64_t molecule_count = std::accumulate(molecules.begin(), molecules.end(), 0,
                [] ( int64_t sum, const molecule_descriptor &cpy) -> int64_t{
                return sum + cpy.get_depth();
        });

//Calculate expected number of molecules after pcr
        int64_t expected_number_after_pcr = std::pow((1 + efficiency), cycles) * molecule_count;

        double drop_ratio = static_cast<double>(number_of_target_reads) / expected_number_after_pcr;

        for( const molecule_descriptor &pcp : molecules){
            string base_id = pcp.get_id();
            vector<int> positions(pcp.size());
            std::iota(positions.begin(), positions.end(), 0);

            for(int cycle = 0; cycle < cycles; ++cycle){
                do_pcr(ost, pcp, cycle, positions, vector<std::pair<int, char>>{}, drop_ratio);
            }
        }
        return 0;
    }
};

class PCR_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {

        string preset_string = "presets (Cha, R. S., & Thilly, W. G. (1993). Specificity, efficiency, and fidelity of PCR. Genome Research, 3(3), S18-S29.)\n";
        for(const auto &p: pcr_presets){
            preset_string += ("- " + p.first + ": " + "error-rate: " + std::to_string(p.second.first) + ", efficiency-rate: " + std::to_string(p.second.second) + "\n");
        }

        // clang-format off
        options.add_options("main")
            (
                "i,input",
                "input mdf file",
                cxxopts::value<string>()
            )(
                "o,output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "molecule-count",
                "Number of molecules to generate",
                cxxopts::value<size_t>()
            )(
                "cycles",
                "Number of PCR cycles to perform",
                cxxopts::value<int>()
            )(
                "error-rate",
                "Error rate for PCR",
                cxxopts::value<double>()
            )(
                "efficiency",
                "PCR efficiency",
                cxxopts::value<double>()
            )(
                "x,preset",
                preset_string,
                cxxopts::value<string>()
            )
        ;
        // clang-format on
        return options.parse(argc, argv);
    }

    map<string, pair<double, double>> pcr_presets{
        {"Taq-setting1", {2 * std::pow(0.1, 4), 0.88}}, {"Taq-setting2", {7.2 * std::pow(0.1, 5), 0.36}},
        {"Klenow", {1.3 * std::pow(0.1, 4), 0.80}},     {"T7", {3.4 * std::pow(0.1, 5), 0.90}},
        {"T4", {3.0 * std::pow(0.1, 6), 0.56}},         {"Vent", {4.5 * std::pow(0.1, 5), 0.70}},
    };

    cxxopts::ParseResult args;


public:
    impl (int argc, char **argv) : tksm_module{"PCR", "PCR amplification module"}, args(parse(argc, argv)) {}

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "molecule-count", "cycles"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                loge("{} is required!", param);
                ++missing_parameters;
            }
        }
        // Other parameter checks here



        if (args["preset"].count() > 0) {
            string preset = args["preset"].as<string>();
            auto it       = pcr_presets.find(preset);
            if (it == pcr_presets.end()) {
                loge("Preset {} not found", preset);
                ++missing_parameters;
            }
        }
        else{
            if (args["error-rate"].count() == 0) {
                loge("Error rate is required!");
                ++missing_parameters;
            }
            if (args["efficiency"].count() == 0) {
                loge("Efficiency is required!");
                ++missing_parameters;
            }
        }

        if (missing_parameters > 0) {
            fmt::print(stderr, "{}\n", options.help());
            return 1;
        }
        return 0;
    }
    int run() {
        if (process_utility_arguments(args)) {
            return 0;
        }
        if (validate_arguments()) {
            return 1;
        }
        describe_program();
        string input_file = args["input"].as<string>();
        string output_file = args["output"].as<string>();
        
        double error_rate = 0;
        double efficiency = 0;
        size_t molecule_count = args["molecule-count"].as<size_t>();
        int cycles = args["cycles"].as<int>();

        if( args["preset"].count()){
            std::pair<double, double> preset = pcr_presets[args["preset"].as<string>()];
            error_rate = preset.first;
            efficiency = preset.second;
        }
        if( args["error-rate"].count()){
            error_rate = args["error-rate"].as<double>();
        }
        if( args["efficiency"].count()){
            efficiency = args["efficiency"].as<double>();
        }

        ifstream input(input_file);

        vector<molecule_descriptor> molecules = parse_mdf(input);

        if( molecules.size() > 2 * molecule_count){
            std::shuffle(molecules.begin(), molecules.end(), rand_gen);
            molecules.resize( 2 * molecule_count);
        }
        ofstream output (output_file);
        return PCR{cycles, efficiency, error_rate, tksm_module::rand_gen}.perform(molecules, molecule_count, output);
        return 0;
    }

    void describe_program() {
        logi("Running PCR");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        logi("Molecule count: {}", args["molecule-count"].as<size_t>());
        logi("Cycles: {}", args["cycles"].as<int>());
        if (args.count("preset")) {
            logi("Used preset: {}", args["preset"].as<string>());
            std::pair<double, double> preset = pcr_presets[args["preset"].as<string>()];
            
            if(args.count("error-rate")){
                logi("Error rate is overriden from {} to {}", preset.first, args["error-rate"].as<double>());
            }
            else{
                logi("Error rate: {}", preset.first);
            }
            if(args.count("efficiency")){
                logi("Efficiency is overridden from {} to {}", preset.second, args["efficiency"].as<double>());
            }
            else{
                logi("Efficiency: {}", preset.second);
            }
        }
        else{
            if(args.count("error-rate")){
                logi("Error rate: {}", args["error-rate"].as<double>());
            }
            if(args.count("efficiency")){
                logi("Efficiency: {}", args["efficiency"].as<double>());
            }
        }
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(PCR_module);

