
#include <string>
#include <map>
#include <vector>
#include <random>
#include <fstream>
#include <numeric>
#include <functional>
#include <algorithm>

#include "fasta.h"
#include "tree.h"
#include "interval.h"
#include "reverse_complement.h"
#include <cxxopts/cxxopts.hpp>
#include <tuple>
#include "mdf.h"

using std::string;
using std::map;
using std::vector;
using std::tuple;

using std::ostream;
using std::ofstream;
using std::ifstream;


//Random number generator, seed is set in the main function
std::mt19937 rand_gen{std::random_device{}()};

class PCRRunner{
    int     cycles;
    double  efficiency;
    double  error_rate;
    std::uniform_real_distribution<> z1dist  {0,1};
    std::uniform_int_distribution<> basedist {0,3};
    char bases[4] = {'A','C','T','G'};

    public:
    PCRRunner( int cycles, double efficiency, double error_rate) :
        cycles{cycles},
        efficiency{efficiency},
        error_rate{(4*error_rate)/3} //Our error definition is a base to be assigned to random base not changing to other bases so we are adjusting the probability
    {}

    void do_pcr (ostream &ost, const molecule_descriptor& md, int step, const vector<int> &positions, vector<std::pair<int, char>> mutations, double drop_ratio) {
            if( z1dist(rand_gen) > efficiency) { //Molecule is not captured by PCR
                return;
            }
            double expected_mutation_count = error_rate * positions.size();

            vector<int> mutation_pos;
            std::sample(positions.begin(), positions.end(), std::back_inserter(mutation_pos), (int)expected_mutation_count, rand_gen);
            

            for(int pos : mutation_pos){ //TODO: We can make this a normal distribution    
                mutations.push_back(std::make_pair(pos, bases[basedist(rand_gen)]));
            }


            molecule_descriptor mdc{md};
            mdc.update_errors(mutations);
            mdc.id(mdc._id + "." + std::to_string(step));
            if ( z1dist(rand_gen) > drop_ratio){ // Molecule is captured by the sequencing
                ost << mdc;
            }
            for(int cycle = step+1; cycle < cycles; ++cycle){
                do_pcr(ost, mdc, cycle, positions, mutations, drop_ratio);
            }
    } 


    int run(const vector<molecule_descriptor> &molecules, ostream &ost, int number_of_target_reads){
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
int main(int argc, char **argv){

    cxxopts::Options options("tksm PCR module", "PCR PCR PCR PTR");
    
    map<string, tuple<double, double>> pcr_presets {
        { "Taq-setting1", {2*std::pow(0.1,4),0.88}},
        { "Taq-setting2", {7.2*std::pow(0.1,5),0.36}},
        { "Klenow", {1.3*std::pow(0.1,4),0.80}},
        { "T7", {3.4*std::pow(0.1,5),0.90}},
        { "T4", {3.0*std::pow(0.1,6),0.56}},
        { "Vent", {4.5*std::pow(0.1,5),0.70}},
    };

    string preset_string = "presets (Cha, R. S., & Thilly, W. G. (1993). Specificity, efficiency, and fidelity of PCR. Genome Research, 3(3), S18-S29.)\n";
    preset_string.reserve(500);
    for(const auto &p : pcr_presets){
        preset_string+= ("- " +p.first + ": " + "efficiency: " + std::to_string(std::get<1>(p.second)) +", error-rate:" + std::to_string(std::get<0>(p.second)) + "\n");
    }

    options.add_options()
        ("i,input", "Molecule description file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("molecule-count", "Number of molecules to simulate", cxxopts::value<size_t>()->default_value("100000"))
        ("cycles", "Number of pcr cycles to simulate", cxxopts::value<int>()->default_value("6"))
        ("efficiency", "Probability of a molecule being duplicated during each pcr cycle", cxxopts::value<double>()->default_value("0.75"))
        ("error-rate", "Probability of substition errors for each base during each pcr cycle", cxxopts::value<double>()->default_value("0.000001"))
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
        ("x,preset", preset_string, cxxopts::value<string>())
    ;

    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }

    std::vector<string> mandatory =  {"input", "output"};
    if(args.count("fastq") > 0 && args["fastq"].as<bool>()){
        mandatory.push_back("references");
    }
    int missing_parameters = 0;
    for( string &param : mandatory){
        std::cout << param << "\n";
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }

    if(missing_parameters  > 0){
        std::cerr << options.help() << std::endl;
        return 1;
    }
    
    int cycles = args["cycles"].as<int>();
    double pcr_efficiency = args["efficiency"].as<double>();
    double error_rate = args["error-rate"].as<double>();
    int number_of_target_reads = args["molecule-count"].as<size_t>();

    if(args["preset"].count() > 0){
        if (pcr_presets.find(args["preset"].as<string>()) == pcr_presets.end()){
            std::cerr << "Preset doesn't exist!\n";
            return -1;
        }
        tuple<double, double> setting = pcr_presets[args["preset"].as<string>()];

        error_rate = std::get<0>(setting);
        pcr_efficiency = std::get<1>(setting);
    }

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);

    ifstream md_file(args["input"].as<string>());
    vector<molecule_descriptor> molecules = parse_mdf(md_file);

    if( molecules.size() > 2 * args["molecule-count"].as<size_t>()){
        std::shuffle(molecules.begin(), molecules.end(), rand_gen);
        molecules.resize( 2 * args["molecule-count"].as<size_t>());
    }


    ofstream out_stream(args["output"].as<string>());
    return PCRRunner{cycles, pcr_efficiency, error_rate}.run(molecules, out_stream, number_of_target_reads);
}


