#include "cut.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"
#include "stats.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

#include "pimpl.h"
inline std::pair<molecule_descriptor, molecule_descriptor>
split_molecule(const molecule_descriptor &md, int split_pos) {
    int bases_so_far = 0;
    molecule_descriptor md1;
    molecule_descriptor md2;
    for (const einterval &g : md.cget_segments()) {
        if (g.size() + bases_so_far < split_pos) {
            md1.append_segment(g);
        }
        else {
            if (bases_so_far < split_pos) {  // Transition interval
                einterval cpy = g;

                einterval cpy2 = g;
                if (g.plus_strand) {
                    cpy.truncate(0, split_pos - bases_so_far);
                    cpy2.truncate(split_pos - bases_so_far, cpy2.size());
                }
                else {
                    cpy.truncate(split_pos - bases_so_far, cpy.size());
                    cpy2.truncate(0, split_pos - bases_so_far);
                }
                md1.append_segment(cpy);
                md2.append_segment(cpy2);
            }
            else {
                md2.append_segment(g);
            }
        }
        bases_so_far += g.size();
    }
    return {md1, md2};
}
class Cut_module::impl : public tksm_module {

    static vector<molecule_descriptor> split_molecule_n(const molecule_descriptor &md, const vector<double> &cuts){
        vector<molecule_descriptor> result;
        std::pair<molecule_descriptor, molecule_descriptor> mdp = {{}, md};
        int index = 0;
        for( double d : cuts){
            mdp = split_molecule(mdp.second, d);
            if(mdp.first.size() == 0 ){
                logd("Molecule of size zero generated, skipping");
                continue;
            }
            result.push_back(mdp.first);
            result.back().id(fmt::format("{}_C{}",md.get_id(),  index++))->depth(1);
        }
        return result;
    }
    static vector<molecule_descriptor> cut_molecule(const molecule_descriptor &md, double probability, auto &randgen) {
        std::binomial_distribution<size_t> number_of_cuts_distribution{md.size(), probability};

        size_t number_of_cuts = number_of_cuts_distribution(randgen);
        if (number_of_cuts == 0) {
            return {md};
        }
        vector<int> dirichelet_weights(number_of_cuts + 1, 1);  // TODO: Explore non uniform dirichelet
        dirichelet_distribution<> dirichelet{
            dirichelet_weights, static_cast<double>(md.size())};  // TODO: Maybe add a constructor that accepts number
                                                                  // of params and sets them to one instead
        vector<double> cuts = dirichelet(randgen);
        if (md.has_comment("circular")) {
            std::uniform_int_distribution<unsigned long> uniform_dist(0, md.size());
            size_t pos = uniform_dist(randgen);
            auto mdp = split_molecule(md, pos);

            // Find start position.
            for (const einterval &ee : mdp.first.get_segments()) {
                mdp.second.append_segment(ee);
            }
            mdp.second.id(md.get_id());
            cuts.pop_back();
            return split_molecule_n(mdp.second, cuts); 
        }
        return split_molecule_n(md, cuts);
    }

    cxxopts::ParseResult parse(int argc, char **argv) {
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
                "p,probability",
                "Per base probability of fragmentation",
                cxxopts::value<double>()
             )
            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"<MODULE>", "<MODULE> description"}, args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "probability"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                loge("{} is required!", param);
                ++missing_parameters;
            }
        }
        // Other parameter checks here

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

        double probability = args["probability"].as<double>();

        ifstream input(input_file);

        ofstream output(output_file);
        std::uniform_int_distribution<long> lendist(0, LONG_MAX);
        for (auto &md : stream_mdf(input, true)) {
            auto mds = cut_molecule(md, probability, rand_gen);
            for(const auto &mm : mds){
                output << mm;
            }
        }
        return 0;
    }

    void describe_program() {
        logi("Running [MODULE]");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Cut_module);
