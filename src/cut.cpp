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

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

#include "pimpl.h"
inline std::pair<molecule_descriptor, molecule_descriptor> split_molecule(molecule_descriptor &md, int split_pos){
    
    int bases_so_far = 0;
    molecule_descriptor md1;
    molecule_descriptor md2;
    for(const einterval &g : md.get_segments()){
        std::cerr << fmt::format("Cutting at {}; Current base {}; Segment size {};", split_pos, bases_so_far, g.size()) << "\n";
        if(g.size() + bases_so_far < split_pos){
            md1.append_segment(g);
        }
        else{
            if(bases_so_far  < split_pos) { // Transition interval
                einterval cpy = g;

                einterval cpy2 = g;
               if(g.plus_strand){
                    cpy.truncate(0, split_pos - bases_so_far);
                    cpy2.truncate(split_pos - bases_so_far, cpy2.size());
               }
               else{
                   cpy.truncate(split_pos - bases_so_far, cpy.size());
                   cpy2.truncate(0, split_pos - bases_so_far);
               }
                md1.append_segment(cpy);
                md2.append_segment(cpy2);
            }
            else{
                md2.append_segment(g);
            }
        }
        bases_so_far += g.size();
    }
    return {md1,md2};
}
class Cut_module::impl : public tksm_module {
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
        std::vector<string> mandatory = {"input", "output"};
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

        ifstream input(input_file);

        ofstream output(output_file);
        std::uniform_int_distribution<long> lendist(0, LONG_MAX);
        for (auto &md : stream_mdf(input,true)) {
            auto iter = md.meta.find("circular");
            long pos = lendist(rand_gen) % md.size();
            auto mdp = split_molecule(md, pos);
            if( iter != md.meta.end()){ // Cut circular
                //Find start position.
                std::cerr<<"Circular " << pos << "\n";
                for(const einterval &ee : mdp.first.get_segments()){
                    mdp.second.append_segment(ee);
                }
                mdp.second.id(md.get_id());
                output << mdp.second;
            }
            else{
                mdp.first.id(md.get_id() + "-0");

                mdp.second.id(md.get_id() + "-1");
                output << mdp.first << mdp.second;
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
