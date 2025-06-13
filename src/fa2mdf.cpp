#include "fa2mdf.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <variant>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::ifstream;
using std::ofstream;
using std::set;
using std::string;
using std::vector;

#include "pimpl.h"

class FA2MDF_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char ** argv){
        // clang-format off
        options.add_options("main")
        (
            "r,reference",
            "Whole genome reference file",
            cxxopts::value<string> ()
        )(
            "frag-len-dist",
            "Fragment length distribution",
            cxxopts::value<string> ()
        )(
            "o,output",
            "output mdf file",
            cxxopts::value<string> ()
        )(
            "circular",
            "Contigs are considered to be circular",
            cxxopts::value<bool> ()->default_value("false")->implicit_value("true")
          )
        ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char ** argv) : tksm_module{"<FA2MDF>", "<FA2MDF> description"}, args(parse(argc, argv)){ }

    ~impl() = default;

    int validate_arguments(){
        std::vector<string> mandatory = { "reference", "output"};
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
    } // validate_arguments

    int run(){
        if (process_utility_arguments(args)) {
            return 0;
        }
        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        string reference_file = args["reference"].as<string>();

        // Read fasta index file .fai
        string fai_file = reference_file + ".fai";
        ifstream fai(fai_file);
        string line;
        int64_t ref_length = 0;

        vector<string> ref_names_index;

        vector<int64_t> ref_lens_so_far;
        vector<int64_t> ref_lens;

        string output_file = args["output"].as<string>();

        ofstream output (output_file);
        while (std::getline(fai, line)) {
            std::istringstream iss(line);
            string ref_name;
            int64_t ref_len;
            int64_t ref_offset;
            int64_t ref_line_bases;
            int64_t ref_line_width;
            iss >> ref_name >> ref_len >> ref_offset >> ref_line_bases >> ref_line_width;
            ref_length += ref_len;
            ref_names_index.push_back(ref_name);
            ref_lens_so_far.push_back(ref_length);
            ref_lens.push_back(ref_len);

            molecule_descriptor mol{ ref_name, true};
            mol.append_segment({ ref_name, 0, ref_len, true });
            if ( args["circular"].as<bool>()){
                mol.add_comment("circular", "true");  
            }
            output << mol;
        }

        return 0;
    } // run

    void describe_program(){
        logi("Running [FA2MDF]");
        logi("Reference file: {}", args["reference"].as<string>());
        logi("Output file: {}", args["output"].as<string>());

        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(FA2MDF_module);
