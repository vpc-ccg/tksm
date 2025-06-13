#include "plasmids.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"


#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
 
#include "pimpl.h"


class plasmids_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "f,fasta",
                "Comma separated list of plasmid fasta files as the database",
                cxxopts::value<vector<string>>()
            )(
                "o,output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "output-ref",
                "output fasta reference file",
                cxxopts::value<string>()
            )(
                "length-threshold",
                "Length threshold to differentiate small and large plasmids",
                cxxopts::value<int>()->default_value("20000")
            )(
                "few-copy-distribution",
                "Length distribution of few copy plasmids(µ,σ)",
                cxxopts::value<vector<double>>()->default_value("2,3")
            )(
                "many-copy-distribution",
                "Length distribution of many copy plasmids (µ,σ)",
                cxxopts::value<vector<double>>()->default_value("80,40")
            )(
                "large-plasmid-many-copy-probability",
                "Large plasmid can be either few or many copy. This parameter decides if a small plasmid is many copy",
                cxxopts::value<double>()->default_value("0.01")
            )(
                "small-plasmid-many-copy-probability",
                "Small plasmid can be either few or many copy. This parameter decides if a small plasmid is many copy",
                cxxopts::value<double>()->default_value("0.66")
            )
            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"plasmids module", "Generates plasmid sequences and molecules"}, args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"fasta", "output"};
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

        vector<string> input_files = args["fasta"].as<vector<string>>();


        string output_file = args["output"].as<string>();
        string reference_output_file = args["output-ref"].as<string>();

        int length_threshold = args["length-threshold"].as<int>();

        vector<double> few_copy_distribution = args["few-copy-distribution"].as<vector<double>>();
        vector<double> many_copy_distribution = args["many-copy-distribution"].as<vector<double>>();

        double large_plasmid_many_copy_probability = args["large-plasmid-many-copy-probability"].as<double>();
        double small_plasmid_many_copy_probability = args["small-plasmid-many-copy-probability"].as<double>();

        std::uniform_real_distribution<> pick_type(0,1);

        std::normal_distribution<>       pick_coverage_few( few_copy_distribution[0], few_copy_distribution[1]);
        std::normal_distribution<>       pick_coverage_many( many_copy_distribution[0], many_copy_distribution[1]);
        auto generate_copy_count = [&] (int length) -> int{
            if(length > length_threshold ){
                if(pick_type(rand_gen) < large_plasmid_many_copy_probability)
                    return pick_coverage_many(rand_gen);
                else
                    return pick_coverage_few(rand_gen);
            }
            else{
                if(pick_type(rand_gen) < small_plasmid_many_copy_probability)
                    return pick_coverage_many(rand_gen);
                else
                    return pick_coverage_few(rand_gen);
            
            }
        };

        ofstream output(output_file);
        ofstream ref_output(reference_output_file);
        int index = 0;
        for(const string &ifile : input_files){
            gzFile fp = gzopen(ifile.c_str(), "r");
            if(!fp){
                loge("Cannot open: {}", ifile);
                return -1;
            }

            kseq_t *seq = kseq_init(fp);

            while(kseq_read(seq) >= 0){
                int copy_count = std::max( 1, generate_copy_count(seq->seq.l));
                molecule_descriptor md{fmt::format("{}_{}",seq->name.s,index), true};
                md.append_segment({string{seq->name.s}, 0, (int) seq->seq.l, true})->depth(copy_count);
                output << md;
                ref_output << fmt::format(">{} {}\n{}\n", seq->name.s, seq->comment.s, seq->seq.s);
                ++index;
            }
            kseq_destroy(seq);
            gzclose(fp);
            //Read fasta, generate contigs
            //Assign copy number considering its size
            //Append plasmid sequence to the output fasta
            //Append plasmid molecule to the output md
        }
        return 0;
    }

    void describe_program() {
        logi("Running [plasmids]");

        for(const string &ifile : args["fasta"].as<vector<string>>()){
            logi("Input file: {}", ifile);
        }
        logi("Output file: {}", args["output"].as<string>());
        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(plasmids_module);
