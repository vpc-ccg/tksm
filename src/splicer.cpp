#include "splicer.h"
#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "pimpl.h"
#include "gtf.h"
#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"
#include "fusion.cpp"


using std::string;
using std::vector;

class Splicer_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char** argv) {
        // clang-format off
        options.add_options("main")
        (
            "g,gtf",
            "Path to gtf annotation file",
            cxxopts::value<string>()
        )(
            "a,abundance",
            "Path to tab separated abundance table (formatted as transcript_id\\tpm\\cell-barcode)",
            cxxopts::value<string>()
        )(
            "use-whole-id",
            "Do not trim the transcript version",
            cxxopts::value<bool>()->default_value("false")->implicit_value("true")
        )(
            "molecule-count",
            "Number of molecules to simulate",
            cxxopts::value<int>()
        )(
            "o,output",
            "Output path",
            cxxopts::value<string>()
        )(
            "non-coding",
            "Process non-coding genes/transcripts as well",
            cxxopts::value<bool>()->default_value("false")->implicit_value("true")
        )(
            "default-depth",
            "Default depth for transcripts that are not in expression table",
            cxxopts::value<int>()->default_value("0")
        );
        // clang-format on
        return options.parse(argc, argv);
    }

    Fusion_submodule fusion_submodule;
    cxxopts::ParseResult args;

    friend class Fusion_submodule;
public:
    impl(int argc, char** argv) : tksm_module{"Splicer", "RNA Splicing module"},fusion_submodule(options, rand_gen), args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"gtf", "abundance", "output", "molecule-count"};
        int missing_parameters        = 0;
        for (string& param : mandatory) {
            if (args.count(param) == 0) {
                loge("Missing mandatory parameter {}", param);
                ++missing_parameters;
            }
        }
        if(fusion_submodule.receive_arguments(args)==Fusion_submodule::submodule_status::ERROR){
            ++missing_parameters;
        }
        if (missing_parameters > 0) {
            fmt::print("{}\n", options.help());
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


        std::string gtf_file               = args["gtf"].as<string>();
        std::string abundance_file         = args["abundance"].as<string>();
        std::string output_file            = args["output"].as<string>();
        int molecule_count                 = args["molecule-count"].as<int>();
        [[maybe_unused]] bool use_whole_id = args["use-whole-id"].as<bool>();
        [[maybe_unused]] bool non_coding   = args["non-coding"].as<bool>();
        int default_depth                  = args["default-depth"].as<int>();

        std::uniform_real_distribution<> dist(0, 1);
        std::ofstream outfile{output_file};

        logi("Reading GTF file {}", gtf_file);
        std::map<string, molecule_descriptor> isoforms = read_gtf_transcripts(gtf_file, default_depth);

        std::ifstream abundance_file_stream{abundance_file};
        std::string buffer;
        logi("Reading abundance file {} and printing simulated molecules to {}!", abundance_file, output_file);
        if(!abundance_file_stream.is_open()){
            loge("Could not open abundance file {}!", abundance_file);
            return 1;
        }
        std::getline(abundance_file_stream, buffer);  // skip header
        vector<std::tuple<string, double, string>> abundances;
        while (std::getline(abundance_file_stream, buffer)) {
            string tid     = "BEG";
            double tpm     = 0;
            string comment = "";
            std::istringstream(buffer) >> tid >> tpm >> comment;
            abundances.emplace_back(tid, tpm, comment);
        }

        auto fusion_status = fusion_submodule.receive_arguments(args);
        if(fusion_status == Fusion_submodule::submodule_status::ERROR){
            return 1;
        }
        else if (fusion_status == Fusion_submodule::submodule_status::RUN){
            auto fusions = fusion_submodule.run(this, abundances);
            for( auto &[gene, transcripts] : fusions){
                for( auto &t : transcripts){
                    abundances.emplace_back(t.info.at("transcript_id"), t.get_abundance(), t.get_comment());
                    isoforms.emplace(t.info.at("transcript_id"), t);
                }
            }
        }
    
        size_t index = 0;
        for(auto& [tid, tpm, comment] : abundances){
            format_annot_id(tid, !args["use-whole-id"].as<bool>());
            auto md_ptr = isoforms.find(tid);
            if (md_ptr == isoforms.end()) {
                logw("Isoform {} is not found in the input GTF {}!", tid, gtf_file);
                continue;
            }
            molecule_descriptor molecule = isoforms[tid];
            double count = tpm * molecule_count / 1'000'000;
            double carry = count - int(count);

            if (dist(rand_gen) < carry) {
                count++;
            }
            if (int(count) == 0) {
                continue;
            }

            molecule.
                add_comment("tid", tid)->
                add_comment("CB", comment)->
                depth(count)->
                id("Molecule_" + std::to_string(index++));
            outfile << molecule;
        }
        return 0;
    }

    void describe_program() {
        logi("Splicer module");
        logi("Input GTF file: {}", args["gtf"].as<string>());
        logi("Input abundance file: {}", args["abundance"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        logi("Molecule count: {}", args["molecule-count"].as<int>());
        logi("Use whole transcript id: {}", args["use-whole-id"].as<bool>());
        logi("Process non-coding genes/transcripts: {}", args["non-coding"].as<bool>());
        logi("Default depth: {}", args["default-depth"].as<int>());
        fusion_submodule.describe_program(args);
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Splicer_module);
