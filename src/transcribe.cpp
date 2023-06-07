#include "transcribe.h"

#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "fusion.cpp"
#include "gtf.h"
#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "pimpl.h"
#include "util.h"

using std::string;
using std::vector;

class Splicer_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char** argv) {
        // clang-format off
        options.add_options("main")
        (
            "g,gtf",
            "Path to gtf annotation files",
            cxxopts::value<vector<string>>()
        )(
            "a,abundance",
            "Path to tab separated abundances tables (formatted as transcript_id\\tpm\\cell-barcode)",
            cxxopts::value<vector<string>>()
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
        )(
            "molecule-prefix",
            "Prefix for molecule names",
            cxxopts::value<string>()->default_value("M")
        )(
            "w,weights",
            "Comma separated weights of each provided abundance file",
            cxxopts::value<vector<double>>()->default_value("1")
         )
        ;
        // clang-format on
        return options.parse(argc, argv);
    }

    auto process_file_weights(const cxxopts::ParseResult&args){


        auto W  = args["weights"].as<vector<double>>();
        if(W.size() == 1){
            return vector<double>(W.size(), W[0]/args["abundance"].as<vector<string>>().size());
        }
        assert(v_strings.size() == args["abundance"].as<vector<string>>().size());
    

        double sum = std::accumulate(W.begin(), W.end(), 0.0);
        for(auto& w : W){
            w /= sum;
        }
        return W;
    }


    Fusion_submodule fusion_submodule;
    cxxopts::ParseResult args;

    friend class Fusion_submodule;
public:
    impl(int argc, char** argv)
        : tksm_module{"Splicer", "RNA Splicing module"}, fusion_submodule(options, rand_gen), args(parse(argc, argv)) {}

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
        if (fusion_submodule.receive_arguments(args) == Fusion_submodule::submodule_status::ERROR) {
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

        vector<string> gtf_files               = args["gtf"].as<vector<string>>();
        vector<string> abundance_files        = args["abundance"].as<vector<string>>();
        string output_file            = args["output"].as<string>();
        int molecule_count                 = args["molecule-count"].as<int>();
        [[maybe_unused]] bool use_whole_id = args["use-whole-id"].as<bool>();
        [[maybe_unused]] bool non_coding   = args["non-coding"].as<bool>();
        int default_depth                  = args["default-depth"].as<int>();
        string molecule_prefix             = args["molecule-prefix"].as<string>();
        vector<double> file_weights                =  process_file_weights(args);

        std::uniform_real_distribution<> dist(0, 1);
        std::ofstream outfile{output_file};

        logi("Reading GTF files {}", fmt::join(gtf_files, ", "));
        std::unordered_map<string, transcript> isoforms;
        for(auto& gtf_file : gtf_files) {
            isoforms.merge(read_gtf_transcripts_deep(gtf_file, default_depth));
        }
        auto w_iter = file_weights.begin();
        for(const string &abundance_file : abundance_files) {
            double file_W = *w_iter;
            ++w_iter;
            std::ifstream abundance_file_stream{abundance_file};
            std::string buffer;
            logi("Reading abundance file {} and printing simulated molecules to {}!", abundance_file, output_file);
            if (!abundance_file_stream.is_open()) {
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
                format_annot_id(tid, !use_whole_id);
                abundances.emplace_back(tid, tpm, comment);
            }

            auto fusion_status = fusion_submodule.receive_arguments(args);
            if (fusion_status == Fusion_submodule::submodule_status::ERROR) {
                return 1;
            }
            else if (fusion_status == Fusion_submodule::submodule_status::RUN) {
                fusion_submodule.run(this, abundances, isoforms);
            }

            size_t index = 0;

            double sum_tpm = std::accumulate(abundances.begin(), abundances.end(), 0.0, [](double sum, auto& t) {
                return sum + std::get<1>(t);
            });

            for (auto& [tid, tpm, comment] : abundances) {
                format_annot_id(tid, !args["use-whole-id"].as<bool>());
                auto md_ptr = isoforms.find(tid);
                if (md_ptr == isoforms.end()) {
                    logw("Isoform {} is not found in the input GTFs!", tid);
                    continue;
                }
                molecule_descriptor molecule = md_ptr->second;
                double count                 = file_W * tpm * molecule_count / sum_tpm;
                double carry                 = count - int(count);

                if (dist(rand_gen) < carry) {
                    count++;
                }
                logd("Isoform {} has abundance {} and will be simulated {} times", tid, tpm, count);
                if (int(count) == 0) {
                    continue;
                }

                molecule.add_comment("tid", tid)
                    ->add_comment("CB", comment)
                    ->depth(count)
                    ->id(molecule_prefix + std::to_string(index++));
                outfile << molecule;
            }
        }
        return 0;
    }

    void describe_program() {
        logi("Splicer module");
        string gtf_files = fmt::format("{}",fmt::join(args["gtf"].as<vector<string>>(), ", "));
        logi("Input GTF files: {}", gtf_files);
        string abundance_files = fmt::format("{}",fmt::join(args["abundance"].as<vector<string>>(), ", "));
        logi("Input abundance files: {}", abundance_files);
        logi("Output file: {}", args["output"].as<string>());
        logi("Molecule count: {}", args["molecule-count"].as<int>());
        logi("Use whole transcript id: {}", args["use-whole-id"].as<bool>());
        logi("Process non-coding genes/transcripts: {}", args["non-coding"].as<bool>());
        logi("Default depth: {}", args["default-depth"].as<int>());
        logi("Molecule prefix: {}", args["molecule-prefix"].as<string>());

        fusion_submodule.describe_program(args);
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Splicer_module);
