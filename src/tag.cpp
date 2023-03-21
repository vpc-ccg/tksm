#include "tag.h"

#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "pimpl.h"
#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::string;
using std::vector;


class TAG_module::impl : public tksm_module {
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
                "f,umi-fasta",
                "output umi fasta file",
                cxxopts::value<string>()
            )(
                "5,format5",
                "5' TAG format",
                cxxopts::value<string>()->default_value("")
            )(
                "3,format3",
                "3' TAG format",
                cxxopts::value<string>()->default_value("")
            )(
                "contig-prefix",
                "Prefix of the umi contigs in the mdf and umi-fasta",
                cxxopts::value<string>()->default_value("tksm_umi_ctg")
            )(
                "skip-tag-hashing",
                "Skip hashing of the TAGs for saving space",
                cxxopts::value<bool>()->default_value("false")->implicit_value("true")
             );
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"umi", "TAGging module"}, args(parse(argc, argv)) {}
    ~impl() = default;
    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "umi-fasta"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                loge("{} is required!", param);
                ++missing_parameters;
            }
        }

        if (args["format5"].as<string>().empty() && args["format3"].as<string>().empty()) {
            loge("At least one of the TAG formats must be provided");
            ++missing_parameters;
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

        string format5 = args["format5"].as<string>();
        string format3 = args["format3"].as<string>();

        string mdf_file_path{args["input"].as<string>()};
        std::ifstream mdf_file{mdf_file_path};



        logi("Adding tags");

        if (isdigit(format5[0])) {
            int len = std::stoi(format5);
            format5 = std::string(len, 'N');
        }
        if (isdigit(format3[0])) {
            int len = std::stoi(format3);
            format3 = std::string(len, 'N');
        }

        fmt2seq make_seq5(format5);
        fmt2seq make_seq3(format3);

        string umi_ctg_prefix = args["contig-prefix"].as<string>();

        string outfile_name = args["output"].as<string>();
        logi("Adding TAGs and printing to: {}", outfile_name);
        fmtlog::poll(true);
        std::ofstream outfile{outfile_name};
        std::string umi_ref_file = args["umi-fasta"].as<string>();
        std::ofstream umifile{umi_ref_file};



        if(args["skip-tag-hashing"].as<bool>()){
            int index = 0;
            for(auto &md : stream_mdf(mdf_file, true)) {
                string umi_seq5        = make_seq5[rand_gen];
                string umi_seq3        = make_seq3[rand_gen];
                string umi_ctg_name    = fmt::format("{}_{}", umi_ctg_prefix, index);

                umifile << fmt::format(">{}\n", umi_ctg_name);
                umifile << umi_seq5 << umi_seq3 << "\n";

                int len5 = static_cast<int>(umi_seq5.size());
                int len3 = static_cast<int>(umi_seq3.size());
                if (len5 > 0) {
                    md.prepend_segment(ginterval{umi_ctg_name, 0, len5, true});
                }
                if (len3 > 0) {
                    md.append_segment(ginterval{umi_ctg_name, len5, len5 + len3, true});
                }
                ++index;
                outfile << md;
            }
        }
        else{
            int index = 0;
            std::unordered_map<string, int> umi_ctg_map;
            auto find_ctg = [&umi_ctg_map, &index, umi_ctg_prefix] (string &umi_seq) {
                auto it = umi_ctg_map.find(umi_seq);
                if(it == umi_ctg_map.end()){
                    umi_ctg_map[umi_seq] = index;
                    string umi_ctg_name    = fmt::format("{}_{}", umi_ctg_prefix, index);
                    ++index;
                    return umi_ctg_name;
                }
                else{
                    return fmt::format("{}_{}", umi_ctg_prefix, it->second);
                }
            };

            for(auto &md : stream_mdf(mdf_file, true)) {
                string umi_seq5        = make_seq5[rand_gen];

                string umi_seq3        = make_seq3[rand_gen];

                int len5 = static_cast<int>(umi_seq5.size());
                int len3 = static_cast<int>(umi_seq3.size());
                if (len5 > 0) {
                    string umi_ctg_name5 = find_ctg(umi_seq5);
                    md.prepend_segment(ginterval{umi_ctg_name5, 0, len5, true});
                }
                if (len3 > 0) {
                    string umi_ctg_name3 = find_ctg(umi_seq3);
                    md.append_segment(ginterval{umi_ctg_name3, len5, len5 + len3, true});
                }
                outfile << md;
            }

            for(auto [umi_seq, index] : umi_ctg_map){
                string umi_ctg_name    = fmt::format("{}_{}", umi_ctg_prefix, index);
                umifile << fmt::format(">{}\n", umi_ctg_name);
                umifile << umi_seq << "\n";
            }
        }
        return 0;
    }

    void describe_program() {
        logi("Running TAGging module");
        logi("Input MDF: {}", args["input"].as<string>());
        logi("Output MDF: {}", args["output"].as<string>());
        logi("Output TAG FASTA: {}", args["umi-fasta"].as<string>());
        logi("5' tag format: {}", args["format5"].as<string>());
        logi("3' tag format: {}", args["format3"].as<string>());
        logi("Contig prefix: {}", args["contig-prefix"].as<string>());
        logi("Seed {}", args["seed"].as<int>());
        fmtlog::poll(true);
    }
};


MODULE_IMPLEMENT_PIMPL_CLASS(TAG_module)
