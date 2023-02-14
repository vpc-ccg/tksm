#ifndef _UMI_H_
#define _UMI_H_
#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

class UMI_module : public tksm_module {
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
                "5' UMI format",
                cxxopts::value<string>()->default_value("")
            )(
                "3,format3",
                "3' UMI format",
                cxxopts::value<string>()->default_value("")
            )(
                "contig-prefix",
                "Prefix of the umi contigs in the mdf and umi-fasta",
                cxxopts::value<string>()->default_value("tksm_umi_ctg")
            );
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    UMI_module(int argc, char **argv) : tksm_module{"umi", "UMI tagging module"}, args(parse(argc, argv)) {}

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
            loge("At least one of the UMI formats must be provided");
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

        auto streamer = stream_mdf(mdf_file);

        logi("Adding UMI tags");

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
        logi("Adding UMIs and printing to: {}", outfile_name);
        fmtlog::poll(true);
        std::ofstream outfile{outfile_name};
        std::string umi_ref_file = args["umi-fasta"].as<string>();
        std::ofstream umifile{umi_ref_file};

        int index = 0;
        while (streamer) {
            molecule_descriptor md = streamer();
            string umi_seq5        = make_seq5[rand_gen];
            string umi_seq3        = make_seq3[rand_gen];
            string umi_ctg_name    = fmt::format("{}_{}", umi_ctg_prefix, index);

            umifile << fmt::format(">{}\n", umi_ctg_name);
            umifile << umi_seq5 << umi_seq3 << "\n";

            int len5 = static_cast<int>(umi_seq5.size());
            int len3 = static_cast<int>(umi_seq3.size());
            if (len5 > 0) {
                md.prepend_segment(ginterval{umi_ctg_name, 0, len5, "+"});
            }
            if (len3 > 0) {
                md.append_segment(ginterval{umi_ctg_name, len5, len5 + len3, "+"});
            }
            ++index;
            outfile << md;
        }

        return 0;
    }

    void describe_program() {
        logi("Running UMI tagging module");
        logi("Input MDF: {}", args["input"].as<string>());
        logi("Output MDF: {}", args["output"].as<string>());
        logi("Output UMI FASTA: {}", args["umi-fasta"].as<string>());
        logi("5' UMI tag format: {}", args["format5"].as<string>());
        logi("3' UMI tag format: {}", args["format3"].as<string>());
        logi("Random seed: {}", args["seed"].as<int>());
        fmtlog::poll(true);
    }
};

#endif
