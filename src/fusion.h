#ifndef _FUSION_H_
#define _FUSION_H_
#include <random>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

class Fusion_module : public tksm_module {
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
                "g,gtf",
                "Path to GTF annotation file",
                cxxopts::value<string>()
            )(
                "fusion-file",
                "Path to tab separated fusion file",
                cxxopts::value<string>()
            )(
                "fusion-count",
                "Number of random fusions to generate",
                cxxopts::value<int>()->default_value("0")
            )(
                "disable-deletions",
                "Disables deletions (from fusions) that removes expression on the overlapping genes",
                cxxopts::value<bool>()->default_value("false")->implicit_value("true")
            )(
                "translocation-ratio",
                "Ratio of translocated fusions",
                cxxopts::value<double>()->default_value("0")
            )
            ;
        // clang-format on

        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    Fusion_module(int argc, char **argv) : tksm_module{"fusion", "Fusion module"}, args(parse(argc, argv)) {}

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "gtf"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                report_missing_parameter(param);
                ++missing_parameters;
            }
        }

        if (missing_parameters > 0) {
            std::cerr << options.help() << std::endl;
            return 1;
        }

        if (args["fusion-count"].as<int>() == 0 && args["fusion-file"].as<string>().empty()) {
            loge("Either fusion-file or fusion-count must be specified");
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

        string mdf_file_path{args["input"].as<string>()};
        string output_file_path{args["output"].as<string>()};
        string gtf_file_path{args["gtf"].as<string>()};

        std::ifstream mdf_file{mdf_file_path};
        auto streamer = stream_mdf(mdf_file);

        vector<gtf> transcripts = read_gtf(gtf_file_path);

        int fusion_count           = args["fusion-count"].as<int>();
        double translocation_ratio = args["translocation-ratio"].as<double>();
        bool disable_deletions     = args["disable-deletions"].as<bool>();

        vector<molecule_descriptor> fusions;  // Using md for now, but should be changed to a new class maybe
        if (args["fusion-file"].count() > 0) {
            string fusion_file_path{args["fusion-file"].as<string>()};
            std::ifstream fusion_file{fusion_file_path};
            fusions = read_fusions(fusion_file);
        }
        generate_fusions(fusions, fusion_count, transcripts, translocation_ratio);

        map<string, double> genes_to_be_removed;
        if (!disable_deletions) {
            genes_to_be_removed = get_genes_to_be_removed(fusions, transcripts);
        }

        std::ofstream output_file{output_file_path};

        std::uniform_real_distribution<double> dist(0.0, 1.0);
        while (streamer) {
            auto md = streamer();
            auto it = genes_to_be_removed.find(md.gene);
            if (it == genes_to_be_removed.end() || it->second < dist(rand_gen)) {
                output_file << md;
            }
        }

        // Generate md for the fusion transcripts

        return 0;
    }

    void describe_program() {
        logi("Running Fusion module");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        logi("GTF file: {}", args["gtf"].as<string>());
        logi("Fusion file: {}", args["fusion-file"].as<string>());
        logi("Fusion count: {}", args["fusion-count"].as<int>());
        if (args["disable-deletions"].as<bool>()) {
            logi("Deletions are disabled");
        }
        logi("Translocation ratio: {}", args["translocation-ratio"].as<double>());
        
        fmtlog::poll(true);
    }
};

#endif
