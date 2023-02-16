#include "single-cell-barcoder.h"

#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "pimpl_impl.h"
#include "util.h"

class SingleCellBarcoder_module::impl : public tksm_module {
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
                "f,barcode-fasta",
                "output barcode fasta file",
                cxxopts::value<string>()
            )(
                "keep-meta-barcodes",
                "Keep the barcodes in the mdf metadata",
                cxxopts::value<bool>()
            );
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv)
        : tksm_module{"single-cell-barcoder", "Single cell barcode module"}, args(parse(argc, argv)) {}

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "barcode-fasta"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                std::cerr << param << " is required!\n";
                ++missing_parameters;
            }
        }

        if (missing_parameters > 0) {
            std::cerr << options.help() << std::endl;
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
        std::ifstream mdf_file{mdf_file_path};
        auto streamer = stream_mdf(mdf_file);

        string outfile_name = args["output"].as<string>();
        logi("Adding barcodess and printing to: {}", outfile_name);
        fmtlog::poll(true);
        std::ofstream outfile{outfile_name};
        std::string fasta_file_name = args["barcode-fasta"].as<string>();
        std::ofstream fastafile{fasta_file_name};
        std::unordered_set<string> used_barcodes;

        while (streamer) {
            molecule_descriptor md      = streamer();
            const string &barcode_str   = md.get_comment("CB")[0];
            const string barcode_ctg_id = "CB_" + barcode_str;
            if (barcode_str != ".") {
                md.prepend_segment(ginterval{barcode_ctg_id, 0, (int)barcode_str.size(), "+"});
                if (used_barcodes.count(barcode_str) == 0) {
                    fastafile << ">" << barcode_ctg_id << "\n" << barcode_str << "\n";
                    used_barcodes.insert(barcode_str);
                }
            }
            if (!args["keep-meta-barcodes"].as<bool>()) {
                md.drop_comment("CB");
            }
            outfile << md;
        }
        return 0;
    }

    void describe_program() {
        logi("Running Single-cell barcoding module");
        logi("Input MDF: {}", args["input"].as<string>());
        logi("Output MDF: {}", args["output"].as<string>());
        logi("Output barcode fasta: {}", args["barcode-fasta"].as<string>());
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPLE_CLASS(SingleCellBarcoder_module);
