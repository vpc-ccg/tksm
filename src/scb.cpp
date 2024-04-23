#include "scb.h"

#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "pimpl.h"
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
                "keep-meta-barcodes",
                "Keep the barcodes in the mdf metadata",
                cxxopts::value<bool>()
            );
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"scb", "Single cell barcode module"}, args(parse(argc, argv)) {}
    ~impl() {}

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output"};
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

        string outfile_name = args["output"].as<string>();
        logi("Adding barcodess and printing to: {}", outfile_name);
        fmtlog::poll(true);
        std::ofstream outfile{outfile_name};

        for (auto &md : stream_mdf(mdf_file_path, true)) {
            const string &barcode_str = md.get_comment("CB")[0];
            if (barcode_str != ".") {
                md.append_segment(ginterval{barcode_str, 0, (int)barcode_str.size(), true});
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

        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(SingleCellBarcoder_module);
