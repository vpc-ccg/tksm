#include "random_wgs.h"

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

class RWGS_module::impl : public tksm_module {
    ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "r,reference",
                "Whole genome reference file",
                cxxopts::value<string>()
            )(
                "frag-len-dist",
                "Fragment length distribution",
                cxxopts::value<string>()
            )(
                "o,output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "base-count",
                "Number of bases to generate",
                cxxopts::value<int64_t>()
             )(
                 "depth",
                 "Depth of coverage",
                 cxxopts::value<double>()
              )
            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"<RWGS>", "<RWGS> description"}, args(parse(argc, argv)) {}

    ~impl() = default;

    set<string> implemented_dists = {"normal", "uniform", "lognormal", "exponential"};

    auto parse_dist_str(const string &st) const -> std::tuple<string, int64_t, int64_t> {
        string frag_len_dist;
        int64_t frag_len_dist_mean;
        int64_t frag_len_dist_std{0};

        std::istringstream parse_dist(args["frag-len-dist"].as<string>());
        parse_dist >> frag_len_dist >> frag_len_dist_mean;
        if (frag_len_dist != "exponential") {
            parse_dist >> frag_len_dist_std;
        }
        return {frag_len_dist, frag_len_dist_mean, frag_len_dist_std};
    }
    auto get_dist(string dist_name, int64_t a, int64_t b) const
        -> std::variant<std::normal_distribution<double>, std::uniform_real_distribution<double>,
                        std::lognormal_distribution<double>, std::exponential_distribution<double>> {
        if (dist_name == "normal") {
            return std::normal_distribution<double>(a, b);
        }
        else if (dist_name == "uniform") {
            return std::uniform_real_distribution<double>(a, b);
        }
        else if (dist_name == "lognormal") {
            return std::lognormal_distribution<double>(a, b);
        }
        else if (dist_name == "exponential") {
            return std::exponential_distribution<double>(a);
        }
        else {
            throw std::invalid_argument("Invalid distribution name");
        }
    }

    int validate_arguments() {
        std::vector<string> mandatory = {"reference", "output", "frag-len-dist"};
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

        if (args.count("base-count") == 0 && args.count("depth") == 0) {
            loge("Either base-count or depth is required!");
            return 1;
        }

        auto [frag_len_dist, frag_len_dist_mean, frag_len_dist_std] =
            parse_dist_str(args["frag-len-dist"].as<string>());
        if (implemented_dists.find(frag_len_dist) == implemented_dists.end()) {
            loge("Invalid fragment length distribution");
            return 1;
        }
        if (frag_len_dist_mean <= 0 || frag_len_dist_std < 0) {
            loge("Invalid fragment length distribution parameters");
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

        string reference_file = args["reference"].as<string>();

        // Read fasta index file .fai
        string fai_file = reference_file + ".fai";
        ifstream fai(fai_file);
        string line;
        int64_t ref_length = 0;

        vector<string> ref_names_index;

        vector<int64_t> ref_lens_so_far;
        vector<int64_t> ref_lens;
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
        }

        logi("Reference length: {}", ref_length);

        fmtlog::poll(true);

        string output_file = args["output"].as<string>();

        int64_t base_count = 0;
        if (args.count("base-count") > 0) {
            base_count = args["base-count"].as<int64_t>();
        }
        else if (args.count("depth") > 0) {
            double depth = args["depth"].as<double>();
            base_count   = depth * ref_length;
        }

        auto [frag_len_dist, frag_len_dist_mean, frag_len_dist_std] =
            parse_dist_str(args["frag-len-dist"].as<string>());
        auto frag_length_dist = get_dist(frag_len_dist, frag_len_dist_mean, frag_len_dist_std);
        auto position_dist    = std::uniform_int_distribution<int64_t>(0, ref_length - 1);
        auto strand_dist      = std::uniform_int_distribution<int64_t>(0, 1);

        ofstream output(output_file);

        int64_t generated_bases = 0;
        int64_t index           = 0;
        while (generated_bases < base_count) {
            int64_t pos       = position_dist(rand_gen);
            int64_t ref_index = 0;
            while (pos > ref_lens_so_far[ref_index]) {
                ++ref_index;
            }
            int ref_pos  = pos - ref_lens_so_far[ref_index] + ref_lens[ref_index];
            int frag_len = std::visit([&](auto &&arg) { return arg(rand_gen); }, frag_length_dist);
            if (frag_len > ref_lens[ref_index] - ref_pos) {
                frag_len = ref_lens[ref_index] - ref_pos;
            }
            bool plus_strand = strand_dist(rand_gen) == 0;
            molecule_descriptor mol{fmt::format("{}_{}:{}-{}{}", index, ref_names_index[ref_index], ref_pos,
                                                ref_pos + frag_len, plus_strand ? "+" : "-"),
                                    plus_strand};
            mol.append_segment({ref_names_index[ref_index], ref_pos, ref_pos + frag_len, plus_strand});
            output << mol;
            generated_bases += frag_len;
            ++index;
        }

        return 0;
    }

    void describe_program() {
        logi("Running [RWGS]");
        logi("Reference file: {}", args["reference"].as<string>());
        logi("Output file: {}", args["output"].as<string>());

        if (args.count("base-count") > 0) {
            logi("Base count: {}", args["base-count"].as<int64_t>());
        }
        else if (args.count("depth") > 0) {
            logi("Depth: {}", args["depth"].as<double>());
        }

        logi("Fragment length distribution: {}", args["frag-len-dist"].as<string>());

        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(RWGS_module);
