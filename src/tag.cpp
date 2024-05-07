#include "tag.h"

#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "pimpl.h"
#include "util.h"

using std::string;
using std::vector;

class fmt2seq {
public:
    std::array<char, 4> table{{'A', 'T', 'C', 'G'}};
    std::array<string, 128> lookup;
    string fmt;
    fmt2seq(const string &fmt) : fmt(fmt) {
#define set_l(A, B) lookup[A[0]] = B
#define set_sl(A, B) set_l(#A, #B)
        set_sl(A, A);
        set_sl(G, G);
        set_sl(T, T);
        set_sl(C, C);
        set_sl(U, U);
        set_sl(R, GA);
        set_sl(Y, TC);
        set_sl(K, GT);
        set_sl(M, AC);
        set_sl(S, GC);
        set_sl(W, AT);
        set_sl(B, GTC);
        set_sl(D, GAT);
        set_sl(H, ACT);
        set_sl(V, GCA);
        set_sl(N, AGCT);
#undef set_l
#undef set_sl
    }
    template <class RANDGEN>
    string operator[](RANDGEN &gen) const {
        std::stringstream st;

        for (char c : fmt) {
            string buffer;

            std::sample(lookup[c].begin(), lookup[c].end(), std::back_inserter(buffer), 1, gen);
            st << buffer;
        }
        return st.str();
    }
};


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
                "5,format5",
                "5' TAG format",
                cxxopts::value<string>()->default_value("")
            )(
                "3,format3",
                "3' TAG format",
                cxxopts::value<string>()->default_value("")
            );
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"umi", "TAGging module"}, args(parse(argc, argv)) {}
    ~impl() = default;
    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output"};
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

        string outfile_name = args["output"].as<string>();
        logi("Adding TAGs and printing to: {}", outfile_name);
        fmtlog::poll(true);
        std::ofstream outfile{outfile_name};

        for (auto &md : stream_mdf(mdf_file_path, true)) {
            string umi_seq5 = make_seq5[rand_gen];

            string umi_seq3 = make_seq3[rand_gen];

            int len5 = static_cast<int>(umi_seq5.size());
            int len3 = static_cast<int>(umi_seq3.size());
            if (len5 > 0) {
                md.prepend_segment(ginterval{umi_seq5, 0, len5, true});
            }
            if (len3 > 0) {
                md.append_segment(ginterval{umi_seq3, 0, len3, true});
            }
            outfile << md;
        }

        return 0;
    }

    void describe_program() {
        logi("Running TAGging module");
        logi("Input MDF: {}", args["input"].as<string>());
        logi("Output MDF: {}", args["output"].as<string>());

        logi("5' tag format: {}", args["format5"].as<string>());
        logi("3' tag format: {}", args["format3"].as<string>());

        logi("Seed {}", args["seed"].as<int>());
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(TAG_module)
