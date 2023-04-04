#include "polyA.h"

#include <cxxopts.hpp>
#include <variant>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "pimpl.h"
#include "util.h"


using namespace std::string_literals;
#ifndef POLYA_REF_PREFIX
#define POLYA_REF_PREFIX ""s
#endif

class PolyA_module::impl : public tksm_module {
private:
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options()
            (
                "i,input",
                "Input file",
                cxxopts::value<std::string>()
            )(
                "o,output",
                "Output file",
                cxxopts::value<std::string>()
            )(
                "gamma",
                "Use Gamma distribution [α,β]",
                cxxopts::value<vector<double>>()
            )(
                "poisson",
                "Use Poisson distribution [λ]",
                cxxopts::value<vector<double>>()
            )(
                "weibull",
                "Use Weibull distribution [α,β]",
                cxxopts::value<vector<double>>()
            )(
                "normal",
                "Use Normal distribution [μ,σ]",
                cxxopts::value<vector<double>>()
            )(
                "min-length",
                "Minimum length of polyA",
                cxxopts::value<int>()->default_value("0")
            )(
                "max-length",
                "Maximum length of polyA",
                cxxopts::value<int>()->default_value("5000")
            );
        // clang-format on
        auto result = options.parse(argc, argv);
        return result;
    }
    cxxopts::ParseResult args;

    int validate_arguments() {
        vector<string> mandatory = {"input", "output"};
        int missing_arguments = 0;
        for (auto &arg : mandatory) {
            if (!args.count(arg)) {
                report_missing_parameter(arg);
                missing_arguments++;
            }
        }
        if (!args.count("gamma") && !args.count("poisson") && !args.count("weibull") && !args.count("normal")) {
            loge("No distribution specified");
            missing_arguments++;
        }
        if (args.count("gamma") + args.count("poisson") + args.count("weibull") + args.count("normal") > 1) {
            loge("Multiple distributions specified");
            missing_arguments++;
        }
        if (args.count("gamma")) {
            auto vec = args["gamma"].as<vector<double>>();
            if (vec.size() != 2) {
                loge("Gamma distribution requires two parameters");
                missing_arguments++;
            }
        }
        if (args.count("poisson")) {
            auto vec = args["poisson"].as<vector<double>>();
            if (vec.size() != 1) {
                loge("Poisson distribution requires one parameter");
                missing_arguments++;
            }
        }
        if (args.count("weibull")) {
            auto vec = args["weibull"].as<vector<double>>();
            if (vec.size() != 2) {
                loge("Weibull distribution requires two parameters");
                missing_arguments++;
            }
        }
        if (args.count("normal")) {
            auto vec = args["normal"].as<vector<double>>();
            if (vec.size() != 2) {
                loge("Normal distribution requires two parameters");
                missing_arguments++;
            }
        }
        if (args["min-length"].as<int>() < 0) {
            loge("Minimum length of polyA cannot be negative");
            missing_arguments++;
        }
        if (args["max-length"].as<int>() < 0) {
            loge("Maximum length of polyA cannot be negative");
            missing_arguments++;
        }
        if (args["min-length"].as<int>() > args["max-length"].as<int>()) {
            loge("Minimum length of polyA cannot be greater than maximum length of polyA");
            missing_arguments++;
        }
        if(missing_arguments > 0){
            std::cerr << options.help() <<  std::endl;
        }
        fmtlog::poll(true);
        return missing_arguments;
    }
    string ctg;
public:
    impl(int argc, char **argv)
        : tksm_module("polyA module", "Adds polyA tails to molecules with given size distribution"),
          args(parse(argc, argv)), ctg(args["max-length"].as<int>(), 'A'){

    }

    template <class Distribution>
    static molecule_descriptor add_polyA(const molecule_descriptor &md, Distribution &dist, int min_polya_len, int max_polya_len, const string &ctg,auto &rand_gen) {
        int poly_a_len = dist(rand_gen);
        if (poly_a_len < min_polya_len) {
            poly_a_len = min_polya_len;
        }
        if (poly_a_len > max_polya_len) {
            poly_a_len = max_polya_len;
        }

        molecule_descriptor new_md = md;
        if(poly_a_len > 0) {
            new_md.append_segment(ginterval{ctg.substr(0,poly_a_len), 0, poly_a_len, true});
        }
        return new_md;
    }

    static auto polya_transformer(int min_length, int max_length, const string &ctg, auto &dist, auto &rand_gen) {
        return std::ranges::views::transform([&, min_length, max_length, dist](const auto &md) {
            return std::visit([&](auto dist) { return add_polyA(md, dist, min_length, max_length, ctg, rand_gen); }, dist);
        });
    }
    auto get_dist() -> std::variant<std::gamma_distribution<double>, std::poisson_distribution<int>,
                                    std::weibull_distribution<double>, std::normal_distribution<double>> {
        if (args.count("gamma")) {
            auto vec = args["gamma"].as<vector<double>>();
            return std::gamma_distribution<double>(vec[0], vec[1]);
        }
        if (args.count("poisson")) {
            auto vec = args["poisson"].as<vector<double>>();
            return std::poisson_distribution<int>(vec[0]);
        }
        if (args.count("weibull")) {
            auto vec = args["weibull"].as<vector<double>>();
            return std::weibull_distribution<double>(vec[0], vec[1]);
        }
        if (args.count("normal")) {
            auto vec = args["normal"].as<vector<double>>();
            return std::normal_distribution<double>(vec[0], vec[1]);
        }
        return std::gamma_distribution<double>(1, 1);
    }
    auto operator()() {
        auto dist =  get_dist();
        return polya_transformer(args["min-length"].as<int>(), args["max-length"].as<int>(), ctg, dist, rand_gen);
    }
    int run() {
        if (process_utility_arguments(args)) {
            return 0;
        }

        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        int min_length          = args["min-length"].as<int>();
        int max_length          = args["max-length"].as<int>();
        std::string input_file  = args["input"].as<std::string>();
        std::string output_file = args["output"].as<std::string>();

        std::ofstream output(output_file);
        if (!output.is_open()) {
            loge("Cannot open output file: {}", output_file);
            return 1;
        }

        auto dist = get_dist();
        for (const auto &md : stream_mdf(input_file, true) | polya_transformer(min_length, max_length, ctg, dist, rand_gen)) {
            output << md;
        }

        return 0;
    }

    void describe_program() {  // Use logi
        logi("Running polyA module");
        logi("Input file: {}", args["input"].as<std::string>());
        logi("Output file: {}", args["output"].as<std::string>());
        logi("Minimum length of polyA: {}", std::to_string(args["min-length"].as<int>()));
        logi("Maximum length of polyA: {}", std::to_string(args["max-length"].as<int>()));

        logi("Distribution: ");
        if (args.count("gamma")) {
            auto vec = args["gamma"].as<vector<double>>();
            logi("Gamma distribution [α,β]: {},{}", vec[0], vec[1]);
        }
        else if (args.count("poisson")) {
            auto vec = args["poisson"].as<vector<double>>();
            logi("Poisson distribution [λ]: {}", vec[0]);
        }
        else if (args.count("weibull")) {
            auto vec = args["weibull"].as<vector<double>>();
            logi("Weibull distribution [α,β]: {},{}", vec[0], vec[1]);
        }
        else if (args.count("normal")) {
            auto vec = args["normal"].as<vector<double>>();
            logi("Normal distribution [μ,σ]: {},{}", vec[0], vec[1]);
        }
        else {
            logi("No distribution specified");
        }
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(PolyA_module);
