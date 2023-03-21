#include "filter.h"
#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>
#include <fstream>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;

#include "pimpl.h"

class FilterCondition{
    public:
        std::function<bool(const molecule_descriptor &)> cond_func;
    FilterCondition(const string &condition) {
        const vector<string> fields = rsplit(condition, " ");
        if (fields.size() != 2) {
            throw std::runtime_error("Invalid condition: " + condition);
        }
        const string condition_type_str = fields[0];
        const string condition_expression_str = fields[1];
        if (condition_type_str == "info"){
            cond_func = [condition_expression_str](const molecule_descriptor &md) {
                auto iter = md.meta.find(condition_expression_str);
                if( iter == md.meta.end()) {
                    return false;
                }
                if(iter->second.empty()){
                    return false;
                }
                if(iter->second[0]== "."){
                    return false;
                }
                return true;
            };
        }
    }
    operator bool() const {
        return true;
    }
};

class Filter_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "i,input",
                "input mdf file",
                cxxopts::value<string>()
            )(
                "t,true-output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "f,false-output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "c,condition",
                "Comma separated conditions to filter (and)",
                cxxopts::value<vector<string>>()
             )
            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"Filter", "Splits the input to 2 files w.r.t. condition"}, args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "true-output", "condition"};
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

        string input_file  = args["input"].as<string>();

        vector<ofstream> output_files;
        output_files.emplace_back(args["true-output"].as<string>());
        if(args.count("false-output")) {
            output_files.emplace_back(args["false-output"].as<string>());
        }
        
        vector<FilterCondition> conditions;
        for(const auto &condition : args["condition"].as<vector<string>>()) {
            conditions.emplace_back(condition);
        }

        for(auto &md : stream_mdf(input_file)) {
            bool removed = false;
            for(const auto &condition : conditions) {

                if(!condition.cond_func(md)) {
                    output_files[1] << md;
                    removed = true;
                }
            }
            if(!removed) {
                output_files[0] << md;
            }
        }
        for(auto &output_file : output_files) {
            output_file.close();
        }
        return 0;
    }

    void describe_program() {
        logi("Running filter module");
        logi("Input file: {}", args["input"].as<string>());
        logi("True output file: {}", args["true-output"].as<string>());
        if(args.count("false-output")) {
            logi("False output file: {}", args["false-output"].as<string>());
        }
        logi("Conditions: {}", fmt::join(args["condition"].as<vector<string>>(), " and "));
        //Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Filter_module);

