#include "size_selection.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

#include "pimpl.h"

template<class T=double>
class logistic_function{
    using real_type = T;
    real_type location;
    real_type scale;
    static constexpr real_type sech(real_type x) {
        return real_type{1} / std::cosh(x);
    }
    public:
    logistic_function(real_type location, real_type scale) : location{location}, scale{scale} {}
    logistic_function() : location{0}, scale{1} {}
    
    real_type operator() (real_type x) {
        real_type intermediate = (x - location) / (scale * 2);
        return sech(intermediate) * sech(intermediate) * (real_type{.25}/scale);
    }
};

template<class T=double>
class sigmoid_function{
    using real_type = T;
    real_type location;
    real_type scale;

    public:
    sigmoid_function(real_type location, real_type scale) : location{location}, scale{scale} {}
    sigmoid_function() : location{0}, scale{1} {}
    
    real_type operator() (real_type x) {
        real_type _x = (x-location)/scale;
        real_type _exp = std::exp(_x);
        return _exp / (_exp + 1);
    }
};

class Size_selection_module::impl : public tksm_module {
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
                "l,location",
                "",
                cxxopts::value<int>()
             )(
                 "s,scale",
                 "",
                 cxxopts::value<int>()
              )
             (
                "f,function",
                "",
                cxxopts::value<string>()->default_value("sigmoid")
             )
            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"size-select", "size-select module selects molecules by size smoothed out with a sigmoid function"}, args(parse(argc, argv)) {}

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
        // Other parameter checks here

        if (missing_parameters > 0) {
            fmt::print(stderr, "{}\n", options.help());
            return 1;
        }
        return 0;
    }
    int run() {
        sigmoid_function<> foo{0,10};
        

        return 0;
        if (process_utility_arguments(args)) {
            return 0;
        }
        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        string input_file = args["input"].as<string>();

        string output_file = args["output"].as<string>();

        ifstream input(input_file);
        


        ofstream output(output_file);

        for (auto &md : stream_mdf(input)) {
            output << md;
        }
        return 0;
    }

    void describe_program() {
        logi("Running [MODULE]");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Size_selection_module);
