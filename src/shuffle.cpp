#include "shuffle.h"
#include <cxxopts.hpp>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <limits>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;

#include "pimpl.h"


class Shuffle_module::impl : public tksm_module {

    void shuffle_stream( const string &input_file, ostream &output, size_t buffer_size = std::numeric_limits<size_t>::max()) {
        std::uniform_int_distribution<size_t> disko(0, buffer_size - 1);
        vector<molecule_descriptor> buffer;

        for(auto &md : stream_mdf(input_file, true)) {
            if(buffer_size > buffer.size()){
                buffer.push_back(md);
            }
            else{
                size_t pos = disko(rand_gen);
                output << buffer[pos];
                buffer[pos] = md;
            }
        }
        if(buffer.size() > 0) {
            std::ranges::shuffle(buffer, rand_gen);
            for(auto &md : buffer) {
                output << md;
            }
        }
    }
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
                "buffer-size",
                "Sets the buffer size for shuffling, if not set, the whole file is read into memory",
                cxxopts::value<size_t>()
            )

            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"shuffle", "shuffles mdfs"}, args(parse(argc, argv)) {}

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
        if (process_utility_arguments(args)) {
            return 0;
        }
        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        string input_file  = args["input"].as<string>();

        string output_file = args["output"].as<string>();
        ofstream output(output_file);

        size_t buffer_size = [&] () ->size_t{
            if(args["buffer-size"].count() > 0 ){
                return args["buffer-size"].as<size_t>();
            }
            else{
                return std::numeric_limits<size_t>::max();
            }
        }
        ();

        shuffle_stream( input_file, output, buffer_size);
        return 0;
    }

    void describe_program() {
        logi("Running Shuffle Module");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        //Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Shuffle_module);

