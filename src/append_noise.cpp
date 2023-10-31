#include "append_noise.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <string>
#include <variant>
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
class DistVisitor {
    
    std::variant<std::normal_distribution<>, std::lognormal_distribution<>> dist;
    public:

    auto get_dist(const vector<string> &dist_setting)
        -> std::variant<std::normal_distribution<>, std::lognormal_distribution<>> {
        if (dist_setting[0] == "normal") {
            auto [mean, stdev] = [&dist_setting] () -> std::pair<double, double>{
                return {stod(dist_setting[1]), stod(dist_setting[2])};
                        }();
            return std::normal_distribution<>{mean, stdev};
        }
        else if (dist_setting[0]== "lognormal") {
                // return
            auto [mean, stdev] = [&dist_setting] () -> std::pair<double, double>{
                return {stod(dist_setting[1]), stod(dist_setting[2])};
                        }();
            return std::lognormal_distribution<>{mean, stdev};
        }
        else{
            loge("Distribution not implemented! At: {}:{} {}", __FILE_NAME__, __PRETTY_FUNCTION__, __LINE__);
            exit(1);
        }
    }
        DistVisitor( const vector<string> &dist_settings): dist{get_dist(dist_settings)}{}

        int operator () (auto &rand_gen) { //Generate random val using the settings
            return std::visit([&](auto &dist) {return dist(rand_gen);}, dist);
        }
};

class NoiseAdder{
    bool palindromic;
    string alphabet;
    DistVisitor length_dist_picker;
    std::uniform_int_distribution<size_t> character_picker;
    std::uniform_real_distribution<double> error_dist;
    double error_rate;
    public:
        NoiseAdder(
                bool palindromic,
                const string &alphabet,
                const vector<string> &length_dist_str,
                double error_rate //if palindormic
                ) :
            palindromic{ palindromic},
            alphabet{ alphabet},
            length_dist_picker{ length_dist_str},
            character_picker{0, alphabet.size()-1},
            error_dist{0, 1},
            error_rate {error_rate} {}

        string generate_sequence(int length, auto &rand_gen){
            string seq;
            seq.reserve(length);
            for(int i = 0; i < length;++i){
                seq.push_back(alphabet[character_picker(rand_gen)]);
            }
            return seq;
        }

        int operator () (molecule_descriptor &md, auto &rand_gen, bool prepend = false) {
            int noise_length = length_dist_picker(rand_gen);
            if(noise_length <= 0){
                return 0;
            }
            else if (palindromic){
                // TO DO
                int pal_len_so_far = 0;
                const auto &mdseg = md.cget_segments();
                std::remove_const<std::remove_reference<decltype(mdseg)>::type>::type new_segments;
                for(auto it = mdseg.rbegin(); it != mdseg.rend(); ++it){
                    pal_len_so_far += it->size();
                    new_segments.push_back(*it);
                    new_segments.back().plus_strand = !new_segments.back().plus_strand; 
                    if(pal_len_so_far > noise_length){
                        int extra_len =  pal_len_so_far - noise_length;
                        if(it->plus_strand){
                            new_segments.back().end -= extra_len; 
                        }
                        else{
                            new_segments.back().start += extra_len; 
                        }
                        break;
                    }
                }
                for(auto &seg : new_segments){
                    
                    for(int i =0; i < seg.size(); ++i){
                        if(error_dist(rand_gen) < error_rate){
                            seg.add_error(i, alphabet[character_picker(rand_gen)]);
                        }
                    }
                    md.append_segment(seg);
                }
            }
            else{
                string noise_seq = generate_sequence(noise_length, rand_gen);
                if(prepend){
                    md.prepend_segment({noise_seq, 0, static_cast<int>(noise_seq.size()), true});
                }
                else{
                    md.append_segment({noise_seq, 0, static_cast<int>(noise_seq.size()), true}); // Fix this difference :)
                }
            }
            return noise_length;
        }
};

class AppendNoise_module::impl : public tksm_module {
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
                "alphabet",
                "Alphabet of the noise sequences (Not used for palindromic noise), same character can be used more than once to change occurance rate (AAAGTC, 50\% of the characters will be A)",
                cxxopts::value<string>()->default_value("AGTC")

             )(
                "palindromic",
                "Generate palindromic noise (Instead of random)",
                cxxopts::value<bool>()->default_value("false")
             )(
                "error-rate",
                "Error rate of the palindromic sequences",
                cxxopts::value<double>()->default_value("0.5")
             )(
                "length-dist",
                "Noise length distribution comma separated list of distribution_name,args,... (normal,0,0.5 for example for normal distribution with 0 mean and 0.5 stddev)",
                cxxopts::value<vector<string>>()
              )

            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"Append Noise module", "Append Noise module description"}, args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "length-dist"};
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

        string input_file = args["input"].as<string>();

        string output_file = args["output"].as<string>();

        NoiseAdder noiser {
            args["palindromic"].as<bool>(),
            args["alphabet"].as<string>(),
            args["length-dist"].as<vector<string>>(),
            args["error-rate"].as<double>()
        };

        ifstream input(input_file);

        ofstream output(output_file);

        for (auto &md : stream_mdf(input)) {
            noiser(md, rand_gen);
            output << md;
        }
        return 0;
    }

    void describe_program() {
        logi("Running Append Noise module");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(AppendNoise_module);
