#include "mutate.h"

#include <cxxopts.hpp>
#include <fstream>
#include <random>
#include <string>
#include <variant>
#include <vector>
#include <sstream>
#include <algorithm>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::map;
#include "pimpl.h"

class Mutate_module::impl : public tksm_module {

    using position_type = int;
    static bool is_number(const std::string& s)
    {
        return !s.empty() && std::find_if(s.begin(), 
            s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
    }
    struct Mod{
        using einterval_or_null = std::optional<einterval>;
        struct snv {
            char replacement;
            std::pair<std::nullopt_t, std::nullopt_t> operator () (einterval &i, position_type pos) const {
                i.add_error(pos - i.start, replacement); // Needs strand logic
                return {std::nullopt, std::nullopt};
            }
        };
        struct ins { //TODO DISTRIBUTE THE EXISTING MODS
            string insertion;
            std::pair<einterval, einterval> operator() (einterval &i, position_type pos) const {
                position_type end = i.end;
                i.end = pos;
                if(insertion[0] != '.'){
                    
                }
                return std::make_pair(einterval{insertion.substr(1),0,(int)insertion.size()-1, "+"}, einterval{i.chr, pos, end, "+"}); // Make sure end is correct
            }
        };//TODO DISTRIBUTE THE EXISTING MODS
        struct del { // Add 2 mods per deletion execute both (Only one should work).
            position_type pos1;
            position_type pos2;
            std::pair<einterval, std::nullopt_t> operator() (einterval &i, position_type pos) const {
                position_type from, to;
                if( pos == pos1){
                    from = pos;
                    to = pos2;
                
                }else if (pos == pos2){
                    
                }

                return {i, std::nullopt};
            }
            del reversed(){
                return {pos2, pos1};
            }
        };
        std::variant<snv, ins, del> vari;

        public:
        position_type pos;
            Mod(position_type pos, const string &mod_string) : pos{pos}{
                if(is_number(mod_string)){
                    vari = del{pos, stoi(mod_string)};
                }
                else if (mod_string.size() == 1){
                    vari = snv{mod_string[0]};
                }
                else{
                    vari = ins{mod_string};
                }
            }
        
            auto operator () (einterval &i) const -> std::pair<einterval_or_null, einterval_or_null> {
                return std::visit( [ &i, this] (auto var) -> decltype(this->operator()(i)){
                        return var(i,pos);
                        }, this->vari);
            }

            auto reversed() const {
                Mod other = *this;
                other.vari = std::get<del>(other.vari).reversed();
                return other;
            }

    };
    static void apply_mods(molecule_descriptor &md, const map<string, vector<Mod>> &mods){
        for(einterval &e : md.get_segments()){
            if(!mods.contains(e.chr)){
                continue;
            }
            const vector<Mod> &mods_in_chr = mods.at(e.chr);
            auto iter = std::ranges::lower_bound(mods_in_chr, e.start, std::ranges::less{}, &Mod::pos);
            while(iter != mods_in_chr.end() && iter->pos < e.end){
                (*iter)(e);
            }
        }
    }
    auto read_modifications(const string &path) -> map<string, vector< Mod>>{
        
        map<string, vector<Mod>> modification_forest;
        
        ifstream file (path);
        string buffer;
        buffer.reserve(100);
        while(std::getline(file,buffer)){
            string chr;
            position_type pos;
            string modification_string;
            
            std::istringstream{buffer} >> chr >> pos >> modification_string;
            modification_forest[chr].emplace_back(pos, modification_string);
            auto &back =  modification_forest[chr].back();
            if ( std::holds_alternative<Mod::del> (back.vari)){
                modification_forest[chr].push_back(back.reversed());
            }
        }
        for(auto &chr_tree : modification_forest){
            std::ranges::sort(chr_tree.second, [](const auto &A, const auto &B){ return A.pos < B.pos;});
            //std::ranges::sort(chr_tree.second, [](const auto &A, const auto &B){ return A.pos < B.pos;});
        }
        return modification_forest;
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
                "t,tsv",
                "Mutations in a tabular file of 'chr\tpos\tmodification'",
                //
                // SNV : chr1 100 A
                //      MDF: chr1 75 125 + 25A     
                // DEL : chr1 100 150
                //      MDF: chr1 75 125 + 25A     
                // INS : chr1 100 .AGT
                cxxopts::value<string>()
             )

            ;
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"<Mutate>", "<Mutate> description"}, args(parse(argc, argv)) {}

    ~impl() = default;

    int validate_arguments() {
        vector<string> mandatory = {"input", "output"};
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

        string mutation_file = args["tsv"].as<string>();

        auto mutation_forest = read_modifications(mutation_file);

        ifstream input(input_file);

        ofstream output(output_file);

        for (auto &md : stream_mdf(input)) {
            apply_mods(md, mutation_forest);
            output << md;
        }
        return 0;
    }

    void describe_program() {
        logi("Running [Mutate]");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        // Other parameters logs are here
        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Mutate_module);
