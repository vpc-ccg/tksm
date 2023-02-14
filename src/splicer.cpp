#include <exception>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <fstream>
#include <set>


#include <climits>

#include <cgranges/IITree.h>
#include <cxxopts.hpp>

#include "interval.h"
#include "tree.h"
#include "gtf.h"
#include "mdf.h"

#include "umi.h"
#include <cstring>
#include <utility>


#include <fmt/core.h>
#define FMTLOG_HEADER_ONLY
#include <fmtlog/fmtlog.h>

using std::pair;
using std::ofstream;
using std::vector;
using std::string;
using std::map;

std::mt19937 rand_gen{std::random_device{}()};


void describe_program(const cxxopts::ParseResult &args){
    logi("Running Splicer Module:");
    logi("gtf annotation: {}", args["gtf"].as<string>());
    logi("abundance: {}", args["abundance"].as<string>());
    logi("output: {}", args["output"].as<string>());
    logi("Non-coding transcripts are {}",  args["non-coding"].as<bool>() ? "not skipped" : "skipped");
    logi("seed: {}", args["seed"].as<int>());
    fmtlog::poll(true);
}

int main(int argc, char **argv){

    cxxopts::Options options("tksm Splicer", "Splicer module of tksm");

    options.add_options()
        ("g,gtf",  "Path to gtf annotation file", cxxopts::value<string>())
        ("a,abundance",  "Path to tab separated abundance table (Formatted as transcript_id\\tcount\\tpm)", cxxopts::value<string>())
        ("use-whole-id", "Use whole transcript id instead of first 15 characters (ENSEMBL ids)", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("molecule-count", "Number of molecules to simulate", cxxopts::value<int>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("non-coding", "Process non-coding genes/transcripts as well", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("default-depth", "Default depth for transcripts that are not in expression table", cxxopts::value<int>()->default_value("0"))
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("verbosity", "choose verbosity among [DEBUG, INFO, WARN, ERROR, OFF]", cxxopts::value<string>()->default_value("WARN"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    fmtlog::setLogLevel(LogLevels::parse_loglevel(args["verbosity"].as<string>()));

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"gtf", "output", "abundance", "molecule-count"};
    
    int missing_parameters = 0;
    for( string &param : mandatory){
        if(args.count(param) == 0){
            loge("{} is required!", param);
            ++missing_parameters;
        }
    }
    if(missing_parameters  > 0){
        fmt::print(stderr, "{}", options.help());
        return 1;
    }

    describe_program(args);


    //parameters will be moved to argument parser when done
    string path_to_gtf   {args["gtf"].as<string>()};
    string out_path {args["output"].as<string>()};

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);
   
    logi("Reading gtf annotation: {}...", path_to_gtf);
    fmtlog::poll(true);
    auto isoforms = read_gtf_transcripts(path_to_gtf, args["default-depth"].as<int>());


    std::ifstream table_file(args["abundance"].as<string>());
    ofstream outfile( out_path);
    if(!table_file){
        loge("Cannot open {}! Terminating!", args["abundance"].as<string>());
        return 1;
    }
    double molecule_count = args["molecule-count"].as<int>();

    string tid = "BEG";
    double tpm;
    string comment;

    std::uniform_real_distribution<> dist(0, 1);

    std::string buffer;
    std::getline(table_file, buffer); //Read header

    logi("Simulating molecule descriptions: {}...", out_path);
    fmtlog::poll(true);
    while(std::getline(table_file, buffer)){
        std::istringstream(buffer) >> tid >> tpm >> comment;
        format_annot_id(tid, !args["use-whole-id"].as<bool>());
        auto md_ptr = isoforms.find(tid);
        if(md_ptr == isoforms.end()){
            logw("Isoform {} is not found in the input GTF {}!", tid, path_to_gtf);
            continue;
        }
        molecule_descriptor molecule = isoforms[tid];
        comment.append(";");
        comment = "CB=" + comment;
        double count = tpm*molecule_count/1'000'000;
        double carry = count - int(count);

        if( dist(rand_gen) < carry){
            count++;
        }
        if( int(count) == 0){
            continue;
        }
        molecule.comment(comment)->depth(count); //rounded down to int

        outfile << molecule;
    }
    logi("Splicing Done...");
    return 0;
}


