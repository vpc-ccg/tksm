#include <exception>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <fstream>
#include <set>

#include <climits>

#include "extern/IITree.h"
#include "extern/cxxopts.hpp"
#include "interval.h"
#include "tree.h"
#include "gtf.h"
#include "mdf.h"

using std::ofstream;
using std::vector;
using std::string;
using std::map;

std::mt19937 rand_gen{std::random_device{}()};

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser Splicer", "Splicer module of RNAInfuser");

    options.add_options()
        ("g,gtf",  "Path to gtf annotation file", cxxopts::value<string>())
        ("e,expression-table",  "Path to tab separated expression table (Formatted as transcript_id\\tcount)", cxxopts::value<string>())
        ("a,abundance-table",  "Path to tab separated abundance table (Formatted as transcript_id\\tcount\\tpm)", cxxopts::value<string>())
        ("use-whole-id", "Use whole transcript id instead of first 15 characters (ENSEMBL ids)", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("molecule-count", "Number of molecules to simulate, (requires abundance table", cxxopts::value<int>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("non-coding", "Process non-coding genes/transcripts as well", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("default-depth", "Default depth for transcripts that are not in expression table", cxxopts::value<int>()->default_value("0"))
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory = {"gtf", "output"};
    
    int missing_parameters = 0;
    for( string &param : mandatory){
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }
    if(missing_parameters  > 0){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    //parameters will be moved to argument parser when done
    string path_to_gtf   {args["gtf"].as<string>()};
    string out_path {args["output"].as<string>()};

    bool coding_only = args["non-coding"].as<bool>();

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);
    


    string buffer;

    bool coding_gene {true};
    bool coding_transcript {true};
    gtf *current_gene = nullptr;
    gtf *current_transcript = nullptr;

    int default_depth = args["default-depth"].as<int>();

    vector<pcr_copy> isoforms;


    map<string, int> transcript2count;
    map<string, double> transcript2tpm;


    std::cerr << "Reading expression\n";
    if(args.count("expression-table") > 0){
        ifstream table_file(args["expression-table"].as<string>());
        string tid = "BEG";
        int count;
        while(!table_file.eof()){
            table_file >> tid >> count;
            if(!args["use-whole-id"].as<bool>()){
                tid = tid.substr(0,15);
            }
            transcript2count[tid] = count;
        }
    }
    if(args.count("abundance-table") > 0){
        ifstream table_file(args["abundance-table"].as<string>());
        string tid = "BEG";
        double count;
        double tpm;
        table_file >> tid >> tid >> tid;
        while(!table_file.eof()){
            table_file >> tid >> count >> tpm;
            if(!args["use-whole-id"].as<bool>()){
                tid = tid.substr(0,15);
            }
            transcript2count[tid] = count;
            transcript2tpm[tid] = tpm;
        }
    
    }

    ifstream gtfile( path_to_gtf);


    std::uniform_real_distribution<> carry_over_random(0, 1.0);

    std::cerr << "Simulating\n";
    while(std::getline(gtfile, buffer)){
        if(buffer[0] == '#'){
            continue;
        }
        gtf *entry = new gtf(buffer);

        switch(entry->type){
            case gtf::entry_type::gene:
                coding_gene = (entry->info["gene_biotype"] == "protein_coding");
                if(coding_only && !coding_gene){
                    continue;
                }
                delete current_gene;
                current_gene = entry;
                break;
            case gtf::entry_type::transcript:
                coding_transcript = (entry->info["transcript_biotype"] == "protein_coding");
                if(coding_only && coding_transcript){
                    continue;
                }
                delete current_transcript;
                current_transcript = entry;
                isoforms.emplace_back(current_transcript->info["transcript_id"]);
                if(args["molecule-count"].count() == 0){//Empty block so I can define iter;
                    auto iter = transcript2count.find(isoforms.back().id);
                    if( iter!= transcript2count.end()){
                        isoforms.back().depth = iter->second;
                        isoforms.back().reversed = !current_transcript->plus_strand;
                    }else{
                        isoforms.back().depth = default_depth;
                        isoforms.back().reversed = !current_transcript->plus_strand;
                    }
                }
                else{
                    auto iter = transcript2tpm.find(isoforms.back().id);
                    double mc = args["molecule-count"].as<int>();
                    double tpm_multiplier = mc / 1000000;
                    if( iter!= transcript2tpm.end()){
                        isoforms.back().depth = iter->second * tpm_multiplier;
                        double real_carry = isoforms.back().depth - iter->second * tpm_multiplier;
                        if( real_carry > 0){
                            if(carry_over_random(rand_gen) < real_carry){
                                isoforms.back().depth ++ ;
                            }
                        }
                        isoforms.back().reversed = !current_transcript->plus_strand;
                    }else{
                        isoforms.back().depth = default_depth;
                        isoforms.back().reversed = !current_transcript->plus_strand;
                    }

                }
                break;
            case gtf::entry_type::exon:
                if(coding_only && !(coding_transcript && coding_gene)){
                    continue;
                }
                isoforms.back().segments.push_back(*entry);
                delete entry;
                break;
            default:
                break;
        }
    }
    
    ofstream outfile( out_path);
    for( const pcr_copy &pcp : isoforms){
        if(pcp.depth == 0){
            continue;
        }
        print_mdf(outfile, pcp.id, pcp, vector< vector<std::pair<int, char>>>{pcp.segments.size()});
    }

    return 0;
}
