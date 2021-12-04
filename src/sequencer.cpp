
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include <algorithm>
#include <thread>

#include <vector>


#include "extern/cxxopts.hpp"
#include "fasta.h"
#include "util.h"

using std::ifstream;
using std::vector;
using std::string;

int main(int argc, char **argv){
    cxxopts::Options options("RNAInfuser sequencer module", "Sequencer");
    const string software_name = "rnainfuser";
    options.add_options()
        ("m,molecule-description", "Molecule description file", cxxopts::value<string>())
        ("r,references",  "List of comma separated references, (required for fastq output)", cxxopts::value<vector<string>>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("temp", "Temp folder", cxxopts::value<string>()->default_value("."))
        ("n,name", "Project/Sample name", cxxopts::value<string>()->default_value("simulation"))
        ("badread", "path to badread-runner.py (default checks $PATH)", cxxopts::value<string>()->default_value("badread-runner.py"))
        ("t,threads", "Number of threads/batches", cxxopts::value<int>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("fasta", "Input is fasta", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("h,help", "Help screen")
    ;

    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }

    std::vector<string> mandatory =  {"molecule-description", "output"};
    if(args.count("fasta") == 0){
        mandatory.push_back("references");
    }

    int missing_parameters = 0;
    for( string &param : mandatory){
        std::cout << param << "\n";
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }

    if(missing_parameters  > 0){
        std::cerr << options.help() << std::endl;
        return 1;
    }
  

    vector<string> batch_files;
    vector<int64_t> actual_throughputs;
    string input_filename {args["m"].as<string>()};
    if(args.count("fasta") > 0 ){ //Input is fasta
        map<string, string> seqs = read_fasta_fast(input_filename);
        int64_t total_throughput = std::accumulate(seqs.begin(),seqs.end(),0, [] (int64_t sum, const std::pair<string, string> &sp) { return sum + sp.second.size();});
        int64_t throughput_per_batch = total_throughput / args["threads"].as<int>() + 1;
        int64_t throughput_so_far = throughput_per_batch;
        int current_batch = -1;
        std::ofstream ost;
        for( const auto &sp : seqs){
            int64_t seq_len = sp.second.size();
            if(throughput_so_far + seq_len > throughput_per_batch && current_batch + 1< args["t"].as<int>()){
                if(current_batch >= 0){
                    actual_throughputs.push_back(throughput_so_far);

                }
                ++current_batch;

                throughput_so_far = 0;
                string batch_file = args["temp"].as<string>() + "/" + software_name + "." +
                    args["name"].as<string>() + ".batch-" + std::to_string(current_batch) + ".fasta";
                batch_files.push_back(batch_file);
                ost.close();
                ost.open(batch_file);
            }
            ost << ">" << sp.first << "\n" << sp.second << "\n";
            throughput_so_far += seq_len;
        }
        actual_throughputs.push_back(throughput_so_far);
    }
    else{//Generate fasta from mdf

    }
    
    auto throughput_iter = actual_throughputs.begin();

    vector<std::thread> threads;
    for(const string &bf : batch_files){
        string batch_out_name = bf.substr(0,bf.find_last_of(".")) + ".fastq";
        string command = "PYTHONHASHSEED=0 " + args["badread"].as<string>() + " simulate --reference=" + bf + " --length 1000000,0 --seed 42" +
            " --quantity=" + std::to_string(*throughput_iter) + " --glitches=0,0,0 --junk_reads=0 --random_reads=0 --chimeras=0" + " > " + batch_out_name + " 2> " + batch_out_name +".log";
        std::cout << command << "\n";
        ++throughput_iter;
        threads.push_back(std::thread{[](const std::string &command){
                std::system(command.c_str());
        },command});
    }

    for( std::thread &t : threads){
        t.join();
    }
    
    string cat_command = "cat ";
    for(const string &bf : batch_files){
        string batch_out_name = bf.substr(0,bf.find_last_of(".")) + ".fastq";
        cat_command += (batch_out_name + " ");
    }
    cat_command += " > " + args["output"].as<string>();
    std::system(cat_command.c_str());

    for(const string &bf : batch_files){

        string batch_out_name = bf.substr(0,bf.find_last_of(".")) + ".fastq";
        string rm_command = "rm " + batch_out_name;
        system(rm_command.c_str());
    }
    return 0;
}
