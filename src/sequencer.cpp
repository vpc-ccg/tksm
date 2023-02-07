
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include <algorithm>
#include <thread>

#include <vector>


#include <cxxopts/cxxopts.hpp>
#include "fasta.h"
#include "util.h"
#include "interval.h"
#include "fasta.h"
#include "mdf.h"
#include "reverse_complement.h"

#define FMT_HEADER_ONLY 1
#include<fmt/format.h>

#include <future>

using std::ifstream;
using std::vector;
using std::string;



//PYTHONHASHSEED0 ../Badread/badread-runner.py simulate --reference=long-reads_temp/rnainfuser.simulation.batch-31.fasta --length 1000000,0 --seed 42 --quantity=4412980 --glitches=0,0,0 --junk_reads=0 --random_reads=0 --chimeras=0 > long-reads_temp/rnainfuser.simulation.batch-31.fastq 2> long-reads_temp/rnainfuser.simulation.batch-31.fastq.log=
int main(int argc, char **argv){
    cxxopts::Options options("RNAInfuser sequencer module", "Sequencer");
    const string software_name = "rnainfuser";
    options.add_options("RNAInfuser")
        ("i,input", "Molecule description file", cxxopts::value<string>())
        ("r,references",  "List of comma separated references, (required for fastq output)", cxxopts::value<vector<string>>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("temp", "Temp folder", cxxopts::value<string>()->default_value("."))
        ("n,name", "Project/Sample name", cxxopts::value<string>()->default_value("simulation"))
        ("badread", "path to badread-runner.py (default checks $PATH)", cxxopts::value<string>()->default_value("badread-runner.py"))
        ("t,threads", "Number of threads/batches", cxxopts::value<int>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("keep-temp", "Keep generated temp files", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("h,help", "Help screen")
    ;
    options.add_options("Badread")
        ("length", "Badread length distribution (Don't use it, truncate module should be used to modify read lengths)", cxxopts::value<string>()->default_value("1000000,0"))
        ("glitches", "Badread glitches parameter", cxxopts::value<string>()->default_value("0,0,0"))
        ("junk_reads", "Badread junk reads parameter", cxxopts::value<string>()->default_value("0"))
        ("random_reads", "Badread random reads parameter", cxxopts::value<string>()->default_value("0"))
        ("chimeras", "Badread chimeras parameter", cxxopts::value<string>()->default_value("0"))
        ("other-br", "Custom parameters to be passed to Badread", cxxopts::value<string>()->default_value(""))
    ;
    auto args = options.parse(argc, argv);

    vector<string> help_groups = {"RNAInfuser", "Badread"};

    if(args.count("help") > 0){
        std::cout << options.help(help_groups) << std::endl;
        return 0;
    }

    std::vector<string> mandatory =  {"input", "output", "references"};


    int missing_parameters = 0;
    for( string &param : mandatory){
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }

    if(missing_parameters  > 0){
        std::cerr << options.help(help_groups) << std::endl;
        return 1;
    }
  

    vector<string> batch_files;
    vector<int64_t> actual_throughputs;
    string input_filename {args["input"].as<string>()};

    std::ifstream infile(input_filename);
    vector<molecule_descriptor> molecules = parse_mdf(infile);
    map<string, string> refs;
    for(const string &path : args["references"].as<vector<string>>()){
        refs.merge(read_fasta_fast(path));
    }
    int64_t total_throughput = std::accumulate(molecules.begin(),molecules.end(),0L,
    [] (int64_t sum, const molecule_descriptor &pc) {
        for(const auto &ival : pc.cget_segments()){
            sum+=(ival.end - ival.start);
        }
        return sum;
    });
    int64_t throughput_per_batch = total_throughput / args["threads"].as<int>() + 1;
    int64_t throughput_so_far = throughput_per_batch;
    int current_batch = -1;
    std::ofstream ost;

    std::cerr << "Throughput per thread :" << throughput_per_batch << "\n";

    for( const molecule_descriptor& pc : molecules){
        pcr_molecule pcm{pc};
        for(const auto& seg : pc.cget_segments()){
            throughput_so_far += (seg.end - seg.start);
        }
        if(throughput_so_far > throughput_per_batch && current_batch + 1< args["t"].as<int>()){
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
        ost << ">" << pc.get_id()<< " depth=1\n";
        
        size_t seq_index = 0;
        for(const auto& seg : pc.cget_segments()){
            string seq = refs.at(seg.chr).substr(seg.start, seg.end - seg.start);
            for(const auto &pr : pc.errors_so_far){
                if((unsigned) pr.first > seq_index && (unsigned) pr.first < seq_index + seq.size()){
                    seq[ pr.first - seq_index] = pr.second;
                }
            }
            if(!seg.plus_strand){
                reverse_complement::complement_inplace(seq);
            }
            ost << seq << "\n";
            seq_index += seq.size();
        }
    }
    actual_throughputs.push_back(throughput_so_far);

    
    auto throughput_iter = actual_throughputs.begin();
    string other = args["other-br"].as<string>();
    std::cerr << other << "!\n";
    vector<string> br_cmd {{"junk_reads", "length", "glitches", "random_reads", "chimeras"}};

    vector<std::future<int>> promises;
    vector<int> return_values(args["t"].as<int>(),0);
    for(int i = 0; i< args["t"].as<int>(); ++i){
//    for(const string &bf : batch_files){
        const string &bf = batch_files[i];
        string batch_out_name = bf.substr(0,bf.find_last_of(".")) + ".fastq";
      
        for(string cmd : br_cmd){
            if (other.find(cmd) == string::npos){
                other = fmt::format("{} --{}={}", other, cmd, args[cmd].as<string>());
            }
        }
        string command = fmt::format( "PYTHONHASHSEED=0 {} simulate --reference {} --quantity {} --seed {} {} > {} 2> {}", 
                args["badread"].as<string>(),
                bf,
                std::to_string(*throughput_iter),
                args["seed"].as<int>(),
                other,
                batch_out_name,
                batch_out_name + ".log");
//        string command = "PYTHONHASHSEED=0 " + args["badread"].as<string>() + " simulate --reference=" + bf + 
//            " --length " + args["length"].as<string>() + " --seed=" + std::to_string(args["seed"].as<int>()) +
//            " --quantity=" + std::to_string(*throughput_iter) + " --glitches=" + args["glitches"].as<string>() + " --junk_reads="  + 
//            args["junk-reads"].as<string>() +  " --random_reads=" + args["random-reads"].as<string>() +
//            " --chimeras=" + args["chimeras"].as<string>() + " > " + batch_out_name + " 2> " + batch_out_name +".log";
        std::cout << command << "\n\n";
        ++throughput_iter;

        promises.push_back(std::async([i](const std::string &command) -> int{
            int ret = std::system(command.c_str());

            if(WEXITSTATUS(ret)){ //Didn't exit normally
                return i;
            }
            return 0;
        },command));
    }

    for(int i = 0; i< args["t"].as<int>(); ++i){
    //for( std::thread &t : threads){
        int ret = promises[i].get();

        if(ret){
            std::cerr << "Error: Please check following log file\n";
            std::cerr <<  batch_files[i].substr(0,batch_files[i].find_last_of(".")) + ".fastq.log\n";
            return -1;
        }
    }
    
    string cat_command = "cat ";
    for(const string &bf : batch_files){
        string batch_out_name = bf.substr(0,bf.find_last_of(".")) + ".fastq";
        cat_command += (batch_out_name + " ");
    }
    cat_command += " > " + args["output"].as<string>();
    std::system(cat_command.c_str());

    if(!args["keep-temp"].as<bool>()){
        for(const string &bf : batch_files){

            string batch_out_name = bf.substr(0,bf.find_last_of(".")) + ".fastq";
            string rm_command = "rm " + batch_out_name;
            system(rm_command.c_str());
        }
    }
    return 0;
}
