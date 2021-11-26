
#include <string>
#include <map>
#include <vector>
#include <random>
#include <fstream>
#include <numeric>
#include <functional>
#include <algorithm>

#include "fasta.h"
#include "tree.h"
#include "interval.h"
#include "reverse_complement.h"
#include "extern/cxxopts.hpp"

using std::string;
using std::map;
using std::vector;

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cerr;

//Random number generator, seed is set in the main function
std::mt19937 rand_gen{std::random_device{}()};

void run_pcr( 
        const vector<isoform> &isoforms,
        const map<string, string> &sequences,
        std::ofstream &ost,
        int cycles,
        double pcr_duplication_rate,
        double error_rate,
        int number_of_target_reads,
        double random_pairing_rate_per_cycle){

    vector< vector< pcr_molecule>> saved_for_rp(cycles);

    using mut_tree = tree<int, vector<std::pair< int, char >>>;

//Calculate number of molecules to be produced with pcr
    int64_t molecule_count = std::accumulate(isoforms.begin(), isoforms.end(), 0, [] ( int64_t sum, const isoform &vec) -> int64_t{
            return sum + vec.depth;
            });

//Calculate expected number of molecules after pcr
    int64_t expected_number_after_pcr = std::pow((1 + pcr_duplication_rate), cycles) * molecule_count;

    double drop_ratio = static_cast<double>(number_of_target_reads) / expected_number_after_pcr;
    
//Convert error_rate to random_error_rate where characters can mutate into same chars (A->A or T->T etc.)
    error_rate = error_rate * 4.0 / 3.0; 

    std::uniform_real_distribution<> z1dist(0,1);
    std::uniform_int_distribution<> basedist(0,3);
    char bases[] = {'A','C','T','G'};

    //PCR lambda function
    std::function<void(const pcr_molecule&, mut_tree &,int, int, vector<int>&)> do_pcr = 
        [&do_pcr, &z1dist, pcr_duplication_rate, &bases, &basedist, random_pairing_rate_per_cycle, &saved_for_rp, cycles]
        (const pcr_molecule &iso, mut_tree &mt, int current_cycle, int expected_mutation_count, vector<int> &positions) mutable -> void{
        for(int cyc = current_cycle + 1; cyc <= cycles; ++cyc){
            if(z1dist(rand_gen) < pcr_duplication_rate){ // Molecule is caught in this cycle of pcr
                                                         // GC bias can be introduced here
                std::shuffle(positions.begin(),positions.end(), rand_gen);
                vector<std::pair<int, char>> mutations;
                for(int j = 0; j < expected_mutation_count; ++j){ //TODO: We can make this a normal distribution    
                    mutations.push_back(std::make_pair(positions[j], bases[basedist(rand_gen)]));
                }
                if(z1dist(rand_gen) < random_pairing_rate_per_cycle){
                    vector<std::pair<int,char>> mutations_so_far;
                    mut_tree *current = &mt;
                    while(current !=nullptr){
                        for( auto pr : current->data){
                            mutations_so_far.push_back(pr);
                        }
                        current = current->parent;
                    }
                    pcr_molecule other(iso);
                    other.paired[0].errors_so_far = mutations_so_far;
                    saved_for_rp[current_cycle].push_back(other);
                }
                else{
                    mt.add_child(cyc, mutations);
                    do_pcr(iso, mt[cyc], cyc, expected_mutation_count, positions);
                }
            }
        }
    };

    //Makes a printer lambda function to be run on the tree
    auto make_pcr_printer = [&z1dist, drop_ratio] (ostream &ost, string *_base_id, int copy_number, const map<string, string> *sequences, const pcr_molecule *_pm){
        return  [&ost, &z1dist, drop_ratio, _base_id, copy_number, sequences, _pm] (int depth, const mut_tree *mt) mutable -> void{
            string base_id{*_base_id};
            const pcr_molecule &pm {*_pm};
            if( z1dist(rand_gen) < drop_ratio) { //Molecule is captured by sequencing
                ost << ">" << base_id << "_" << copy_number << " ";
                mut_tree *current = mt->parent;
                ost << "lineage:";
                ost << mt->identity;
                while( current != nullptr){ //Print pcr Lineage
                    ost << "_" << current->identity;
                    current = current->parent;
                }
                ost << "\n";

                for( const pcr_copy &pcp : pm.paired){
                    unsigned seq_index = 0;
                    const isoform &iso = pcp.iso;
                    for(const exon &e : iso.segments){
                        if( sequences->find(e.chr) == sequences->end()){
                            cerr << e << "!\n";
                            throw std::runtime_error("Chr not in fasta");
                        }
                        if( (unsigned) e.end > sequences->at(e.chr).size()){
                            cerr << e << "!\n";
                            throw std::runtime_error("Exon is outside of the chr");
                        }
                        string seq = sequences->at(e.chr).substr( e.start, e.end - e.start); // Maybe string view in the future?
                       
                        for( const std::pair<int, char> &pr : pcp.errors_so_far){//Errors so far, only used for random_pairings
                            if((unsigned) pr.first > seq_index && (unsigned) pr.first < seq_index + seq.size()){
                                seq[ pr.first - seq_index] = pr.second;
                            }
                        }
                        std::function<void(const mut_tree *, string&)> apply_mutation = [&apply_mutation, seq_index] ( const mut_tree *mt, string &seq) mutable -> void{
                            if(mt->parent != nullptr){
                                apply_mutation(mt->parent, seq);
                            }
                            for(const auto &pr: mt->data){
                                if((unsigned) pr.first > seq_index && (unsigned) pr.first < seq_index + seq.size()){
                                    seq[ pr.first - seq_index] = pr.second;
                                }
                            }
                        };

                        apply_mutation(mt, seq);
                        if(! e.plus_strand){
                            reverse_complement::complement_inplace(seq);
                        }

#ifdef DEBUG
                        ost << e << "\t";
#endif
                        ost << seq << "\n";
                        seq_index += seq.size();
                    }
                }
            }
        };
    };
    for( const isoform &iso : isoforms){
        string base_id = iso.gene;
        int molecule_size = 0;
        for(const exon &e : iso.segments){
            molecule_size += (e.end - e.start);
        }

        for( int i = 0; i < iso.depth; ++i){
            pcr_molecule pcm{iso};
            mut_tree mutation_tree;
            vector<int> positions(molecule_size);
            std::iota(positions.begin(), positions.end(), 0);
            int expected_mutation_count = error_rate * molecule_size; 

            do_pcr(pcm, mutation_tree, 0, expected_mutation_count, positions);
            auto printer = make_pcr_printer(ost, &base_id, i, &sequences, &pcm);
            mutation_tree.df_execute2(printer);
        }

    }
    int cycle = 1;

    vector<pcr_molecule> random_pairings;
    for( auto vec_iter = saved_for_rp.begin(); vec_iter !=saved_for_rp.end(); ++vec_iter){
        vector<pcr_molecule> &pc_vec = *vec_iter; 
        //After shuffle, pair adjacent elements
        std::shuffle(pc_vec.begin(), pc_vec.end(), rand_gen);
        for(auto iter = pc_vec.begin(); iter + 1 != pc_vec.end() && iter !=pc_vec.end(); iter= std::next(iter,2)){
            auto next_iter = std::next(iter);
            pcr_molecule rp(*iter, *next_iter);
            int molecule_size = 0;
            for( const pcr_copy &cp :rp.paired){
                for(const exon &e : cp.iso.segments){
                    molecule_size += (e.end - e.start);
                }
            }
            mut_tree mutation_tree;
            vector<int> positions(molecule_size);
            std::iota(positions.begin(), positions.end(), 0);
            int expected_mutation_count = error_rate * molecule_size; 

            do_pcr(rp, mutation_tree, cycle, expected_mutation_count, positions);
            string base_id = "RP";
            for( const pcr_copy &pc : rp.paired){
                base_id += "_";
                base_id += pc.iso.gene;
            }
            auto printer = make_pcr_printer(ost, &base_id, cycle, &sequences, &rp);
            mutation_tree.df_execute2(printer);
        }
        
        //Continue PCR
        ++cycle;
    }
}


int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser PCR module", "PCR PCR PCR PTR");

    options.add_options()
        ("m,molecule-description", "Molecule description file", cxxopts::value<string>())
        ("r,references",  "List of comma separated references", cxxopts::value<vector<string>>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("read-count", "Number of reads to simulate", cxxopts::value<int>()->default_value("100000"))

        ("pcr", "Number of pcr cycles to simulate", cxxopts::value<int>()->default_value("6"))
        ("pcr-efficiency", "Probability of a molecule being duplicated during each pcr cycle", cxxopts::value<double>()->default_value("0.75"))
        ("pcr-error-rate", "Probability of substition errors for each base during each pcr cycle", cxxopts::value<double>()->default_value("0.01"))
        ("pcr-random-pairing-rate", "Probability of molecules randomly attaching to each other during each pcr cycle", cxxopts::value<double>()->default_value("0.0001"))

        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
    ;

    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory {{"molecule-description", "references", "output"}};

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
    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);

    ifstream md_file(args["molecule-description"].as<string>());
    string buffer;
    buffer.reserve(1000);

    vector<isoform> isoforms;
    while(std::getline(md_file, buffer)){
        auto fields = rsplit(buffer, "\t");
        string id {fields[0]};
        string gene{fields[1]};
        int depth{stoi(fields[2])};
        int exon_count{stoi(fields[3])};
        vector<exon> exons;
        for(int i = 0; i< exon_count; ++i){
            std::getline(md_file, buffer);
            auto fields = rsplit(buffer, "\t");
            string chr{fields[0]};
            int start {stoi(fields[1])};
            int end   {stoi(fields[2])};
            string strand{fields[3]};
            exons.emplace_back(chr, start, end, strand, gene + "_" + std::to_string(i), nullptr);
        }
        isoforms.emplace_back(exons, depth, gene);
    }

    map<string, string> sequences;
    for(const string &p : args["references"].as<vector<string>>()){
        sequences.merge(read_fasta_fast(p));
    }

    ofstream out_stream(args["output"].as<string>());

    int cycles = args["pcr"].as<int>();
    double pcr_duplication_rate = args["pcr-efficiency"].as<double>();
    double error_rate = args["pcr-error-rate"].as<double>();
    int number_of_target_reads = args["read-count"].as<int>();
    double random_pairing_rate_per_cycle = args["pcr-random-pairing-rate"].as<double>();

    run_pcr( 
        isoforms,
        sequences,
        out_stream,
        cycles,
        pcr_duplication_rate,
        error_rate,
        number_of_target_reads,
        random_pairing_rate_per_cycle);



    return 0;
}
