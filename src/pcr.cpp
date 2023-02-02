
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
#include <cxxopts/cxxopts.hpp>
#include <tuple>
#include "mdf.h"

using std::string;
using std::map;
using std::vector;
using std::tuple;

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cerr;


//Random number generator, seed is set in the main function
std::mt19937 rand_gen{std::random_device{}()};
/*
void run_pcr( 
        const vector<molecule_descriptor> &molecules,
        const map<string, string> &sequences,
        std::ofstream &ost,
        int cycles,
        double pcr_duplication_rate,
        double error_rate,
        int number_of_target_reads,
        double random_pairing_rate_per_cycle,
        bool fastq_output){

    vector< vector< pcr_molecule>> saved_for_rp(cycles);

    using mut_tree = tree<int, vector<std::pair< int, char >>>;

//Calculate number of molecules to be produced with pcr
    int64_t molecule_count = std::accumulate(molecules.begin(), molecules.end(), 0, [] ( int64_t sum, const molecule_descriptor &cpy) -> int64_t{
            return sum + cpy.get_depth();
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
        (const pcr_molecule &iso, mut_tree &mt, int current_cycle, double expected_mutation_count, vector<int> &positions) mutable -> void{
        for(int cyc = current_cycle + 1; cyc <= cycles; ++cyc){
            if(z1dist(rand_gen) < pcr_duplication_rate){ // Molecule is caught in this cycle of pcr
                                                         // GC bias can be introduced here
                std::shuffle(positions.begin(),positions.end(), rand_gen);
                vector<std::pair<int, char>> mutations;
                int actual_mutation_count = int(expected_mutation_count);

                double mutation_int_diff = expected_mutation_count - int(expected_mutation_count);
                if(mutation_int_diff < 1){
                    actual_mutation_count += (z1dist(rand_gen) < mutation_int_diff?1:0);
                }

                for(int j = 0; j < actual_mutation_count; ++j){ //TODO: We can make this a normal distribution    
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

    auto make_pcr_mdbed_printer = [&z1dist, drop_ratio] (ostream &ost, string *_base_id, int copy_number, const pcr_molecule *_pm){

        return  [&ost, &z1dist, drop_ratio, _base_id, _pm] (int depth, const mut_tree *mt) mutable -> void{
            string base_id{*_base_id};
            const pcr_molecule &pm {*_pm};
            size_t molecule_count = std::accumulate(pm.paired.begin(), pm.paired.end(), 0, [] ( size_t sum, const molecule_descriptor &cpy) -> int64_t{
                    return sum + cpy.cget_segments().size();
                    });
            if( z1dist(rand_gen) < drop_ratio) { //Molecule is captured by sequencing
                vector<vector<std::pair<int, char>>> errors_per_segment {molecule_count};
                const mut_tree *current = mt;
                string lineage = std::to_string(mt->identity);
                while( current->parent != nullptr){ //Print pcr Lineage
                    current = current->parent;
                    lineage += ( "_" + std::to_string(current->identity));
                }

                interval so_far {0, 0};
                int copy_end = 0;
                for( const molecule_descriptor &pcp : pm.paired){
                    int interval_counter = 0;
                    for(const ginterval &ival : pcp.cget_segments()){
                        so_far = interval{so_far.end, so_far.end +(ival.end - ival.start)};
                        for( const std::pair<int, char> &pr : pcp.errors_so_far){//Errors so far, only used for random_pairings
                            if(so_far.contains(pr.first)){
                                errors_per_segment[interval_counter].push_back(pr);
                            }
                        }
                        interval_counter++;
                    }
                    so_far = interval{copy_end,copy_end};
                    interval_counter = 0;
                    for(const ginterval &ival : pcp.cget_segments()){
                        so_far = interval{so_far.end, so_far.end +(ival.end - ival.start)};
                        std::function<void(const mut_tree *)> apply_mutation = [&errors_per_segment, &apply_mutation, interval_counter, so_far] ( const mut_tree *mt) mutable -> void{
                            if(mt->parent != nullptr){
                                apply_mutation(mt->parent);
                            }
                            for(const auto &pr: mt->data){
                                if(so_far.contains(pr.first)){
                                    errors_per_segment[interval_counter].push_back(pr);
                                }
                            }
                        };
                        apply_mutation(mt);
                        copy_end = ival.end;
                        interval_counter++;
                    }

                }

                print_mdf(ost, base_id+"_"+lineage, pm, errors_per_segment);
            }
        };
    };
    //Makes a prnter lambda function to be run on the tree
    auto make_pcr_fastq_printer = [&z1dist, drop_ratio] (ostream &ost, string *_base_id, int copy_number, const map<string, string> *sequences, const pcr_molecule *_pm){
        return  [&ost, &z1dist, drop_ratio, _base_id, copy_number, sequences, _pm] (int depth, const mut_tree *mt) mutable -> void{
            string base_id{*_base_id};
            const pcr_molecule &pm {*_pm};
            if( z1dist(rand_gen) < drop_ratio) { //Molecule is captured by sequencing
                ost << ">" << base_id << "_" << copy_number << "_";
                mut_tree *current = mt->parent;
//                ost << "lineage:";
                ost << mt->identity;
                while( current != nullptr){ //Print pcr Lineage
                    ost << "_" << current->identity;
                    current = current->parent;
                }
                ost << "\n";

                for( const molecule_descriptor &pcp : pm.paired){
                    unsigned seq_index = 0;
                    for(const ginterval &ival : pcp.cget_segments()){
                        if( sequences->find(ival.chr) == sequences->end()){
                            cerr << ival << "!\n";
                            throw std::runtime_error("Chr not in fasta");
                        }
                        if( (unsigned) ival.end > sequences->at(ival.chr).size()){
                            cerr << ival << "!\n";
                            throw std::runtime_error("Exon is outside of the chr");
                        }
                        string seq = sequences->at(ival.chr).substr( ival.start, ival.end - ival.start); // Maybe string view in the future?
                       
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
                        if(! ival.plus_strand){
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
    for( const molecule_descriptor &pcp : molecules){
        string base_id = pcp.get_id();
        int molecule_size = 0;
        for(const ginterval &ival : pcp.cget_segments()){
            molecule_size += (ival.end - ival.start);
        }

        for( int i = 0; i < pcp.get_depth(); ++i){
            pcr_molecule pcm{pcp};
            pcm.paired.front().depth (1) ;
            mut_tree mutation_tree;
            vector<int> positions(molecule_size);
            std::iota(positions.begin(), positions.end(), 0);
            double expected_mutation_count = error_rate * molecule_size; 


            do_pcr(pcm, mutation_tree, 0, expected_mutation_count, positions);
            if(fastq_output){
                auto printer = make_pcr_fastq_printer(ost, &base_id, i, &sequences, &pcm);
                mutation_tree.df_execute2(printer);
            }
            else{
                auto printer = make_pcr_mdbed_printer(ost, &base_id, i, &pcm);
                mutation_tree.df_execute2(printer);
            }
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
            for( const molecule_descriptor &cp :rp.paired){
                for(const ginterval &ival : cp.cget_segments()){
                    molecule_size += (ival.end - ival.start);
                }
            }
            mut_tree mutation_tree;
            vector<int> positions(molecule_size);
            std::iota(positions.begin(), positions.end(), 0);
            double expected_mutation_count = error_rate * molecule_size; 

            do_pcr(rp, mutation_tree, cycle, expected_mutation_count, positions);
            string base_id = "RP";
            for( const molecule_descriptor &pc : rp.paired){
                base_id += "_";
                base_id += pc.get_id();
            }

            if(fastq_output){
                auto printer = make_pcr_fastq_printer(ost, &base_id, cycle, &sequences, &rp);
                mutation_tree.df_execute2(printer);
            }
            else{
                auto printer = make_pcr_mdbed_printer(ost, &base_id, cycle, &rp);
                mutation_tree.df_execute2(printer);
            }
        }
        
        //Continue PCR
        ++cycle;
    }
}
*/

class PCRRunner{
    int     cycles;
    double  efficiency;
    double  error_rate;
    using mut_tree = tree<int, vector<std::pair< int, char >>>;


    std::uniform_real_distribution<> z1dist  {0,1};
    std::uniform_int_distribution<> basedist {0,3};
    char bases[4] = {'A','C','T','G'};

    public:
    PCRRunner( int cycles, double efficiency, double error_rate) :
        cycles{cycles},
        efficiency{efficiency},
        error_rate{(4*error_rate)/3} //Our error definition is a base to be assigned to random base not changing to other bases so we are adjusting the probability
    {}

    void do_pcr (ostream &ost, const molecule_descriptor& md, int step, const vector<int> &positions, vector<std::pair<int, char>> mutations, double drop_ratio) {
            if( z1dist(rand_gen) > efficiency) { //Molecule is not captured by PCR
                return;
            }
            double expected_mutation_count = error_rate * positions.size();

            vector<int> mutation_pos;
            std::sample(positions.begin(), positions.end(), std::back_inserter(mutation_pos), (int)expected_mutation_count, rand_gen);
            

            for(int pos : mutation_pos){ //TODO: We can make this a normal distribution    
                mutations.push_back(std::make_pair(pos, bases[basedist(rand_gen)]));
            }


            molecule_descriptor mdc{md};
            mdc.update_errors(mutations);
            mdc.id(mdc._id + "." + std::to_string(step));
            if ( z1dist(rand_gen) > drop_ratio){ // Molecule is captured by the sequencing
                ost << mdc;
            }
            for(int cycle = step+1; cycle < cycles; ++cycle){
                do_pcr(ost, mdc, cycle, positions, mutations, drop_ratio);
            }
    } 


    int run(const vector<molecule_descriptor> &molecules, ostream &ost, int number_of_target_reads){
//Calculate number of molecules to be produced with pcr
        int64_t molecule_count = std::accumulate(molecules.begin(), molecules.end(), 0,
                [] ( int64_t sum, const molecule_descriptor &cpy) -> int64_t{
                return sum + cpy.get_depth();
        });

//Calculate expected number of molecules after pcr
        int64_t expected_number_after_pcr = std::pow((1 + efficiency), cycles) * molecule_count;

        double drop_ratio = static_cast<double>(number_of_target_reads) / expected_number_after_pcr;

        for( const molecule_descriptor &pcp : molecules){
            string base_id = pcp.get_id();
            vector<int> positions(pcp.size());
            std::iota(positions.begin(), positions.end(), 0);

            for( int i = 0; i < pcp.get_depth(); ++i){
                mut_tree mutation_tree;
                molecule_descriptor pc_c{pcp};
                pc_c.id(pc_c.get_id() + "_" + std::to_string(i));
                for(int cycle = 0; cycle < cycles; ++cycle){
                    do_pcr(ost, pc_c, cycle, positions, vector<std::pair<int, char>>{}, drop_ratio);
                }
            }

        }
        return 0;
    }
};
int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser PCR module", "PCR PCR PCR PTR");
    
    map<string, tuple<double, double>> pcr_presets {
        { "Taq-setting1", {2*std::pow(0.1,4),0.88}},
        { "Taq-setting2", {7.2*std::pow(0.1,5),0.36}},
        { "Klenow", {1.3*std::pow(0.1,4),0.80}},
        { "T7", {3.4*std::pow(0.1,5),0.90}},
        { "T4", {3.0*std::pow(0.1,6),0.56}},
        { "Vent", {4.5*std::pow(0.1,5),0.70}},
    };

    string preset_string = "presets (Cha, R. S., & Thilly, W. G. (1993). Specificity, efficiency, and fidelity of PCR. Genome Research, 3(3), S18-S29.)\n";
    preset_string.reserve(500);
    for(const auto &p : pcr_presets){
        preset_string+= ("- " +p.first + ": " + "efficiency: " + std::to_string(std::get<1>(p.second)) +", error-rate:" + std::to_string(std::get<0>(p.second)) + "\n");
    }

    options.add_options()
        ("i,molecule-description", "Molecule description file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("read-count", "Number of reads to simulate", cxxopts::value<size_t>()->default_value("100000"))
        ("cycles", "Number of pcr cycles to simulate", cxxopts::value<int>()->default_value("6"))
        ("efficiency", "Probability of a molecule being duplicated during each pcr cycle", cxxopts::value<double>()->default_value("0.75"))
        ("error-rate", "Probability of substition errors for each base during each pcr cycle", cxxopts::value<double>()->default_value("0.000001"))
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
        ("x,preset", preset_string, cxxopts::value<string>())
    ;

    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }

    std::vector<string> mandatory =  {"molecule-description", "output"};
    if(args.count("fastq") > 0 && args["fastq"].as<bool>()){
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
    
    int cycles = args["cycles"].as<int>();
    double pcr_efficiency = args["efficiency"].as<double>();
    double error_rate = args["error-rate"].as<double>();
    int number_of_target_reads = args["read-count"].as<size_t>();

    if(args["preset"].count() > 0){
        if (pcr_presets.find(args["preset"].as<string>()) == pcr_presets.end()){
            std::cerr << "Preset doesn't exist!\n";
            return -1;
        }
        tuple<double, double> setting = pcr_presets[args["preset"].as<string>()];

        error_rate = std::get<0>(setting);
        pcr_efficiency = std::get<1>(setting);
    }

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);

    ifstream md_file(args["molecule-description"].as<string>());
    vector<molecule_descriptor> molecules = parse_mdf(md_file);

    if( molecules.size() > 2 * args["read-count"].as<size_t>()){
        std::shuffle(molecules.begin(), molecules.end(), rand_gen);
        molecules.resize( 2 * args["read-count"].as<size_t>());
    }


    ofstream out_stream(args["output"].as<string>());
    return PCRRunner{cycles, pcr_efficiency, error_rate}.run(molecules, out_stream, number_of_target_reads);
}


